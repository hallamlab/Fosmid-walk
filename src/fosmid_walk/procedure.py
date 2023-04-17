from Bio import SeqIO
import uuid
from pathlib import Path
import json

import subprocess
import os

from numpy import record

# Avery, Kat, and Tony, 2023

def estimate(
        sample: Path, out_dir: Path,
        backbone: Path=Path('/app/pcc1.fasta'), temp_dir: Path=Path('/ws'),
        threads: int|None=None,
    ):

    BACKBONE = "backbone"
    BACKBONE_HEADER = BACKBONE
    sample = sample.absolute()
    backbone = backbone.absolute()
    sample_name = ".".join(sample.name.split(".")[:-1])
    workspace = temp_dir.joinpath(f"{sample_name}-{uuid.uuid4().hex}")

    os.makedirs(workspace, exist_ok=True)
    os.chdir(workspace)

    try:
        # Copies fosmid backbone senquence to new directory
        BACKBONE_SIGNATURE = None # the first 7 nucleotides of backbone
        with open(backbone) as f:
            _ = f.readline() # original header
            with open(workspace.joinpath(f"{BACKBONE}.fasta"), 'w') as new_f:
                new_f.write(f">{BACKBONE_HEADER}\n") # make header predictable
                l = f.readline()
                BACKBONE_SIGNATURE = l[:7]
                new_f.write(l)
                for l in f: new_f.write(l)
        subprocess.call([str('bwa index '+f'{BACKBONE}.fasta')], shell = True)

        # Aligns reads to vector backbone ***PATH TO FILES UPDATED HERE***
        os.symlink(sample, f"{sample_name}.fastq")
        subprocess.call([str('bwa mem '+f'{BACKBONE}.fasta '+sample_name+'.fastq > '+f'{BACKBONE}_aln.sam')], shell = True)
        subprocess.call([str('samtools sort '+f'{BACKBONE}_aln.sam > '+f'{BACKBONE}_aln.bam')], shell = True)
        subprocess.call([str('samtools index '+f'{BACKBONE}_aln.bam')], shell = True)
        subprocess.call([str('samtools view '+f'{BACKBONE}_aln.bam {BACKBONE_HEADER}:0-20 -o '+f'{BACKBONE}-5-020.bam')], shell = True)
        subprocess.call([str('samtools view -F 0x10 '+f'{BACKBONE}-5-020.bam -o '+'020-fwd.bam')], shell = True)
        subprocess.call([str('samtools view -f 0x10 '+f'{BACKBONE}-5-020.bam -o '+'020-rev.bam')], shell = True)
        subprocess.call([str('samtools fasta '+'020-fwd.bam > '+'020-fwd.fasta')], shell = True)
        subprocess.call([str('samtools fasta '+'020-rev.bam > '+'020-rev.fasta')], shell = True)

        fwd = SeqIO.parse(str('020-fwd.fasta'), "fasta")
        recs = []
        for f in fwd:
            recs.append((True, f))

        rev = SeqIO.parse(str('020-rev.fasta'), "fasta")
        for r in rev:
            rev_r = r.seq.reverse_complement()
            r.seq = rev_r
            recs.append((False, r))

        kept_indicies = []
        final_recs = []
        hit_lens = []
        CUT = 100
        for i, (is_fwd, rec) in enumerate(recs):
            if BACKBONE_SIGNATURE not in rec.seq: continue
            new_seq = rec.seq[:rec.seq.index(BACKBONE_SIGNATURE)]
            hit_lens.append(len(new_seq))
            if len(new_seq) >= CUT:
                trim_seq = new_seq[len(new_seq)-CUT:len(new_seq)]
                final_recs.append((is_fwd, trim_seq))
                kept_indicies.append(i)
        assert len(final_recs)>0, "no hits to vector backbone"
        with open(f'{BACKBONE}-5-020.fasta', 'w') as f:
            for i, (is_fwd, seq) in enumerate(final_recs):
                f.write(f">{i}\n{seq}\n")

        with open(out_dir.joinpath("stats.json"), "w") as j:
            json.dump({
                "hit_count": len(hit_lens),
                "hits_kept": len(final_recs),
                "hit_lengths": hit_lens,
            }, j)

        uclust_cmd = [
            "usearch",
            "-cluster_fast", f"{BACKBONE}-5-020.fasta",
            "-uc", f"{BACKBONE}-5-020.uc",
            "-id", str(0.9),
            "-consout", f"{BACKBONE}-5-020_clusters.fasta",
            "-sizeout",
            "-centroids", f"{BACKBONE}-5-020_centroids.fasta",
            "-minsize", str(2),
        ]
        if threads is not None:
            print(f"\nusing {threads} threads for usearch")
            uclust_cmd+= ["-threads", threads]
        subprocess.call(' '.join(uclust_cmd), shell = True)

        os.makedirs(out_dir, exist_ok=True)
        final_output =  open(out_dir.joinpath(f'{sample_name}_full.fasta'), 'w')
        hits =          open(out_dir.joinpath(f'{sample_name}_hits.fasta'), 'w')
        try:
            centroids = SeqIO.parse(str(f'{BACKBONE}-5-020_centroids.fasta'), "fasta")
            size_ones, total = 0, 0
            for c in centroids:
                id, size, _ = c.id.split(";")
                id = int(id)
                rev, original = recs[kept_indicies[id]]
                header = f">{original.id};{'forward' if not rev else 'reverse_compliment'};{size}"
                final_output.write(f"{header}\n{original.seq}\n")
                hits.write(f"{header}\n{c.seq}\n")

                total += 1
                if size == "size=1":
                    size_ones += 1

            total = size_ones + total
            print(f"""\
                {total} clusters resolved and {size_ones} had only 1 member.
                This leaves {total-size_ones} results
            """.replace("  ", ""))
        finally:
            final_output.close()
            hits.close()

    finally:
        # shutil.rmtree(Path("./").absolute())
        pass
