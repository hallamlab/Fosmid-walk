import shutil
from Bio import SeqIO
import uuid
from pathlib import Path

import subprocess
import os
import re
import csv

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
        rev = SeqIO.parse(str('020-rev.fasta'), "fasta")
        recs = []
        for f in fwd:
            recs.append(f)

        for r in rev:
            rev_r = r.seq.reverse_complement()
            r.seq = rev_r
            recs.append(r)

        final_recs = []
        for rec in recs:
            if BACKBONE_SIGNATURE in rec.seq: # first n nt in vector backbone
                new_seq = rec.seq[:str(rec.seq).index(BACKBONE_SIGNATURE)]
                if len(new_seq) >= 100:
                    trim_seq = new_seq[len(new_seq)-100:len(new_seq)]
                    rec.seq = trim_seq
                    final_recs.append(rec)
        assert len(final_recs)>0, "no hits to vector backbone"
        SeqIO.write(final_recs, str(f'{BACKBONE}-5-020.fasta'), "fasta")

        uclust_cmd = ["usearch",
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

        final_output = []
        clusters = SeqIO.parse(str(f'{BACKBONE}-5-020_clusters.fasta'), "fasta")
        size_ones = 0
        for c in clusters:
            if "size=1;" not in c.id:
                cluster = re.findall(r'\w+', c.id)[0]
                size = re.findall(r'(?<=size\=)\w+', c.id)[0]
                output_row = [cluster, size, c.seq]
                final_output.append(output_row)
            else:
                size_ones += 1

        total = size_ones + len(final_output)
        print(f"""\
            {total} clusters resolved, but {size_ones} had only 1 member.
            This leaves {total-size_ones} results
        """.replace("  ", ""))

        os.makedirs(out_dir, exist_ok=True)
        with open(out_dir.joinpath(f'{sample_name}_foswalk.csv'), 'w') as writeFile:
            writer = csv.writer(writeFile)
            writer.writerow(['Cluster', 'Size', 'Sequence'])
            writer.writerows(final_output)
    finally:
        # shutil.rmtree(Path("./").absolute())
        pass
