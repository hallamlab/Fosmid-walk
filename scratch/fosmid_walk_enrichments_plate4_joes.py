
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

import subprocess
import os
import re
import csv

# Distinguishing portion of .fastq file name
# sample_names = ["52639.4.404742.AAGAGCCA-AAGAGCCA", "52639.4.404742.AGAAGCGT-AGAAGCGT", "52639.4.404742.ATGGTCCA-ATGGTCCA", "52639.4.404742.CCAAGCAA-CCAAGCAA",
# "52639.4.404742.CCATACGT-CCATACGT", "52639.4.404742.GTCTAGGT-GTCTAGGT", "52639.4.404742.TAACCGGT-TAACCGGT", "52639.4.404742.TATTCGCC-TATTCGCC", "52639.4.404742.TGTACCGT-TGTACCGT", "52639.4.404742.TTGAGGCA-TTGAGGCA", "52639.4.404742.TTGGACGT-TTGGACGT"]

os.chdir("./cache")
sample_names = ["10k_1"]

final_output = [['Sample', 'Cluster', 'Size', 'Sequence']]

for s in sample_names:
# Creates directory for each sample
    subprocess.call([str(f'rm -r {s} && mkdir -p '+s)], shell = True)
# Copies fosmid backbone senquence to new directory
    subprocess.call([str('cp pcc1fos_kat.fasta ' +s+'/pcc1fos.fasta')], shell = True)
    subprocess.call([str('bwa index '+s+'/pcc1fos.fasta')], shell = True)
# Aligns reads to pCC1 vector backbone ***PATH TO FILES UPDATED HERE***
    subprocess.call([str('bwa mem '+s+'/pcc1fos.fasta '+s+'.fastq > '+s+'/pcc1fos_aln.sam')], shell = True)
    subprocess.call([str('samtools sort '+s+'/pcc1fos_aln.sam > '+s+'/pcc1fos_aln.bam')], shell = True)
    subprocess.call([str('samtools index '+s+'/pcc1fos_aln.bam')], shell = True)
    subprocess.call([str('samtools view '+s+'/pcc1fos_aln.bam pCC1fos:0-20 -o '+s+'/pCC1-5-020.bam')], shell = True)
    subprocess.call([str('samtools view -F 0x10 '+s+'/pCC1-5-020.bam -o '+s+'/020-fwd.bam')], shell = True)
    subprocess.call([str('samtools view -f 0x10 '+s+'/pCC1-5-020.bam -o '+s+'/020-rev.bam')], shell = True)
    subprocess.call([str('samtools fasta '+s+'/020-fwd.bam > '+s+'/020-fwd.fasta')], shell = True)
    subprocess.call([str('samtools fasta '+s+'/020-rev.bam > '+s+'/020-rev.fasta')], shell = True)

    fwd = SeqIO.parse(str(s+'/020-fwd.fasta'), "fasta")

    rev = SeqIO.parse(str(s+'/020-rev.fasta'), "fasta")

    recs = []
    for f in fwd:
        recs.append(f)

    for r in rev:
        rev_r = r.seq.reverse_complement()
        r.seq = rev_r
        recs.append(r)

    final_recs = []
    for rec in recs:
        if "GTGGGAT" in rec.seq: # first n nt in vector backbone
            new_seq = rec.seq[:str(rec.seq).index('GTGGGAT')]
            if len(new_seq) >= 100:
                trim_seq = new_seq[len(new_seq)-100:len(new_seq)]
                rec.seq = trim_seq
                final_recs.append(rec)

    SeqIO.write(final_recs, str(s+'/pCC1-5-020.fasta'), "fasta")

    uclust_cmd = ["./usearch11.0.667_i86linux32",
                  "-cluster_fast", s + os.sep + "pCC1-5-020.fasta",
                  "-uc", s + os.sep + "pCC1-5-020.uc",
                  "-id", str(0.9),
                  "-consout", s + os.sep + "pCC1-5-020_clusters.fasta",
                  "-sizeout",
                  "-centroids", s + os.sep + "pCC1-5-020_centroids.fasta",
                  "-minsize", str(2)]
    subprocess.call(' '.join(uclust_cmd), shell = True)

    clusters = SeqIO.parse(str(s+'/pCC1-5-020_clusters.fasta'), "fasta")
    for c in clusters:
        if "size=1;" not in c.id:
            cluster = re.findall('\w+', c.id)[0]
            size = re.findall('(?<=size\=)\w+', c.id)[0]
            output_row = [s, cluster, size, c.seq]
            final_output.append(output_row)

with open('output.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    for fo in final_output:
        writer.writerow(fo)
