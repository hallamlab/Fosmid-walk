import os
os.chdir('./cache')

# Avery, 2023

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

import subprocess
import os
import re
import csv

# sample_names = ['CD1', 'CD2', 'CH1', 'CH2', 'SD1_1', 'SD1_2', 'SD2_1', 'SD2_2', 'SH1_1', 'SH1_2', 'SH2_1', 'SH2_2', 'SWT1_1', 'SWT1_2', 'SWT2_1', 'SWT2_2']
sample_names = ['10k_1', '10k_2']
# data_path = 'S5'
data_path = './cache'
final_output = [['Sample', 'Cluster', 'Size', 'Sequence', 'Reverse Complement']]
fasta_name = 'pcc1fos_kat'
fasta_name = 'pCC1fos'

# Indexing vector backbone
bwa_index = ['bwa', 'index',
             str('data/' + fasta_name + '.fasta')]
subprocess.call(bwa_index)

for s in sample_names:
    command_list = []

    # # Generate fastq from .bam files
    # fastq_gen = ['samtools', 'fastq', '-0',
    #              str('data/' + s + '.fastq'),
    #              str(data_path + '/' + s + '.bam')]
    # command_list.append(fastq_gen)

    mkdir = ['mkdir',s]
    command_list.append(mkdir)

    # Align reads to vector backbone
    bwa_mem = ['bwa', 'mem',
               str('data/' + fasta_name + '.fasta'),
               str('data/'+s+'.fastq'),
               '>',
               str(s+'/'+fasta_name+'_aln.sam')]
    command_list.append(bwa_mem)

    # Generate .bam file from .sam file
    sam_sort = ['samtools', 'sort',
                str(s+'/'+fasta_name+'_aln.sam'),
                '>',
                str(s+'/'+fasta_name+'_aln.bam')]
    command_list.append(sam_sort)

    # Index .bam file
    sam_index = ['samtools', 'index',
                 str(s+'/'+fasta_name+'_aln.bam')]
    command_list.append(sam_index)

    # Output .bam with reads aligned to pCC1fos:0-40
    sam_view1 = ['samtools', 'view',
                 str(s+'/'+fasta_name+'_aln.bam'),
                 str(fasta_name+':0-40'),
                 '-o',
                 str(s+'/'+fasta_name+'_5_040.bam')]
    command_list.append(sam_view1)

    # Output .bam with fwd reads aligned
    sam_view2 = ['samtools', 'view', '-F', '0x10',
                 str(s+'/'+fasta_name+'_5_040.bam'),
                 '-o',
                 str(s+'/040-fwd.bam')]
    command_list.append(sam_view2)

    # Output .bam with rev reads aligned
    sam_view3 = ['samtools', 'view', '-f', '0x10',
                 str(s+'/'+fasta_name+'_5_040.bam'),
                 '-o',
                 str(s+'/040-rev.bam')]
    command_list.append(sam_view3)

    # Generate .fasta file with fwd reads
    sam_fasta1 = ['samtools', 'fasta',
                  str(s+'/040-fwd.bam'),
                  '>',
                  str(s+'/040-fwd.fasta')]
    command_list.append(sam_fasta1)

    # Generate .fasta file with rev reads
    sam_fasta2 = ['samtools', 'fasta',
                  str(s+'/040-rev.bam'),
                  '>',
                  str(s+'/040-rev.fasta')]
    command_list.append(sam_fasta2)

    # Run commands in terminal
    for command in command_list:
        subprocess.call(' '.join(command), shell = True)

    # Open fwd and rev .fasta files
    fwd = SeqIO.parse(str(s+'/040-fwd.fasta'), "fasta")

    rev = SeqIO.parse(str(s+'/040-rev.fasta'), "fasta")

    # Append fwd reads and reverse_complement of rev reads to recs = []
    recs = []
    for f in fwd:
        recs.append(f)

    for r in rev:
        rev_r = r.seq.reverse_complement()
        r.seq = rev_r
        recs.append(r)

    # Trim sequences to 100bp from defined site
    final_recs = []
    for rec in recs:
        # Hard-coded sequence to trim to (GTGGGAT)
        if "GTGGGAT" in rec.seq:
            new_seq = rec.seq[:str(rec.seq).index('GTGGGAT')]
            if len(new_seq) >= 100:
                trim_seq = new_seq[len(new_seq)-100:len(new_seq)]
                rec.seq = trim_seq
                final_recs.append(rec)

    # Generate .fasta file
    SeqIO.write(final_recs, str(s+'/'+s+'-5-040.fasta'), "fasta")

    # Cluster trimmed sequences in .fasta
    # uclust_cmd = ["usearch",
    uclust_cmd = ["./usearch11.0.667_i86linux32",
                  "-cluster_fast", s + os.sep + s + '-5-040.fasta',
                  "-uc", s + os.sep + s + '-5-040.uc',
                  "-id", str(0.9),
                  "-consout", s + os.sep + s + '-5-040_clusters.fasta',
                  "-sizeout",
                  "-centroids", s + os.sep + s + '-5-040_centroids.fasta']
    subprocess.call(' '.join(uclust_cmd), shell = True)

    # Remove singletons
    clusters = SeqIO.parse(str(s+'/' + s + '-5-040_clusters.fasta'), "fasta")
    for c in clusters:
        if "size=1;" not in c.id:
            cluster = re.findall('\w+', c.id)[0]
            size = re.findall('(?<=size\=)\w+', c.id)[0]
            output_row = [s, cluster, size, c.seq, c.seq.reverse_complement()]
            final_output.append(output_row)

with open('output.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    for fo in final_output:
        writer.writerow(fo)
