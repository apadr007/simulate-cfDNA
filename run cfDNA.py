import simulate_cfDNA
import os
from Bio import SeqIO


input_name = 'GRCh38_cDNA.fa'
bps = [180, 280, 320]
prob = [0.3, 0.5, 0.2]
output_name = 'human_frag_sizes.fa'

# simulate cfDNA
simulate_cfDNA.sequence_parser(input_name, bps, prob, output_name)


# calculate size distribution of each sequence in the file
myfasta = 'human_frag_sizes.fa'
with open(myfasta, 'r') as fasta:
   fasta_length = []
   for fasta_sequences in SeqIO.parse(fasta, "fasta"):
      fasta_length.append(len(fasta_sequences.seq))
print(fasta_length)
