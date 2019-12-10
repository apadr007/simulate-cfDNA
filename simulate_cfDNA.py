import os
os.chdir("/Users/ganglion/Google drive/Coursera/better coding/cfDNA simulation")
print("Current Working Directory " , os.getcwd())

def sample_frag_size(sequence, bps, prob):
        """ (str, list of int, list of int) -> tuple

        sequence: nucleotide sequence
        bps: list of size values
        prob: weights for the probability of picking a size range value

        Return index for the sequence being presented by:
                1) checking that the size of the sequence is at least 180 bp long.
                2) Randomly selecting either 180, 280, 320bp fragement using biased
                   weights.


        """
        import numpy as np
        import random

        frag_size = np.random.choice(bps, p=prob)

        if len(sequence) >= frag_size:
                start = random.randint(0, frag_size)
                stop = start + frag_size
                if stop > len(sequence):
                        start = 0
                        stop = start + frag_size


                if stop <= len(sequence):
                        #print('Frag size:', frag_size, "seq length:", len(sequence))
                        return( tuple((int(start), int(stop))) )
                else:
                        pass
        else:
                pass

def sequence_parser(input_name, bps, prob, output_name):
        """ (str of filename, list of int, list of int, str) -> fasta file

        sequence_parser takes a fasta file name with the fragment size (bp),
        and weights of selecting those sizes (prob), and writes
        them into a file (outputname) of size (size).

        The final size of the output file (output_name) will be the length of the
        input fasta file (input_name)

        check to make the output file
        >>> sequence_parser('head_GRCh38_cDNA.fa', [180, 280, 320], [0.3, 0.5, 0.2], 'test_human_frag_sizes.fa')
        >>> os.path.isfile('test_human_frag_sizes.fa') == True
        True

        check that output file is fasta
        >>> print(is_fasta('test_human_frag_sizes.fa'))
        True

        """
        from Bio import SeqIO

        with open(input_name, 'r') as fasta, open(output_name, 'w') as output_file:
                for fasta_sequences in SeqIO.parse(fasta, "fasta"):
                        name = fasta_sequences.id
                        fasta_seq = fasta_sequences.seq

                        myindex = sample_frag_size(fasta_seq, bps, prob)
                        if myindex != None:
                                output_file.write(">" + name + '\n' + str(fasta_seq[myindex[0]:myindex[1]]) + '\n')
                        else:
                                pass


def is_fasta(filename):
        """ (str (filename) ) -> bool

        check if a file is in fasta format. This function was made to use for
        the sequence_parser function

        >>> myfile = is_fasta('test.fa')
        >>> myfile == True
        True
        """

        from Bio import SeqIO
        with open(filename, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file




input_name = 'GRCh38_cDNA.fa'
bps = [180, 280, 320]
prob = [0.3, 0.5, 0.2]
output_name = 'human_frag_sizes.fa'



sequence_parser(input_name, bps, prob, output_name)



# run doct test if this is main
if __name__ == '__main__':
    import doctest
    doctest.testmod()
