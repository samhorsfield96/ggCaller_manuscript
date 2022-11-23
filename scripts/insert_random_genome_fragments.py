from Bio import SeqIO
import argparse
import numpy as np

# seed the pseudorandom number generator
from random import seed
from random import randint

def get_options():
	description = "Adds random fragment to fasta file."
	parser = argparse.ArgumentParser(description=description,
                                     prog='python insert_random_genome_fragments.py')

	IO = parser.add_argument_group('Input/options.out')
	IO.add_argument('--frag_file',
       	             help='File to generate fragments from (FASTA format)')
	IO.add_argument('--out',
                    help='File to append to (FASTA format)')
	IO.add_argument('--frag_size',
                    type=int,
		    help='Size of fragment to add')
	IO.add_argument('--seed',
                    type=int,
		    default=None,
                    help='Seed for random number generation')
	return parser.parse_args()

def insert_random(infile, outfile, frag_size, seed_no):
	# determine if to pass on file
	s = np.random.poisson(1)
	if s == 0:
		return
	
	# read in fasta file to fragment
	seq_dict = {}
	fasta_index = 0
	for fasta in SeqIO.parse(open(infile),'fasta'):
		seq_dict[fasta_index] = str(fasta.seq)
		fasta_index += 1
	
	#set seed
	seed(seed_no)
	
	# for each fragment to be added, generate index to pull sequence from
	for i in range(0, s):
		finished = False
		while not finished:
			entry_index = randint(0, fasta_index - 1)
			entry_seq = seq_dict[entry_index]
			if len(entry_seq) > frag_size:
				finished = True
	
		# generate index for sequence to slice, ensuring it is > fragsize from end
		finished = False
		while not finished:
			seq_index = randint(0, len(entry_seq))
			if (len(entry_seq) - 1) - seq_index > frag_size:
				finished = True
	
		# slice sequence
		random_seq = entry_seq[seq_index:seq_index + frag_size]
	
		# insert into outfile
		with open(outfile, "a") as f:
			f.write(">Contaminant_" + str(i) + "\n" + random_seq + "\n")

def main():
	options = get_options()
	insert_random(options.frag_file, options.out, options.frag_size, options.seed)
	return 0

if __name__ == '__main__':
    main()
