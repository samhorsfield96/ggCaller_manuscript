from Bio import SeqIO
import argparse
import numpy as np
import os
import math
import random

def get_options():
	description = "Randomly fragments fasta file."
	parser = argparse.ArgumentParser(description=description,
                                     prog='python fragment_fasta.py')

	IO = parser.add_argument_group('Input/options.out')
	IO.add_argument('--ref_dir',
                     help='Directory of fasta files to get empirical lengths from')
	IO.add_argument('--sample_dir',
                     help='Directory of fasta files to fragment in place')
	IO.add_argument('--min_length',
                     type=int,
                     help='Minimum length of sequence allowed')
	IO.add_argument('--frag_file',
       	             help='File to generate fragments from (FASTA format)')
	IO.add_argument('--frag_size',
                    type=int,
		    help='Mean size of fragment')
	IO.add_argument('--std_dev', 
			type=int,
			help='Std deviation of fragment size')
	return parser.parse_args()

def fragment_random(infile, frag_size, std_dev):
	# read in fasta file to fragment
	seq_dict = {}
	fasta_index = 0
	for fasta in SeqIO.parse(open(infile),'fasta'):
		seq_dict[fasta_index] = str(fasta.seq)
		fasta_index += 1
	
	# create list for fragmenting
	fragment_list = []
		
	# iterate over all the contigs, fragmenting at position determined by normal distribution
	for index, contig in seq_dict.items():
		frag_start = 0
		frag_end = 0
		contig_end = len(contig)
		finished = False
		while not finished:
			# draw from normal distribution
			s = round(np.random.normal(loc=frag_size, scale=std_dev, size=1)[0])
			frag_end = frag_start + s
			if frag_end < contig_end:
				fragment_list.append((index, frag_start, frag_end))
				frag_start = frag_end
			else:
				fragment_list.append((index, frag_start, contig_end))
				finished = True
	
	# fragment original fasta file
	with open(infile, "w") as f:
		contig_id = 0
		for fragment in fragment_list:
			index, frag_start, frag_end = fragment
			f.write(">" + str(index) + "_" + str(frag_start) + "_" + str(frag_end) + "\n" + seq_dict[index][frag_start : frag_end] + "\n")

def count_fragments(dir):
	dir_list = []
	for file in os.listdir(dir):
		if file.endswith(".fasta"):
			total_length = 0
			fasta_list = []
			for fasta in SeqIO.parse(open(dir + "/" + file),'fasta'):
				fasta_list.append(len(fasta.seq))
				total_length += len(fasta.seq)
			fasta_list.sort()
			prop_list = [x / total_length for x in fasta_list]
			dir_list.append(prop_list)
	return dir_list

def fragment_empirical(ref_dir, sample_dir, min_length):
	dir_list = count_fragments(ref_dir)
	for file in os.listdir(sample_dir):
		if file.endswith(".fasta"):
			seq_list = []
			prop_list = random.sample(dir_list, 1)[0]
			prop_list = random.sample(prop_list, len(prop_list))
			for fasta in SeqIO.parse(open(sample_dir + "/" + file),'fasta'):
				total_len = len(fasta.seq)
				index1 = 0
				prev_prop = 0
				total_slice = 0
				for i, prop in enumerate(prop_list):
					if (i < len(prop_list) - 2) and (total_len - total_slice > min_length):
						index2 = index1 + math.ceil((total_len * prev_prop)) + math.ceil((total_len * prop))
						slice = str(fasta.seq)[index1:index2]
						if len(slice) < min_length:
							prev_prop += prop
							continue
						prev_prop = 0
						total_slice += len(slice)
						seq_list.append(slice)
						index1 = index2
					else:
						slice = str(fasta.seq)[index1:]
						seq_list.append(slice)
						break
			with open(sample_dir + "/" + file, "w") as f:
				contig_id = 0
				for fragment in seq_list:
					f.write(">" + str(contig_id) + "\n" + fragment + "\n")
					contig_id += 1

def main():
	options = get_options()
	fragment_empirical(options.ref_dir, options.sample_dir, options.min_length)
	return 0

if __name__ == '__main__':
    main()
