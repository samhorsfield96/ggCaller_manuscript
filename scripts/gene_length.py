from Bio import SeqIO
import argparse

def get_options():
	description = 'Generates file with gene lengths.'
	parser = argparse.ArgumentParser(description=description,
					prog='python gene_length.py')

	IO = parser.add_argument_group('Input/Output options')
	IO.add_argument('--infile',
			default=None,
			help='Input file ')
	IO.add_argument('--outfile',
			default=None,
			help='Output file ')
	return parser.parse_args()

if __name__ == "__main__":
	
	options = get_options()
	infile = options.infile
	outfile = options.outfile

	len_dist = []

	fasta_sequences = SeqIO.parse(open(infile),'fasta')
	for fasta in fasta_sequences:
		len_dist.append(len(fasta.seq))

	with open(outfile, "w") as f:
		for entry in len_dist:
			f.write(str(entry) + "\n")

