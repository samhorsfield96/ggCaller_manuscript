import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def get_options():
	description = 'Fragments sequence at each gene location. Can also be used to translate CDS sequences.'
	parser = argparse.ArgumentParser(description=description,
									 prog='python fragment_at_gene.py')
	IO = parser.add_argument_group('Input/Output options')
	IO.add_argument('--fastafile',
					default=None,
					help='Input protein fasta file. ')
	IO.add_argument('--panaroodata',
					default=None,
					help='Input panaroo data file. ')
	IO.add_argument('--targets',
					default=None,
					help='Targets to search for. Can be comma separated. ')
	IO.add_argument('--prokkamap',
					default=None,
					help='Mapping from prokka names to COG annotation.'
						 'Generate using ```grep "<target> *.gff > prokkamap.txt"``` ')
	IO.add_argument('--outpref',
					default="results",
					help='Output prefix ')
	IO.add_argument('--ignorepseudo',
					default=False,
					action="store_true",
					help='Ignore pseudogenes.'
						 'Default = False')
	return parser.parse_args()

def parse_fasta(infile, tag, outpref, ignore_pseudogenes):
	fasta_sequences = SeqIO.parse(open(infile),'fasta')
	output_gen = []
	for fasta in fasta_sequences:
		name, description, sequence = fasta.id, fasta.description, fasta.seq
		if ignore_pseudogenes and "psuedogene" in description:
			continue
		if any(gene in description for gene in tag):
			output_gen.append(SeqRecord(sequence, description=tag, id=name))
	SeqIO.write(output_gen, outpref + "_" + tag + ".faa", "fasta")

def parse_panaroo(gene_set, data, outpref, ignore_pseudogenes):

	output_gen = []

	with open(data, "r") as f:
		lines = f.readlines()[1:]
		for line in lines:
			splitted = line.split(",")
			name = splitted[3]
			if ignore_pseudogenes and ("_len" in name or "_stop" in name):
				continue
			if name in gene_set:
				sequence = splitted[4]
				#print(sequence)
				output_gen.append(SeqRecord(Seq(sequence), description="", id=name))

	SeqIO.write(output_gen, outpref + "panaroo_"+ tag + ".faa", "fasta")

if __name__ == "__main__":
	options = get_options()
	fasta_infile = options.fastafile
	data = options.panaroodata
	targets = options.targets.split(",")
	prokkamap = options.prokkamap

	outpref = options.outpref

	ignore_pseudogenes = options.ignorepseudo

	if fasta_infile is not None:
		if targets is None:
			print("Please specific list of targets.")
			sys.exit(1)
		parse_fasta(fasta_infile, targets, outpref, ignore_pseudogenes)
	elif data is not None:
		gene_set = set()
		with open(prokkamap, "r") as f:
			for line in f:
				split = line.split("locus_tag=")[-1].split(";")[0]
				gene_set.add(split)

		parse_panaroo(gene_set, data, outpref, ignore_pseudogenes)
	else:
		print("Please either specify a panaroo data file or fasta file")
		sys.exit(1)


