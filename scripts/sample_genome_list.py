import random
import gzip
import shutil
import os
import argparse
from Bio import SeqIO

def get_options():
	description = 'Randomly samples genomes from list of file paths.'
	parser = argparse.ArgumentParser(description=description,
									 prog='python sample_genome_list.py')
	IO = parser.add_argument_group('Input/Output options')
	IO.add_argument('--infile',
					required=True,
					help='Infile containing paths to genomes to sample from, one per line')
	IO.add_argument('--outpref',
					default="result",
					help='Output prefix ')
	IO.add_argument('--genome_outdir',
					default=".",
					help='Directory to save unzipped genomes.')
	IO.add_argument('--sample_sizes',
					required=True,
					help='Sample sizes in form "a,b,c".')
	IO.add_argument('--include_Ns',
					default=False,
					action="store_true",
					help='Include assemblies with Ns. Default = False')
	return parser.parse_args()

def main():
	options = get_options()

	infile = options.infile
	outpref = options.outpref
	genome_outdir = options.outdir
	sample_sizes = [int(a) for a in options.samples.split(",")]
	avoid_Ns = not options.include_Ns

	total_sample = []

	with open(infile, "r") as f:
		for line in f:
			line = line.strip()
			total_sample.append(line)

	prev_size = 0
	ggc_samples = set()

	prev_sampled = set()
	for size in sample_sizes:
		sample_size = 0
		prok_samples = set()
		while sample_size < (size - prev_size):
			print("Sample size: " + str(sample_size))
			print("Total length: " + str(len(total_sample)))
			curr_sample = set(random.sample(total_sample, (size - prev_size) - sample_size))
			prev_sampled.update(curr_sample)

			to_remove = set()
			for file in curr_sample:
				outname = os.path.splitext(os.path.basename(file))[0]
				if ".gzip" in file:
					with gzip.open(file, 'rt') as f_in:
						if avoid_Ns:
							N_present = False
							for record in SeqIO.parse(f_in, "fasta"):
								if record.seq.count('N') > 0:
									N_present = True
									break
							if N_present:
								to_remove.add(file)
								continue
					with gzip.open(file, 'rb') as f_in:
						with open(genome_outdir + "/" + outname, 'wb') as f_out:
							shutil.copyfileobj(f_in, f_out)
				else:
					with open(file, 'r') as f_in:
						if avoid_Ns:
							N_present = False
							for record in SeqIO.parse(f_in, "fasta"):
								if record.seq.count('N') > 0:
									N_present = True
									break
							if N_present:
								to_remove.add(file)
								continue
					shutil.copy(file, genome_outdir + "/" + outname)

			# remove items with Ns
			curr_sample = list(set(curr_sample) - to_remove)
			ggc_samples.update(curr_sample)
			prok_samples.update(curr_sample)
			# remove previosly sampled items
			total_sample = set(total_sample)
			total_sample = total_sample.difference(prev_sampled)
			total_sample = list(total_sample)
			sample_size += len(curr_sample)

		prev_size = size

		# write prokka sample
		with open(outpref + "_prokka_N" + str(size) + ".txt", "w") as f:
			for entry in prok_samples:
				outname = os.path.splitext(os.path.basename(entry))[0]
				f.write(genome_outdir + "/" + outname + "\n")

		# write ggCaller file
		with open(outpref + "_ggc_N" + str(size) + ".txt", "w") as f:
			for entry in ggc_samples:
				outname = os.path.splitext(os.path.basename(entry))[0]
				f.write(genome_outdir + "/" + outname + "\n")


if __name__ == "__main__":
	main()



