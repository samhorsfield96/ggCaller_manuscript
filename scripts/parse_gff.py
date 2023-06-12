import argparse
import numpy as np
import scipy
from scipy import stats
from Bio import SeqIO
import json

def get_options():
	parser = argparse.ArgumentParser(description='Print cluster size information from PEPPAN gff file and roary clusters file.', prog='python parse_gff.py')

	# input options
	parser.add_argument('--peppan',
						required=True,
						help='PEPPAN .gff file to analysis, or pair of PEPPAN .genes and .encode.csv in format "x.genes,x.encode.csv"')
	parser.add_argument('--roary',
						default=None,
						help='Roary clusters matrix file to analysis')
	parser.add_argument('--len_dict',
						default=None,
						help='Length dictionary .json. Required if --roary specified')
	parser.add_argument('--ignore_singletons',
						default=False,
						action="store_true",
						help='Ignore singletons in cluster size calculations')
	parser.add_argument('--ignore_pseudogenes',
						default=False,
						action="store_true",
						help='Ignore psuedogenes.')
	parser.add_argument('--verbose',
						default=False,
						action="store_true",
						help='Print average and stddev values to console.')
	parser.add_argument('--outpref',
						default="parsed_graph",
						help='Output prefix. Default = parsed_graph')

	return parser.parse_args()

def compare_peppan(infile, ignore_pseudogenes):
	len_dict = {}

	# iterate over paired files
	cluster_dict = {}

	singleton_id = 0

	with open(infile, "r") as gff:
		for line in gff:
			if line[0] != "#":
				split_line = line.rstrip().split("\t")

				# get isolate
				#iso = split_line[0].split(":")[0]

				type = split_line[1]
				if type == "misc_feature":
					continue

				#iso_set.add(iso)

				length = int(split_line[4]) - int(split_line[3])

				# attempt to parse name
				if "old_locus_tag" in split_line[-1]:
					old_loci = split_line[-1].split("old_locus_tag=")[-1].split(";")[0].split(",")
					for old_locus in old_loci:
						name = old_locus.split(":")[0]
						old_length_range = old_locus.split(":")[1].split("-")
						old_length = int(old_length_range[1]) - int(old_length_range[0])
						len_dict[name] = old_length

				# ignore psuedogene
				if ignore_pseudogenes and type == "pseudogene":
					continue

				if "ortholog_group" not in split_line[-1]:
					# add singleton
					cluster_id = "single_" + str(singleton_id)
					singleton_id += 1
				else:
					info = split_line[-1].split("ortholog_group:")[-1].split(";")[0].split("/")[0].split(":")
					cluster_id = "_".join((info[0], info[1]))

				if cluster_id not in cluster_dict:
					cluster_dict[cluster_id] = []
				cluster_dict[cluster_id].append(length)

	return cluster_dict, len_dict

def compare_roary(infile, len_dict):
	# iterate over paired files
	cluster_dict = {}

	with open(infile, "r") as mat:
		for line in mat:
			split_line = line.rstrip().split("\t")
			cluster_id = split_line[0]

			for gene in split_line[1::]:
				# get isolate
				#iso = split_line[0].split(":")[0]

				#iso_set.add(iso)

				# if gene not in len_dict:
				# 	continue
				length = len_dict[gene]

				if cluster_id not in cluster_dict:
					cluster_dict[cluster_id] = []
				cluster_dict[cluster_id].append(length)

	return cluster_dict

def parse_peppan(peppan_genes, peppan_encode):
	# get encodings for peppan
	id_dict = {}
	with open(peppan_encode) as f:
		for line in f:
			#print(line)
			if ":" not in line:
				continue
			line.rstrip()
			split_line = line.split(",")
			id_dict[int(split_line[1])] = split_line[0].split(":")[1]

	fasta_sequences = SeqIO.parse(open(peppan_genes),'fasta')
	
	len_dict = {}
	for fasta in fasta_sequences:
		name, sequence = int(fasta.id), fasta.seq
		len_dict[id_dict[name]] = len(sequence)

	return {}, len_dict

if __name__ == "__main__":
	options = get_options()
	peppan_infile = options.peppan

	roary_infile = options.roary
	len_dict_file = options.len_dict
	outpref = options.outpref
	ignore_singletons = options.ignore_singletons
	ignore_pseudogenes = options.ignore_pseudogenes
	verbose = options.verbose

	if "," in peppan_infile:
		split_peppan = peppan_infile.split(",")
		peppan_genes = split_peppan[0]
		peppan_encode = split_peppan[1]
		peppan_dict_list, len_dict = parse_peppan(peppan_genes, peppan_encode)
	else:
		peppan_dict_list, len_dict = compare_peppan(peppan_infile, ignore_pseudogenes)
	
	# iterate over peppan inputs and calculate statistics
	stat_list_IQR = []
	stat_list_prop_IQR = []
	stat_list_MAD = []
	stat_list_prop_MAD = []
	stat_list_prop_stdev = []
	stat_list_stdev = []
	stat_size_cluster = []

	# get stats per cluster size
	stat_list_IQR_size = []
	stat_list_MAD_size = []
	stat_list_stdev_size = []

	for cluster_id, ORF_sizes in peppan_dict_list.items():
		ORF_sizes = np.array(ORF_sizes).astype(int)
		cluster_size = ORF_sizes.size
		if cluster_size == 1 and ignore_singletons:
			continue

		stat_size_cluster.append(cluster_size)
		q75, q50, q25 = np.percentile(ORF_sizes, [75, 50, 25])

		if cluster_size == 1:
			stat_list_MAD.append(0.0)
			stat_list_prop_MAD.append(0.0)
			stat_list_MAD_size.append(0.0)
		else:
			stat_list_MAD.append(scipy.stats.median_abs_deviation(ORF_sizes))
			stat_list_MAD_size.append(scipy.stats.median_abs_deviation(ORF_sizes) / cluster_size)
			stat_list_prop_MAD.append(scipy.stats.median_abs_deviation(ORF_sizes) / q50)

		stat_list_IQR.append(float(q75 - q25))
		stat_list_IQR_size.append(float((q75 - q25) / cluster_size))
		stat_list_prop_IQR.append(float((q75 - q25) / q50))

		stat_list_prop_stdev.append(float(np.std(ORF_sizes)) / float(np.mean(ORF_sizes)))
		stat_list_stdev.append(float(np.std(ORF_sizes)))
		stat_list_stdev_size.append(float(np.std(ORF_sizes) / cluster_size))

	# convert to numpy arrays
	stat_list_IQR = np.array(stat_list_IQR)
	stat_list_IQR_size = np.array(stat_list_IQR_size)
	stat_list_prop_IQR = np.array(stat_list_prop_IQR)

	stat_list_MAD = np.array(stat_list_MAD)
	stat_list_MAD_size = np.array(stat_list_MAD_size)
	stat_list_prop_MAD = np.array(stat_list_prop_MAD)

	stat_list_stdev = np.array(stat_list_stdev)
	stat_list_stdev_size = np.array(stat_list_stdev_size)
	stat_list_prop_stdev = np.array(stat_list_prop_stdev)

	stat_size_cluster = np.array(stat_size_cluster)

	combined_stdev = np.column_stack((stat_list_prop_stdev, stat_list_stdev, stat_list_stdev_size))
	combined_IQR = np.column_stack((stat_list_prop_IQR, stat_list_IQR, stat_list_IQR_size))
	combined_MAD = np.column_stack((stat_list_prop_MAD, stat_list_MAD, stat_list_MAD_size))

	np.savetxt(outpref + "_peppan_stddev.txt", combined_stdev, delimiter=",")
	np.savetxt(outpref + "_peppan_IQR.txt", combined_IQR, delimiter=",")
	np.savetxt(outpref + "_peppan_MAD.txt", combined_MAD, delimiter=",")
	np.savetxt(outpref + "_peppan_sizes.txt", stat_size_cluster, delimiter=",")

	if verbose:
		print("Input: " + peppan_infile)
		print("Average prop. standard deviation: {}".format(np.mean(stat_list_prop_stdev)))
		print("Stdev prop. standard deviation: {}".format(np.std(stat_list_prop_stdev)))
		print("Average standard deviation: {}".format(np.mean(stat_list_stdev)))
		print("Stdev standard deviation: {}".format(np.std(stat_list_stdev)))
		print("Average standard deviation (size adjusted): {}".format(np.mean(stat_list_stdev_size)))
		print("Stdev standard deviation (size adjusted: {}".format(np.std(stat_list_stdev_size)))

		print("Average prop. IQR: {}".format(np.mean(stat_list_prop_IQR)))
		print("Stdev prop. IQR: {}".format(np.std(stat_list_prop_IQR)))
		print("Average IQR: {}".format(np.mean(stat_list_IQR)))
		print("Stdev IQR: {}".format(np.std(stat_list_IQR)))
		print("Average IQR (size adjusted): {}".format(np.mean(stat_list_IQR_size)))
		print("Stdev IQR (size adjusted): {}".format(np.std(stat_list_IQR_size)))

		print("Average prop. MAD: {}".format(np.mean(stat_list_prop_MAD)))
		print("Stdev prop. MAD: {}".format(np.std(stat_list_prop_MAD)))
		print("Average MAD: {}".format(np.mean(stat_list_MAD)))
		print("Stdev MAD: {}".format(np.std(stat_list_MAD)))
		print("Average MAD (size adjusted): {}".format(np.mean(stat_list_MAD_size)))
		print("Stdev MAD (size adjusted): {}".format(np.std(stat_list_MAD_size)))

		print("Average cluster size: {}".format(np.mean(stat_size_cluster)))
		print("Stdev cluster size: {}".format(np.std(stat_size_cluster)))

	if roary_infile != None:
		with open(len_dict_file, "r") as f:
			len_dict = json.load(f)

		roary_dict_list = compare_roary(roary_infile, len_dict)

		# iterate over peppan inputs and calculate statistics
		stat_list_IQR = []
		stat_list_prop_IQR = []
		stat_list_MAD = []
		stat_list_prop_MAD = []
		stat_list_prop_stdev = []
		stat_list_stdev = []
		stat_size_cluster = []

		# get stats per cluster size
		stat_list_IQR_size = []
		stat_list_MAD_size = []
		stat_list_stdev_size = []

		for cluster_id, ORF_sizes in roary_dict_list.items():
			ORF_sizes = np.array(ORF_sizes).astype(int)
			cluster_size = ORF_sizes.size
			if cluster_size == 1 and ignore_singletons:
				continue
			stat_size_cluster.append(cluster_size)

			q75, q50, q25 = np.percentile(ORF_sizes, [75, 50, 25])

			if cluster_size == 1:
				stat_list_MAD.append(0.0)
				stat_list_MAD_size.append(0.0)
				stat_list_prop_MAD.append(0.0)
			else:
				stat_list_MAD.append(scipy.stats.median_abs_deviation(ORF_sizes))
				stat_list_MAD_size.append(scipy.stats.median_abs_deviation(ORF_sizes) / cluster_size)
				stat_list_prop_MAD.append(scipy.stats.median_abs_deviation(ORF_sizes) / q50)

			stat_list_IQR.append(float(q75 - q25))
			stat_list_IQR_size.append(float((q75 - q25) / cluster_size))
			stat_list_prop_IQR.append(float((q75 - q25) / q50))

			stat_list_prop_stdev.append(float(np.std(ORF_sizes)) / float(np.mean(ORF_sizes)))
			stat_list_stdev.append(float(np.std(ORF_sizes)))
			stat_list_stdev_size.append(float(np.std(ORF_sizes) / cluster_size))

		# convert to numpy arrays
		stat_list_IQR = np.array(stat_list_IQR)
		stat_list_IQR_size = np.array(stat_list_IQR_size)
		stat_list_prop_IQR = np.array(stat_list_prop_IQR)

		stat_list_MAD = np.array(stat_list_MAD)
		stat_list_MAD_size = np.array(stat_list_MAD_size)
		stat_list_prop_MAD = np.array(stat_list_prop_MAD)

		stat_list_stdev = np.array(stat_list_stdev)
		stat_list_stdev_size = np.array(stat_list_stdev_size)
		stat_list_prop_stdev = np.array(stat_list_prop_stdev)

		stat_size_cluster = np.array(stat_size_cluster)

		combined_stdev = np.column_stack((stat_list_prop_stdev, stat_list_stdev, stat_list_stdev_size))
		combined_IQR = np.column_stack((stat_list_prop_IQR, stat_list_IQR, stat_list_IQR_size))
		combined_MAD = np.column_stack((stat_list_prop_MAD, stat_list_MAD, stat_list_MAD_size))

		np.savetxt(outpref + "_roary_stddev.txt", combined_stdev, delimiter=",")
		np.savetxt(outpref + "_roary_IQR.txt", combined_IQR, delimiter=",")
		np.savetxt(outpref + "_roary_MAD.txt", combined_MAD, delimiter=",")
		np.savetxt(outpref + "_roary_sizes.txt", stat_size_cluster, delimiter=",")

		if verbose:
			print("Input: " + roary_infile)
			print("Average prop. standard deviation: {}".format(np.mean(stat_list_prop_stdev)))
			print("Stdev prop. standard deviation: {}".format(np.std(stat_list_prop_stdev)))
			print("Average standard deviation: {}".format(np.mean(stat_list_stdev)))
			print("Stdev standard deviation: {}".format(np.std(stat_list_stdev)))
			print("Average standard deviation (size adjusted): {}".format(np.mean(stat_list_stdev_size)))
			print("Stdev standard deviation (size adjusted: {}".format(np.std(stat_list_stdev_size)))

			print("Average prop. IQR: {}".format(np.mean(stat_list_prop_IQR)))
			print("Stdev prop. IQR: {}".format(np.std(stat_list_prop_IQR)))
			print("Average IQR: {}".format(np.mean(stat_list_IQR)))
			print("Stdev IQR: {}".format(np.std(stat_list_IQR)))
			print("Average IQR (size adjusted): {}".format(np.mean(stat_list_IQR_size)))
			print("Stdev IQR (size adjusted): {}".format(np.std(stat_list_IQR_size)))

			print("Average prop. MAD: {}".format(np.mean(stat_list_prop_MAD)))
			print("Stdev prop. MAD: {}".format(np.std(stat_list_prop_MAD)))
			print("Average MAD: {}".format(np.mean(stat_list_MAD)))
			print("Stdev MAD: {}".format(np.std(stat_list_MAD)))
			print("Average MAD (size adjusted): {}".format(np.mean(stat_list_MAD_size)))
			print("Stdev MAD (size adjusted): {}".format(np.std(stat_list_MAD_size)))

			print("Average cluster size: {}".format(np.mean(stat_size_cluster)))
			print("Stdev cluster size: {}".format(np.std(stat_size_cluster)))






