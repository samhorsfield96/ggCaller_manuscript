import networkx as nx
import numpy as np
import scipy
from scipy import stats
import argparse


def get_options():
	parser = argparse.ArgumentParser(description='Print cluster size information from gml file', prog='python parse_gml.py')

	# input options
	parser.add_argument('--infile',
						required=True,
						help='.gml file to analysis')
	parser.add_argument('--ignore_singletons',
						default=False,
						action="store_true",
						help='Ignore singletons in cluster size calculations')
	parser.add_argument('--ignore_refound',
						default=False,
						action="store_true",
						help='Ignore refound genes.')
	parser.add_argument('--verbose',
						default=False,
						action="store_true",
						help='Print average and stddev values to console.')
	parser.add_argument('--outpref',
						default="parsed_graph",
						help='Output prefix. Default = parsed_graph')

	return parser.parse_args()

def get_node_stats(G, ignore_singletons, verbose, ignore_refound):
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

	for n in G.nodes():
		node_info = G.nodes[n]
		ORF_sizes = node_info['lengths']
		if ignore_refound:
			gene_ids = node_info['geneIDs'].split(";")
			to_remove = []
			for i in range(len(gene_ids)):
				gene_nums = gene_ids[i].split("_")
				if "refound" in gene_nums:
					to_remove.append(i)

			# reverse to ensure indexing still supported
			to_remove.reverse()
			for index in to_remove:
				try:
					del ORF_sizes[index]
				except IndexError:
					pass

		ORF_sizes = np.array(ORF_sizes).astype(int)
		cluster_size = ORF_sizes.size
		if cluster_size == 0 or (cluster_size == 1 and ignore_singletons):
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

	if verbose:
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

	combined_stdev = np.column_stack((stat_list_prop_stdev, stat_list_stdev, stat_list_stdev_size))
	combined_IQR = np.column_stack((stat_list_prop_IQR, stat_list_IQR, stat_list_IQR_size))
	combined_MAD = np.column_stack((stat_list_prop_MAD, stat_list_MAD, stat_list_MAD_size))

	return combined_stdev, combined_IQR, combined_MAD, stat_size_cluster

if __name__ == "__main__":
	options = get_options()
	infile = options.infile
	outpref = options.outpref
	verbose = options.verbose
	ignore_singletons = options.ignore_singletons
	ignore_refound = options.ignore_refound

	G = nx.read_gml(infile)
	if verbose:
		print("Input: " + infile)

	combined_stdev, combined_IQR, combined_MAD, stat_size_cluster = get_node_stats(G, ignore_singletons, verbose, ignore_refound)

	np.savetxt(outpref + "_stddev.txt", combined_stdev, delimiter=",")
	np.savetxt(outpref + "_IQR.txt", combined_IQR, delimiter=",")
	np.savetxt(outpref + "_MAD.txt", combined_MAD, delimiter=",")
	np.savetxt(outpref + "_sizes.txt", stat_size_cluster, delimiter=",")
