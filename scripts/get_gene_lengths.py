import json
import glob
import os
import argparse
from multiprocessing import Pool

def get_options():
    description = 'Generates json with all gene lengths from directory.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python sample_genome_list.py')
    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--indir',
                    required=True,
                    help='Directory containing gffs.')
    IO.add_argument('--outpref',
                    default="lengths",
                    help='Output prefix for file ')
    IO.add_argument('--threads',
                    default=1,
                    type=int,
                    help='Number of threads. Default = 1 ')
    return parser.parse_args()

def map_analysis(infile):
    gff_dict = {}

    with open(infile, "r") as gff:
        for line in gff:
            if line[0] != "#":
                split_line = line.rstrip().split("\t")
                if split_line[2] == "gene":
                    length = int(split_line[4]) - int(split_line[3])
                    locus_tag = split_line[-1].split("locus_tag=")[-1]
                    gff_dict[locus_tag] = length

    return gff_dict

def main():
    options = get_options()
    indir = options.indir
    outpref = options.outpref
    threads = options.threads

    all_files = glob.glob(os.path.join(indir, "*.gff"))

    gff_dict = {}

    with Pool(processes=threads) as pool:
        for gff_dict_temp in pool.map(map_analysis, all_files):
            gff_dict.update(gff_dict_temp)

    with open(outpref + ".json", "w") as f:
        json.dump(gff_dict, f)

if __name__ == "__main__":
    main()