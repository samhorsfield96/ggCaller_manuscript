from Bio import SeqIO
import argparse
import numpy as np
from pylab import *
import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO

def get_options():
    description = 'Compares start and end positions of genes based on alignment fasta.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python gene_end_comparison.py')
    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--indir',
                    required=True,
                    help='Input directory containing ggCaller and panaroo alignment files. '
                         'Files to be analysised must be in format <tool>_<gene>.aln')
    IO.add_argument('--id-file',
                    default=None,
                    help='Read in previously generated pairwise average amino acid identity matrix')
    IO.add_argument('--outpref',
                    default="result",
                    help='Output prefix ')
    return parser.parse_args()

def count_gaps(infile, tool_dict):
    # list of tuples, [0] = start, [1] = end, [2] = within
    gap_count = []

    file_pref = os.path.splitext(os.path.basename(infile))[0]
    split_file_pref = file_pref.split("_")
    tool = split_file_pref[0]
    gene = split_file_pref[1]

    if tool_dict != None:
        tool = tool_dict[tool]

    align_len = 0

    for record in SeqIO.parse(infile, "fasta"):
        sequence = str(record.seq)
        total_gaps = sequence.count('-')
        seq_len = len(sequence) - total_gaps

        if not align_len:
            align_len = len(sequence)

        # iterate over string forwards to count forward gaps
        start_gaps = 0
        for letter in sequence:
            if letter != "-":
                break
            start_gaps += 1

        end_gaps = 0
        for letter in reversed(sequence):
            if letter != "-":
                break
            end_gaps += 1

        gap_count.append((start_gaps, end_gaps, total_gaps - (start_gaps + end_gaps), seq_len))

    reference = gap_count[0]

    diff_start = np.zeros(len(gap_count) - 1)
    diff_end = np.zeros(len(gap_count) - 1)
    diff_within = np.zeros(len(gap_count) - 1)

    # taken as number of gaps in reference - number of gaps in entry, over the reference sequence length
    # negative number means reference starts before entry
    # positive number means entry starts before reference
    for index in range(1, len(gap_count)):
        entry = gap_count[index]
        index = index - 1
        diff_start[index] = (reference[0] - entry[0]) / align_len
        diff_end[index] = (reference[1] - entry[1]) / align_len
        diff_within[index] = (entry[2] - reference[2]) / align_len

    # create pandas dataframe
    li = []
    for name, array in zip(("start", "end", "within"), (diff_start, diff_end, diff_within)):
        data = np.column_stack(([tool] * array.size, [gene] * array.size, [name] * array.size, array))
        df = pd.DataFrame(data, columns=["Tool", "Gene", "Type", "Diff"])

        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=True)

    return frame, tool, gene

def read_files(in_dir, aai=False, tool_dict=None, prefix="", ext="txt"):
    all_files = glob.glob(os.path.join(in_dir, prefix + "*." + ext))

    li = []
    aai_li = []
    for filename in all_files:
        df, tool, gene = count_gaps(filename, tool_dict)

        # commented out as takes long time
        if aai:
            aai_df = get_aai(filename, tool, gene)
        else:
            aai_df = pd.DataFrame()
        aai_li.append(aai_df)

        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=True)
    aai_frame = pd.concat(aai_li, axis=0, ignore_index=True)

    frame['Diff'] = frame['Diff'].astype(float)
    return frame, aai_frame

def get_aai(filename, tool, gene):
    align = AlignIO.read(filename, "fasta")
    align_len = len(align)
    id_list = []

    col_start = 0

    for i1, align1 in enumerate(align):
        align1_arr = np.array(list(str(align1.seq)))
        align1_arr_gaps = align1_arr == "-"
        for i2 in range(col_start, align_len):
            if i1 == i2:
                continue
            align2_arr = np.array(list(str(align[i2].seq)))
            align2_arr_gaps = align2_arr == "-"
            num_match = np.count_nonzero(align1_arr == align2_arr)
            num_match_gaps = np.count_nonzero(logical_and(align1_arr_gaps, align2_arr_gaps))
            id_list.append((num_match - num_match_gaps) / (align1_arr.size - num_match_gaps))
        col_start += 1

    df = pd.DataFrame(id_list, columns=['perc_id'])

    df['Tool'] = tool
    df['Gene'] = gene

    df['perc_id'] = df['perc_id'].astype(float)

    return df

def main():
    options = get_options()
    indir = options.indir
    id_file = options.id_file
    outpref = options.outpref

    # determine whether to generate aai matrix
    aai = False
    if id_file is None:
        aai = True

    tool_dict = {"GGC": "ggCaller", "PAN": "Prokka + Panaroo"}
    data_full, aai_full = read_files(indir, aai=aai, tool_dict=tool_dict, ext="aln")

    # save aai file
    if aai:
        aai_full.to_csv(outpref + "_aai_mat.csv", index=False)

    sns.set(font_scale=1.5)
    sns.set(style='whitegrid')

    # set colour palette
    colors = ["#FF0B04", "#5A5A5A"]
    sns.set_palette(sns.color_palette(colors))

    #plot average amino acid identity
    if not aai:
        aai_full = pd.read_csv(id_file)

    plot = sns.FacetGrid(aai_full, col="Gene", hue="Tool", sharey=False, legend_out=True)

    plot.map(sns.histplot, "perc_id", stat='probability', binwidth=0.025, binrange=(0, 1.0)).add_legend()
    plot._legend.set_title("Workflow")

    plot.set(xlabel='Average amino acid identity', ylabel='Proportion')

    plt.savefig(outpref + '_aai_hist.png')

    plt.clf()

    # plot starts
    data = data_full[(data_full['Type'] == "start")]

    if tool_dict != None:
        args = dict(x="Tool", y="Diff", data=data, hue="Tool", hue_order=[tool_dict["GGC"], tool_dict["PAN"]],
                    order=[tool_dict["GGC"], tool_dict["PAN"]])
    else:
        args = dict(x="Tool", y="Diff", data=data, hue="Tool", hue_order=['GGC', 'PAN'], order=['GGC', 'PAN'])

    plot = sns.catplot(
        data=data, x='Tool', y='Diff',
        col='Gene', hue="Tool", kind='box',
        sym="", dodge=False, hue_order=args["hue_order"], order=args["order"]
    )

    plot.map(sns.stripplot, args["x"], args["y"], args["hue"], hue_order=args["hue_order"], order=args["order"],
          palette=colors, dodge=False, alpha=0.4, ec='k', linewidth=1)

    for ax_n in plot.axes:
        for ax in ax_n:
            ax.axhline(0, linewidth=2, color='gray', linestyle="--", alpha=0.6)

    plot.set(xlabel='Workflow', ylabel='Start difference (prop. alignment length)')

    plt.savefig(outpref + '_boxcompare_start.png')

    plt.clf()

    # plot end
    data = data_full[(data_full['Type'] == "end")]

    plot = sns.catplot(
        data=data, x='Tool', y='Diff',
        col='Gene', hue="Tool", kind='box',
        sym="", dodge=False, hue_order=args["hue_order"], order=args["order"]
    )

    plot.map(sns.stripplot, args["x"], args["y"], args["hue"], hue_order=args["hue_order"], order=args["order"],
          palette=colors, dodge=False, alpha=0.4, ec='k', linewidth=1)

    plot.set(xlabel='Workflow', ylabel='Stop difference (prop. alignment length)')

    for ax_n in plot.axes:
        for ax in ax_n:
            ax.axhline(0, linewidth=2, color='gray', linestyle="--", alpha=0.6)

    plt.savefig(outpref + '_boxcompare_end.png')

    plt.clf()

    # plot within
    data = data_full[(data_full['Type'] == "within")]

    plot = sns.catplot(
        data=data, x='Tool', y='Diff',
        col='Gene', hue="Tool", kind='box',
        sym="", dodge=False, hue_order=args["hue_order"], order=args["order"]
    )

    plot.map(sns.stripplot, args["x"], args["y"], args["hue"], hue_order=args["hue_order"], order=args["order"],
          palette=colors, dodge=False, alpha=0.4, ec='k', linewidth=1)

    for ax_n in plot.axes:
        for ax in ax_n:
            ax.axhline(0, linewidth=2, color='gray', linestyle="--", alpha=0.6)

    plot.set(xlabel='Workflow', ylabel='Difference in gaps from reference')

    plt.savefig(outpref + '_boxcompare_within.png')

if __name__ == "__main__":
    main()

