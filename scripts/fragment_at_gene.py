import os

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import listdir
from os.path import isfile, join
import random

def get_options():
    description = 'Fragments sequence at each gene location. Can also be used to translate CDS sequences.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python fragment_at_gene.py')
    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='Input fasta file. ')
    IO.add_argument('--CDS',
                    required=True,
                    help='Fasta file containing CDS sequences. Each sequence header must be in format'
                         '<contig_id>_<CDS_ID>, with contig_ID matching that in infile fasta')
    IO.add_argument('--translate',
                    action='store_true',
                    default=False,
                    help='Translate sequences in infile')
    return parser.parse_args()

def make_contigs(CDS_file, infile):
    CDS_list = []
    # get all known CDS sequences
    for seq_record in SeqIO.parse(CDS_file, "fasta"):
        ID = seq_record.id
        ID = ID.split("_")[0]
        CDS_list.append((ID, str(seq_record.seq), len(seq_record.seq)))

    #read in ref files

    ref_contigs = []
    for seq_record in SeqIO.parse(infile, "fasta"):
        ID = seq_record.id
        ID = ID.split("_")[0]
        seq = str(seq_record.seq)

        # get coordinates of CDS
        contig_index = 0
        for CDS in CDS_list:
            CDS_ID, CDS_seq, CDS_len = CDS
            if CDS_ID == ID:
                coord = seq.find(CDS_seq)
                if coord == -1:
                    coord = seq.find(str(Seq(CDS_seq).reverse_complement()))
                if coord == -1:
                    continue
                # generate a random break within the CDS
                break_coord = random.randint(coord, coord + CDS_len)
                ref_contigs.append(
                    SeqRecord(Seq(seq[0:break_coord + 1]), id=ID + "_" + str(contig_index),
                              description=""))
                contig_index += 1
                seq = seq[break_coord + 1:]

        # append the last segment
        ref_contigs.append(
            SeqRecord(Seq(seq), id=ID + "_" + str(contig_index),
                      description=""))

        outfile = os.path.splitext(infile)[0] + "_fragmented.fasta"
        SeqIO.write(ref_contigs, outfile, "fasta")

def translate(infile):
    CDS_list = []
    outfile = os.path.splitext(infile)[0] + "translated.fasta"
    # get all known CDS sequences
    for seq_record in SeqIO.parse(infile, "fasta"):
        CDS_list.append(SeqRecord(seq_record.seq.translate(), id=seq_record.id,
                  description=seq_record.description))
    SeqIO.write(CDS_list, outfile, "fasta")

def main():
    options = get_options()
    infile = options.infile
    CDS_file = options.CDS
    if options.translate:
        translate(infile)
    else:
        make_contigs(CDS_file, infile)

if __name__ == "__main__":
    main()
