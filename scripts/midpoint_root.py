import dendropy
import argparse

def get_options():
    description = 'Generates midpoint rooted tree'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python midpoint_root.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='Input nwk tree ')
    IO.add_argument('--outfile',
                    default="midpoint_root.nwk",
                    help='Output filename ')
    return parser.parse_args()

def main():
    # parse command line arguments for ggCaller
    options = get_options()
    infile = options.infile
    outfile = options.outfile

    tree = dendropy.Tree.get(path=infile, schema="newick")
    tree.reroot_at_midpoint(update_bipartitions=False)
    tree.write(path=outfile, schema="newick")

if __name__ == '__main__':
    main()