#!/usr/bin/env python3

__description__='''

Create .csv states file for PastML.

'''


if __name__ == '__main__':

    import sys
    import argparse
    from readFilesFunctions import *

    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("-msa",type=str,help="MSA in fasta format.")
    app.add_argument("-tree",type=str,help="Tree in .nwk format.")
    app.add_argument("-out",type=str,help="Output file.")
    args=app.parse_args()

    # Get MSA
    msa=readAlign(args.msa,"fasta")
    msadic=align2dic(msa)
    
    # Get tree names
    tree_names=getTreenames(args.tree,"newick")

    # Write output
    out=open(args.out,"w")
    # Header
    out.write("ID")
    for name in msadic:
        seq=msadic[name]
        break
    for col in range(len(seq)):
        out.write(",col"+str(col))
    out.write("\n")
    # Gap/residue as binary 0/1
    for name in tree_names:
        seq=msadic[name]
        out.write(name)
        for col in seq:
            if col=="-":
                state=0
            else:
                state=1
            out.write(","+str(state))
        out.write("\n")
    out.close()
