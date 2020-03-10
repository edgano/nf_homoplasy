#!/usr/bin/env python3

__description__='''

Python version of the Cedric's script for the analysis of homoplasy.

'''

def csv2pandas(fil):
    import pandas as pd
    df=pd.read_csv(fil,header=[0,1],index_col=0)
    return(df)

def pandas2headers(df):
    # Get headers: families, aligners, trees
    families=list(df.index)
    nseqs=[x[0] for x in df["nseq"].values.tolist()]
    aligners_rep=[x[0] for x in list(df.columns)][1:]  # Note that they are repeated (just as in the csv file)
    trees_rep=[x[1] for x in list(df.columns)][1:]
    aligners,trees=[],[]
    for x in aligners_rep:
        if x not in aligners:
            aligners.append(x)
    for x in trees_rep:
        if x not in trees:
            trees.append(x)
    return(families,nseqs,aligners,trees)

def csv2analyse(metricdf,scoredf,families,nseqs,aligners,trees):
    total={"fam":[],"aln":[],"global":0}
    unused={"fam":0,"aln":0,"global":0}
    for aligner in aligners:
        total["aln"]+=1
        for family in families:
            for tr1 in trees:
                for tr2 in trees:
                    # Metrics
                    h1=metricdf[aligner][tr1][family]
                    h2=metricdf[aligner][tr2][family]
                    D1=h1-h2
                    # Score
                    s1=scoredf[aligner][tr1][family]
                    s2=scoredf[aligner][tr2][family]
                    D2=s1-s2




if __name__ == '__main__':

    import sys
    import argparse

    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("-dir",type=str,help="Directory where the input (tc,sp,homo,whomo,whomo2,ngap,ngap2).csv files are stored.")
    app.add_argument("-score",type=str,choices=["tc","sp"])
    app.add_argument("-metrics",type=str,choices=["homo","whomo","whomo2","ngap","ngap2"])
    app.add_argument("-norm",type=str,choices=["PerLen","PerSeq","PerLenSeq","ByLen","BySeq","ByLenSeq"],default=None,help="Data normalization. \
                    Divided by length, number of sequence, or both. Or multiplied by Length, number of sequence, or both.")
    app.add_argument("-mrdelta",type=float,default=0)
    app.add_argument("-minseq",type=int,default=0)
    app.add_argument("-maxseq",type=int,default=100000)
    args=app.parse_args()

    # Read csv to panda Data Frames
    metricfile=args.dir+"/"+args.metrics+".norm"+args.norm+".csv"
    scorefile==args.dir+"/"+args.score+".csv"
    metricdf=csv2pandas(metricfile)
    scoredf=csv2pandas(scorefile)
    # Headers & index
    families,nseqs,aligners,trees=pandas2headers(metricdf)




