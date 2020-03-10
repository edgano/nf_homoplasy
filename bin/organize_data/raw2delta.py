
def raw2delta(fil):

    import pandas as pd
    import numpy as np

    incsv=pd.read_csv(fil,header=[0,1],index_col=0)
    families,aligners,trees,mixtrees,nseqs=get_headers(incsv)
    head=create_outheader(aligners,mixtrees)

    delt=[]
    delt.append(families)
    delt.append(nseqs)
    for aln in aligners:
        for tr in mixtrees:
            tr1=tr.split(".")[0]
            tr2=tr.split(".")[1]
            val1=incsv[aln].loc[:,tr1]
            val2=incsv[aln].loc[:,tr2]
            delta=val1-val2
            delt.append(delta)
    delt=map(list,zip(*delt))  # Transpose
    #=head+delt  # List of list with header and delta values
    df=pd.DataFrame(delt)  # Pandas data frame

    return(head,df)  


def get_headers(incsv):
    
    # Get headers: families, aligners, trees
    families=list(incsv.index)
    nseqs=[x[0] for x in incsv["nseq"].values.tolist()]
    aligners_rep=[x[0] for x in list(incsv.columns)][1:]  # Note that they are repeated (just as in the csv file)
    trees_rep=[x[1] for x in list(incsv.columns)][1:]
    aligners,trees=[],[]
    for x in aligners_rep:
        if x not in aligners:
            aligners.append(x)
    for x in trees_rep:
        if x not in trees:
            trees.append(x)

    # Get mixed trees: tree1.tree2
    mixtrees=[]
    for tr1 in trees:
        for tr2 in trees:
            if tr1==tr2:
                continue
            tr12=tr1+"."+tr2
            tr21=tr2+"."+tr1
            if tr12 not in mixtrees and tr21 not in mixtrees:   
                mixtrees.append(tr12)
    
    return(families,aligners,trees,mixtrees,nseqs)


def create_outheader(aligners,mixtrees):
    
    head=[["Family","nseq"]]
    for aln in aligners:
        for tr in mixtrees:
            head[-1].append(aln)
    head.append(["Family","nseq"])
    for aln in aligners:
        for tr in mixtrees:
            head[-1].append(tr)
    
    return(head)


if __name__ == '__main__':

    import sys
    import argparse
    import pandas as pd

    #fil="test.csv"
    #out="delta-test.csv"

    app = argparse.ArgumentParser()
    app.add_argument('input',type=str,help="Input raw.csv filename.")
    app.add_argument('output',type=str,help="Output delta.csv filename.")
    args = app.parse_args()

    fil=args.input
    out=args.output

    head,df=raw2delta(fil)
    o=open(out,"w")
    for i in head:
        o.write(",".join(i)+"\n")
    o.close()
    df.to_csv(out,header=False,mode="a",index=False,na_rep="NA",float_format="%g")

