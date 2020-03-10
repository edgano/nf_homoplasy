#!/usr/bin/env python3

__description__='''

It counts the number of state changes in the tree

'''

def readPastML(fil):
    cols={}
    with open(fil) as f:
        for line in f:
            line=line.strip("\n")
            fields=line.split("\t")
            # Initialize
            if fields[0]=="node":
                ncol=len(fields)-1
                for i in range(ncol):
                    cols[i]=[]
            # Add values
            else:
                for i in range(ncol):
                    cols[i].append(fields[i+1])
    return(cols)


def cols2transition(cols):
    trans={}
    for n in range(len(cols)):
        col=cols[n]
        trans[n]=0
        for i in range(1,len(col)):
            prev,current=col[i-1],col[i]
            if current!=prev:
                trans[n]+=1
    return(trans)


def transition2table(trans,out):
    o=open(out,"w")
    summ=0
    for i in trans:
        o.write("col"+str(i)+"\t"+str(trans[i])+"\n")
        summ+=trans[i]
    o.write("Total\t"+str(summ)+"\n")
    o.close()



if __name__ == '__main__':

    import sys
    import argparse

    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("i",type=str,help="Input file: the output from PastML.")
    app.add_argument("o",type=str,help="Output file.")
    args=app.parse_args()

    cols=readPastML(args.i)
    trans=cols2transition(cols)
    transition2table(trans,args.o)
