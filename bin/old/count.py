
fname="test_seatoxin_CO/alignments/seatoxin.dpa_1000.CLUSTALO.with.codnd.tree.aln"


def get_aln_len(fname):
    f=open(fname)
    line=f.readline()
    line=f.readline().strip("\n")
    alnlen=len(line)  
    return(alnlen)


def get_nseq(fname):
    nseq=0
    with open(fname) as f:
        for line in f:
            if line[0]==">":
                nseq+=1
    return(nseq)
    

def msa2numpy(fname):
    ''' Read the msa into a numpy matrix '''
    import numpy as np
    nseq=get_nseq(fname)
    alnlen=get_aln_len(fname)
    
    msa=np.zeros( (nseq,alnlen),dtype=str )
    
    i=0
    with open(fname) as f:
        for line in f:
            line=line.strip("\n")
            if line[0]==">":
                continue
            for j in range(alnlen):
                msa[i,j]=line[j]
            i+=1
    return(msa)
    

def count_gap(msa):
    ''' 
    Count gaps:
    ngap:  the number of gaps in the whole alignment
    ngap2: sum of the square of the number of gaps found in each column 
    '''
    ngap=0
    ngap2=0
    dim=msa.shape
    for col in range(dim[1]):
        tmp=0
        for lin in range(dim[0]):
            if msa[lin,col]=='-':
                tmp+=1
        ngap+=tmp
        ngap2+=tmp*tmp 
    return(ngap,ngap2)


msa=msa2numpy(fname)
alnlen=get_aln_len(fname) 
ngap,ngap2=count_gap(msa)


print("LEN: "+str(alnlen-2))   # substract start and end of msa data structure   --- WHY??
print("NGAP: "+str(ngap))
print("NGAP2: "+str(ngap2))

