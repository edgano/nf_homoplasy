import os,sys
import argparse
from IPython.display import display, HTML
import pandas as pd

bucket=sys.argv[1]  # Bucket size: 100, 1000, etc
datatyp=sys.argv[2]   # "raw" | "delta"
folder="data/"+datatyp+"/bucket"+bucket  # Where (homo,whomo,whomo2,ngap,ngap2,tc,sp).csv files are stored
excelfolder="excels/"  # Output will be stored here

# Metrics
scores=["tc","sp"]
metrics=["homo","whomo","whomo2","ngap","ngap2"]
otherparam=["col","len"]

# Load CSV to pandas
csvdic={}
csvdic["tc"]=pd.read_csv(folder+"/tc.csv",header=[0,1],index_col=0)
csvdic["sp"]=pd.read_csv(folder+"/sp.csv",header=[0,1],index_col=0)
csvdic["col"]=pd.read_csv(folder+"/col.csv",header=[0,1],index_col=0)
csvdic["len"]=pd.read_csv(folder+"/len.csv",header=[0,1],index_col=0)

csvdic["homo"]=pd.read_csv(folder+"/homo.csv",header=[0,1],index_col=0)
csvdic["whomo"]=pd.read_csv(folder+"/whomo.csv",header=[0,1],index_col=0)
csvdic["whomo2"]=pd.read_csv(folder+"/whomo2.csv",header=[0,1],index_col=0)
csvdic["ngap"]=pd.read_csv(folder+"/ngap.csv",header=[0,1],index_col=0)
csvdic["ngap2"]=pd.read_csv(folder+"/ngap2.csv",header=[0,1],index_col=0)

# Create nseq df
tmp=csvdic["tc"].copy()
n=tmp.shape[1]
for i in range(1,n):
    tmp[[i]]=tmp.iloc[:,0]
csvdic["nseq"]=tmp


# Normalize data: per len & by len
csvdic_normPerLen={}
csvdic_normByLen={}
for i in metrics:
    csv=csvdic[i]
    length=csvdic["len"]
    # Per len
    tmp=csv/length
    tmp["nseq"]=csv["nseq"]
    csvdic_normPerLen[i]=tmp
    # By len
    tmp=csv*length
    tmp["nseq"]=csv["nseq"]
    csvdic_normByLen[i]=tmp
for i in scores+otherparam:
    csvdic_normPerLen[i]=csvdic[i]
    csvdic_normByLen[i]=csvdic[i]
    

# Normalize data: per nseq & by nseq
csvdic_normPerSeq={}
csvdic_normBySeq={}
for i in metrics:
    csv=csvdic[i]
    nseq=csvdic["nseq"]
    # Per nseq
    tmp=csv/nseq
    tmp["nseq"]=csv["nseq"]
    csvdic_normPerSeq[i]=tmp
    # By nseq
    tmp=csv*nseq
    tmp["nseq"]=csv["nseq"]
    csvdic_normBySeq[i]=tmp
for i in scores+otherparam:
    csvdic_normPerSeq[i]=csvdic[i]
    csvdic_normBySeq[i]=csvdic[i]
    

# Normalize data: per/by len & nseq
csvdic_normPerLenSeq={}
csvdic_normByLenSeq={}
for i in metrics:
    csv=csvdic[i]
    length=csvdic["len"]
    nseq=csvdic["nseq"]
    # Per len & nseq
    tmp=csv/length/nseq
    tmp["nseq"]=csv["nseq"]
    csvdic_normPerLenSeq[i]=tmp
    # By len & nseq
    tmp=csv*length*nseq
    tmp["nseq"]=csv["nseq"]
    csvdic_normByLenSeq[i]=tmp
for i in scores+otherparam:
    csvdic_normPerLenSeq[i]=csvdic[i]
    csvdic_normByLenSeq[i]=csvdic[i]
    

# Change to excel directory
os.chdir("excels/"+datatyp+"/bucket"+str(bucket))

dics={"original":csvdic,
      "normPerLen":csvdic_normPerLen,"normByLen":csvdic_normByLen,
      "normPerSeq":csvdic_normPerSeq,"normBySeq":csvdic_normBySeq,
      "normPerLenSeq":csvdic_normPerLenSeq,"normByLenSeq":csvdic_normByLenSeq}

# Write excels
for typ,d in dics.items():
    out="homoplasy."+typ+".xlsx"
    with pd.ExcelWriter(out) as writer: 
        for i in scores+metrics+otherparam:
            df=d[i]
            pd.DataFrame([("Family",)+x for x in zip(*df.columns)]).to_excel(writer,sheet_name=i,header=False,index=False)
            df.to_excel(writer,sheet_name=i,header=False,startrow=1,na_rep="NA")


# Change to csv directory
os.chdir("../../../data/"+datatyp+"/bucket"+str(bucket))
   
# Write csv
for typ,d in dics.items():
    for i in scores+metrics+otherparam:
        df=d[i]
        out=i+"."+typ+".csv"
        with open(out,"w") as f:
            f.write('\n'.join([','.join(h) for h in [("Family",) + x for x in zip(*df.columns)] ]) + '\n')
        df.to_csv(out,mode="a",na_rep="NA",header=False)
        

