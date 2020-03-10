import sys
import argparse
import pandas as pd 

app = argparse.ArgumentParser()
app.add_argument('folder',type=str,help="Folder where the input (homo,whomo,whomo2,ngap,ngap2,tc,sp).csv files are stored.")
app.add_argument('typ',type=str,choices=['','original','normPerLen','normPerSeq','normPerLenSeq','normByLen','normBySeq','normByLenSeq'])
app.add_argument('out',type=str,help="Output .xlsx file.")
args = app.parse_args()

folder=args.folder
typ=args.typ
out=args.out

# Files path
files={}
if typ=="":
    metrics=["tc","sp","col","len","homo","whomo","whomo2","ngap","ngap2"]
    for i in metrics:
        files[i]=folder+"/"+i+".csv"
else:
    metrics=["tc","sp","homo","whomo","whomo2","ngap","ngap2"]
    for i in metrics:
        files[i]=folder+"/"+i+"."+typ+".csv"


# Read data frame
dfs={}
for name,fil in files.items():
    dfs[name]=pd.read_csv(fil,header=[0,1],index_col=0)

# Write excel
with pd.ExcelWriter(out) as writer: 
    for i in metrics:
        df=dfs[i]
        pd.DataFrame([("Family",)+x for x in zip(*df.columns)]).to_excel(writer,sheet_name=i,header=False,index=False)
        df.to_excel(writer,sheet_name=i,header=False,startrow=1,na_rep="NA")

