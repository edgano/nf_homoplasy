import sys
import pandas as pd

folder=sys.argv[1]  
typ=sys.argv[2]

out=folder+"/delta."+typ+".csv"
o=open(out,"w")

# Files path
metrics=["tc","sp","homo","whomo","whomo2","ngap","ngap2"]
files={"sp":folder+"/sp."+typ+".csv",
        "tc":folder+"/tc."+typ+".csv",
        #"col":folder+"/col."+typ+".csv",
        "homo":folder+"/homo."+typ+".csv",
        "whomo":folder+"/whomo."+typ+".csv",
        "whomo2":folder+"/whomo2."+typ+".csv",
        #"len":folder+"/len."+typ+".csv",
        "ngap":folder+"/ngap."+typ+".csv",
        "ngap2":folder+"/ngap2."+typ+".csv"}

# Read data frame
dfs={}
for name,fil in files.items():
    dfs[name]=pd.read_csv(fil,header=[0,1],index_col=0)
    dfs[name]=dfs[name].fillna("NA")

# Get headers
families=list(dfs["sp"].index)
aligners_rep=[x[0] for x in dfs["sp"].columns][1:]
trees_rep=[x[1] for x in dfs["sp"].columns][1:]
aligners,trees=[],[]
for x in aligners_rep:
    if x not in aligners:
        aligners.append(x)
for x in trees_rep:
    if x not in trees:
        trees.append(x)

# # Prepare data frame
# o.write("family aligner tree1 tree2 tc sp homo whomo whomo2 ngap ngap2\n")
# for fam in families:
#     for aligner in aligners:
#         for tree in trees:
#             tree1=tree.split(".")[0]
#             tree2=tree.split(".")[1]
#             o.write("{} {} {} {}".format(fam,aligner,tree1,tree2))
#             for metric in metrics:
#                 df=dfs[metric]
#                 value=df[aligner].loc[fam,tree]
#                 o.write(" "+str(value))
#             o.write("\n")
            

# Prepare data frame
df={}
df["family"]=[item for item in families for i in range(len(trees_rep))]
df["aligner"]=aligners_rep*len(families)
tree_col=trees_rep*len(families)
df["tree1"]=[x.split(".")[0] for x in tree_col]
df["tree2"]=[x.split(".")[1] for x in tree_col]
for metric in metrics:
    df[metric]=[]
    for fam in families:
        df[metric]+=list(dfs[metric].loc[fam,:])[1:]
l=[]
for x in ["family","aligner","tree1","tree2"]+metrics:
    l.append(df[x])
lt=map(list, zip(*l))

o.write(" ".join(["family","aligner","tree1","tree2"]+metrics) +"\n")
for i in lt:
    tmp=" ".join([str(x) for x in i])
    o.write(tmp+"\n")
o.close()

