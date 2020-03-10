
##
metrics=c("tc","homo","w_homo","w_homo2","ngap","ngap2")


# Import data
#fil=read.csv("/users/cn/sjin/projects/homoplasy/nf_homoplasty/data/data_14tree_above1000/delta_data.csv",header=TRUE,sep=" ",row.names=NULL)
fil=read.csv("/users/cn/sjin/projects/homoplasy/nf_homoplasty/data/data_14tree_above1000/delta_symmetric_data.csv",header=TRUE,sep=" ",row.names=NULL)

fil$tree=paste(fil$tree1,fil$tree2,sep="_")



# Remove guide trees
rm_tree=function(df,treename)
{
    pos=c(which(df$tree1==treename))
    pos=c(pos,which(df$tree2==treename))
    df2=df[-pos,]
    return(df2)
}

# to_be_rm=c("parttreednd2","parttreednd2size")
# for (tr in to_be_rm)
# {
#     fil=rm_tree(fil,tr)
# }




# Function
count_quadrants=function(df,metric,score)
{
    pos=which(abs(get(metric,df))<=0.001 & abs(get(score,df))<=0.001) # & is.na(get(metric,df)) & is.na(get(score,df)))
    unused=length(pos)
    df=df[-pos,]
    q=c()
    ## Quadrant I (upper right): +,+
    pos=which(get(metric,df)>0 & get(score,df)>0)
    q=c(q,length(pos))
    ## Quadrant II (upper left): -,+
    pos=which(get(metric,df)<0 & get(score,df)>0)
    q=c(q,length(pos))
    ## Quadrant III (bottom left): -,-
    pos=which(get(metric,df)<0 & get(score,df)<0)
    q=c(q,length(pos))
    ## Quadrant IV (bottom right): +,-
    pos=which(get(metric,df)>0 & get(score,df)<0)
    q=c(q,length(pos))
    q=c(q,unused)
    return(q)
}




# Count by FAMILY
fams=levels(fil$family)
score="tc"
countList=list()
n=0
for (metric in metrics[2:length(metrics)])
{
    q1=c();q2=c();q3=c();q4=c()
    for (fam in fams)
    {   
        df=fil[which(fil$family==fam),]
        q=count_quadrants(df,metric,score)
        q1=c(q1,q[1]); q2=c(q2,q[2]); q3=c(q3,q[3]); q4=c(q4,q[4])
    }
    positive=q1+q3
    negative=q2+q4
    df2=data.frame(family=fams,q1=q1,q2=q2,q3=q3,q4=q4,positive=positive,negative=negative,diff=positive-negative,np_ratio=negative/positive)
    n=n+1
    countList[[n]]=df2
}
names(countList)=metrics[2:length(metrics)]

# for (metric in metrics[2:length(metrics)])
# {
#     write.table(get(metric,countList),paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/tables/tables_14trees_above1000/quadrant_counts_by_fam.",metric,".csv",sep=""),sep=" ",row.names=FALSE,quote=FALSE)
# }



# Count by ALIGNER
aligners=levels(fil$aligner)
score="tc"
countList=list()
n=0
for (metric in metrics[2:length(metrics)])
{
    q1=c();q2=c();q3=c();q4=c();unused=c()
    for (aln in aligners)
    {   
        df=fil[which(fil$aligner==aln),]
        q=count_quadrants(df,metric,score)
        q1=c(q1,q[1]); q2=c(q2,q[2]); q3=c(q3,q[3]); q4=c(q4,q[4]); unused=c(unused,q[5])
    }
    positive=q1+q3
    negative=q2+q4
    df2=data.frame(aligner=aligners,q1=q1,q2=q2,q3=q3,q4=q4,positive=positive,negative=negative,diff=positive-negative,np_ratio=negative/positive,nperc=negative/(q1+q2+q3+q4),unused=unused/nrow(df)
    n=n+1
    countList[[n]]=df2
}
names(countList)=metrics[2:length(metrics)]





# Count by FAMILY & ALIGNER
Sys.setlocale(, "C")
fams=sort(levels(fil$family))
aligners=levels(fil$aligner)
score="sp"
countList=list()
n=0
for (metric in metrics[2:length(metrics)])
{
    q1=c();q2=c();q3=c();q4=c();al=c();unused=c()
    for (aln in aligners)
    {
        for (fam in fams)
        {
        df=fil[which(fil$family==fam & fil$aligner==aln),]
        q=count_quadrants(df,metric,score)
        q1=c(q1,q[1]); q2=c(q2,q[2]); q3=c(q3,q[3]); q4=c(q4,q[4]); unused=c(unused,q[5])
        }
        al=c(al,rep(aln,length(fams)))
    }
    positive=q1+q3
    negative=q2+q4
    df2=data.frame(family=rep(fams,length(aligners)),aligner=al,q1=q1,q2=q2,q3=q3,q4=q4,positive=positive,negative=negative,diff=positive-negative,np_ratio=negative/positive,nperc=negative/(q1+q2+q3+q4),used=(positive+negative),unused=unused)
    n=n+1
    countList[[n]]=df2
}
names(countList)=metrics[2:length(metrics)]

# for (metric in metrics[2:length(metrics)])
# {
#     write.table(get(metric,countList),paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/tables/tables_14trees_above1000/quadrant_counts_by_fam+aligner.",metric,".csv",sep=""),sep=" ",row.names=FALSE,quote=FALSE)
# }



# Count Positive and Negative quadrants

