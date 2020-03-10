

##
metrics=c("tc","homo","w_homo","w_homo2","ngap","ngap2")


# Import data
fil=read.csv("delta_data.csv",header=TRUE,sep=" ",row.names=NULL)

#fil=read.csv("delta_viceversa_data.csv",header=TRUE,sep=" ",row.names=NULL)
fil$tree=paste(fil$tree1,fil$tree2,sep="_")




# Remove guide trees
rm_tree=function(df,treename)
{
    pos=c(which(df$tree1==treename))
    pos=c(pos,which(df$tree2==treename))
    df2=df[-pos,]
    return(df2)
}

to_be_rm=c("parttreednd2","parttreednd2size")
for (tr in to_be_rm)
{
    fil=rm_tree(fil,tr)
}




# Slope function
get_slope_interc=function(df,metric)
{
    slop=lm(df$tc ~ get(metric,df),na.action=na.exclude)$coefficients[2]
    interc=lm(df$tc ~ get(metric,df),na.action=na.exclude)$coefficients[1]
    return(c(slop,interc))
}



## Get the slope by GROUP (fam, aln, tree)
# Use symmetric viceversa data 

# By FAMILY
fams=levels(fil$family)

slopeArray_by_fam=array(0.0,dim=c(length(fams),length(metrics)-1))
intercArray_by_fam=array(0.0,dim=c(length(fams),length(metrics)-1))
colnames(slopeArray_by_fam)=metrics[2:length(metrics)]
colnames(intercArray_by_fam)=metrics[2:length(metrics)]
for (f in 1:length(fams))
{
    fam=fams[f]
    for (i in 2:length(metrics))
    {
        metric=metrics[i]
        df2=fil[which(fil$family==fam),]
        tmp=get_slope_interc(df2,metric)
        slopeArray_by_fam[f,(i-1)]=tmp[1]
        intercArray_by_fam[f,(i-1)]=tmp[2]
    }
}
slopeDf_by_fam=cbind(data.frame(family=fams),as.data.frame(slopeArray_by_fam))
intercDf_by_fam=cbind(data.frame(family=fams),as.data.frame(intercArray_by_fam))

write.table(slopeDf_by_fam,'/users/cn/sjin/projects/homoplasy/nf_homoplasty/tables/slope_by_fam.csv',sep=" ",row.names=FALSE,quote=FALSE)
write.table(intercDf_by_fam,'/users/cn/sjin/projects/homoplasy/nf_homoplasty/tables/interc_by_fam.csv',sep=" ",row.names=FALSE,quote=FALSE)



# By FAMILY & ALIGNER
fams=levels(fil$family)
aligners=levels(fil$aligner)

slopeDf=data.frame()
intercDf=data.frame()
colnames(slopeArray_by_fam)=c("family","aligner",metrics[2:length(metrics)])
colnames(intercArray_by_fam)=c("family","aligner",metrics[2:length(metrics)])
for (a in 1:length(aligners))
{
    slope_tmp=array(0.0,dim=c(length(fams),length(metrics)-1))
    interc_tmp=array(0.0,dim=c(length(fams),length(metrics)-1))
    colnames(slope_tmp)=metrics[2:length(metrics)]
    colnames(interc_tmp)=metrics[2:length(metrics)]
    aln=aligners[a]
    for (f in 1:length(fams))
    {
        fam=fams[f]
        for (i in 2:length(metrics))
        {
            metric=metrics[i]
            df2=fil[which(fil$family==fam & fil$aligner==aln),]
            tmp=get_slope_interc(df2,metric)
            slope_tmp[f,(i-1)]=tmp[1]
            interc_tmp[f,(i-1)]=tmp[2]
        }
    }
    slope_tmp=cbind(data.frame(family=fams,aligner=rep(aln,length(fams))),as.data.frame(slope_tmp))
    interc_tmp=cbind(data.frame(family=fams,aligner=rep(aln,length(fams))),as.data.frame(interc_tmp))
    slopeDf=rbind(slopeDf,slope_tmp)
    intercDf=rbind(intercDf,interc_tmp)
}
write.table(slopeDf,'/users/cn/sjin/projects/homoplasy/nf_homoplasty/tables/slope_by_fam+aligner.csv',sep=" ",row.names=FALSE,quote=FALSE)
write.table(intercDf,'/users/cn/sjin/projects/homoplasy/nf_homoplasty/tables/interc_by_fam+aligner.csv',sep=" ",row.names=FALSE,quote=FALSE)


# By FAMILY & ALIGNER & TREE
# Use symmetric viceversa data
fams=levels(fil$family)
aligners=levels(fil$aligner)
trees=levels(fil$tree1)

slopeDf=data.frame()
intercDf=data.frame()
colnames(slopeArray_by_fam)=c("family","aligner","tree1",metrics[2:length(metrics)])
colnames(intercArray_by_fam)=c("family","aligner","tree1",metrics[2:length(metrics)])
for (a in 1:length(aligners))
{
    aln=aligners[a]
    for (f in 1:length(fams))
    {
        fam=fams[f]
        for (t in 1:length(trees))
        {
            tr=trees[t]
            slope_tmp=array(0.0,dim=c(1,length(metrics)-1))
            interc_tmp=array(0.0,dim=c(1,length(metrics)-1))
            colnames(slope_tmp)=metrics[2:length(metrics)]
            colnames(interc_tmp)=metrics[2:length(metrics)]
            for (i in 2:length(metrics))
            {
                metric=metrics[i]
                df2=fil[which(fil$family==fam & fil$aligner==aln & fil$tree1==tr),]
                tmp=get_slope_interc(df2,metric)
                slope_tmp[1,(i-1)]=tmp[1]
                interc_tmp[1,(i-1)]=tmp[2]
            }
            slope_tmp=cbind(data.frame(family=fam,aligner=aln,tree1=tr),as.data.frame(slope_tmp))
            interc_tmp=cbind(data.frame(family=fam,aligner=aln,tree1=tr),as.data.frame(interc_tmp))
            slopeDf=rbind(slopeDf,slope_tmp)
            intercDf=rbind(intercDf,interc_tmp)
        }
    }
}
write.table(slopeDf,'/users/cn/sjin/projects/homoplasy/nf_homoplasty/tables/slope_by_fam+aligner+tree.csv',sep=" ",row.names=FALSE,quote=FALSE)
write.table(intercDf,'/users/cn/sjin/projects/homoplasy/nf_homoplasty/tables/interc_by_fam+aligner+tree.csv',sep=" ",row.names=FALSE,quote=FALSE)



