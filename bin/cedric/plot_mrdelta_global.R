library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)

ori="/users/cn/sjin/projects/homoplasy/nf_homoplasty"
ylim_min=-10
ylim_max=15
xlim_min=0
xlim_max=1
bucket1=100
bucket2=1000
foldname=paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/global/bucket",bucket1,"_",bucket2,sep="")  # Output folder

# Metrics
nmet=5
metrics=c("homo","whomo","whomo2","ngap","ngap2",
            "homoPerLen","whomoPerLen","whomo2PerLen","ngapPerLen","ngap2PerLen",
            "homoByLen","whomoByLen","whomo2ByLen","ngapByLen","ngap2ByLen",
            "homoPerSeq","whomoPerSeq","whomo2PerSeq","ngapPerSeq","ngap2PerSeq",
            "homoBySeq","whomoBySeq","whomo2BySeq","ngapBySeq","ngap2BySeq",
            "homoPerLenSeq","whomoPerLenSeq","whomo2PerLenSeq","ngapPerLenSeq","ngap2PerLenSeq",
            "homoByLenSeq","whomoByLenSeq","whomo2ByLenSeq","ngapByLenSeq","ngap2ByLenSeq")
typs=c("original","perlen","bylen","perseq","byseq","perlenseq","bylenseq")
met=list()
for (i in 1:length(typs))
{
    met[[i]]=metrics[(nmet*(i-1)+1):(nmet*i)]
    #print(c((nmet*(i-1)+1),nmet*i))
}
names(met)=typs



# Read file
f_global1=read.table(paste(ori,"/data/mrdelta/out.bucket",bucket1,".global.tsv",sep=""),header=TRUE)
f_global2=read.table(paste(ori,"/data/mrdelta/out.bucket",bucket2,".global.tsv",sep=""),header=TRUE)
f_global1$metric=factor(f_global1$metric,levels=metrics,ordered=TRUE)
f_global2$metric=factor(f_global2$metric,levels=metrics,ordered=TRUE)


# Functions
get_color=function(metrics)
{
    colorset=hue_pal()(length(metrics))
    names(colorset)=metrics
    return(colorset)
}

plot_vsmrdelta=function(df,typ,metrics,ylim_min,ylim_max)
{
    df2=df[which(df$metric %in% metrics),]
    colorset=get_color(metrics)
    g_deltaacc=ggplot(df2[which(df2$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+
               geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    g_Nratio=ggplot(df2[which(df2$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+
               geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    g_usedRatio=ggplot(df2[which(df2$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+
               geom_point()+geom_line()+ylab("usedRatio")+ylim(0,1)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    g=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,nrow=3)
    ggsave(paste(foldname,"/mrdelta.",typ,".global.png",sep=""),plot=g,width=17.5,height=10.5)
}

plot_vsunused=function(df,typ,metrics,ylim_min,ylim_max)
{
    df2=df[which(df$metric %in% metrics),]
    colorset=get_color(metrics)
    g_deltaacc=ggplot(df2[which(df2$typ=="deltaacc"),],aes(x=unusedRatio,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+
               geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    g_Nratio=ggplot(df2[which(df2$typ=="Nratio"),],aes(x=unusedRatio,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+
               geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    g=ggarrange(g_deltaacc,g_Nratio,common.legend=TRUE,nrow=3)
    ggsave(paste(foldname,"/mrdelta.unused.",typ,".global.png",sep=""),plot=g,width=17.5,height=10.5)
}



## =========================================
## Plot Nratio/deltaacc/usedRatio vs mrdelta

df_global1=data.frame(value=c(f_global1$Nratio,f_global1$deltaacc,1-f_global1$unusedRatio),
                         mrdelta=rep(f_global1$mrdelta,3),
                         score=rep(as.character(f_global1$score),3),
                         metric=rep(as.character(f_global1$metric),3),
                         bucket=rep("bucket100",3*nrow(f_global1)),
                         typ=c(rep("Nratio",nrow(f_global1)),rep("deltaacc",nrow(f_global1)),rep("usedRatio",nrow(f_global1)))
                        )
df_global2=data.frame(value=c(f_global2$Nratio,f_global2$deltaacc,1-f_global2$unusedRatio),
                         mrdelta=rep(f_global2$mrdelta,3),
                         score=rep(as.character(f_global2$score),3),
                         metric=rep(as.character(f_global2$metric),3),
                         bucket=rep("bucket1000",3*nrow(f_global2)),
                         typ=c(rep("Nratio",nrow(f_global2)),rep("deltaacc",nrow(f_global2)),rep("usedRatio",nrow(f_global2)))
                        )
df_global=rbind(df_global1,df_global2)

for (typ in typs)
{
    m=get(typ,met)
    plot_vsmrdelta(df_global,typ,m,ylim_min,ylim_max)
}





## ==============================
## Nratio/deltaacc vs unusedRatio

df_global1=data.frame(value=c(f_global1$Nratio,f_global1$deltaacc),
                         mrdelta=rep(f_global1$mrdelta,2),
                         score=rep(as.character(f_global1$score),2),
                         metric=rep(as.character(f_global1$metric),2),
                         unusedRatio=rep(f_global1$unusedRatio,2),
                         bucket=rep("bucket100",2*nrow(f_global1)),
                         typ=c(rep("Nratio",nrow(f_global1)),rep("deltaacc",nrow(f_global1)))
                        )
df_global2=data.frame(value=c(f_global2$Nratio,f_global2$deltaacc),
                         mrdelta=rep(f_global2$mrdelta,2),
                         score=rep(as.character(f_global2$score),2),
                         metric=rep(as.character(f_global1$metric),2),
                         unusedRatio=rep(f_global2$unusedRatio,2),
                         bucket=rep("bucket1000",2*nrow(f_global2)),
                         typ=c(rep("Nratio",nrow(f_global2)),rep("deltaacc",nrow(f_global2)))
                        )
df_global=rbind(df_global1,df_global2)

for (typ in typs)
{
    m=get(typ,met)
    plot_vsunused(df_global,typ,m,ylim_min,ylim_max)
}





