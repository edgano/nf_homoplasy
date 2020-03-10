
library(ggplot2)
library(scales)

bucket=1000

foldname="deltaacc"
x_axis="mrdelta"
y_axis="deltaacc"
ylim_min=-1
ylim_max=15


metrics=c("homo","whomo","whomo2","ngap","ngap2",
            "homoPerLen","whomoPerLen","whomo2PerLen","ngapPerLen","ngap2PerLen",
            "homoByLen","whomoByLen","whomo2ByLen","ngapByLen","ngap2ByLen",
            "homoPerSeq","whomoPerSeq","whomo2PerSeq","ngapPerSeq","ngap2PerSeq",
            "homoBySeq","whomoBySeq","whomo2BySeq","ngapBySeq","ngap2BySeq",
            "homoPerLenSeq","whomoPerLenSeq","whomo2PerLenSeq","ngapPerLenSeq","ngap2PerLenSeq",
            "homoByLenSeq","whomoByLenSeq","whomo2ByLenSeq","ngapByLenSeq","ngap2ByLenSeq")
metrics_original=c("homo","whomo","whomo2","ngap","ngap2")
metrics_perlen=c("homoPerLen","whomoPerLen","whomo2PerLen","ngapPerLen","ngap2PerLen")
metrics_bylen=c("homoByLen","whomoByLen","whomo2ByLen","ngapByLen","ngap2ByLen")
metrics_perseq=c("homoPerSeq","whomoPerSeq","whomo2PerSeq","ngapPerSeq","ngap2PerSeq")
metrics_byseq=c("homoBySeq","whomoBySeq","whomo2BySeq","ngapBySeq","ngap2BySeq")
metrics_perlenseq=c("homoPerLenSeq","whomoPerLenSeq","whomo2PerLenSeq","ngapPerLenSeq","ngap2PerLenSeq")
metrics_bylenseq=c("homoByLenSeq","whomoByLenSeq","whomo2ByLenSeq","ngapByLenSeq","ngap2ByLenSeq")


f_fam=read.table(paste("out.bucket",bucket,".fam.tsv",sep=""),header=TRUE)
f_aln=read.table(paste("out.bucket",bucket,".aln.tsv",sep=""),header=TRUE)
f_global=read.table(paste("out.bucket",bucket,".global.tsv",sep=""),header=TRUE)


# Data frames
f_aln_sp=f_aln[which(f_aln$score=="sp"),]
f_aln_tc=f_aln[which(f_aln$score=="tc"),]
f_global_sp=f_global[which(f_global$score=="sp"),]
f_global_tc=f_global[which(f_global$score=="tc"),]

df_aln_sp=data.frame(mrdelta=c(f_global_sp$mrdelta,f_aln_sp$mrdelta),
                    metric=c(as.character(f_global_sp$metric),as.character(f_aln_sp$metric)),
                    aligner=c(as.character(f_global_sp$aligner),as.character(f_aln_sp$aligner)),
                    Nratio=c(f_global_sp$Nratio,f_aln_sp$Nratio),
                    deltaacc=c(f_global_sp$deltaacc,f_aln_sp$deltaacc),
                    minacc=c(f_global_sp$minacc,f_aln_sp$minacc),
                    maxacc=c(f_global_sp$maxacc,f_aln_sp$maxacc),
                    unusedRatio=c(f_global_sp$unusedRatio,f_aln_sp$unusedRatio))
df_aln_tc=data.frame(mrdelta=c(f_global_tc$mrdelta,f_aln_tc$mrdelta),
                    metric=c(as.character(f_global_tc$metric),as.character(f_aln_tc$metric)),
                    aligner=c(as.character(f_global_tc$aligner),as.character(f_aln_tc$aligner)),
                    Nratio=c(f_global_tc$Nratio,f_aln_tc$Nratio),
                    deltaacc=c(f_global_tc$deltaacc,f_aln_tc$deltaacc),
                    minacc=c(f_global_tc$minacc,f_aln_tc$minacc),
                    maxacc=c(f_global_tc$maxacc,f_aln_tc$maxacc),
                    unusedRatio=c(f_global_tc$unusedRatio,f_aln_tc$unusedRatio))
df_aln_sp$metric=factor(df_aln_sp$metric,levels=metrics)
df_aln_tc$metric=factor(df_aln_tc$metric,levels=metrics)
# Original metrics
df_aln_sp_original=df_aln_sp[which(df_aln_sp$metric %in% metrics_original),]
df_aln_tc_original=df_aln_tc[which(df_aln_tc$metric %in% metrics_original),]
# Metrics normalized per length
df_aln_sp_perlen=df_aln_sp[which(df_aln_sp$metric %in% metrics_perlen),]
df_aln_tc_perlen=df_aln_tc[which(df_aln_tc$metric %in% metrics_perlen),]
# Metrics normalized by length
df_aln_sp_bylen=df_aln_sp[which(df_aln_sp$metric %in% metrics_bylen),]
df_aln_tc_bylen=df_aln_tc[which(df_aln_tc$metric %in% metrics_bylen),]
# Metrics normalized per nseq
df_aln_sp_perseq=df_aln_sp[which(df_aln_sp$metric %in% metrics_perseq),]
df_aln_tc_perseq=df_aln_tc[which(df_aln_tc$metric %in% metrics_perseq),]
# Metrics normalized by nseq
df_aln_sp_byseq=df_aln_sp[which(df_aln_sp$metric %in% metrics_byseq),]
df_aln_tc_byseq=df_aln_tc[which(df_aln_tc$metric %in% metrics_byseq),]
# Metrics normalized per length & nseq
df_aln_sp_perlenseq=df_aln_sp[which(df_aln_sp$metric %in% metrics_perlenseq),]
df_aln_tc_perlenseq=df_aln_tc[which(df_aln_tc$metric %in% metrics_perlenseq),]
# Metrics normalized by length & nseq
df_aln_sp_bylenseq=df_aln_sp[which(df_aln_sp$metric %in% metrics_bylenseq),]
df_aln_tc_bylenseq=df_aln_tc[which(df_aln_tc$metric %in% metrics_bylenseq),]



# Plot by aligner

color_axis="aligner"

# # Using all metrics
# g_aln_sp=ggplot(df_aln_sp,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric,ncol=5)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
# g_aln_tc=ggplot(df_aln_tc,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric,ncol=5)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
# g_aln_sp
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.allmetrics.aln.png",sep=""))
# g_aln_tc
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.allmetrics.aln.png",sep=""))

# Using original metrics
g_aln_sp_original=ggplot(df_aln_sp_original,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_original=ggplot(df_aln_tc_original,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_original
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.originalmetrics.aln.png",sep=""))
g_aln_tc_original
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.originalmetrics.aln.png",sep=""))
# Using metrics normalized per length
g_aln_sp_perlen=ggplot(df_aln_sp_perlen,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_perlen=ggplot(df_aln_tc_perlen,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_perlen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricperlen.aln.png",sep=""))
g_aln_tc_perlen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricperlen.aln.png",sep=""))
# Using metrics multiplied by length
g_aln_sp_bylen=ggplot(df_aln_sp_bylen,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_bylen=ggplot(df_aln_tc_bylen,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_bylen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricbylen.aln.png",sep=""))
g_aln_tc_bylen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricbylen.aln.png",sep=""))
# Using metrics normalized per nseq
g_aln_sp_perseq=ggplot(df_aln_sp_perseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_perseq=ggplot(df_aln_tc_perseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_perseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricperseq.aln.png",sep=""))
g_aln_tc_perseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricperseq.aln.png",sep=""))
# Using metrics multiplied by nseq
g_aln_sp_byseq=ggplot(df_aln_sp_byseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_byseq=ggplot(df_aln_tc_byseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_byseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricbyseq.aln.png",sep=""))
g_aln_tc_byseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricbyseq.aln.png",sep=""))
# Using metrics normalized per length & nseq
g_aln_sp_perlenseq=ggplot(df_aln_sp_perlenseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_perlenseq=ggplot(df_aln_tc_perlenseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_perlenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricperlenseq.aln.png",sep=""))
g_aln_tc_perlenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricperlenseq.aln.png",sep=""))
# Using metrics multiplied by length & nseq
g_aln_sp_bylenseq=ggplot(df_aln_sp_bylenseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_bylenseq=ggplot(df_aln_tc_bylenseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_bylenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricbylenseq.aln.png",sep=""))
g_aln_tc_bylenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricbylenseq.aln.png",sep=""))




# Plot by metric

color_axis="metric"

# # All metric
# g_aln_sp=ggplot(df_aln_sp,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
# g_aln_tc=ggplot(df_aln_tc,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
# g_aln_sp
# ggsave(paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.allmetrics.bymetric.aln.png",sep=""))
# g_aln_tc
# ggsave(paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.allmetrics.bymetric.aln.png",sep=""))

# Original metrics
g_aln_sp_original=ggplot(df_aln_sp_original,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_original=ggplot(df_aln_tc_original,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_original
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.originalmetrics.bymetric.aln.png",sep=""))
g_aln_tc_original
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.originalmetrics.bymetric.aln.png",sep=""))

# Metrics normalized per length
g_aln_sp_perlen=ggplot(df_aln_sp_perlen,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_perlen=ggplot(df_aln_tc_perlen,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_perlen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricperlen.bymetric.aln.png",sep=""))
g_aln_tc_perlen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricperlen.bymetric.aln.png",sep=""))

# Metrics multiplied by length
g_aln_sp_bylen=ggplot(df_aln_sp_bylen,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_bylen=ggplot(df_aln_tc_bylen,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_bylen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricbylen.bymetric.aln.png",sep=""))
g_aln_tc_bylen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricbylen.bymetric.aln.png",sep=""))

# Metrics normalized per nseq
g_aln_sp_perseq=ggplot(df_aln_sp_perseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_perseq=ggplot(df_aln_tc_perseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_perseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricperseq.bymetric.aln.png",sep=""))
g_aln_tc_perseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricperseq.bymetric.aln.png",sep=""))

# Metrics multiplied by nseq
g_aln_sp_byseq=ggplot(df_aln_sp_byseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_byseq=ggplot(df_aln_tc_byseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_byseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricbyseq.bymetric.aln.png",sep=""))
g_aln_tc_byseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricbyseq.bymetric.aln.png",sep=""))

# Metrics normalized per length & nseq
g_aln_sp_perlenseq=ggplot(df_aln_sp_perlenseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_perlenseq=ggplot(df_aln_tc_perlenseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_perlenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricperlenseq.bymetric.aln.png",sep=""))
g_aln_tc_perlenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricperlenseq.bymetric.aln.png",sep=""))

# Metrics multiplied by length & nseq
g_aln_sp_bylenseq=ggplot(df_aln_sp_bylenseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_tc_bylenseq=ggplot(df_aln_tc_bylenseq,aes_string(x=x_axis,y=y_axis,color=color_axis))+facet_wrap(~aligner)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
g_aln_sp_bylenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.sp.metricbylenseq.bymetric.aln.png",sep=""))
g_aln_tc_bylenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.tc.metricbylenseq.bymetric.aln.png",sep=""))








# By aligner. Nratio & unusedRatio
rat=y_axis
x_axis="mrdelta"
y_axis="ratio"
color_axis="aligner"
shape_axis="type"



# Data frame
if (rat=="deltaacc")
{
    df_ratio_sp=data.frame(mrdelta=rep(df_aln_sp$mrdelta,2),
                        aligner=rep(df_aln_sp$aligner,2),
                        metric=rep(df_aln_sp$metric,2),
                        ratio=c(get(rat,df_aln_sp),rescale(1-df_aln_sp$unusedRatio,to=c(ylim_min,ylim_max),from=c(0,1))),
                        type=c( rep(rat,nrow(df_aln_sp)),rep("usedRatio",nrow(df_aln_sp)) ))
    df_ratio_tc=data.frame(mrdelta=rep(df_aln_tc$mrdelta,2),
                        aligner=rep(df_aln_tc$aligner,2),
                        metric=rep(df_aln_tc$metric,2),
                        ratio=c(get(rat,df_aln_tc),rescale(1-df_aln_sp$unusedRatio,to=c(ylim_min,ylim_max),from=c(0,1))),
                        type=c( rep(rat,nrow(df_aln_tc)),rep("usedRatio",nrow(df_aln_tc)) ))
}else{
    df_ratio_sp=data.frame(mrdelta=rep(df_aln_sp$mrdelta,2),
                        aligner=rep(df_aln_sp$aligner,2),
                        metric=rep(df_aln_sp$metric,2),
                        ratio=c(get(rat,df_aln_sp),1-df_aln_sp$unusedRatio),
                        type=c( rep(rat,nrow(df_aln_sp)),rep("usedRatio",nrow(df_aln_sp)) ))
    df_ratio_tc=data.frame(mrdelta=rep(df_aln_tc$mrdelta,2),
                        aligner=rep(df_aln_tc$aligner,2),
                        metric=rep(df_aln_tc$metric,2),
                        ratio=c(get(rat,df_aln_tc),1-df_aln_tc$unusedRatio),
                        type=c( rep(rat,nrow(df_aln_tc)),rep("usedRatio",nrow(df_aln_tc)) ))
}
df_ratio_sp$metric=factor(df_ratio_sp$metric,levels=metrics)
df_ratio_tc$metric=factor(df_ratio_tc$metric,levels=metrics)
# Original metrics
df_ratio_sp_original=df_ratio_sp[which(df_ratio_sp$metric %in% metrics_original),]
df_ratio_tc_original=df_ratio_tc[which(df_ratio_tc$metric %in% metrics_original),]
# Metrics normalized per length
df_ratio_sp_perlen=df_ratio_sp[which(df_ratio_sp$metric %in% metrics_perlen),]
df_ratio_tc_perlen=df_ratio_tc[which(df_ratio_tc$metric %in% metrics_perlen),]
# Metrics normalized by length
df_ratio_sp_bylen=df_ratio_sp[which(df_ratio_sp$metric %in% metrics_bylen),]
df_ratio_tc_bylen=df_ratio_tc[which(df_ratio_tc$metric %in% metrics_bylen),]
# Metrics normalized per nseq
df_ratio_sp_perseq=df_ratio_sp[which(df_ratio_sp$metric %in% metrics_perseq),]
df_ratio_tc_perseq=df_ratio_tc[which(df_ratio_tc$metric %in% metrics_perseq),]
# Metrics normalized by nseq
df_ratio_sp_byseq=df_ratio_sp[which(df_ratio_sp$metric %in% metrics_byseq),]
df_ratio_tc_byseq=df_ratio_tc[which(df_ratio_tc$metric %in% metrics_byseq),]
# Metrics normalized per length & nseq
df_ratio_sp_perlenseq=df_ratio_sp[which(df_ratio_sp$metric %in% metrics_perlenseq),]
df_ratio_tc_perlenseq=df_ratio_tc[which(df_ratio_tc$metric %in% metrics_perlenseq),]
# Metrics normalized by length & nseq
df_ratio_sp_bylenseq=df_ratio_sp[which(df_ratio_sp$metric %in% metrics_bylenseq),]
df_ratio_tc_bylenseq=df_ratio_tc[which(df_aln_tc$metric %in% metrics_bylenseq),]

# # Plot
# # All metric
# g_aln_sp=ggplot(df_ratio_sp,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric,ncol=5)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
# g_aln_tc=ggplot(df_ratio_tc,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric,ncol=5)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+ylim(ylim_min,ylim_max)
# g_aln_sp
# ggsave(paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.sp.allmetrics.bymetric.aln.png",sep=""))
# g_aln_tc
# ggsave(paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.tc.allmetrics.bymetric.aln.png",sep=""))

g_aln_sp_original=ggplot(df_ratio_sp_original,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))

# Original metrics
g_aln_sp_original=ggplot(df_ratio_sp_original,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_tc_original=ggplot(df_ratio_tc_original,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_sp_original
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.sp.originalmetrics.bymetric.aln.png",sep=""))
g_aln_tc_original
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.tc.originalmetrics.bymetric.aln.png",sep=""))

# Metrics normalized per length
g_aln_sp_perlen=ggplot(df_ratio_sp_perlen,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_tc_perlen=ggplot(df_ratio_tc_perlen,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_sp_perlen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.sp.metricperlen.bymetric.aln.png",sep=""))
g_aln_tc_perlen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.tc.metricperlen.bymetric.aln.png",sep=""))

# Metrics multiplied by length
g_aln_sp_bylen=ggplot(df_ratio_sp_bylen,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_tc_bylen=ggplot(df_ratio_tc_bylen,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_sp_bylen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.sp.metricbylen.bymetric.aln.png",sep=""))+scale_y_continuous(rat,sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_tc_bylen
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.tc.metricbylen.bymetric.aln.png",sep=""))+scale_y_continuous(rat,sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))

# Metrics normalized per nseq
g_aln_sp_perseq=ggplot(df_ratio_sp_perseq,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_tc_perseq=ggplot(df_ratio_tc_perseq,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_sp_perseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.sp.metricperseq.bymetric.aln.png",sep=""))
g_aln_tc_perseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.tc.metricperseq.bymetric.aln.png",sep=""))

# Metrics multiplied by nseq
g_aln_sp_byseq=ggplot(df_ratio_sp_byseq,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_tc_byseq=ggplot(df_ratio_tc_byseq,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_sp_byseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.sp.metricbyseq.bymetric.aln.png",sep=""))
g_aln_tc_byseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.tc.metricbyseq.bymetric.aln.png",sep=""))

# Metrics normalized per length & nseq
g_aln_sp_perlenseq=ggplot(df_ratio_sp_perlenseq,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_tc_perlenseq=ggplot(df_ratio_tc_perlenseq,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_sp_perlenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.sp.metricperlenseq.bymetric.aln.png",sep=""))
g_aln_tc_perlenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.tc.metricperlenseq.bymetric.aln.png",sep=""))

# Metrics multiplied by length & nseq
g_aln_sp_bylenseq=ggplot(df_ratio_sp_bylenseq,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("sp - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_tc_bylenseq=ggplot(df_ratio_tc_bylenseq,aes_string(x=x_axis,y=y_axis,color=color_axis,shape=shape_axis))+facet_wrap(~metric)+geom_point()+geom_line()+labs(title=paste("tc - bucket",bucket,sep=""))+scale_y_continuous(rat,limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="UsedRatio"))
g_aln_sp_bylenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.sp.metricbylenseq.bymetric.aln.png",sep=""))
g_aln_tc_bylenseq
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/bucket",bucket,"/bucket",bucket,".mrdelta.used.tc.metricbylenseq.bymetric.aln.png",sep=""))












