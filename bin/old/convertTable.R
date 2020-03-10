# https://stackoverflow.com/questions/37756389/r-convert-multiple-columns-into-single-column

  
## Libraries & packages
#install.packages("data.table")
library(data.table)

setwd("~/homoplasy/data")

metrics=c("tc","hom","w_hom","w_hom2","ngap","ngap2")
files=c(paste(aligner,"_TC.csv",sep=""),
        paste(aligner,"_HOM.csv",sep=""),
        paste(aligner,"_W_HOM.csv",sep=""),
        paste(aligner,"_W_HOM2.csv",sep=""),
        paste(aligner,"_NGAP.csv",sep=""),
        paste(aligner,"_NGAP2.csv",sep=""))

aligner="clustalo"


## Import data
import_data=function(aligner,metrics,files){
    "
    Input  -  aligner : name of the aligner (clustalo, )
              metrics : vector with the name of the metrics (tc, hom, w_hom, w_hom2, ngap, ngap2)
              files   : vector with the name of the .csv files to be imported
    Output -  dataList: list of imported data 
    "
    dataList=list()
    n=length(metrics)
    for (i in 1:n){
        metric=metrics[i]
        fil=files[i]
        dataList[[i]]=read.csv(file=fil,header=TRUE,sep=",")
    }
    names(dataList)=metrics
    return(dataList)
}



## Convert to 1 column dataframe
to1coldataList=function(dataList,metrics){

    df1colList=list()
    n=length(metrics)
    for (i in 1:n){
        data=dataList[[i]]
        nam=metrics[i]
        print(nam)
        setDT(data)
        df1colList[[i]]=melt(data, id = 1, value.name = nam, variable.name = "tree")
    }
    names(df1colList)=metrics
    return(df1colList)
}


#combine both datasets
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
tc = merge(fullTC, fullHOM, by.x=c("family", "tree"), by.y=c("family", "tree"))
total = merge(tc, fullGAP, by.x=c("family", "tree"), by.y=c("family", "tree"))
