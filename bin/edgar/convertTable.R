# https://stackoverflow.com/questions/37756389/r-convert-multiple-columns-into-single-column

#install.packages("data.table")
library(data.table)

setwd("~/Downloads")

## import data
#TOP20
tc_homo_top20 <- read.csv(file="TOP20clustalo_TC.csv", header=TRUE, sep=",")
hom_homo_top20 <- read.csv(file="TOP20clustalo_HOM.csv", header=TRUE, sep=",")
gap_homo_top20 <- read.csv(file="TOP20clustalo_GAP.csv", header=TRUE, sep=",")

#OVER1000
tc_homo_over1000 <- read.csv(file="OVERclustalo_TC.csv", header=TRUE, sep=",")
hom_homo_over1000 <- read.csv(file="OVERclustalo_HOM.csv", header=TRUE, sep=",")
gap_homo_over1000 <- read.csv(file="OVERclustalo_GAP.csv", header=TRUE, sep=",")

#ALL
tc_homo_all <- read.csv(file="clustalo_TC.csv", header=TRUE, sep=",")
hom_homo_all <- read.csv(file="clustalo_HOM.csv", header=TRUE, sep=",")
gap_homo_all <- read.csv(file="clustalo_GAP.csv", header=TRUE, sep=",")

####
#### TOP20
####

## convert to 1 column dataframe
setDT(tc_homo_top20)
fullTC_top20=melt(tc_homo_top20, id = 1, value.name = "tc", variable.name = "tree")

setDT(hom_homo_top20)
fullHOM_top20=melt(hom_homo_top20, id = 1, value.name = "hom", variable.name = "tree")

setDT(gap_homo_top20)
fullGAP_top20=melt(gap_homo_top20, id = 1, value.name = "gap", variable.name = "tree")

#combine both datasdts
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
tc_hom_top20 = merge(fullTC_top20, fullHOM_top20, by.x=c("family", "tree"), by.y=c("family", "tree"))
total_top20 = merge(tc_hom_top20, fullGAP_top20, by.x=c("family", "tree"), by.y=c("family", "tree"))

## analisys
#install.packages("ggplot2")
library(ggplot2)
#install.packages("plotly")
library(plotly)

p_tree_gap_top20 <- ggplot(data = total_top20, aes(x = gap, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("top20") +
  xlab("Avg Gap") + ylab("") +
  ylim(0,100)

p_tree_hom_top20 <- ggplot(data = total_top20, aes(x = hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("top20") +
  xlab("Homoplasty") + ylab("") +
  ylim(0,100) +
  theme(legend.position="none")

ggplotly(p_tree_gap_top20)
ggplotly(p_tree_hom_top20)

####
#### OVER1000
####

## convert to 1 column dataframe
setDT(tc_homo_over1000)
fullTC_over1000=melt(tc_homo_over1000, id = 1, value.name = "tc", variable.name = "tree")

setDT(hom_homo_over1000)
fullHOM_over1000=melt(hom_homo_over1000, id = 1, value.name = "hom", variable.name = "tree")

setDT(gap_homo_over1000)
fullGAP_over1000=melt(gap_homo_over1000, id = 1, value.name = "gap", variable.name = "tree")

#combine both datasdts
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
tc_hom_over1000 = merge(fullTC_over1000, fullHOM_over1000, by.x=c("family", "tree"), by.y=c("family", "tree"))
total_over1000 = merge(tc_hom_over1000, fullGAP_over1000, by.x=c("family", "tree"), by.y=c("family", "tree"))

## analisys
#install.packages("ggplot2")
#library(ggplot2)
#install.packages("plotly")
#library(plotly)

p_tree_gap_over1000 <- ggplot(data = total_over1000, aes(x = gap, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("Over 1000") +
  xlab("Avg Gap") + ylab("") +
  ylim(0,100) +
  theme(legend.position="none")

p_tree_hom_over1000 <- ggplot(data = total_over1000, aes(x = hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("Over 1000") +
  xlab("Homoplasy")+ ylab("") +
  ylim(0,100) +
  theme(legend.position="none")

ggplotly(p_tree_gap_over1000)
ggplotly(p_tree_hom_over1000)

####
#### ALL
####

## convert to 1 column dataframe
setDT(tc_homo_all)
fullTC_all=melt(tc_homo_all, id = 1, value.name = "tc", variable.name = "tree")

setDT(hom_homo_all)
fullHOM_all=melt(hom_homo_all, id = 1, value.name = "hom", variable.name = "tree")

setDT(gap_homo_all)
fullGAP_all=melt(gap_homo_all, id = 1, value.name = "gap", variable.name = "tree")

#combine both datasdts
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
tc_hom_all = merge(fullTC_all, fullHOM_all, by.x=c("family", "tree"), by.y=c("family", "tree"))
total_all = merge(tc_hom_all, fullGAP_all, by.x=c("family", "tree"), by.y=c("family", "tree"))

## analisys
#install.packages("ggplot2")
#library(ggplot2)
#install.packages("plotly")
#library(plotly)

p_tree_gap_all <- ggplot(data = total_all, aes(x = gap, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("full dataset") +
  xlab("Avg Gap")+ ylab("") +
  ylim(0,100) +
  theme(legend.position="none")

p_tree_hom_all <- ggplot(data = total_all, aes(x = hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("full dataset") +
  xlab("Homoplasy")+ ylab("") +
  ylim(0,100) +
  theme(legend.position="none")

ggplotly(p_tree_gap_all)
ggplotly(p_tree_hom_all)

#######
### GRID WITH ALL PLOTS
#######
#install.packages("gridExtra")
#library(gridExtra)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p_tree_gap_top20)
p_tree_gap_top20 <- p_tree_gap_top20 + theme(legend.position="none")

#grid.arrange(p_tree_gap_top20, p_tree_hom_top20,p_tree_gap_over1000, p_tree_hom_over1000,p_tree_gap_all, p_tree_hom_all,legend, nrow = 3)


#install.packages("ggpubr")
library(ggpubr)

figure <-ggarrange(p_tree_gap_top20, p_tree_hom_top20,p_tree_gap_over1000, p_tree_hom_over1000,p_tree_gap_all, p_tree_hom_all, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")


annotate_figure(figure,
                top = text_grob("Alignment with CO", size = 14),
                left = text_grob("TC Score", rot = 90)
)

####################################################
#remove RND trees
mBed_famsa =subset(total, delta=="mBed.FAMSA")
mBed_pt =subset(total,  delta=="mBed.PartTree" )
pt_famsa =subset(total, delta=="PartTree.FAMSA")
homoplasy_nonRND = rbind(mBed_famsa, mBed_pt, pt_famsa)

## count how many points in each sector

## upp - left
length(which(homoplasy_nonRND$hom<0 & homoplasy_nonRND$tc>0))
## upp - right
length(which(homoplasy_nonRND$hom>0 & homoplasy_nonRND$tc>0))
## down - left
length(which(homoplasy_nonRND$hom<0 & homoplasy_nonRND$tc<0))
## down - right
length(which(homoplasy_nonRND$hom>0 & homoplasy_nonRND$tc<0))



