# https://stackoverflow.com/questions/37756389/r-convert-multiple-columns-into-single-column

#install.packages("data.table")
library(data.table)

setwd("~/Downloads")

## import data
#TOP20
tc_homo_top20 <- read.csv(file="homoplasty\ -\ top20DELTA_clustalo_TC.csv", header=TRUE, sep=",")
hom_homo_top20 <- read.csv(file="homoplasty\ -\ top20DELTA_clustalo_HOM.csv", header=TRUE, sep=",")
w_homo_top20 <- read.csv(file="homoplasty\ -\ top20DELTA_clustalo_W_HOM.csv", header=TRUE, sep=",")
gap_homo_top20 <- read.csv(file="homoplasty\ -\ top20DELTA_clustalo_GAP.csv", header=TRUE, sep=",")

#OVER1000
tc_homo_over1000 <- read.csv(file="homoplasty\ -\ overDELTA_clustalo_TC.csv", header=TRUE, sep=",")
hom_homo_over1000 <- read.csv(file="homoplasty\ -\ overDELTA_clustalo_HOM.csv", header=TRUE, sep=",")
w_homo_over1000 <- read.csv(file="homoplasty\ -\ overDELTA_clustalo_W_HOM.csv", header=TRUE, sep=",")
gap_homo_over1000 <- read.csv(file="homoplasty\ -\ overDELTA_clustalo_GAP.csv", header=TRUE, sep=",")

#ALL
tc_homo_all <- read.csv(file="homoplasty\ -\ DELTA_clustalo_TC.csv", header=TRUE, sep=",")
hom_homo_all <- read.csv(file="homoplasty\ -\ DELTA_clustalo_HOM.csv", header=TRUE, sep=",")
w_homo_all <- read.csv(file="homoplasty\ -\ DELTA_clustalo_W_HOM.csv", header=TRUE, sep=",")
gap_homo_all <- read.csv(file="homoplasty\ -\ DELTA_clustalo_GAP.csv", header=TRUE, sep=",")

####
#### TOP20
####

## convert to 1 column dataframe
setDT(tc_homo_top20)
fullTC_top20=melt(tc_homo_top20, id = 1, value.name = "tc", variable.name = "tree")

setDT(hom_homo_top20) 
fullHOM_top20=melt(hom_homo_top20, id = 1, value.name = "hom", variable.name = "tree")

setDT(w_homo_top20) 
fullwHOM_top20=melt(w_homo_top20, id = 1, value.name = "w_hom", variable.name = "tree")

setDT(gap_homo_top20)
fullGAP_top20=melt(gap_homo_top20, id = 1, value.name = "gap", variable.name = "tree")

#combine both datasdts
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
tc_hom_top20 = merge(fullTC_top20, fullHOM_top20, by.x=c("family", "tree"), by.y=c("family", "tree"))
total_top20 = merge(tc_hom_top20, fullGAP_top20, by.x=c("family", "tree"), by.y=c("family", "tree"))
full_top20 = merge(total_top20, fullwHOM_top20, by.x=c("family", "tree"), by.y=c("family", "tree"))

## analisys
#install.packages("ggplot2")
library(ggplot2)
#install.packages("plotly")
library(plotly)

p_tree_gap_top20 <- ggplot(data = full_top20, aes(x = gap, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("top20") +
  xlab("Avg Gap") + ylab("") +
  ylim(-100,100)
  
p_tree_hom_top20 <- ggplot(data = full_top20, aes(x = hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("top20") +
  xlab("Homoplasty") + ylab("") +
  ylim(-100,100) +theme(legend.position="none")

p_tree_w_hom_top20 <- ggplot(data = full_top20, aes(x = w_hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("top20") +
  xlab("W_Homoplasty") + ylab("") +
  ylim(-100,100) +theme(legend.position="none")

#ggplotly(p_tree_gap_top20)
#ggplotly(p_tree_hom_top20)
#ggplotly(p_tree_w_hom_top20)

####
#### OVER1000
####

## convert to 1 column dataframe
setDT(tc_homo_over1000)
fullTC_over1000=melt(tc_homo_over1000, id = 1, value.name = "tc", variable.name = "tree")

setDT(hom_homo_over1000)
fullHOM_over1000=melt(hom_homo_over1000, id = 1, value.name = "hom", variable.name = "tree")

setDT(w_homo_over1000) 
fullwHOM_over1000=melt(w_homo_over1000, id = 1, value.name = "w_hom", variable.name = "tree")

setDT(gap_homo_over1000)
fullGAP_over1000=melt(gap_homo_over1000, id = 1, value.name = "gap", variable.name = "tree")

#combine both datasdts
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
tc_hom_over1000 = merge(fullTC_over1000, fullHOM_over1000, by.x=c("family", "tree"), by.y=c("family", "tree"))
total_over1000 = merge(tc_hom_over1000, fullGAP_over1000, by.x=c("family", "tree"), by.y=c("family", "tree"))
full_over1000 = merge(total_over1000, fullwHOM_over1000, by.x=c("family", "tree"), by.y=c("family", "tree"))

## analisys
#install.packages("ggplot2")
#library(ggplot2)
#install.packages("plotly")
#library(plotly)

p_tree_gap_over1000 <- ggplot(data = full_over1000, aes(x = gap, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("Over 1000") +
  xlab("Avg Gap") + ylab("") +
  ylim(-100,100) +theme(legend.position="none")


p_tree_hom_over1000 <- ggplot(data = full_over1000, aes(x = hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("Over 1000") +
  xlab("Homoplasy")+ ylab("") +
  ylim(-100,100) +theme(legend.position="none")

p_tree_w_hom_over1000 <- ggplot(data = full_over1000, aes(x = w_hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("Over 1000") +
  xlab("W_Homoplasty") + ylab("") +
  ylim(-100,100) +theme(legend.position="none")

#ggplotly(p_tree_gap_over1000)
#ggplotly(p_tree_hom_over1000)
#ggplotly(p_tree_w_hom_over1000)

####
#### ALL
####

## convert to 1 column dataframe
setDT(tc_homo_all)
fullTC_all=melt(tc_homo_all, id = 1, value.name = "tc", variable.name = "tree")

setDT(hom_homo_all)
fullHOM_all=melt(hom_homo_all, id = 1, value.name = "hom", variable.name = "tree")

setDT(w_homo_all) 
fullwHOM_all=melt(w_homo_all, id = 1, value.name = "w_hom", variable.name = "tree")

setDT(gap_homo_all)
fullGAP_all=melt(gap_homo_all, id = 1, value.name = "gap", variable.name = "tree")

#combine both datasdts
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
tc_hom_all = merge(fullTC_all, fullHOM_all, by.x=c("family", "tree"), by.y=c("family", "tree"))
total_all = merge(tc_hom_all, fullGAP_all, by.x=c("family", "tree"), by.y=c("family", "tree"))
full_all = merge(total_all, fullwHOM_all, by.x=c("family", "tree"), by.y=c("family", "tree"))

## analisys
#install.packages("ggplot2")
#library(ggplot2)
#install.packages("plotly")
#library(plotly)

p_tree_gap_all <- ggplot(data = full_all, aes(x = gap, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("full dataset") +
  xlab("Avg Gap")+ ylab("") +
  ylim(-100,100) #+theme(legend.position="none")

p_tree_hom_all <- ggplot(data = full_all, aes(x = hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("full dataset") +
  xlab("Homoplasy")+ ylab("") +
  ylim(-100,100) #+theme(legend.position="none")

p_tree_w_hom_all <- ggplot(data = full_all, aes(x = w_hom, y = tc, color = tree)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  ggtitle("full dataset") +
  xlab("W_Homoplasty") + ylab("") +
  ylim(-100,100) #+theme(legend.position="none")

#ggplotly(p_tree_gap_all)
#ggplotly(p_tree_hom_all)
#ggplotly(p_tree_w_hom_all)

#######
### GRID WITH ALL PLOTS
#######

#install.packages("ggpubr")
library(ggpubr)

figure_complete <-ggarrange(p_tree_gap_top20, p_tree_hom_top20,p_tree_gap_over1000, p_tree_hom_over1000,p_tree_gap_all, p_tree_hom_all, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")

annotate_figure(figure_complete,
                top = text_grob("DELTAS - CO align", size = 14),
                left = text_grob("TC Score", rot = 90)
)

figure_top20 <-ggarrange(p_tree_gap_top20, p_tree_hom_top20,p_tree_w_hom_top20, ncol=1, nrow=3, common.legend = TRUE, legend="right")

annotate_figure(figure_top20,
                top = text_grob("DELTAS - CO align", size = 14),
                left = text_grob("TC Score", rot = 90)
)

figure_over <-ggarrange(p_tree_gap_over1000, p_tree_hom_over1000,p_tree_w_hom_over1000, ncol=1, nrow=3, common.legend = TRUE, legend="right")

annotate_figure(figure_over,
                top = text_grob("DELTAS - CO align", size = 14),
                left = text_grob("TC Score", rot = 90)
)

figure_all <-ggarrange(p_tree_gap_all, p_tree_hom_all,p_tree_w_hom_all, ncol=1, nrow=3, common.legend = TRUE, legend="right")

annotate_figure(figure_all,
                top = text_grob("DELTAS - CO align", size = 14),
                left = text_grob("TC Score", rot = 90)
)
####################################################
#remove RND trees

## > ALL ++> 202 over 282
mBed_famsa =subset(full_all, tree=="mBed.FAMSA")
mBed_pt =subset(full_all,  tree=="mBed.PartTree" )
pt_famsa =subset(full_all, tree=="PartTree.FAMSA")

homoplasy_nonRND_full = rbind(mBed_famsa, mBed_pt, pt_famsa)

## count how many points in each sector
## upp - left
length(which(homoplasy_nonRND_full$hom<0 & homoplasy_nonRND_full$tc>0))
## upp - right
length(which(homoplasy_nonRND_full$hom>0 & homoplasy_nonRND_full$tc>0))
## down - left
length(which(homoplasy_nonRND_full$hom<0 & homoplasy_nonRND_full$tc<0))
## down - right
length(which(homoplasy_nonRND_full$hom>0 & homoplasy_nonRND_full$tc<0))

## > TOP20 ++> 116 over 120
mBed_famsa =subset(full_top20, tree=="mBed.FAMSA")
mBed_pt =subset(full_top20,  tree=="mBed.PartTree" )
pt_famsa =subset(full_top20, tree=="PartTree.FAMSA")

homoplasy_nonRND_top20 = rbind(mBed_famsa, mBed_pt, pt_famsa)

## count how many points in each sector
## upp - left
length(which(homoplasy_nonRND_top20$hom<0 & homoplasy_nonRND_top20$tc>0))
## upp - right
length(which(homoplasy_nonRND_top20$hom>0 & homoplasy_nonRND_top20$tc>0))
## down - left
length(which(homoplasy_nonRND_top20$hom<0 & homoplasy_nonRND_top20$tc<0))
## down - right
length(which(homoplasy_nonRND_top20$hom>0 & homoplasy_nonRND_top20$tc<0))


## > OVER1000 ++> 202 over 225
mBed_famsa =subset(full_over1000, tree=="mBed.FAMSA")
mBed_pt =subset(full_over1000,  tree=="mBed.PartTree" )
pt_famsa =subset(full_over1000, tree=="PartTree.FAMSA")

homoplasy_nonRND_over1000 = rbind(mBed_famsa, mBed_pt, pt_famsa)

## count how many points in each sector
## upp - left
length(which(homoplasy_nonRND_over1000$hom<0 & homoplasy_nonRND_over1000$tc>0))
## upp - right
length(which(homoplasy_nonRND_over1000$hom>0 & homoplasy_nonRND_over1000$tc>0))
## down - left
length(which(homoplasy_nonRND_over1000$hom<0 & homoplasy_nonRND_over1000$tc<0))
## down - right
length(which(homoplasy_nonRND_over1000$hom>0 & homoplasy_nonRND_over1000$tc<0))



