###
### Script to plot correlation btwn deltaTC vs deltaHomoplasy
###

## set directory
setwd("~/CBCRG/tc_homo")

###
## Install packages
###
install.packages("ggplot2")
library(ggplot2)
install.packages("plotly")
library(plotly)

###
##  Function to perform CORRELATION by GROUP
###
require(plyr)
func <- function(xx)
{
  return(data.frame(COR = cor(xx$TC, xx$HOM)))
}
###

homoplasy <- read.csv(file="co_top20_reformat.csv", header=TRUE, sep=",")

## CALL correlation by group
ddply(homoplasy, .(tree.compa), func)


p_family <- ggplot(data = homoplasy, aes(x = HOM, y = TC, color = family)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)

p_family + theme(axis.ticks.y = element_blank())

p_trees <- ggplot(data = homoplasy, aes(x = HOM, y = TC, color = tree.compa)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

p_trees = p_trees + scale_y_continuous(breaks=seq(-50,50,10))

ggplotly(p_family)
ggplotly(p_trees)


## upp - left
length(which(homoplasy$HOM<0 & homoplasy$TC>0))
## upp - right
length(which(homoplasy$HOM>0 & homoplasy$TC>0))

## down - left
length(which(homoplasy$HOM<0 & homoplasy$TC<0))
## down - right
length(which(homoplasy$HOM>0 & homoplasy$TC<0))



mBed_famsa =subset(homoplasy, tree.compa=="mBed-FAMSA")
mBed_pt =subset(homoplasy,  tree.compa=="mBed-PartTree" )
pt_famsa =subset(homoplasy, tree.compa=="PartTree-FAMSA")


homoplasy_nonRND = rbind(mBed_famsa, mBed_pt, pt_famsa)
## upp - left
length(which(homoplasy_nonRND$HOM<0 & homoplasy_nonRND$TC>0))
## upp - right
length(which(homoplasy_nonRND$HOM>0 & homoplasy_nonRND$TC>0))

## down - left
length(which(homoplasy_nonRND$HOM<0 & homoplasy_nonRND$TC<0))
## down - right
length(which(homoplasy_nonRND$HOM>0 & homoplasy_nonRND$TC<0))


p_trees <- ggplot(data = homoplasy_nonRND, aes(x = HOM, y = TC, color = tree.compa)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

ggplotly(p_trees)
