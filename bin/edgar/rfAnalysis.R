#install.packages("phangorn")
library(phangorn)
 
setwd("~/CBCRG/homo_Tree_eval")

treeCO <- read.tree("seatoxin.codnd.dnd")
treeFFTNS1 <- read.tree("seatoxin.fftns1dnd.dnd")
RF.dist(treeCO, treeFFTNS1, normalize = TRUE)

finalTree <-read.tree("seatoxin.result")
result <-RF.dist(finalTree, normalize = TRUE)

##treesNames<- c("blength","codnd","cwdnd","cwqdnd","dpparttree","dpparttreednd","dpparttreednd0","dpparttreednd1","dpparttreednd2","dpparttreednd2size","fastaparttreednd","fastparttree","fftns1dnd","fftns1dndmem","fftns2dnd","fftns2dndmem","kmdnd","mafftdnd","nj","parttreednd","parttreednd0","parttreednd1","parttreednd2","parttreednd2size","swldnd")
treesNames<- c("codnd","dpparttreednd1","dpparttreednd2","dpparttreednd2size","fastaparttreednd","fftns1dnd","fftns1dndmem","fftns2dnd","fftns2dndmem","mafftdnd","parttreednd0","parttreednd1","parttreednd2","parttreednd2size")

resultMatrix<-as.matrix(result)

rownames(resultMatrix) <- treesNames
colnames(resultMatrix) <- treesNames

resultMatrix
upper.tri.remove(resultMatrix)


write.table(resultMatrix, file = "seatoxin_Analysis.csv", sep=",")
