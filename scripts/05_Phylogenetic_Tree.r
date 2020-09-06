library(readxl)
library(phangorn)
library("ggtree")
library("colorspace")
library(ggplot2)
#153.xls contains columnwise samples with snp information given based on presence and absence indicated by 1 and 0 respectively
file <- read_excel("153.xlsx")
data1 <- as.phyDat(t(file), type="USER", levels = c(0, 1))
stree = NJ(dist.gene(t(file)))
treeRatchet1 <- pratchet(data1, start=stree)
treeRatchet2 <- acctran(treeRatchet1, data1)

#Saving the tree 
saveRDS(treeRatchet2, file = "153.Rds")

treeRatchet2 <- readRDS("153.Rds")
#Create list for groups
cls <- list(PT=c("153PT1","153PT2","153PT3","153PT4","153PT5","153PT6","153PT7","153PT8","153PT9","153PT10","153PT11","153PT12","153PT13","153PT14","153PT15","153PT16","153PT17","153PT18"),
            LMT=c("153LMT1","153LMT2","153LMT3","153LMT4","153LMT5","153LMT6","153LMT7","153LMT8","153LMT9","153LMT10","153LMT11","153LMT12","153LMT13","153LMT14","153LMT15","153LMT16","153LMT17","153LMT18","153LMT19","153LMT20"))
treeRatchet2 <- groupOTU(treeRatchet2, cls)

#Plot the tree 
pdf("153.pdf") 
ggtree(treeRatchet2, layout = "unrooted",aes(color=group)) + geom_tiplab(aes(label=label)) + scale_color_manual(values=c("deeppink3","dodgerblue"), breaks=1:2,labels=c("PT", "LMT")) + theme(legend.position="right")
dev.off()
