library(GEOquery)
GSE9893=getGEO("GSE9893",destdir=".")
GSE9893_a=GSE9893[[1]]
GSE9893_dat=exprs(GSE9893_a)
GSE9893_dat=as.data.frame(GSE9893_dat)

library(org.Hs.eg.db)
library(clusterProfiler)



gpl5049=read.table("GPL5049.txt",header=T,sep='\t',stringsAsFactors = F)
rownames(gpl5049)=gpl5049$ID
overlapid=intersect(rownames(GSE9893_dat),gpl5049$ID)
gpl5049=gpl5049[overlapid,]
gpl5049_filter=gpl5049[!duplicated(gpl5049$Gene_Symbol),]
gpl5049_filter=gpl5049_filter[!is.na(gpl5049_filter$Gene_Symbol),]
rownames(gpl5049_filter)=gpl5049_filter$Gene_Symbol
gene_ids=bitr(rownames(gpl5049_filter), fromType="SYMBOL", toType=c("SYMBOL", "GENENAME","ENTREZID"), OrgDb="org.Hs.eg.db")
gpl5049_filter=gpl5049_filter[gene_ids$SYMBOL,]
data_filter=GSE9893_dat[gpl5049_filter$ID,]

library(genefu)

dannot=data.frame(row.names =rownames(gpl5049_filter), probe=gpl5049_filter$ID,"Gene.Symbol"=rownames(gpl5049_filter),"EntrezGene.ID"=gene_ids$ENTREZID)
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(data_filter),annot=dannot,do.mapping=TRUE)

subtype=as.data.frame(PAM50Preds$subtype)

write.table(file="GSE9893_subtype.txt",subtype,row.names = T,col.names = F,quote = F,sep='\t')
