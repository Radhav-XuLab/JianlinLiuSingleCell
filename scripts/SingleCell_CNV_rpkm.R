library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(circlize)
#Create rpkm values first exome wise then add them genewise to get an expression matrix data for all samples separately 
#Provide these files as input files 
sample_list <- read.table("sample.list")
sample_list <- t(sample_list)
filename <- sample_list[,1]
sample1<- read.table(paste(filename,".txt",sep = ""),sep = '\t',header = TRUE)
a <- data.frame(sample1)

for (i in 2:length(sample_list)){
  filename <- sample_list[,i]
  sample2<- read.table(paste(filename,".txt",sep = ""),sep = '\t',header = TRUE)
  b <- data.frame(sample2)
  a <- merge(a,b)
  rm(b)
  }
a$chr <- paste("chr",sep = "",a$chr)
head(a)
exprSet=a[,4:ncol(a)]
row.names(exprSet) <- exprSet$gene
exprSet <- exprSet[,-1]

#Filter the expression matrix
dim(exprSet)
exprSet=exprSet[apply(exprSet,1,function(x) all(x>0)),]
dim(exprSet)
boxplot(exprSet,las=2)
normExprSet=apply(exprSet, 2, function(x) log(1e6*x/sum(x)+1))
boxplot(normExprSet,las=2)
normExprSet=as.data.frame(normExprSet)

#Combine the gene coordinate information to the expression matrix
normExprSet$gene=rownames(normExprSet)
anno=a[,1:4]
res=merge(anno,normExprSet,by='gene')

#Coordinate the genes
bed=res[,1:4]
head(bed)
normExprSet=res[,5:ncol(res)]
table(bed$chr)
bed$chr= factor(bed$chr,levels =paste0('chr',c(1:19,'X')) )
sort_bed=bed[order(bed$chr,bed$start),]
sort_normExprSet=normExprSet[order(bed$chr,bed$start),]
table(sort_bed$chr)
head(sort_bed)
dim(normExprSet)
res=cbind(sort_bed,sort_normExprSet)
head(res)
dim(res)

#Calculate CNV
all_cnv <- lapply(split(res,res$chr), function(x){
  anno=x[,1:4]
  dat=x[,5:103]
  cnv <- lapply(51:(nrow(x)-50), function(i){
    this_cnv <- unlist( lapply(5:103, function(j){
      sum(x[(i-50):(i+50),j])/101
    }))
    return(this_cnv)
  })
  cnv=do.call(rbind,cnv)
  cnv=cbind(x[51:(nrow(x)-50),1:4],cnv)
})
all_cnv=do.call(rbind,all_cnv) 
head(all_cnv)
dim(all_cnv)
D=t(scale(all_cnv[,5:103] ))
D[abs(D) < 1] <- 0
dim(D)

F <- as.data.frame(t(D))
colnames(F) <- sample_list
F <-F %>% rownames_to_column('ID') 
DF <- F %>%  mutate_if(is.numeric, round, digits = 1)
DF <-DF %>% column_to_rownames('ID')
G <- DF[abs(DF$WTB) <= 1,]
K <- t(G)
K <- K[-1,]
row.names(G) -> A
Z <- unlist(strsplit(A,"\\."))
chr_loca <- rle(Z[seq(1, length(Z), 2)])

#Plot the heatmap
cols = colorRamp2(c(-4, -3, -2, -1, 0, 1, 2, 3, 4), c("#2166ac","#4393c3" ,"#92c5de" ,"#d1e5f0" 
                                                      ,"white","#fddbc7","#f4a582","#d6604d","#b2182b"))

h1 = Heatmap(K,name = "heatmap", col = cols,cluster_rows = FALSE,cluster_columns = FALSE,row_names_gp = gpar(fontsize = 3),column_names_gp = gpar(fontsize =0.01 ),heatmap_legend_param = list(title = "Scale"))
png(filename="SC_CNV.png", width = 20000, height = 10000, res = 2000)
h1
loca = 0
for (l in 1:(length(chr_loca$lengths)-1)){
  loca = loca + chr_loca$lengths[l]
  decorate_heatmap_body("heatmap", {
    x = loca/ncol(K)
    grid.lines(c(x, x), c(0, 1), gp = gpar(fill="black", col = "black",lty = 2))
  })
}
decorate_heatmap_body("heatmap", {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  
})
#Adding gridlines to indicate groups
decorate_heatmap_body("heatmap", {
  y = 67/nrow(K)
  grid.lines(c(0,1), c(y, y), gp = gpar(fill="black", col = "black",lwd = 1))
})
decorate_heatmap_body("heatmap", {
  y = 58/nrow(K)
  grid.lines(c(0,1), c(y, y), gp = gpar(fill="black", col = "black",lwd = 1))
})
decorate_heatmap_body("heatmap", {
  y = 20/nrow(K)
  grid.lines(c(0,1), c(y, y), gp = gpar(fill="black", col = "black",lwd = 1))
})
decorate_heatmap_body("heatmap", {
  y = 91/nrow(K)
  grid.lines(c(0,1), c(y, y), gp = gpar(fill="black", col = "black",lwd = 1))
})
decorate_heatmap_body("heatmap", {
  y = 74/nrow(K)
  grid.lines(c(0,1), c(y, y), gp = gpar(fill="black", col = "black",lwd = 2))
})
decorate_heatmap_body("heatmap", {
  y = 38/nrow(K)
  grid.lines(c(0,1), c(y, y), gp = gpar(fill="black", col = "black",lwd = 2))
})

dev.off()

