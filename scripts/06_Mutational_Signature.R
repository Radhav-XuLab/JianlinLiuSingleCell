library(deconstructSigs)
#Input file contains columns for chromosome number, variant position/coordinate, reference and alternate bases.
df_table<-read.table("all.txt", sep='\t')

#Add sample name Bulk
df<-data.frame(sample=c(rep('Bulk', times=(length(df_table$V3)))),
               chr=df_table$V1, pos=df_table$V2, ref=df_table$V3,
               alt=df_table$V4)

#Create input to whichSignatures
sigs.input <- mut.to.sigs.input(mut.ref = df, 
                                sample.id = "sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

#Create mutational signature chart
wchSig = whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = 'Bulk', 
                         contexts.needed = TRUE,
                         tri.counts.method = 'default')

pdf('Bulk_all.sig.pdf', width = 10, height = 10)
chart<-plotSignatures(wchSig)
dev.off()
