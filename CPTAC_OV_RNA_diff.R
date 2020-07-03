setwd("G:/CPTAC-analysis/")
alter <- read.table("OV_altered_sample_rna_count.tsv",sep='\t',header = T,stringsAsFactors = T)
unalter <- read.table("OV_unaltered_sample_rna_count.tsv",sep='\t',header = T,stringsAsFactors = T)
ncol(alter)#20-1=19
ncol(unalter)#108-1-4=103
dup = which((colnames(unalter) %in% colnames(alter))[-1] == TRUE)
unalter <- unalter[,-dup]
#altered_samples_list from Cbioportal has 4 dups with unaltered_samples_list due to some unknown errors
#so I removed them from unaltered samples
#it will be fixed in next try.
which(alter[,1] != unalter[,1])#0 
rnaseq=cbind(alter,unalter[,-1])#|rowname[,1]|alter[,2:20]|unalter[21:123]|
rnaseq=rnaseq[-((nrow(rnaseq)-4):nrow(rnaseq)),]#drop meaningless rows
remove(alter,unalter)
row.names(rnaseq) <- rnaseq[,1]
rnaseq <- rnaseq[,-1]
colnames(rnaseq) <- gsub('\\.','-',colnames(rnaseq))

#######################################

library(limma)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
grouplist <- c(rep("BRCA",19),rep("non_BRCA",103))
design <- model.matrix(~0+grouplist)
colnames(design) <- c("BRCA","non_BRCA")
rownames(design) <- colnames(rnaseq)

dge <- DGEList(counts=rnaseq)
density_df <- data.frame(dge$counts[,1:40])#20alter:20unalter for visualization
cols = colorRampPalette(brewer.pal(9, "Set1"))(40)
plot(density(log(density_df[,1]+1)),col='red')
for (i in 2:ncol(density_df)) {
  lines(density(log(density_df[,i]+1)),col='red')
}

#filter
keep <- filterByExpr(dge, design)
dge2 <- dge[keep,,keep.lib.sizes=FALSE]

density_df2 <- data.frame(dge2$counts[,1:40])
lines(density(log(density_df2[,1]+1)),col='blue')
for (i in 2:ncol(density_df2)) {
  lines(density(log(density_df2[,i]+1)),col='blue')
}


dge <- calcNormFactors(dge2)
logCPM <- cpm(dge,log=T,prior.count = 1)#edgeR advise to use this matrix for plotting
density_df3 <- data.frame(logCPM[,1:40])
lines(density(density_df3[,1]),col='purple')
for (i in 2:ncol(density_df3)) {
  lines(density(density_df3[,i]),col='purple')
}
#plotSA()
v <- voom(dge,design,normalize.method = "quantile")

fit <- lmFit(v,design)#linear modelling
plotSA(fit)
contrasts <- "BRCA-non_BRCA"
cont.matrix <- makeContrasts(contrasts=contrasts,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix)#calc the diff
plotSA(fit2)
fit2 <- eBayes(fit2)#use Bayes test
plotSA(fit2)
DEG <- topTable(fit2,coef=contrasts,n=Inf)
DEG <- na.omit(DEG)

row_names <- rownames(DEG)[which(DEG$P.Value<0.05)]
Result <- cbind(rnaseq[row_names,],P.V=DEG[row_names,]$P.Value,logFC=DEG[row_names,]$logFC)
write.table(Result,'DIFF_RNA_RESULT.csv',sep=',',quote = F,row.names = T,fileEncoding = 'UTF-8')
