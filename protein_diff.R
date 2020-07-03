library(readxl)
protein <- as.data.frame(read_excel("G:/CPTAC-analysis/CPTAC_protein.xlsx"))
rownames(protein) <- protein$refseq_peptide
colnames(protein)[3:ncol(protein)]=substr(colnames(protein)[3:ncol(protein)],5,17)
group <- read.table('G:/CPTAC-analysis/group.txt',sep='\t',header = T)
group$Sample <- gsub("\\.","-",substr(group$Sample,1,12))
non_BRCA <- which((colnames(protein)[3:ncol(protein)] %in% group$Sample[which(group$Group==0)]) == TRUE)
BRCA <- which((colnames(protein)[3:ncol(protein)] %in% group$Sample[which(group$Group==1)]) == TRUE)

protein <- protein[,c(BRCA,non_BRCA)]
grouplist <- c(rep("BRCA",length(BRCA)),rep("non_BRCA",length(non_BRCA)))
design <- model.matrix(~0+grouplist)
colnames(design) <- c("BRCA","non_BRCA")
rownames(design) <- colnames(protein)
#########################
library(limma)
fit <- lmFit(protein,design)
contrasts <- "BRCA-non_BRCA"
cont.matrix <- makeContrasts(contrasts=contrasts,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit3 <- eBayes(fit2)
DEG <- topTable(fit3,sort.by = 'B',number = nrow(protein))
DEG <- na.omit(DEG)
