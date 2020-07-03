setwd("G:\\CPTAC")
library(readxl)
CNA <- read_excel("CPTAC_CNA.xlsx")
CNA <- CNA[-1,]
colnames(CNA)[6:ncol(CNA)] <- substr(colnames(CNA)[6:ncol(CNA)],1,16)
genelist <- CNA$`Hybridization REF`
maf <- read.table('alter_MAF.csv',sep=',',header = T,stringsAsFactors = F)
colnames(maf)[2:ncol(maf)] <- gsub('\\.','-',colnames(maf)[2:ncol(maf)])
genelist2 <- maf$GENESYMBOL
gene <- intersect(genelist,genelist2)
cohort <- intersect(colnames(CNA)[6:ncol(CNA)],colnames(maf)[2:ncol(maf)])


