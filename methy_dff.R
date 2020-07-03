setwd("G://CPTAC-analysis/")
alter <- read.table("CPTAC_alter_methy.tsv",sep='\t',header = T,row.names = 1,stringsAsFactors = F)
unalter <- read.table("CPTAC_unalter_methy.tsv",sep='\t',header = T,row.names = 1,stringsAsFactors = F)
dup = intersect(colnames(alter),colnames(unalter))
d_num <- c()
for (i in 1:length(dup)) {
  d_num <- append(d_num,which(colnames(unalter)==dup[i]))
}
unalter <- unalter[,-d_num]
data <- cbind(alter,unalter)
n1 <- ncol(alter)
n2 <- ncol(unalter)
remove(alter,unalter,d_num,dup)
#write.table(data.frame(Sample=colnames(data),Group=c(rep(1,n1),rep(0,n2))),"group.txt",sep='\t',quote = F,fileEncoding = 'UTF-8',row.names = F)
#write.table(data,"data.txt",sep='\t',quote = F,fileEncoding = 'UTF-8')#colnames[1]absent


grouplist <- c(rep("BRCA",n1),rep("non_BRCA",n2))
design <- model.matrix(~0+grouplist)
colnames(design) <- c("BRCA","non_BRCA")
rownames(design) <- colnames(data)

#####################################################
library(limma)
fit <- lmFit(data,design)
contrasts <- "BRCA-non_BRCA"
cont.matrix <- makeContrasts(contrasts=contrasts,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit3 <- eBayes(fit2)
DEG <- topTable(fit3,sort.by = 'B',number = nrow(data))
DEG <- na.omit(DEG)

