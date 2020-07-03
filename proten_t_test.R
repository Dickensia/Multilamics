library(readxl)
protein <- as.data.frame(read_excel("G:/CPTAC-analysis/CPTAC_protein.xlsx"))
rownames(protein) <- protein$refseq_peptide
collist=lapply(strsplit(colnames(protein)[3:ncol(protein)],split = '-'), function(x) paste(x[2:4],collapse = '-'))
colvec=c()
for (i in 1:length(collist)) {
  colvec <- append(colvec,collist[[i]])
}
colnames(protein)[3:ncol(protein)]=colvec

group <- read.table('G:/CPTAC-analysis/group.txt',sep='\t',header = T)
group$Sample <- gsub("\\.","-",substr(group$Sample,1,12))
non_BRCA <- which((colnames(protein)[3:ncol(protein)] %in% group$Sample[which(group$Group==0)]) == TRUE)
BRCA <- which((colnames(protein)[3:ncol(protein)] %in% group$Sample[which(group$Group==1)]) == TRUE)
protein <- protein[3:ncol(protein)]


dup = which(lapply(strsplit(colnames(protein),split = '\\.'),function(x) length(x))!=1)
protein <- protein[,-dup]
rm <- c()
for (i in 1:nrow(protein)){
  if (length(which(is.na(protein[i,])))>30){rm <- append(rm,i)}
  
}
protein <- protein[-rm,]
FC <- apply(protein,1,function(x) mean(as.numeric(x[1:19]),na.rm = T)/mean(as.numeric(x[20:174]),na.rm = T))

#protein <- as.matrix(protein)
alter <- protein[,1:19]
unalter <- protein[,20:ncol(protein)]
v1 <- apply(alter,1,shapiro.test)
v2 <- apply(unalter,1,shapiro.test)
re1 <- data.frame(entrez=c(0),P.V=c(0))
for (i in 1:length(v1)) {
  re1 <- rbind(re1,c(names(v1[i]),v1[[i]]$p.value))
}
re2 <- data.frame(entrez=c(0),P.V=c(0))
for (i in 1:length(v2)) {
  re2 <- rbind(re2,c(names(v2[i]),v2[[i]]$p.value))
}
re1 <- re1[-1,]
re2 <- re2[-1,]
row.names(re1) <- re1$entrez
row.names(re2) <- re2$entrez

norm_row <- protein[intersect(which(re1$P.V>0.05),which(re2$P.V>0.05)),]
p=apply(norm_row,1,function(x) t.test(x[1:19],x[20:174])$p.value)

length(p[p<0.05])

wilcox_row <- protein[-intersect(which(re1$P.V>0.05),which(re2$P.V>0.05)),]
p_=apply(wilcox_row,1,function(x) wilcox.test(x[1:19],x[20:174])$p.value)
length(p_[p_<0.05])

p <- c(p[p<0.05],p_[p_<0.05])
result <- cbind(protein[names(p),],P.V=p,FC=FC[names(p)])
write.table(result,'DIFF_PROTEIN_RESULT.csv',sep=',',quote = F,fileEncoding = 'UTF-8',row.names = T)









