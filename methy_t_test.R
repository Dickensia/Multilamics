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
#remove(alter,unalter,d_num,dup)

rm <- c()
for (i in 1:nrow(data)){
  if (length(which(is.na(data[i,])))>30){rm <- append(rm,i)}
  
}
data <- data[-rm,]
FC <- apply(data,1,function(x) mean(as.numeric(x[1:28]),na.rm = T)/mean(as.numeric(x[28:179]),na.rm = T))

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









