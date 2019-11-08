setwd('C:\\Users/Bioinfo/Desktop/')
getwd()
gx = read.table('TCGA-OV_GX_for_nodes.tsv',sep = '\t',stringsAsFactors = F,header = F)
pr = read.table('OV_protein_for_nodes.tsv',sep = '\t',stringsAsFactors = F,header = T)

sample_gx = c(t(gx[1,2:ncol(gx)]))

sample_pr = pr[,1]
for (i in 1:length(sample_pr)) {
  sample_pr[i] = substr(sample_pr[i],1,16)
}
sample_pr = unique(sample_pr)
sample_common = intersect(sample_gx,sample_pr)

gx<-gx[,c(TRUE,gx[1,2:ncol(gx)] %in% sample_common)]

pr_num = c()
for (i in 1:length(sample_pr)) {
pr_num = append(pr_num,substr(sample_pr[i],1,16) %in% sample_common)
}              
pr = pr[pr_num,]

for (i in 1:nrow(pr)) {
  pr$Sample_ID[i] = substr(pr$Sample_ID[i],1,16)
}

################some reshape###########
gx <- t(gx)
colnames(gx) <- gx[1,] 
gx <- as.data.frame(gx,stringsAsFactors = F)
gx<-gx[-1,]
colnames(gx)[1]<-colnames(pr)[1]
pr<-pr[,-c(2,3,4)]
rownames(gx)<-c(1:nrow(gx))
rownames(pr)<-c(1:nrow(pr))
gx<-gx[,c(1,order(colnames(gx)[1:ncol(gx)]))]
pr<-pr[,c(1,order(colnames(pr)[1:ncol(pr)]))]
###################done#############

##output#########
write.table(file = 'OV_TCGA_gene_expression_for_nodes.tsv', x = gx,sep = '\t',row.names = F,fileEncoding = 'UTF-8')
write.table(file = 'OV_TCGA_protein_expression_for_nodes.tsv', x = pr,sep = '\t',row.names = F,fileEncoding = 'UTF-8')
