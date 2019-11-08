#library(org.Hs.eg.db)
setwd("C:\\Users/Bioinfo/Desktop/")

BRCA = read.table("TCGA-BRCA-L4.csv",sep=',',header = T)

OV = read.table("TCGA-OV-L4.csv",sep=',',header = T)

plst = intersect(colnames(BRCA)[5:ncol(BRCA)],colnames(OV)[5:ncol(OV)])

kegg = read.table("KEGG_opticon.txt",sep='\t',header = T,stringsAsFactors = F)

f <- levels(as.factor(c(kegg$Regulator,kegg$Target.gene)))

node  = plst[plst %in% f]
  
BRCA <- cbind(BRCA[,1:4],BRCA[,node])

OV <- cbind(OV[,1:4],OV[,node])

write.table(BRCA,file = "BRCA_protein_for_nodes.tsv",sep = '\t',row.names = F)
write.table(OV,file = "OV_protein_for_nodes.tsv",sep = '\t',row.names = F)





###############################################################################################
###GX symbol name mapping
library(org.Hs.eg.db)
g2s=unique(toTable(org.Hs.egSYMBOL))
head(g2s)
g2e=unique(toTable(org.Hs.egENSEMBL)) 
head(g2e)
s2e=merge(g2e,g2s,by='gene_id')
table(node %in% s2e$symbol)
gs = c()
for (i in 1:length(node))
{
  if (length(s2e$ensembl_id[which(s2e$symbol==node[i])]) != 1) {
    it = s2e$ensembl_id[which(s2e$symbol==node[i])][1]
    for (j in 2:length(s2e$ensembl_id[which(s2e$symbol==node[i])])) {
      it = paste(it,s2e$ensembl_id[which(s2e$symbol==node[i])][j],sep=',')
    }
    gs<-append(gs,it)
  }
  else {
    gs<-append(gs,s2e$ensembl_id[which(s2e$symbol==node[i])])
  }
}  
map = data.frame(node,gs,stringsAsFactors = F)
add = data.frame()
delst <- c()
for (i in 1:nrow(map)) 
  {
    if (length(strsplit(map$gs[i],split = ',')[[1]]) != 1) {
      delst<-append(delst,i)
      for (j in 1:length(strsplit(map$gs[i],split = ',')[[1]])) {
        print(j)
        add=rbind(add,data.frame(map$node[i],strsplit(map$gs[i],split = ',')[[1]][j]))
      }
    }  
}
colnames(add)<- colnames(map)
map <- map[-delst,]
map = rbind(map,add)
map <- map[order(map$node),]
row.names(map)<-c(1:nrow(map))
colnames(map)<- c('Symbol_ID','Ensembl_ID')
write.table(map,"s2emap.tsv",sep='\t',row.names = F,fileEncoding = 'UTF-8')
