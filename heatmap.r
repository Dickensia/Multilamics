library(pheatmap)

mat <- as.matrix(rnaseq)
n_normal=103
n_cancer=19
sample_type=c(rep("tumor",n_cancer),rep("normal",n_normal))   

annotation_col=data.frame(class=sample_type)
row.names(annotation_col)=colnames(rnaseq)
png("heatmap.png",height=8,width=8,units="in",res=300)
pheatmap(log(mat+1), color = colorRampPalette(c("green", 
"white","red"))(dim(mat)[1]),
border_color = "NA",show_rownames=F,show_colnames=F,
annotation_col=annotation_col)
dev.off()

                                                 