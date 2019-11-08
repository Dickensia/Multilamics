#读入临床数据并进行数表处理
phenoData <- read.table( "C:\\Users/Bioinfo/Desktop/nationwidechildrens.org_clinical_patient_ov.txt",
                         header = T,
                         sep = '\t',
                         quote = '',
                         stringsAsFactors = F)
phenoData<-phenoData[,-1]
#colnames(phenoData)<-phenoData[1,]
phenoData<-phenoData[-(1:2),]
rownames( phenoData ) <- phenoData[ , 1]
colnames( phenoData )[1] <- "Tumor_Sample_Barcode"
phenoData[1:5, 1:5]


maf <- data.table::as.data.table(read.csv(file = "C:\\Users/Bioinfo/Desktop/TCGA.OV.mutect.somatic.maf", 
                                          header = TRUE, sep = '\t', 
                                          stringsAsFactors = FALSE, comment.char = "#"))
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, 1, 12)

require(maftools) 

#phenoData[phenoData$tumor_grade == ""] <- 'no_reported'
#phenoData[phenoData$tumor_grade == "[Not Available]"] <- 'no_reported'

laml <- read.maf(maf, clinicalData = phenoData)
laml 

laml@data <- laml@data[grepl('PASS', laml@data$FILTER), ]

genestable<-read.table("C:\\Users/Bioinfo/Desktop/OV_out.sig_genes.txt",sep="\t",header=T,fill = T,stringsAsFactors = F)
geneslist<-genestable[1:30,1]

library(RColorBrewer)

#突变种类上色
variantClass <- names(table(laml@data$Variant_Classification))
col = c(RColorBrewer::brewer.pal(n = 4, name = 'Set1'),
        RColorBrewer::brewer.pal(n = 5, name = 'Set2'))
names(col) = variantClass
col
#肿瘤级别上色
laml@clinical.data$tumor_grade <- 
  as.factor(laml@clinical.data$tumor_grade)
laml@clinical.data$tumor_grade
gradecolors = RColorBrewer::brewer.pal(n = 7,name = 'Spectral')
names(gradecolors) = levels(laml@clinical.data$tumor_grade)
gradecolors

# #性别上色
# laml@clinical.data$gender <- 
#   as.factor(laml@clinical.data$gender)
# Gendercolors = c("#b3e2cd", "#fb9a99")
# names(Gendercolors) = levels(laml@clinical.data$gender)
# Gendercolors

#绘图需要的list,list里面的变量为laml@clinical里面的列名
phecolors = list(tumor_grade = gradecolors)
                 #gender = Gendercolors)
#画图
png(paste0('oncoplot_top30_phone', "_OV", '','.png'), res = 600,
    width = 6000, height = 4320)

oncoplot(maf = laml, 
         colors = col, 
         bgCol = "#ebebeb", borderCol = "#ebebeb",
         genes = geneslist, GeneOrderSort = F, keepGeneOrder = T,
         fontSize = 0.5 , legendFontSize = 1,
         annotationFontSize = 1,
         titleFontSize = 1,
         sortByMutation = T,
         showTumorSampleBarcodes = F,
         annotationColor = phecolors,
         clinicalFeatures = c("tumor_grade"))
dev.off()

###

