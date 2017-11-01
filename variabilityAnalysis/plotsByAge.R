geneAge <- readLines("dataFiles/HUMAN_PPODv4_OrthoMCL_dollo_age-depth.protein_list")
geneAge <- geneAge[!grepl("^#",geneAge)]
geneAge <- lapply(geneAge,function(gene){unlist(strsplit(gene,"\t"))})
geneAge <- lapply(geneAge,function(gene){gene[c(1,length(gene))]})
geneAge <- do.call(rbind.data.frame,geneAge)
colnames(geneAge) <- c("ENSEMBL","AGE")
dataT <- read.csv("dataFiles/dataT.csv")
par(mfrow=c(2,4))
nameLevel <- c("Human", "Euarchontoglire", "Amniota", "Euteleostomi", "Bilateria", "Opisthokonta", "Eukaryota","Cellular Organisms")
pdf("plotsByAge/plotByAge.pdf")
for(i in (seq_len(8)-1)){
  age <- dataT[dataT[,"ENSEMBL"]%in%as.vector(geneAge[geneAge[,"AGE"]%in%c(0),1]),c("RNAseq_CV","scRNAseq_CV")]
  noAge <- dataT[dataT[,"ENSEMBL"]%in%as.vector(geneAge[geneAge[,"AGE"]%in%c(i),1]),c("RNAseq_CV","scRNAseq_CV")]
  plot(noAge$scRNAseq_CV,log(1+noAge$RNAseq_CV), pch = 16, col="gray75", xlab = "(CV) scRNA-seq", ylab = "log(CV) RNA-seq", ylim=c(0,1), main = nameLevel[i+1])
  points(age[,2],log(1+age[,1]), pch=16)
}
dev.off()
cor(dataT$RNAseq_CV, dataT$scRNAseq_CV)
cor(age$RNAseq_CV,age$scRNAseq_CV)
cor.test(noAge$RNAseq_CV,noAge$scRNAseq_CV)

