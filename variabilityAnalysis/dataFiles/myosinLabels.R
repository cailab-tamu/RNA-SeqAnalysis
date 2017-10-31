dataT <- read.csv("dataFiles/dataT.csv")

# Linear model
lmM <- lm(formula = (log(1+RNAseq_CV))~scRNAseq_CV, data = dataT)

# Confidence intervals
ci <- predict(object = lmM,newdata = data.frame(scRNAseq_CV= dataT$scRNAseq_CV), interval="predict", level = 0.9)
limits <- cbind(gene= dataT$ENSEMBL,x=dataT$scRNAseq_CV,y=log(1+dataT$RNAseq_CV),ci)

# Split the groups
pdf("genesByGroup/T-Myosin.pdf")
g1 <- limits[limits[,3]>limits[,6],]
plot(x = g1[,2], y = g1[,3], col = "gray30", pch = 16, cex= 0.7, ylim=c(0,1), xlab = "CV (scRNA-seq)", ylab= "log(CV) (RNA-seq)", main = "Naive T Cells")
g3 <- limits[limits[,3]<limits[,5],]
points(x = g3[,2], y = g3[,3], col = "gray30", pch = 16, cex= 0.7)
g2 <- limits[((limits[,3]<limits[,6])&(limits[,3]>limits[,5])),]
points(x = g2[,2], y = g2[,3], col = "gray70", pch = 16, cex= 0.7)

# Adding labels of myosin cluster
myosin <- c("ENSG00000196535", "ENSG00000078814", "ENSG00000133392")
myosin <- dataT[dataT$ENSEMBL %in% myosin,c("RNAseq_CV","scRNAseq_CV")]
points(myosin[,2],myosin[,1],col="red", pch=16,cex=0.7)
text (myosin[,2],myosin[,1]+.015, c("MYO18A","MYH7B","MYH11"), cex = 0.7, col="red")

dev.off()

