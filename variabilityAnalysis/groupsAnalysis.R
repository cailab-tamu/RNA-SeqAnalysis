dataM <- read.csv("dataFiles/dataM.csv")
dataT <- read.csv("dataFiles/dataT.csv")

# Linear model
lmM <- lm(formula = (log(1+RNAseq_CV))~scRNAseq_CV, data = dataM)
plot((log(1+dataM$RNAseq_CV))~dataM$scRNAseq_CV)
abline(lmM, col= "red")

# Confidence intervals
ci <- predict(object = lmM,newdata = data.frame(scRNAseq_CV= dataM$scRNAseq_CV), interval="predict", level = 0.95)
limits <- cbind(gene= dataM$ENSEMBL,x=dataM$scRNAseq_CV,y=log(1+dataM$RNAseq_CV),ci)

# Split the groups
pdf("genesByGroup/M-plotGroups.pdf")
g1 <- limits[limits[,3]>limits[,6],]
plot(x = g1[,2], y = g1[,3], col = "red", pch = 16, cex= 0.7, ylim=c(0,1), xlab = "CV (scRNA-seq)", ylab= "log(CV) (RNA-seq)")
g3 <- limits[limits[,3]<limits[,5],]
points(x = g3[,2], y = g3[,3], col = "blue", pch = 16, cex= 0.7)
g2 <- limits[((limits[,3]<limits[,6])&(limits[,3]>limits[,5])),]
points(x = g2[,2], y = g2[,3], col = "green", pch = 16, cex= 0.7)
dev.off()

# Output files
writeLines(text = as.vector(dataM$ENSEMBL[g1[,1]]),con='genesByGroup/M-G1.txt')
writeLines(text = as.vector(dataM$ENSEMBL[g2[,1]]),con='genesByGroup/M-G2.txt')
writeLines(text = as.vector(dataM$ENSEMBL[g3[,1]]),con='genesByGroup/M-G3.txt')

### T - CELLS
dataT <- read.csv("dataFiles/dataT.csv")

# Linear model
lmM <- lm(formula = (log(1+RNAseq_CV))~scRNAseq_CV, data = dataT)
plot((log(1+dataT$RNAseq_CV))~dataT$scRNAseq_CV)
abline(lmM, col= "red")

# Confidence intervals
ci <- predict(object = lmM,newdata = data.frame(scRNAseq_CV= dataT$scRNAseq_CV), interval="predict", level = 0.95)
limits <- cbind(gene= dataT$ENSEMBL,x=dataT$scRNAseq_CV,y=log(1+dataT$RNAseq_CV),ci)

# Split the groups
pdf("genesByGroup/T-plotGroups.pdf")
g1 <- limits[limits[,3]>limits[,6],]
plot(x = g1[,2], y = g1[,3], col = "red", pch = 16, cex= 0.7, ylim=c(0,1), xlab = "CV (scRNA-seq)", ylab= "log(CV) (RNA-seq)")
g3 <- limits[limits[,3]<limits[,5],]
points(x = g3[,2], y = g3[,3], col = "blue", pch = 16, cex= 0.7)
g2 <- limits[((limits[,3]<limits[,6])&(limits[,3]>limits[,5])),]
points(x = g2[,2], y = g2[,3], col = "green", pch = 16, cex= 0.7)
dev.off()

# Output files
writeLines(text = as.vector(dataT$ENSEMBL[g1[,1]]),con='genesByGroup/T-G1.txt')
writeLines(text = as.vector(dataT$ENSEMBL[g2[,1]]),con='genesByGroup/T-G2.txt')
writeLines(text = as.vector(dataT$ENSEMBL[g3[,1]]),con='genesByGroup/T-G3.txt')

