#install.packages("ggplot2")
#install.packages("readr")
library(ggplot2)
library(parallel)
library(readr)
library(Matrix)
##############################################
# NEW FUNCTIONS
##############################################
qN <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
giniCoefficient <- function(data){sum(abs(apply(expand.grid(data,-1*data),1,sum)))/(2*length(data)*sum(data))}
npMEAN <- function(data, replicates=10000) {
  cluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
  output <- mean(parSapply(cluster,seq_len(replicates), function(x){mean(sample(x = data,size = length(data),replace = TRUE))}))
  stopCluster(cluster)
  return(output)
}
npSD <- function(data, replicates=10000){
  cluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
  output <- mean(parSapply(cluster,seq_len(replicates), function(x){sd(sample(x = data,size = length(data),replace = TRUE))}))
  stopCluster(cluster)
  return(output)
}
npVAR <- function(data, replicates=10000){
  cluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
  output <- mean(parSapply(cluster,seq_len(replicates), function(x){var(sample(x = data,size = length(data),replace = TRUE))}))
  stopCluster(cluster)
  return(output)
}
npGC <- function(data, replicates=10000){
  cluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
  output <- mean(parSapply(makeCluster(4),seq_len(replicates), function(x){giniCoefficient(sample(x = data,size = length(data),replace = TRUE))}))
  stopCluster(cluster)
  return(output)
}
##############################################
# M DATA
##############################################
# download.file(url = "https://raw.githubusercontent.com/cailab-tamu/geDatasets/master/BLUEPRINT/M-expressionMatrix_TPM%2BQN%2BLogT.csv",
#               destfile = "dataFiles/M-expressionMatrix_TPM+QN+LogT.csv"
#               )
M_BLUEPRINT <- read_csv("dataFiles/M-expressionMatrix_TPM+QN+LogT.csv")
bpM_ENSEMBL <- M_BLUEPRINT$ENSEMBL
bpM_MEAN <- apply(M_BLUEPRINT[,3:ncol(M_BLUEPRINT)],1,npMEAN)
bpM_SD <- apply(M_BLUEPRINT[,3:ncol(M_BLUEPRINT)],1,npSD)
bpM_VAR <- apply(M_BLUEPRINT[,3:ncol(M_BLUEPRINT)],1,npVAR)
bpM_CV <- bpM_SD/bpM_MEAN
bpM_FF <- bpM_VAR/bpM_MEAN
bpM_GC <- apply(M_BLUEPRINT[,3:ncol(M_BLUEPRINT)],1,npGC)
output <- cbind(ENSEMBL = bpM_ENSEMBL, MEAN=bpM_MEAN, SD=bpM_SD, VAR=bpM_VAR, CV = bpM_CV, FF=bpM_FF, GC = bpM_GC)
write.csv(output,"dataFiles/BLUEPRINT_M.csv", quote = FALSE, row.names = FALSE)
M_10XG <- log(1+qN(readMM("dataFiles/10XGENOMICS_M/matrix.mtx")))
x10gM_ENSEMBL <- read_delim("dataFiles/10XGENOMICS_M/genes.tsv",delim = "\t",col_names = FALSE)[[1]]
x10gM_MEAN <- apply(M_10XG,1,npMEAN)
x10gM_SD <- apply(M_10XG,1,npSD)
x10gM_VAR <- apply(M_10XG,1,npVAR)
x10gM_CV <- apply(M_10XG,1,npCV)
x10gM_GC <- apply(M_10XG,1,npGC)
output <- cbind(ENSEMBL = x10gM_ENSEMBL, MEAN=x10gM_MEAN, SD=x10gM_SD, VAR=x10gM_VAR, CV = x10gM_CV, FF=x10gM_FF, GC = x10gM_GC)
write.csv(output,"dataFiles/10XGENOMICS_Monocytes.csv", quote = FALSE, row.names = FALSE)
# 
# ##############################################
# # T DATA
# ##############################################
# # download.file(url = "https://raw.githubusercontent.com/cailab-tamu/geDatasets/master/BLUEPRINT/T-expressionMatrix_TPM%2BQN%2BLogT.csv",
# #               destfile = "dataFiles/T-expressionMatrix_TPM+QN+LogT.csv"
# #               )
# T_BLUEPRINT <- read_csv("dataFiles/T-expressionMatrix_TPM+QN+LogT.csv")
# bpT_ENSEMBL <- T_BLUEPRINT$ENSEMBL
# bpT_MEAN <- apply(T_BLUEPRINT[,3:ncol(T_BLUEPRINT)],1,npMEAN)
# bpT_SD <- apply(T_BLUEPRINT[,3:ncol(T_BLUEPRINT)],1,npSD)
# bpT_VAR <- apply(T_BLUEPRINT[,3:ncol(T_BLUEPRINT)],1,npVAR)
# bpT_CV <- bpT_SD/bpT_MEAN
# bpT_FF <- bpT_VAR/bpT_MEAN
# bpT_GC <- apply(T_BLUEPRINT[,3:ncol(T_BLUEPRINT)],1,npGC)
# output <- cbind(ENSEMBL = bpT_ENSEMBL, MEAN=bpT_MEAN, SD=bpT_SD, VAR=bpT_VAR, CV = bpT_CV, FF=bpT_FF, GC= bpT_GC)
# write.csv(output,"dataFiles/BLUEPRINT_T.csv", quote = FALSE, row.names = FALSE)

# ##############################################
# # T PLOTS
# ##############################################
# BLUEPRINT_T <- read_csv("dataFiles/BLUEPRINT_T.csv")
# NaiveT <- read_csv("dataFiles/10XGENOMICS_NaiveT.csv")
# merged <- merge(BLUEPRINT_T,NaiveT,by = "ENSEMBL")
# merged <- merged[complete.cases(merged),]
# colnames(merged) <- c("ENSEMBL","RNAseq_MEAN","RNAseq_SD", "RNAseq_VAR", "RNAseq_CV",
#                       "RNAseq_FF", "RNAseq_GC", "scRNAseq_MEAN",
#                       "scRNAseq_SD", "scRNAseq_VAR", "scRNAseq_CV", "scRNAseq_FF", "scRNAseq_GC")
# write.csv(merged, "dataFiles/dataT.csv", row.names = FALSE, quote = FALSE)
# 
# pdf("allDataAnalysis/T-CVvsCV.pdf")
# x <- merged$scRNAseq_CV
# y <- log(1+merged$RNAseq_CV)
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "CV (scRNA-seq)", 
#        y="log(CV) (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/T-MEANvsMEAN.pdf")
# x <- merged$scRNAseq_MEAN
# y <- merged$RNAseq_MEAN
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "MEAN (scRNA-seq)", 
#        y="MEAN (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/T-SDvsSD.pdf")
# x <- merged$scRNAseq_SD
# y <- merged$RNAseq_SD
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "SD (scRNA-seq)", 
#        y="SD (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/T-VARvsVAR.pdf")
# x <- merged$scRNAseq_VAR
# y <- merged$RNAseq_VAR
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "VAR (scRNA-seq)", 
#        y="VAR (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/T-FFvsFF.pdf")
# x <- merged$scRNAseq_FF
# y <- merged$RNAseq_FF
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "FF (scRNA-seq)", 
#        y="FF (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/T-GCvsGC.pdf")
# x <- merged$scRNAseq_GC
# y <- merged$RNAseq_GC
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "GC (scRNA-seq)", 
#        y="GC (RNA-seq)")
# dev.off()
# ##############################################
# # M PLOTS
# ##############################################
# BLUEPRINT_M <- read_csv("dataFiles/BLUEPRINT_M.csv")
# Monocytes <- read_csv("dataFiles/10XGENOMICS_Monocytes.csv")
# merged <- merge(BLUEPRINT_M,Monocytes,by = "ENSEMBL")
# merged <- merged[complete.cases(merged),]
# colnames(merged) <- c("ENSEMBL","RNAseq_MEAN","RNAseq_SD", "RNAseq_VAR", "RNAseq_CV", "RNAseq_FF",
#                       "RNAseq_GC", "scRNAseq_MEAN", "scRNAseq_SD", "scRNAseq_VAR", "scRNAseq_CV",
#                       "scRNAseq_FF", "scRNAseq_GC")
# write.csv(merged, "dataFiles/dataM.csv", row.names = FALSE, quote = FALSE)
# 
# library(ggplot2)
# pdf("allDataAnalysis/M-CVvsCV.pdf")
# x <- merged$scRNAseq_CV
# y <- log(1+merged$RNAseq_CV)
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "CV (scRNA-seq)", 
#        y="log(CV) (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/M-MEANvsMEAN.pdf")
# x <- merged$scRNAseq_MEAN
# y <- merged$RNAseq_MEAN
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "MEAN (scRNA-seq)", 
#        y="MEAN (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/M-SDvsSD.pdf")
# x <- merged$scRNAseq_SD
# y <- merged$RNAseq_SD
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "SD (scRNA-seq)", 
#        y="SD (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/M-VARvsVAR.pdf")
# x <- merged$scRNAseq_VAR
# y <- merged$RNAseq_VAR
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "VAR (scRNA-seq)", 
#        y="VAR (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/M-FFvsFF.pdf")
# x <- merged$scRNAseq_FF
# y <- merged$RNAseq_FF
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "FF (scRNA-seq)", 
#        y="FF (RNA-seq)")
# dev.off()
# 
# pdf("allDataAnalysis/M-GCvsGC.pdf")
# x <- merged$scRNAseq_GC
# y <- merged$RNAseq_GC
# data <- as.data.frame(cbind(x,y))
# ggplot(data, aes(x=x, y=y)) +
#   geom_point(shape=16) +
#   geom_smooth(method=lm, level = 0.95) + 
#   theme_bw() +
#   labs(x = "GC (scRNA-seq)", 
#        y="GC (RNA-seq)")
# dev.off()