correlationAnalysis <- function(dataFile, genesFile, pathwaysFile, nGenes=15, rLimit = 0.85, pLimit = 0.05, CI= 0.95, outFolder){
  data <- read.csv(dataFile)
  if(!dir.exists(outFolder)){
    dir.create(outFolder)
  }
  genes <- read.table(genesFile)
  colnames(genes) <- c("ENSEMBL", "GENECARD")
  data <- merge(genes,data,by = "ENSEMBL")
  pathways <- strsplit(readLines(pathwaysFile),"\t")
  names(pathways) <- unlist(lapply(pathways,function(pathway){pathway[1]}))
  pathways <- lapply(pathways,function(pathway){pathway[3:length(pathway)]})
  data_pathway <- lapply(pathways, function(pathway){
    data[data[,"GENECARD"] %in%pathway,3:14]
  })
  gNames_pathway <- lapply(pathways, function(pathway){
    data[data[,"GENECARD"] %in%pathway,2]
  })
  data_pathway <- data_pathway[unlist(lapply(data_pathway,function(pathway){dim(pathway)[1]>nGenes}))]
  correlations <- cbind(c(1:6),c(7:12))
  cor_pathway <- lapply(data_pathway,function(pathway){
    sapply(seq_len(nrow(correlations)),function(pair){
      x<-pathway[,correlations[pair,1]]
      y<-pathway[,correlations[pair,2]]
      correlationValue <- cor(x,y)
      return(correlationValue)
    })
  })
  cor_pValue <- lapply(data_pathway,function(pathway){
    sapply(seq_len(nrow(correlations)),function(pair){
      x<-pathway[,correlations[pair,1]]
      y<-pathway[,correlations[pair,2]]
      correlationValue <- cor.test(x,y)$p.value
      return(correlationValue)
    })
  })
  cor_pathway <- t(as.data.frame(cor_pathway))
  cor_pValue <- t(as.data.frame(cor_pValue))
  colnames(cor_pathway) <- sapply(seq_len(nrow(correlations)),function(pair){
    paste0(colnames(data)[correlations[pair,1]+2],"|",colnames(data)[correlations[pair,2]+2])
  })
  colnames(cor_pValue) <- colnames(cor_pathway)
  for(i in rownames(cor_pathway)){
    for(j in colnames(cor_pathway)){
      correlationValue <- cor_pathway[i,j]
      pValue <- cor_pValue[i,j]
      variables <- unlist(strsplit(j,"\\|"))
      xlabel <- unlist(strsplit(variables[2],"\\_"))[c(2,1)]
      xlabel <- gsub("RNAseq","RNA-seq",paste0(xlabel[1]," (",xlabel[2],")"))
      ylabel <- unlist(strsplit(variables[1],"\\_"))[c(2,1)]
      ylabel <- gsub("RNAseq","RNA-seq",paste0(ylabel[1]," (",ylabel[2],")"))
      glabel <- as.vector(gNames_pathway[[i]])
      if(!is.na(correlationValue) & correlationValue > rLimit & pValue < pLimit){
        
        title <- paste0(unlist(strsplit(i,"\\_"))[-1], collapse = " ")
        subtitle <- paste(unlist(strsplit(i,"\\_"))[1], "DATABASE", sep = " ")
        ggplot(data_pathway[[i]], aes(x=data_pathway[[i]][,variables[2]], y=data_pathway[[i]][,variables[1]])) +
          geom_point(shape=16) +
          geom_smooth(method=lm, level = CI) +
          theme_bw() +
          geom_text(aes(label=glabel),hjust=0.5, vjust=-1, cex=1.5) +
          labs(x = xlabel,
               y= ylabel,
               caption = paste("r=", round(correlationValue,4), "p=", pValue, sep=" ")) +ggtitle(label = title, subtitle = subtitle)
        
        ggsave(paste0(outFolder,"/",i,"-(",variables[1],"-",variables[2],")",".pdf"))
      }
    }
    
  }
}
  
correlationAnalysis(dataFile = "dataFiles/dataM.csv",
                    genesFile = "dataFiles/genes.tsv", 
                    pathwaysFile = "dataFiles/pathways.txt",
                    nGenes = 15, 
                    rLimit = 0.9, 
                    pLimit = 0.005,
                    CI = 0.95,
                    outFolder = "pathwayAnalysis/M/")

correlationAnalysis(dataFile = "dataFiles/dataT.csv",
                    genesFile = "dataFiles/genes.tsv", 
                    pathwaysFile = "dataFiles/pathways.txt",
                    nGenes = 15, 
                    rLimit = 0.9, 
                    pLimit = 0.005,
                    CI = 0.95,
                    outFolder = "pathwayAnalysis/T/")