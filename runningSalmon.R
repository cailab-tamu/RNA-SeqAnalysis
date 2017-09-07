# runningSalmon
# Daniel Camilo Osorio
# Biomedical Sciences PhD Program
# Cai Lab - Texas A&M University, College Station, TX.

# Activate the input arguments
args <- commandArgs(trailingOnly = TRUE)

# List the files into the inputFolder
donorFiles <- list.files(args[1])

# List the donors and the associated FASTQ files
fileMap <- read.csv(file = args[2],
                    sep = "\t",
                    header = FALSE,
                    stringsAsFactors = FALSE)

# Identify the unique donors
uniqueDonors <- unique(fileMap[,1]) 

# Flag for init the matrix
initMatrix <- FALSE
donorNumber <- 1

# Run SALMON
for (donor in uniqueDonors){
  filesDonor <- fileMap[fileMap[,1]==donor,3]
  filesDonor <- filesDonor[grepl("fastq",filesDonor)]
  filesDonor <- donorFiles[gsub("^\\_[[:alnum:]]+\\_","",donorFiles) %in% gsub("\\/","_",gsub(".cip$","",filesDonor))]
  if (all(filesDonor%in%donorFiles)){
    # If single end
    if (length(filesDonor)==1){
      command <- paste0("salmon quant -i genecodeIndex -l SR -p 8 -r <(gunzip -c ", args[1], filesDonor, ") " , "-o ", args[3], donor)
      writeLines(command, ".command.sh")
      # system("bash .command.sh")
      # If pair end
    } else if(length(filesDonor)==2){
      command <- paste0("salmon quant -i genecodeIndex -l UI -p 8 -1 <(gunzip -c ", args[1], filesDonor[grepl("pair1",filesDonor)], ") -2 <(gunzip -c ", args[1], filesDonor[grepl("pair2",filesDonor)], ") -o ", args[3], donor)
      writeLines(command, ".command.sh")
      # system("bash .command.sh")
    }
  }

  # Set up the matrix
  if(!initMatrix){
    if(file.exists(paste0(args[3],"/",donor,"/quant.sf"))){
      initMatrix <- TRUE
      allValues <<- read.csv(paste0(args[3],"/",donor,"/quant.sf"), sep= "\t")[,c("Name","NumReads")]
      colnames(allValues)[donorNumber <- donorNumber+1] <- donor
      next()
    }
  }
  # Add values
  if(file.exists(paste0(args[3],"/",donor,"/quant.sf"))){
    values <- read.csv(paste0(args[3],"/",donor,"/quant.sf"), sep= "\t")[,c("Name","NumReads")]
    allValues <<- merge.data.frame(x = allValues, y = values,by = "Name",all = TRUE)
    colnames(allValues)[donorNumber <- donorNumber+1] <- donor
  }
}
# Write the output file
write.table(x = allValues, file = "transcriptExpressionValues.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

# Create a dictionary of the Ensambl Gene IDs
GENEID <- gsub("\\.[[:digit:]]+$","",unlist(lapply(strsplit(as.vector(allValues[,"Name"]), "\\|"), function(id){id[2]})))
TXNAME <- as.vector(allValues[,"Name"])
TX2GENE <- as.data.frame(cbind(TXNAME,GENEID))

# Read the SALMON output files
files <- file.path("expressionValues",uniqueDonors,"quant.sf")
fileExists <- file.exists(files)
files <- files[fileExists]

# Calculate the geneExpression
library(tximport)
geneExpression <- tximport(files = files,type = "salmon",tx2gene = TX2GENE)
colnames(geneExpression$counts) <- uniqueDonors[fileExists]

# Write the output file
write.table(x = geneExpression$counts, file = "geneExpressionValues.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Quantile Normalization
library(preprocessCore)
qN <- normalize.quantiles(as.matrix(geneExpression))
colnames(qN) <- colnames(geneExpression)
rownames(qN) <- rownames(geneExpression)

# Log-Transformation
logT <- log2(qN+1)

# Write the output file
write.table(x = logT, file = "qN-LogT_ExpressionValues.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
