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
  filesDonor <- donorFiles[gsub("^\\_[[:alnum:]]+\\_","",donorFiles) %in% gsub(".cip$","",filesDonor)]
  if (all(filesDonor%in%donorFiles)){
    if (length(filesDonor)==1){
      command <- paste0("salmon quant -i genecodeIndex -l SR -p 8 -r <(gunzip -c ", args[1], filesDonor, ") " , "-o ", args[3], donor)
      writeLines(command, ".command.sh")
      system("bash .command.sh")
    }
  }
  
  # Set up the matrix
  if(!initMatrix){
    initMatrix <- TRUE
    allValues <<- read.csv(paste0(args[3],"/",donor,"/quant.sf"), sep= "\t")[,c("Name","NumReads")]
    colnames(allValues)[donorNumber <- donorNumber+1] <- donor
    next()
  }
  # Add values
  values <- read.csv(paste0(args[3],"/",donor,"/quant.sf"), sep= "\t")[,c("Name","NumReads")]
  allValues <<- merge.data.frame(x = allValues, y = values,by = "Name",all = TRUE)
  colnames(allValues)[donorNumber <- donorNumber+1] <- donor
}
# Writing the output file
write.table(x = allValues, file = "expressionValues.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
