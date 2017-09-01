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

# Running SALMON
for (donor in uniqueDonors){
  donorFiles <- fileMap[fileMap[,1]==donor,3]
  donorFiles <- donorFiles[grepl("fastq",donorFiles)]
  donorFiles <- files[gsub("^\\_[[:alnum:]]+\\_","",files) %in% gsub(".cip$","",donorFiles)]
  if (all(donorFiles%in%files)){
    if (length(donorFiles)==1){
      command <- paste0("salmon quant -i genecodeIndexSalmon -l SR -p 8 -r <(gunzip -c ", args[1], donorFiles, ") " , "-o ", args[3], "/", donor)
      writeLines(command, ".command.sh")
      system("bash .command.sh")
    }
  }
}