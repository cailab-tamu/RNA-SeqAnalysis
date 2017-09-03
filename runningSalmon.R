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
}
