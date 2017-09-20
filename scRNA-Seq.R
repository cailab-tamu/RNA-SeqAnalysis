# scRNA-Seq 10X Data
# Daniel Camilo Osorio
# Biomedical Sciences PhD Program
# Cai Lab - Texas A&M University, College Station, TX.

# Libraries
library(R.utils)

# Activate the input arguments
args <- commandArgs(trailingOnly = TRUE)

# Input files
files <- list.files(args[1])

# gencodeIndex
genCode <- args[2]

# Salmon Index
system(paste0("salmon index -t ", genCode, " -i genecodeIndex"))

# SALMON
files <- list.files(args[1],full.names = TRUE)
for (file in files){
  command <- paste0("salmon quant -i genecodeIndex -l A -p 8 -r <(gunzip -c ", file , ") " , "-o ", basename(file), " --seqBias --gcBias --posBias -g", args[3])
  writeLines(command, ".command.sh")
  system("bash .command.sh")
}

