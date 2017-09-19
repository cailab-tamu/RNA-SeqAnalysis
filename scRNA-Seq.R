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

# UMI
UMI <- gsub("[[:punct:]]","",gsub("^\\-[[:alpha:]]+","",gsub("[[:digit:]]","",gsub("[[:lower:]]","", files))))
write(paste0(">","umi_",UMI,"\n",UMI), file = genCode, append = TRUE)

# Salmon Index
system(paste0("salmon index -t ", genCode, " -i genecodeIndex"))

# SALMON
files <- list.files(args[1],full.names = TRUE)
for (file in files){
  system(paste0("salmon quant -i genecodeIndex -l SR -p 8 -r <(gunzip -c ", file , ") " , "-o ", file ))
}