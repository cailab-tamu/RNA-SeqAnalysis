#!/bin/bash
# transcriptExpressionMatrix
# Daniel Camilo Osorio
# Biomedical Sciences PhD Program
# Cai Lab - Texas A&M University, College Station, TX.
# Summarize the expression values into a matrix


# Argument to associate the folder where the fastq files are placed
export inputData = $1

# Argument to associate the mapFile
export mapFile = $2

# Argument to associate an output folder
export outputFolder = $3

# Download the transcript from GeneCodeGenes Release 19 (GRCh37.p13)
wget -nc ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz

# Unzip the downloaded file
gunzip gencode.v19.pc_transcripts.fa.gz

# Enter into the outputFolder
cd $outputFolder

# Create the index for the transcript
salmon index -t ../gencode.v19.pc_transcripts.fa -i genecodeIndex

# Create a output folder for expressionValues
mkdir expressionValues

# Running SALMON
Rscript --vanilla runningSALMON $1 $2 expressionValues/
