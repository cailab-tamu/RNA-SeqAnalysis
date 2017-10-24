#!/bin/bash
# Daniel Osorio
# Cai Lab | Department of Veterinary Integrative Biosciences
# Texas A&M University, College Station, TX.

# Verify the installation
if ! [ -x "$(command -v STAR)" ]; then
  wget -nc https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
  tar -xzf 2.5.3a.tar.gz
  cd STAR-2.5.3a/
  if [[ "$OSTYPE" == "linux-gnu" ]]; then
    export PATH=$(pwd)/bin/Linux_x86_64/:${PATH}
  fi
  if [[ "$OSTYPE" == "darwin"* ]]; then
    export PATH=$(pwd)/bin/MacOSX_x86_64/:${PATH}
  fi
  cd ../
fi

# Download the genome data
wget -nc ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh38.p10.genome.fa.gz
gunzip GRCh38.p10.genome.fa.gz

# Download the annotation data
wget -nc ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
gunzip gencode.v27.annotation.gtf.gz

# STAR index
mkdir genomeDir
STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir genomeDir \
     --genomeFastaFiles GRCh38.p10.genome.fa \
     --sjdbGTFfile gencode.v27.annotation.gtf \
     --sjdbOverhang 99

# STAR mapping
STAR --runThreadN 20 \
     --runMode alignReads \
     --genomeDir genomeDir \
     --readFilesIn $1 ${2:-} \
     --sjdbGTFfile gencode.v27.annotation.gtf \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix $1
