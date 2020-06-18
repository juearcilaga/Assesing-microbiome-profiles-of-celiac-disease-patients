library(dada2); packageVersion("dada2")
# File parsing
pathF <- "/home/lmarcos/JULIANA/DATOS/FW" # CHANGE to the directory containing demultiplexed forward-read fastqs
pathR <- "/home/lmarcos/JULIANA/DATOS/RV" # CHANGE  to the directory containing demultiplexed reverse-read fastqs
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # Filtered reverse files go into the pathR/filtered/ subdirectory
fastqFs <- sort(list.files(pathF, pattern="_R1_001.fastq.gz")) #CHANGE if different extension
fastqRs <- sort(list.files(pathR, pattern="_R2_001.fastq.gz")) #...

#head(file.path(pathF, fastqFs))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS YOU SHOULD CHECK WITH FASTQC
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
