library(dada2); packageVersion("dada2")
# Filename parsing
path <- "/home/lmarcos/JULIANA/SINGLE/" # Change to the directory containing demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
fns <- list.files(path, pattern="fastq.gz") # CHANGE if different file extensions
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
              truncLen=400, maxEE=2, truncQ=20, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
