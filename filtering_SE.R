library(dada2); packageVersion("dada2")
# Filename parsing
path <- "/home/lmarcos/JULIANA/Originales/SRR55_Originales/merged_crudos/" # CHANGE ME to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
fns <- list.files(path, pattern="fastq.gz") # CHANGE if different file extensions
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
              maxEE=2, truncQ=11, minLen=300, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
