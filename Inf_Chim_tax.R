library(dada2); packageVersion("dada2")
# File parsing
filtpath <- "/home/lmarcos/JULIANA/Originales/SRR55_Originales/merged_crudos/filtered/filtered" #Path to filtered reads
# CHANGE ME to the directory  containing your filtered fastq files 
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) 
# CHANGE if different file extensions 
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1) 
# Assumes filename = sample_XXX.fastq.gz 
names(filts) <- sample.names
# Learn error rates
set.seed(100)
errs <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE) 
save.image(file="infer_variants_SRR55_learnErrors.RData")
# Infer sequence variants
derep= derepFastq(filtpath, n = 1e+06, verbose = TRUE, qualityType = "Auto")
#You can also use filtfq in derep 
filtFq <- sort(list.files(filtpath, pattern=".fastq.gz"))
dadaFs <- dada(derep, err=errs, multithread=TRUE)
# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dadaFs) 
saveRDS(seqtab, "/home/lmarcos/JULIANA/Originales/SRR55_Originales/seqtab.rds") 
# CHANGE ME to where you want sequence table saved 
saveRDS(dadaFs, "/home/lmarcos/JULIANA/Originales/SRR55_Originales/dadaFstab.rds") 
# Remove chimeras
NOTchimeras <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
#NOTchimeras <- removeBimeraDenovo(dadaFs, method="consensus", multithread=TRUE)
saveRDS(NOTchimeras, file= "/home/lmarcos/JULIANA/Originales/SRR55_Originales/NOTchimeras.rds.gz") 
# CHANGE ME to where you want sequence table saved
# Assign taxonomy
tax <- assignTaxonomy(NOTchimeras, "/home/lmarcos/JULIANA/silva_nr_v132_train_set.fa", multithread=TRUE) 
save.image("ChimerasTaxonomy.RData")
# Write to disk
saveRDS(tax, file= "home/lmarcos/JULIANA/Originales/SRR55_Originales/tax_SRR55_final.rds") # CHANGE ME ...
