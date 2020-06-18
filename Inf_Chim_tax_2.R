#This is an alternative to Inf_Chim_tax.R, use it if you have an amount of data that makes the analysis to run out of memory

library(dada2); packageVersion("dada2")
# File parsing
filtpath <- "/home/lmarcos/JULIANA/Originales/SRR55_Originales/merged_crudos/filtered/filtered"
# CHANGE ME to the directory  containing your filtered fastq files
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE)
# CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1)
# Assumes filename = sample_XXX.fastq.gz
names(filts) <- sample.names
# Learn forward error rates
errs <- learnErrors(filts, nbases=1e8, multithread=TRUE)
save.image(file="infer_variants_SRR55_learnErrors.RData")
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=errs, multithread=TRUE)
}
# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "/home/lmarcos/JULIANA/Originales/SRR55_Originales/seqtab.rds")
# CHANGE ME to where you want sequence table saved
saveRDS(dds, "/home/lmarcos/JULIANA/Originales/SRR55_Originales/dadaFstab.rds")
# CHANGE ME to where you want sequence table saved
# Remove chimeras
NOTchimeras <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
#NOTchimeras <- removeBimeraDenovo(dds, method="consensus", multithread=TRUE)
saveRDS(NOTchimeras, file= "/home/lmarcos/JULIANA/Originales/SRR55_Originales/NOTchimeras.rds.gz")
# CHANGE ME to where you want sequence table saved
# Assign taxonomy
tax <- assignTaxonomy(NOTchimeras, "/home/lmarcos/JULIANA/silva_nr_v132_train_set.fa", multithread=TRUE)
save.image("ChimerasTaxonomy.RData")
# Write to disk
saveRDS(tax, file= "home/lmarcos/JULIANA/Originales/SRR55_Originales/tax_SRR55_final.rds") # CHANGE ME ...
