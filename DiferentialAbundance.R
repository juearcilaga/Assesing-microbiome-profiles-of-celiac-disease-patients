############################# ANÁLISIS DUODENO ##############################################

#It will be convenient to make sure that Control is the first level in the treatment factor, so that the default log2 fold changes are #calculated as treatment over control and not the other way around. The function relevel achieves this:
#dds$treatment <- relevel( dds$treatment, "Control" )
#as.data.frame( colData(degs_duo) )

#Calcular diferencias 
degs_duo = DESeq(degs_duo, test="Wald", fitType="local");
sigtab_duo = results(degs_duo, cooksCutoff = FALSE);
colnames(sigtab_duo)
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
#[5] "pvalue"         "padj"          
#[1] 518   6

#Filtrar diferencias significativas
alpha = 0.01;
sigtab = sigtab_duo[which(sigtab_duo$padj < alpha), ];
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy.duo.mas1)[rownames(sigtab), ], "matrix"));
dim(sigtab)
#[1] 133  12 #otus*filas

colnames(sigtab)
 #[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
 #[5] "pvalue"         "padj"           "Kingdom"        "Phylum"        
 #[9] "Class"          "Order"          "Family"         "Genus"  

#Exportar tabla para anexos
write.table(sigtab, "sigtab_duodeno.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);

#Se filtraron las VSA para cuyo análisis de expsigtabión diferencial de la abundancia arrojó un valor p (corregido)  < 0.01  y luego se aplicó un #límite de cambio de abundancia "Fold change" de 3, conservando los VSA con Fold Change >3 y < (-3)

sigtab_minus <-filter(sigtab, log2FoldChange <= "-3" );
sigtab_plus <-filter(sigtab, log2FoldChange >= "3" );
dim(sigtab_minus);
#[1] 90 12
dim(sigtab_plus);
#[1] 1 12
#sigtab 
#Filtrar la tabla en base al foldchange, los 20 mas altos y los 20 mas bajos
sigtab_filtered <-rbind(sigtab_plus,sigtab_minus);

#Exportar tabla para artículo
write.table(sigtab_filtered, "sigtab_filtered_duodeno.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE);

#Graficar
theme_set(theme_light());
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
};
# Phylum order
x = tapply(sigtab_filtered$log2FoldChange, sigtab_filtered$Phylum, function(x) max(x));
x = sort(x, TRUE);
sigtab_filtered$Phylum = factor(as.character(sigtab_filtered$Phylum), levels=names(x));

# Genus order
x = tapply(sigtab_filtered$log2FoldChange, sigtab_filtered$Genus, function(x) max(x));
x = sort(x, TRUE);
sigtab_filtered$Genus = factor(as.character(sigtab_filtered$Genus), levels=names(x));
ggplot(sigtab_filtered, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5));


tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/Deseq_duo_caco.tiff", width = 22, height = 12, units = "cm", res = 600, pointsize = 8)
ggplot(sigtab_filtered, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5));
dev.off()



tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/Deseq_stool_caco.tiff", width = 20, height = 12, units = "cm", res = 600, pointsize = 10)
ggplot(sigtab_stool_filtered, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5));
dev.off()



tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/volcano_duo_caco.tiff", width = 12, height = 12, units = "cm", res = 600, pointsize = 10)
# Make a basic volcano plot
with(sigtab_duo, plot(log2FoldChange, -log10(pvalue), pch=21, main="Abundancia diferencial - duodeno", ylim=c(0,20), xlim=c(-5,5), cex =1.5, col="turquoise4" ))
# Add colored points: tomato if padj<0.01, yellow of log2FC>3, turquoise if both)
with(subset(sigtab_duo, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, cex =1, col= "tomato"))
#with(subset(sigtab_duo, padj>=.01 ), points(log2FoldChange, -log10(pvalue), pch=20, cex= 0.7, col= "darkorchid1"))
with(subset(sigtab_duo, padj<.01 & abs(log2FoldChange)>=3), points(log2FoldChange, -log10(pvalue), pch=20, cex=1, col="darkorchid1"))
dev.off()


# Label points with the textxy function from the calibrate plot
#library(calibrate);
#with(subset(sigtab_duo2, padj<.01 & log2FoldChange<=-2 ), textxy(log2FoldChange, -log10(pvalue), labs= Genus, cex=0))
#with(subset(sigtab_duo2, padj<.01 & log2FoldChange>=4), textxy(log2FoldChange, -log10(pvalue), labs= Genus, cex=0))

##################################################################################

#################################################
##ANALISIS DE RUTAS CON PICRUST
phy_duodenum;
otu_duo<- otu_table(phy_duodenum, taxa_are_row=TRUE);
otu_melt<-psmelt(otu_duo);

# Export taxonomy table as "tax.txt"
tax<-as(tax_table(phy_duodenum),"matrix");
tax_cols <- colnames(tax);
tax<-as.data.frame(tax);
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"));
for(co in tax_cols) tax[co]<-NULL;
write.table(tax, "tax_duo.tsv", quote=FALSE, col.names=TRUE, row.names =TRUE, sep="\t");

# Export feature/OTU table

# As a biom file

library(biomformat);packageVersion("biomformat");
## [1] ‘1.6.0’


#otu<-t(as(OTU_vst_duo,"matrix")); # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
#otu<-as(otu_duo,"matrix");
otu<-t(as(otu_duo,"matrix"));
otu_biom<-make_biom(data=otu);
write_biom(otu_biom,"otu_duodeno_biom.biom");

# As a text file

write.table(t(otu_duo), "otu_duodeno.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
#or from the phyloseq object, 't' to transform if taxa_are_rows=FALSE, no 't' if taxa_are_rows=TRUE
#write.table(t(otu_table(ps), "seqtab.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Export metadata (if you have a properly formatted metadata file that you imported in your phyloseq pipeline, you can skip this step and just use that text file directly in QIIME 2)

write.table(sample_data(phy_duodenum),"sample-metadata_duodeno.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE);

#Para extraer las secuencias de las ASVs
library(dada2); packageVersion("dada2");
featureseq <-otu_melt %>% select(OTU,Abundance)
dim(featureseq)
#[1] 69930     2

colnames(featureseq)<-c("sequence","abundance");
featureseq_uniq <- getUniques(featureseq, collapse = TRUE, silence = FALSE)

featureseq2<-cbind(rownames(featureseq_uniq),featureseq_uniq)
colnames(featureseq2)<-c("sequence","abundance");
featureseq2<-as.data.frame(featureseq2);
uniquesToFasta(featureseq2, fout='duo_rep-seqs.fna', ids=featureseq2$sequence);

#son iguales x[!(x %in% y)]


############################# ANÁLISIS SALIVA ##############################################
library(DESeq2);

#Calcular diferencias 
degs_sal = DESeq(degs_sal, test="Wald", fitType="local");
sigtab_sal = results(degs_sal, cooksCutoff = FALSE);

#Filtrar diferencias significativas
alpha = 0.01;
sigtab_sal2 = sigtab_sal[which(sigtab_sal$padj < alpha), ];
sigtab_sal2 = cbind(as(sigtab_sal2, "data.frame"), as(tax_table(phy.saliva.mas1)[rownames(sigtab_sal2), ], "matrix"));
dim(sigtab_sal2)
#[1] 2  12 #otus*filas
#Simonsiella
#Oceanivirga

#Exportar tabla para anexos
write.table(sigtab_sal2, "sigtab_saliva.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/volcano_sal_caco.tiff", width = 12, height = 12, units = "cm", res = 600, pointsize = 10)
# Make a basic volcano plot
with(sigtab_sal, plot(log2FoldChange, -log10(pvalue), pch=21, main="Abundancia diferencial - saliva", ylim=c(0,10), cex =1.5, col="turquoise4"))
# Add colored points: tomato if padj<0.01, yellow of log2FC>3, turquoise if both)
with(subset(sigtab_sal, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20,cex=1, col= "tomato"))
#with(subset(sigtab_sal, abs(log2FoldChange)>=3), points(log2FoldChange, -log10(pvalue), pch=20, col= "yellow4"))
with(subset(sigtab_sal, padj<.01 & abs(log2FoldChange)>=3), points(log2FoldChange, -log10(pvalue), pch=20, cex=1, col="darkorchid1"))
dev.off()

#####################################################

#################################################
##ANALISIS DE RUTAS CON PICRUST
otu_sal<- otu_table(phy_saliva, taxa_are_row=TRUE);
otu_sal_melt<-psmelt(otu_sal);

# Export taxonomy table as "tax.txt"
tax<-as(tax_table(phy_saliva),"matrix");
tax_cols <- colnames(tax);
tax<-as.data.frame(tax);
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"));
for(co in tax_cols) tax[co]<-NULL;
write.table(tax, "tax_sal.txt", quote=FALSE, col.names=TRUE, sep="\t");

# Export feature/OTU table

# As a biom file

library(biomformat);packageVersion("biomformat");
## [1] ‘1.6.0’

#otu<-t(as(otu_sal,"matrix")); # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
otu<-t(as(otu_sal,"matrix"));
otu_biom<-make_biom(data=otu);
write_biom(otu_biom,"otu_saliva_biom.biom");

# As a text file

write.table(t(otu_sal), "otu_saliva.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
#or from the phyloseq object, 't' to transform if taxa_are_rows=FALSE, no 't' if taxa_are_rows=TRUE
#write.table(t(otu_table(ps), "seqtab.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Export metadata (if you have a properly formatted metadata file that you imported in your phyloseq pipeline, you can skip this step and just use that text file directly in QIIME 2)

write.table(sample_data(phy_saliva),"sample-metadata_saliva.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE);

#Para extraer las secuencias de las ASVs
library(dada2); packageVersion("dada2");
featureseq <-otu_sal_melt %>% select(OTU,Abundance)
dim(featureseq)
#[1] 8505     2

colnames(featureseq)<-c("sequence","abundance");
featureseq_uniq <- getUniques(featureseq, collapse = TRUE, silence = FALSE)
featureseq2<-cbind(rownames(featureseq_uniq),featureseq_uniq)
colnames(featureseq2)<-c("sequence","abundance");
featureseq2<-as.data.frame(featureseq2);
uniquesToFasta(featureseq2, fout='sal_rep-seqs.fna', ids=featureseq2$sequence);


#########################################################################################

############################# ANÁLISIS HECES ##############################################
#Calcular diferencias 
degs_stool = DESeq(degs_stool, test="Wald", fitType="local");
sigtab_stool = results(degs_stool, cooksCutoff = FALSE);

#Filtrar diferencias significativas
alpha = 0.01;
sigtab_stool2 = sigtab_stool[which(sigtab_stool$padj < alpha), ];
sigtab_stool2 = cbind(as(sigtab_stool2, "data.frame"), as(tax_table(phy.stool.mas1)[rownames(sigtab_stool2), ], "matrix"));
dim(sigtab_stool2)
#[1] 88 12 #otus*filas

#Exportar tabla para anexos
write.table(sigtab_stool2, "sigtab_stool.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE);

#Se filtraron las VSA para cuyo análisis de expsigtabión diferencial de la abundancia arrojó un valor p (corregido)  < 0.01  y luego se aplicó un #límite de cambio de abundancia "Fold change" de 3, conservando los VSA con Fold Change >3 y < (-3)

sigtab_stool_minus <-filter(sigtab_stool2, log2FoldChange <= "-3" );
sigtab_stool_plus <-filter(sigtab_stool2, log2FoldChange >= "3" );
dim(sigtab_stool_minus);
#[1] 46 12
dim(sigtab_stool_plus);
#[1] 2 12

#Filtrar la tabla en base al foldchange
sigtab_stool_filtered <-rbind(sigtab_stool_plus,sigtab_stool_minus);
#48 12
#Exportar tabla para artículo
write.table(sigtab_stool_filtered, "sigtab_filtered_stool.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE);

#Graficar
theme_set(theme_light());
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
};
# Phylum order
x = tapply(sigtab_stool_filtered$log2FoldChange, sigtab_stool_filtered$Phylum, function(x) max(x));
x = sort(x, TRUE);
sigtab_stool_filtered$Phylum = factor(as.character(sigtab_stool_filtered$Phylum), levels=names(x));

# Genus order
x = tapply(sigtab_stool_filtered$log2FoldChange, sigtab_stool_filtered$Genus, function(x) max(x));
x = sort(x, TRUE);
sigtab_stool_filtered$Genus = factor(as.character(sigtab_stool_filtered$Genus), levels=names(x));
ggplot(sigtab_stool_filtered, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5));


tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/Deseq_stool_caco.tiff", width = 20, height = 12, units = "cm", res = 600, pointsize = 10)
ggplot(sigtab_stool_filtered, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5));
dev.off()


tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/volcano_stool_caco.tiff", width = 12, height = 12, units = "cm", res = 600, pointsize = 10)

# Make a basic volcano plot
#with(sigtab_duo, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,10)))
with(sigtab_stool, plot(log2FoldChange, -log10(pvalue),  pch=21, main="Abundancia diferencial - heces", cex =1.5, col="turquoise4" ))
# Add colored points: tomato if padj<0.01, yellow of log2FC>3, turquoise if both)
with(subset(sigtab_stool, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, cex =1,col= "tomato"))
#with(subset(sigtab_stool, abs(log2FoldChange)>=3), points(log2FoldChange, -log10(pvalue), pch=20, col= "yellow4"))
with(subset(sigtab_stool, padj<.01 & abs(log2FoldChange)>=3), points(log2FoldChange, -log10(pvalue),pch=20, cex=1, col="darkorchid1"))
dev.off()



#######################################################################################

#################################################
##ANALISIS DE RUTAS CON PICRUST
otu_stool<- otu_table(phy_stool, taxa_are_row=TRUE);
otu_stool_melt<-psmelt(otu_stool);

# Export taxonomy table as "tax.txt"
tax<-as(tax_table(phy_stool),"matrix");
tax_cols <- colnames(tax);
tax<-as.data.frame(tax);
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"));
for(co in tax_cols) tax[co]<-NULL;
write.table(tax, "tax_stool.txt", quote=FALSE, col.names=NA, sep="\t");

# Export feature/OTU table

# As a biom file

library(biomformat);packageVersion("biomformat");
## [1] ‘1.6.0’

#otu<-t(as(OTU_vst_duo,"matrix")); # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
otu<-t(as(otu_stool,"matrix"));
otu_biom<-make_biom(data=otu);
write_biom(otu_biom,"otu_stool_biom.biom");

# As a text file

write.table(t(otu_stool), "otu_stool.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
#or from the phyloseq object, 't' to transform if taxa_are_rows=FALSE, no 't' if taxa_are_rows=TRUE
#write.table(t(otu_table(ps), "seqtab.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Export metadata (if you have a properly formatted metadata file that you imported in your phyloseq pipeline, you can skip this step and just use that text file directly in QIIME 2)

write.table(sample_data(phy_stool),"sample-metadata_stool.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE);

#Para extraer las secuencias de las ASVs
library(dada2); packageVersion("dada2");
featureseq <-otu_stool_melt_vst %>% select(OTU,Abundance)
dim(featureseq)
#[1] 67098    2


colnames(featureseq)<-c("sequence","abundance");
featureseq_uniq <- getUniques(featureseq, collapse = TRUE, silence = FALSE)
featureseq2<-cbind(rownames(featureseq_uniq),featureseq_uniq)
colnames(featureseq2)<-c("sequence","abundance");
featureseq2<-as.data.frame(featureseq2);
uniquesToFasta(featureseq2, fout='stool_rep-seqs.fna', ids=featureseq2$sequence);


#############################################################################################

############################# ANÁLISIS FARINGE ##############################################

#Calcular diferencias 
degs_pharynx = DESeq(degs_pharynx, test="Wald", fitType="parametric");
sigtab_pharynx = results(degs_pharynx, cooksCutoff = FALSE);

#Filtrar diferencias significativas
alpha = 0.01;
sigtab_pharynx2 = sigtab_pharynx[which(sigtab_pharynx$padj < alpha), ];
sigtab_pharynx2 = cbind(as(sigtab_pharynx2, "data.frame"), as(tax_table(phy.pharynx.mas1)[rownames(sigtab_pharynx2), ], "matrix"));
dim(sigtab_pharynx2)
#[1] 3 12 #otus*filas

#Exportar tabla para anexos
write.table(sigtab_pharynx2, "sigtab_pharynx.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE);

sigtab_pharynx_filtered <-sigtab_pharynx2

#Graficar
theme_set(theme_light());
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
};
# Phylum order
x = tapply(sigtab_pharynx_filtered$log2FoldChange, sigtab_pharynx_filtered$Phylum, function(x) max(x));
x = sort(x, TRUE);
sigtab_pharynx_filtered$Phylum = factor(as.character(sigtab_pharynx_filtered$Phylum), levels=names(x));

# Genus order
x = tapply(sigtab_pharynx_filtered$log2FoldChange, sigtab_pharynx_filtered$Genus, function(x) max(x));
x = sort(x, TRUE);
sigtab_pharynx_filtered$Genus = factor(as.character(sigtab_pharynx_filtered$Genus), levels=names(x));
ggplot(sigtab_pharynx_filtered, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5));


tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/Deseq_pharynx_caco.tiff", width = 20, height = 12, units = "cm", res = 600, pointsize = 10)
ggplot(sigtab_pharynx_filtered, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5));
dev.off()


tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/volcano_pharynx_caco.tiff", width = 12, height = 12, units = "cm", res = 600, pointsize = 10)
# Make a basic volcano plot
#with(sigtab_duo, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,10)))
with(sigtab_pharynx, plot(log2FoldChange, -log10(pvalue),  pch=21, main="Abundancia diferencial - faringe", cex =1.5, col="turquoise4" ))
# Add colored points: tomato if padj<0.01, yellow of log2FC>3, turquoise if both)
with(subset(sigtab_pharynx, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, cex =1,col= "tomato"))
#with(subset(sigtab_pharynx, abs(log2FoldChange)>=3), points(log2FoldChange, -log10(pvalue), pch=20, col= "yellow4"))
with(subset(sigtab_pharynx, padj<.01 & abs(log2FoldChange)>=3), points(log2FoldChange, -log10(pvalue), pch=20, cex=1, col="darkorchid1"))
dev.off()


#################################################
##ANALISIS DE RUTAS CON PICRUST
otu_pharynx<- otu_table(phy_pharynx, taxa_are_row=TRUE);
otu_pharynx_melt<-psmelt(otu_pharynx);

# Export taxonomy table as "tax.txt"
tax<-as(tax_table(phy_pharynx),"matrix");
tax_cols <- colnames(tax);
tax<-as.data.frame(tax);
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"));
for(co in tax_cols) tax[co]<-NULL;
write.table(tax, "tax_pharynx.txt", quote=FALSE, col.names=FALSE, sep="\t");

# Export feature/OTU table

# As a biom file

library(biomformat);packageVersion("biomformat");
## [1] ‘1.6.0’

#otu<-t(as(OTU_vst_duo,"matrix")); # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
otu<-t(as(otu_pharynx,"matrix"));
otu_biom<-make_biom(data=otu);
write_biom(otu_biom,"otu_pharynx_biom.biom");

# As a text file

write.table(t(otu_pharynx), "otu_pharynx.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
#or from the phyloseq object, 't' to transform if taxa_are_rows=FALSE, no 't' if taxa_are_rows=TRUE
#write.table(t(otu_table(ps), "seqtab.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Export metadata (if you have a properly formatted metadata file that you imported in your phyloseq pipeline, you can skip this step and just use that text file directly in QIIME 2)

write.table(sample_data(phy_pharynx),"sample-metadata_pharynx.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE);

#Para extraer las secuencias de las ASVs
library(dada2); packageVersion("dada2");
featureseq <-otu_pharynx_melt_vst %>% select(OTU,Abundance)
dim(featureseq)
#[1] 9882     2


colnames(featureseq)<-c("sequence","abundance");
featureseq_uniq <- getUniques(featureseq, collapse = TRUE, silence = FALSE)
featureseq2<-cbind(rownames(featureseq_uniq),featureseq_uniq)
colnames(featureseq2)<-c("sequence","abundance");
featureseq2<-as.data.frame(featureseq2);
uniquesToFasta(featureseq2, fout='pharynx_rep-seqs.fna', ids=featureseq2$sequence);

##############################################################################################3

