#Eliminar variables que hayan quedado abiertas de una sesión anterior
rm(list=ls());

#Cargar las librerias necesarias
library(phyloseq);library(ggplot2);library(vegan);library("dplyr");library(tidyverse);library(car);library(ggpubr);library(readxl);

# Restaurar los objetos phyloseq construidos a partir de los resultados de Dada2 para cada estudio
phy_ERR13<-readRDS(file = "phy_ERR13.rds");
phy_ERR15<-readRDS(file = "phy_ERR15.rds");
phy_ERR21<-readRDS(file = "phy_ERR21.rds");
phy_SRR11<-readRDS(file = "phy_SRR11.rds");
phy_SRR39<-readRDS(file = "phy_SRR39.rds");
phy_SRR52<-readRDS(file = "phy_SRR52.rds");
phy_SRR55<-readRDS(file = "phy_SRR55.rds");
phy_SRR60<-readRDS(file = "phy_SRR60.rds");
phy_SRR75<-readRDS(file = "phy_SRR75.rds");
phy_SRR3932<-readRDS(file = "phy_SRR3932.rds");

#Juntar los objetos en uno solo
merged <-merge_phyloseq(phy_ERR13,phy_ERR15,phy_ERR21,phy_SRR11, phy_SRR39,phy_SRR52,phy_SRR55,phy_SRR60,phy_SRR75,phy_SRR3932);
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 53037 taxa and 464 samples ]
#sample_data() Sample Data:       [ 464 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 53037 taxa by 6 taxonomic ranks ]

# Guardar el objeto a un archivo 
saveRDS(merged, file = "phy_all.rds");

merged_NA <-subset_taxa(merged, !is.na(Genus) & !Genus %in% c("", "NA"));

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 32884 taxa and 464 samples ]
#sample_data() Sample Data:       [ 464 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 32884 taxa by 6 taxonomic ranks ]

##################PREPROCESAMIENTO: FILTRADO DE DATOS ############################

# Juntar los taxa que tienen la misma asignación taxonómica a nivel de genero
dedup_tax <-tax_glom(merged, "Genus")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 756 taxa and 464 samples ]
#sample_data() Sample Data:       [ 464 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 756 taxa by 6 taxonomic ranks ]

#Importar la tabla con metadatos de las muestras
METADATA2 <-read_excel("METADATA.xlsx", sheet = "TableS1_buena");

#Definir como nombre las filas del dataframe a las variables de la primera columna; es decir, el nombre de las muestras
row.names(METADATA2) <-METADATA2$sample_name;

METADATA <-sample_data(METADATA2); #convertir el dataframen en un objeto phyloseq
row.names(METADATA) <-METADATA2$sample_name;
dedup_tax <-merge_phyloseq(otu_table(dedup_tax),tax_table(dedup_tax), METADATA);

#Explorar el numero de VSAs por muestra
samplesumsOTU <-sample_sums(dedup_tax);
samplesumsOTU;

#Filtrar muestras con menos de 2000 VSAs
samplesumsOTU_2000  <- prune_samples(sample_sums(OTU_dedup)>2000, dedup_tax)
samplesumsOTU_2000
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 756 taxa and 452 samples ]
#sample_data() Sample Data:       [ 452 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 756 taxa by 6 taxonomic ranks ]

#Filtrar VSAs que no esten presentes en almenos una muestra
abundant_taxa = prune_taxa(taxa_sums(samplesumsOTU_2000)>5,samplesumsOTU_2000);
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 630 taxa and 452 samples ]
#sample_data() Sample Data:       [ 452 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 630 taxa by 6 taxonomic ranks ]

#Remover las VSAs que no tengan asignación a nivel de genero "NA" (Ninguna)
phy_obj_bacteria_filtered_NA <- subset_taxa(abundant_taxa, !is.na(Genus) & !Genus %in% c("", "NA"));

#Exportar los objetos R unclean para diversidad y clean
saveRDS(abundant_taxa, file = "phy_unclean.rds");
saveRDS(phy_obj_bacteria_filtered_NA, file = "phy_clean.rds");

#Obtener la tabla de conteos de ASVs y de Asignación taxonómica como objetos phyloseq a partir del objeto filtrado
OTU_clean= otu_table(phy_obj_bacteria_filtered_NA, taxa_are_rows = FALSE);
TAX_clean= tax_table(phy_obj_bacteria_filtered_NA);

#Para extraer las secuencias de las ASVs

#Convertir el objeto phyloseq a un dataframe
OTU_clean_melt<-psmelt(OTU_clean);
OTU_clean_melt<-as.data.frame(OTU_clean_melt);

#Explorar
dim(OTU_clean_melt);
#[1]   292444      3
colnames(OTU_clean_melt);
#[1] "OTU"       "Sample"    "Abundance"

library(dada2);
OTU_clean_melt2<-OTU_clean_melt%>% select(OTU,Abundance);
colnames(OTU_clean_melt2)<-c("sequence","abundance");
OTU_clean_melt2 <-getUniques(OTU_clean_melt2, collapse = TRUE, silence = TRUE);
uniquesToFasta(OTU_clean_melt2, fout='celiac_rep-seqs.fna', ids=rownames(OTU_clean_melt2));


#construir un objeto phyloseq con las tablas filtradas y los metadatos
phy_obj2 <-phyloseq(OTU_clean,TAX_clean,METADATA);

library(ape);

#Construir un arbol para los taxa y adicionarlo al objeto phyloseq
random_tree <-rtree(ntaxa(phy_obj2), rooted=TRUE, tip.label=taxa_names(phy_obj2));
phy.tree <-merge_phyloseq(phy_obj2, random_tree);

#phyloseq-class experiment-level object
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 647 taxa and 452 samples ]
##sample_data() Sample Data:      [ 452 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 647 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 647 tips and 646 internal nodes ]

################## EVALUACION DE LA DIVERSIDAD ALFA Y BETA ############################
# Construir un dataframe que contenga la cantidad de secuencias obtenidas para cada muestra
sample_sum_df <- data.frame(sum = sample_sums(phy.tree));


#Normalizar la tabla de otus usando VST transformation

#Añadir pseudo conteo
phy.tree.mas1 = transform_sample_counts(phy.tree,function(x) x + 1);

#Transformar a objeto deseq2
library(DESeq2);
OTU.VST = phyloseq_to_deseq2(phy.tree.mas1, ~sampled_tissue); #OTUS #samples

dim(OTU.VST)
#[1] 647 452

#Abundancia normalizada con VsT
OTU.VST= estimateSizeFactors(OTU.VST);
OTU.VST= estimateDispersions(OTU.VST);
OTU.VST= getVarianceStabilizedData(OTU.VST);

#Convertir numeros negativos en 0
#setzero <- function (x)  {return(ifelse(x < 0, 0, x))}
#otu_vst<-setzero(OTU.VST);

#Convertir numeros negativos en positivos
min(OTU.VST);
#[1] -3.624219

posit <- function (x)  {return(x + 3.624219)};
otu_vst<-posit(OTU.VST);


#Construir phyloseq object con vst
otu_vst <- otu_table(otu_vst, taxa_are_row=TRUE);
phy_vst<-phyloseq(otu_vst,METADATA,TAX_clean);

#Construir un arbol para los taxa y adicionarlo al objeto phyloseq
vst_tree <- rtree(ntaxa(phy_vst), rooted=TRUE, tip.label=taxa_names(phy_vst));
phy_vst<-merge_phyloseq(phy_vst,vst_tree);


#Filtrar VSAs que no esten presentes en almenos una muestra
phy_vst =prune_taxa(taxa_sums(phy_vst)>0,phy_vst)

#Calcular la diversidad alfa mediante los indices de Sannon, Chao1 y Simpson
richness_vst <- estimate_richness(otu_vst, split = TRUE, measures = c("Shannon","Simpson"));
richness_vst ;

#The data you have provided does not have any singletons. This is highly suspicious. Results of richness estimates (for example) are probably unreliable, or wrong, if you have already trimmed low-abundance taxa from the data. No hay Chao1 no hay singletons, ojo pues.


richness <- estimate_richness(samplesumsOTU_2000, split = TRUE, measures = c("Chao1");
richness;

#Construir con dataframe con los metadatos y adicionando las columnas que contienen el calculo de indices de diversidad y el numero de secuencias por muestra
METADATA = sample_data(phy_vst);
METADATA$Shannon <- richness$Shannon;
METADATA$Simpson <- richness$Simpson;
METADATA$Chao1 <- richness$Chao1;
METADATA$Shannon_vst <- richness_vst$Shannon;
METADATA$Simpson_vst <- richness_vst$Simpson;
METADATA$sample_sum<-sample_sum_df$sum;
metadata<- data.frame(METADATA);



##########CURVAS DE RAREFACCION #########################################
#Graficar el histograma de la profundidad de secuenciación por muestra en el estudio
tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo tfm/TFM/phy_obj/histodepth_2.tiff", width = 12, height = 12, units = "cm", res = 600, pointsize = 10)
ggplot(metadata, aes(x =sample_sum)) + theme_classic() + geom_histogram(color = "tomato4" , fill="tomato1",binwidth = 15000) + ggtitle("Distribution of sample read counts") + xlab("read counts")
dev.off()

min(metadata$sample_sum)
#[1] 2109 SRR1107506  PRJNA231837
max(metadata$sample_sum)
#[1] 429693 SRR5514966 PRJNA385740

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo tfm/TFM/phy_obj/rarecurv.tiff", width = 12, height = 12, units = "cm", res= 600, pointsize = 10)
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink","brown","cyan","darkgrey");
lty <- c("solid", "dashed", "longdash", "dotdash");
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE);
out <- with(pars[1:40,], rarecurve(otu_table(phy_independent_rare), step = 20, col = col, lty = lty, label = FALSE))
dev.off()

out2 <- with(pars[1:40,], rarecurve(otu_table(phy_independent_rare), step = 20, col = col, lty = "solid", label = FALSE))
library(ranacapa)
rare1 <- ggrare(phy_independent_rare, step = 10, label = NULL, color = "sampled_tissue",  plot = TRUE, parallel = FALSE, se = TRUE)
rare1 <-rare1  + theme_classic();

rare2 <- ggrare(phy_independent_rare, step = 10, label = NULL, color = "study_accession",  plot = TRUE, parallel = FALSE, se = TRUE) + theme_classic();

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo tfm/TFM/phy_obj/rare_study.tiff", width = 12, height = 12, units = "cm", res = 600, pointsize = 10)
rare2 <-rare2 + theme_classic();
rare2
dev.off()




##########################DIVERSIDAD ALPHA###################################

#Explorar el número de muestras que hay por tejido muestreado
table(metadata$sampled_tissue)
#duodenum  pharynx   saliva    stool 
#     135       56       81      180 


#Para comparar índices entre  casos y controles  sin discriminar partes del cuerpo

met_case <-filter(metadata, case_control == "case"); #223 casos
met_control <-filter(metadata, case_control == "control");#229 controles



norm <- function (x)  {return(sqrt(asin(x)))};
shanon_norm<-norm(shannon_vst);


# Prueba de normalidad : Shapiro-Wilk normality test

shannon.test <- shapiro.test(METADATA$Shannon_vst);
print(shannon.test)
#W = 0.89781, p-value < 2.2e-16

#data:  METADATA$Shannon
#W = 0.9696, p-value = 4.5e-08

Simpson.test <- shapiro.test(METADATA$Simpson_vst);
print(Simpson.test)
#W = 0.45249, p-value < 2.2e-16

#data:  METADATA$Simpson
#W = 0.79618, p-value < 2.2e-16


Chao1.test <- shapiro.test(METADATA$Chao1);
print(Chao1.test)
#W = 0.96933, p-value = 3.992e-08


#Los datos no se ajustan a una distribución normal por lo tanto deben compararse con pruebas no paramètricas

#Homocedasticidad: la varianza debe de ser constante entre todos los grupos.
#Levene's Test for Homogeneity of Variance (center = "median")
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


leveneTest(y = METADATA$Chao1, group = METADATA$case_control, center ="median");
#       Df   F value    Pr(>F)    
#group   1  7.1723 0.007675 **
#      450 

leveneTest(y = METADATA$Simpson_vst, group = METADATA$case_control, center ="median");
#       Df  F value Pr(>F)
#group  1  3.8942 0.04906 *
#      450               

leveneTest(y = METADATA$Shannon_vst, group = METADATA$case_control, center ="median");
#       Df F value Pr(>F)
#group 1  8.3304 0.004086 **
#      450
#################### INDICE DE CHAO1 #########################################################################
# Comparación del indice Chao1 entre tejidos muestreados entre casos y controles

#Wilcox test

compare_means(Chao1 ~ case_control, data = metadata, method ="wilcox.test", paired = FALSE,  group.by = "sampled_tissue")
# sampled_tissue .y.   group1  group2       p p.adj p.format p.signif method    
# stool          Chao1 control case   0.126   0.25  0.1257   ns       Wilcoxon
# duodenum       Chao1 control case   0.00812 0.032 0.0081   **       Wilcoxon
# saliva         Chao1 control case   0.487   0.49  0.4872   ns       Wilcoxon
# pharynx        Chao1 control case   0.0181  0.054 0.0181   *        Wilcoxon

#Calcular n de cada grupo comparado
table(metadata$sampled_tissue)
dim(filter(metadata, sampled_tissue =="duodenum", case_control =="case"))
dim(filter(metadata, sampled_tissue =="saliva", case_control =="case"))
#     duodenum pharynx  saliva   stool 
#	 135     56       81      180
#case     61     36	  40      86
#control  74	 20	  41	  94

#  Graficar el índice
pChao_caco <- ggplot(METADATA, aes(sampled_tissue, Chao1))+ theme_classic() + geom_boxplot(aes(fill=case_control),outlier.shape= 1, outlier.color="black") + stat_compare_means(aes(group = case_control),label.y = 112, label = "p.signif") +theme(axis.text=element_text(size=11), axis.title=element_text(size=12,face="bold"))	 

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/Chao1_tissue_caco.tiff", width = 20, height = 12, units = "cm", res = 600, pointsize = 10)
pChao_caco 
dev.off()
#################### INDICE DE SIMPSON #########################################################################
compare_means(Simpson_vst ~ case_control, data = metadata, method ="wilcox.test", paired = FALSE,  group.by = "sampled_tissue")
#  sampled_tissue .y.       group1 group2       p p.adj p.format p.signif method 
# stool          Simpson_… contr… case   0.156   0.31  0.1559   ns       Wilcox…
# duodenum       Simpson_… contr… case   0.00517 0.021 0.0052   **       Wilcox…
# saliva         Simpson_… contr… case   0.461   0.46  0.4612   ns       Wilcox…
# pharynx        Simpson_… contr… case   0.0205  0.062 0.0205   *        Wilcox…

#Graficar
pSimpson_caco <- ggplot(METADATA, aes(sampled_tissue, Simpson_vst))+ theme_classic(base_size = 12) + geom_boxplot(aes(fill=case_control),outlier.shape= 1, outlier.color="black") + stat_compare_means(aes(group = case_control),label.y = 1, label = "p.signif")+theme(axis.text=element_text(size=11), axis.title=element_text(size=12,face="bold"))	

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/Simpson_tissue_caco.tiff", width = 20, height = 12, units = "cm", res = 600, pointsize = 10)
pSimpson_caco  
dev.off()

#################### INDICE DE SHANNON ###########################################
compare_means(Shannon_vst ~ case_control, data = metadata, method ="wilcox.test", paired = FALSE,  group.by = "sampled_tissue")
#  sampled_tissue .y.       group1  group2      p p.adj p.format p.signif method 
# stool          Shannon_… control case   0.0688  0.21 0.069    ns       Wilcox…
# duodenum       Shannon_… control case   0.0254  0.1  0.025    *        Wilcox…
# saliva         Shannon_… control case   0.183   0.37 0.183    ns       Wilcox…
# pharynx        Shannon_… control case   0.503   0.5  0.503    ns       Wilcox…





#Graficar
pShannon_caco <- ggplot(METADATA, aes(sampled_tissue, Shannon_vst))+ theme_classic() + geom_boxplot(aes(fill=case_control), outlier.shape= 1, outlier.color="black") + stat_compare_means(aes(group = case_control),label.y = 5.2, label = "p.signif") +theme(axis.text=element_text(size=11), axis.title=element_text(size=12,face="bold"))


#ggarrange(pShannon_caco, pSimpson_caco, pChao_caco + rremove("x.text"), labels = c("A", "B", "C"), ncol = 3, nrow = 1)

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/Shannon_tissue_caco.tiff", width = 20, height = 12, units = "cm", res = 600, pointsize = 10)
pShannon_caco
dev.off()
###############################DIVERSIDAD BETHA################################################################
#############################FILTRAR LAS MUESTRAS DEPENDIENTES#######################################################
#Algunos individuos fueron muestreados en mas de un tejido, se conservará solo una de las muestras por individuo para conservar la independencia de las muestras antes de hacer los análisis de ordenación

To_remove1 <-filter(metadata,study_accession=="PRJNA385740" & sampled_tissue=="stool");
To_remove2 <-filter(metadata,study_accession=="PRJNA401920" & sampled_tissue=="stool");
To_remove3 <-filter(metadata,study_accession=="PRJNA371697" & sampled_tissue=="duodenum");

'%notin%' <- Negate('%in%')

metadata_independent <- filter(metadata, sample_name %notin% To_remove1$sample_name);
metadata_independent <- filter(metadata_independent, sample_name %notin% To_remove2$sample_name);
metadata_independent <- filter(metadata_independent, sample_name %notin% To_remove3$sample_name);
phy_independent_vst <- subset_samples(phy_vst,sample_name %in% metadata_independent$sample_name);

#356(muestras independientes)
##########################################

library(GUniFrac); #Para alpha diversity deje todos los que estuvieran almenos en una muestra y para betha diversity quité los que no tuvieran almenos 10 secuencias)

ord_independen_unifrac_vst = ordinate(phy_independent_vst, method="MDS", distance="unifrac",weighted=TRUE);

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/betha_div_bias_circle.tiff", width = 19, height = 12, units = "cm", res = 600, pointsize = 5)

a<-plot_ordination(phy_independent_vst, ord_independen_unifrac_vst, color="sequencing_technology", title = "Tecnología de secuenciación", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1)+ theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12), plot.title = element_text(hjust = 0.5))

b<-plot_ordination(phy_independent_vst, ord_independen_unifrac_vst, color="X16S_region", title = "Región 16S secuenciada",axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1) + theme(legend.position="bottom") +  theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12), plot.title = element_text(hjust = 0.5))

ggarrange(b,a, labels = c("A", "B"), ncol = 2, nrow = 1)

dev.off()


tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/betha_div_tissue.tiff", width = 19, height = 12, units = "cm", res = 600, pointsize = 1)

b <- plot_ordination(phy_independent_vst, ord_independen_unifrac_vst, color="sampled_tissue", title = "Tejido muestreado", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1)+ theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12), plot.title = element_text(hjust = 0.5))

a<-plot_ordination(phy_independent_vst, ord_independen_unifrac_vst, color="sequencing_technology", title = "Tecnología de secuenciación", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1)+ theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12), plot.title = element_text(hjust = 0.5))

ggarrange(b,a, labels = c("A", "B"), ncol = 2, nrow = 1)

dev.off()

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/betha_div_CeD.tiff", width = 19, height = 12, units = "cm", res = 600, pointsize = 1)


a <-plot_ordination(phy_independent_vst, ord_independen_unifrac_vst, color="case_control",title = "Caso vs. control", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1)+ theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12), plot.title = element_text(hjust = 0.5))

b<-plot_ordination(phy_independent_vst, ord_independen_unifrac_vst, color="feeding_habit", title = "Dieta con o sin gluten", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1) + theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12), plot.title = element_text(hjust = 0.5))

ggarrange(b,a, labels = c("A", "B"), ncol = 2, nrow = 1)

dev.off()

############################# ANÁLISIS DUODENO ##############################################
#Objeto phyloseq solo muestras duodeno
phy_duodenum <- subset_samples(phy.tree, sampled_tissue =="duodenum");
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 647 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 647 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 647 tips and 646 internal nodes ]

#Filtrar VSAs que no esten presentes en almenos una muestra
phy_duodenum = prune_taxa(taxa_sums(phy_duodenum)>10,phy_duodenum);
phy_duodenum;

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 518 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 518 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 518 tips and 517 internal nodes ]

library(DESeq2);
#Añadir pseudo conteo
phy.duo.mas1 = transform_sample_counts(phy_duodenum,function(x) x + 1);

#Transformar a objeto deseq2
degs_duo = phyloseq_to_deseq2(phy.duo.mas1, ~case_control); #OTUS #samples
dim(degs_duo)
#[1] 518 135

degs_duo$case_control <- relevel( degs_duo$case_control, "control" )
#Abundancia normalizada con VsT
degs_duo= estimateSizeFactors(degs_duo);
degs_duo = estimateDispersions(degs_duo);
diagvst_duo = getVarianceStabilizedData(degs_duo);

otu_duo_vst <- otu_table(diagvst_duo, taxa_are_row=TRUE);
tax_duo= tax_table(phy_duodenum);
met_duodenum = sample_data(phy_duodenum);
phy_duodenum_vst<-phyloseq(otu_duo_vst,tax_duo,met_duodenum);
tree_duo_vst <- rtree(ntaxa(phy_duodenum_vst), rooted=TRUE, tip.label=taxa_names(phy_duodenum_vst));
write.tree(tree_duo_vst, file = "duo_tree.nwk");
phy_duodenum_vst<-merge_phyloseq(phy_duodenum_vst,tree_duo_vst);
#otu_table()   OTU Table:         [ 518 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 518 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 518 tips and 517 internal nodes ]
 

min(otu_duo_vst);
#[1] -4.145409
phy_duodenum_vst = transform_sample_counts(phy_duodenum_vst,function(x) x +4.145409);

ord_duo_unifrac_vst = ordinate(phy_duodenum_vst, method="MDS", distance="unifrac",weighted=TRUE);

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/betha_div_duodeno_2.tiff", width = 19, height = 12, units = "cm", res = 600, pointsize = 1)

a <-plot_ordination(phy_duodenum_vst, ord_duo_unifrac_vst, color ="case_control", title = "Caso vs. control", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1) + theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12),plot.title = element_text(hjust = 0.5))


b <-plot_ordination(phy_duodenum_vst, ord_duo_unifrac_vst, title = "Nº ID. del estudio", color ="study_accession", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1)+  theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12),plot.title = element_text(hjust = 0.5))

ggarrange(a,b, labels = c("A", "B"), ncol = 2, nrow = 1)
dev.off()




############################# ANÁLISIS SALIVA ##############################################
#Objeto phylosec solo muestras saliva
phy_saliva <- subset_samples(phy.tree, sampled_tissue =="saliva");
phy_saliva;
#otu_table()   OTU Table:         [ 647 taxa and 81 samples ]
#sample_data() Sample Data:       [ 81 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 647 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 647 tips and 646 internal nodes ]

#Filtrar VSAs que no esten presentes en almenos 10 lecturas
phy_saliva = prune_taxa(taxa_sums(phy_saliva)>10,phy_saliva);
phy_saliva;
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 105 taxa and 81 samples ]
#sample_data() Sample Data:       [ 81 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 105 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 105 tips and 104 internal nodes ]


library(DESeq2);
#Añadir pseudo conteo
phy.saliva.mas1 = transform_sample_counts(phy_saliva, function(x) x + 1);

#Transformar a objeto deseq2
degs_sal = phyloseq_to_deseq2(phy.saliva.mas1, ~case_control); #OTUS #samples
dim(degs_sal);
#105 81
degs_sal$case_control <- relevel( degs_sal$case_control, "control");

#Abundancia normalizada con VsT
degs_sal= estimateSizeFactors(degs_sal);
degs_sal = estimateDispersions(degs_sal);
diagvst_sal = getVarianceStabilizedData(degs_sal);

#Construir phyloseq object con datos normalizados
otu_sal_vst <- otu_table(diagvst_sal, taxa_are_row=TRUE);
tax_sal= tax_table(phy_saliva);
met_saliva = sample_data(phy_saliva);
phy_sal_vst<-phyloseq(otu_sal_vst,tax_sal,met_saliva);

library(ape);
tree_sal_vst <- rtree(ntaxa(phy_sal_vst), rooted=TRUE, tip.label=taxa_names(phy_sal_vst));
write.tree(tree_sal_vst, file = "saltree.nwk");
phy_sal_vst<-merge_phyloseq(phy_sal_vst,tree_sal_vst);
phy_sal_vst;

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 105 taxa and 81 samples ]
#sample_data() Sample Data:       [ 81 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 105 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 105 tips and 104 internal nodes ]

min(otu_sal_vst);
#-2.47849

phy_sal_vst = transform_sample_counts(phy_sal_vst,function(x) x + 2.47849);

ord_sal_unifrac_vst = ordinate(phy_sal_vst, method="MDS", distance="unifrac",weighted=TRUE);

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/betha_div_saliva.tiff", width = 19, height = 12, units = "cm", res = 600, pointsize = 1)


a <-plot_ordination(phy_sal_vst, ord_sal_unifrac_vst, color ="case_control", title = "Caso vs. control", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1) + theme(legend.position="bottom") +  theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12),plot.title = element_text(hjust = 0.5))


b <-plot_ordination(phy_sal_vst, ord_sal_unifrac_vst, title = "Nº ID. del estudio", color ="study_accession", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1) +  theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12),plot.title = element_text(hjust = 0.5))

ggarrange(a,b, labels = c("A", "B"), ncol = 2, nrow = 1)

dev.off()

############################# ANÁLISIS FARINGE ##############################################

#Objeto phylosec solo muestras pharinx
phy_pharynx <- subset_samples(phy.tree, sampled_tissue =="pharynx");
phy_pharynx;
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 647 taxa and 56 samples ]
#sample_data() Sample Data:       [ 56 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 647 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 647 tips and 646 internal nodes ]

#Filtrar VSAs que no esten presentes en almenos 10 lecturas
phy_pharynx = prune_taxa(taxa_sums(phy_pharynx)>10,phy_pharynx);
phy_pharynx;

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 71 taxa and 56 samples ]
#sample_data() Sample Data:       [ 56 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 71 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 71 tips and 70 internal nodes ]


library(DESeq2);
#Añadir pseudo conteo
phy.pharynx.mas1 = transform_sample_counts(phy_pharynx, function(x) x + 1);

#Transformar a objeto deseq2
degs_pharynx = phyloseq_to_deseq2(phy.pharynx.mas1, ~case_control); #OTUS #samples
dim(degs_pharynx);
#[1] 71 56

degs_sal$case_control <- relevel( degs_sal$case_control, "control");


#Abundancia normalizada con VsT
degs_pharynx= estimateSizeFactors(degs_pharynx);
degs_pharynx= estimateDispersions(degs_pharynx);
diagvst_pharynx = getVarianceStabilizedData(degs_pharynx);


#Construir phyloseq object con datos normalizados
otu_pharynx_vst <- otu_table(diagvst_pharynx, taxa_are_row=TRUE);
tax_pharynx= tax_table(phy_pharynx);
met_pharynx = sample_data(phy_pharynx);
phy_pharynx_vst<-phyloseq(otu_pharynx_vst,tax_pharynx,met_pharynx);

min(otu_pharynx_vst);
# 2.62753


#No aplica porque min es positivo
#phy_pharynx_vst = transform_sample_counts(phy_pharynx_vst,function(x) x + abs(min(otu_pharynx_vst)));

library(ape);
tree_pharynx_vst <- rtree(ntaxa(phy_pharynx_vst), rooted=TRUE, tip.label=taxa_names(phy_pharynx_vst));
phy_pharynx_vst<-merge_phyloseq(phy_pharynx_vst,tree_pharynx_vst);
phy_pharynx_vst;

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 71 taxa and 56 samples ]
#sample_data() Sample Data:       [ 56 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 71 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 71 tips and 70 internal nodes ]


ord_phar_unifrac_vst = ordinate(phy_pharynx_vst, method="MDS", distance="unifrac",weighted=TRUE);

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/betha_div_pharynx.tiff", width = 19, height = 12, units = "cm", res = 600, pointsize = 1)

a <- plot_ordination(phy_pharynx_vst, ord_phar_unifrac_vst, color="case_control", title = "Caso vs. control", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1, alpha=0.75)+  theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12), plot.title = element_text(hjust = 0.5))

b <- plot_ordination(phy_pharynx_vst, ord_phar_unifrac_vst, color="case_control", title = "Unifrac MDS Analysis", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1, alpha=0.75) +  theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12), plot.title = element_text(hjust = 0.5))

ggarrange(a,b, labels = c("A", "B"), ncol = 2, nrow = 1)

dev.off()

############################# ANÁLISIS STOOL ##############################################
#Objeto phylosec solo muestras heces
phy_stool <- subset_samples(phy.tree, sampled_tissue =="stool");
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 647 taxa and 195 samples ]
#sample_data() Sample Data:       [ 195 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 647 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 647 tips and 646 internal nodes ]

#Filtrar VSAs que no esten presentes en almenos 10 lecturas
phy_stool = prune_taxa(taxa_sums(phy_stool)>10,phy_stool);
phy_stool;

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 419 taxa and 195 samples ]
#sample_data() Sample Data:       [ 195 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 419 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 419 tips and 418 internal nodes ]


library(DESeq2);
#Añadir pseudo conteo
phy.stool.mas1 = transform_sample_counts(phy_stool, function(x) x + 1);

#Transformar a objeto deseq2
degs_stool = phyloseq_to_deseq2(phy.stool.mas1, ~case_control); #OTUS #samples
#[1] 338 180
degs_stool$case_control <- relevel( degs_stool$case_control, "control");

#Abundancia normalizada con VsT
degs_stool= estimateSizeFactors(degs_stool);
degs_stool = estimateDispersions(degs_stool);
diagvst_stool = getVarianceStabilizedData(degs_stool);

#Construir phyloseq object con datos normalizados

otu_stool_vst <- otu_table(diagvst_stool, taxa_are_row=TRUE);
tax_stool= tax_table(phy_stool);
met_stool = sample_data(phy_stool);
phy_stool_vst<-phyloseq(otu_stool_vst,tax_stool,met_stool);


min(otu_stool_vst);
#-4.862319

phy_stool_vst = transform_sample_counts(phy_stool_vst,function(x) x + abs(min(otu_stool_vst)));

library(ape);
tree_stool_vst <- rtree(ntaxa(phy_stool_vst), rooted=TRUE, tip.label=taxa_names(phy_stool_vst));
phy_stool_vst<-merge_phyloseq(phy_stool_vst,tree_stool_vst);
phy_stool_vst;

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 419 taxa and 195 samples ]
#sample_data() Sample Data:       [ 195 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 419 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 419 tips and 418 internal nodes ]


ord_stool_unifrac_vst = ordinate(phy_stool_vst, method="MDS", distance="unifrac",weighted=TRUE);

tiff("/media/juliana/Ubuntu/TFM_ubuntu/nuevo_tfm/TFM/phy_obj/betha_div_stool_2.tiff", width = 19, height = 12, units = "cm", res = 600, pointsize = 1)

a <-plot_ordination(phy_stool_vst, ord_stool_unifrac_vst, color ="case_control", title = "Caso vs. control", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1) +  theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12),plot.title = element_text(hjust = 0.5))

b <-plot_ordination(phy_stool_vst, ord_stool_unifrac_vst, title = "Nº ID. del estudio", color ="study_accession", axes=c(1,2)) + theme_classic() + geom_point(size=3, shape=1) +  theme(legend.position="bottom", legend.title = element_blank(), axis.text=element_text(size=11), axis.title=element_text(size=12),plot.title = element_text(hjust = 0.5))

ggarrange(a,b, labels = c("A", "B"), ncol = 2, nrow = 1)

dev.off()




