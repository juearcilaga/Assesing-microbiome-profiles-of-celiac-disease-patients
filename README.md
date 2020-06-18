# Assessing microbiome profiles of celiac disease patients

Case-control study of the microbiome in celiac disease through the analysis of 16S high-throughput sequencing data available in public databases. Data of nine studies accounting for 465 samples taken from 4 different parts of the body (duodenum, saliva, stool, and pharynx) was analyzed to assess the differences in alpha and beta diversity, differential abundance at genus level and correlation coefficients. The metabolic potential of the microbiome was inferred and differential abundance of genes and pathways, as well as correlation coefficients, were obtained. Special emphasis was stated on enzymes involved in Short-Chain Fatty Acids production and proteases with endoprolyl peptidase activity.

here you can find the code used to perform the analysis.

Biological variants detection with Dada2
Taxonomic Assignment with Dada2 and SILVA database
Diversity Analysis with Phyloseq and Vegan
Differential Abundance with DESeq2 included in phyloseq

The first step in the analysis is to filter the reads by quality. You should find the optimal parameters through exploring filtering your data with different parameters. Two options are available for this step: When the reads are paired-end use first cutadapt.h to remove the primers, then flash.sh to merge the pairs and finally filtering_SE.R. If the reads are single-end do the same but skipping the flash.sh step. After having the reads filtered the second step is to infer sequences variants and perform taxonomic assignment using DADA2 and SILVA database. To do this use the Inf_Chim_tax.R script. However, if your data make the analysis run out of memory use the InfChim_tax_2.R script instead.

# Evaluación del perfil taxonómico y funcional del microbioma de pacientes con Enfermedad Celíaca.

Estudio de casos y controles del microbioma en la Enfermedad Celíaca (EC) mediante el análisis de datos de secuenciación de alto rendimiento del gen 16S disponibles en bases de datos públicas. Se analizaron datos de nueve estudios que representaron 465 muestras (casos=168, controles = 192) tomadas de 4 partes diferentes del cuerpo (duodeno =120, saliva =33, heces =195 y faringe =54) para evaluar las diferencias en la diversidad alfa y beta, la abundancia diferencial a nivel de género entre casos y controles y los coeficientes de correlación de los taxa con EC. Se infirió el potencial metabólico del microbioma y se obtuvieron abundancias diferenciales de todos los  genes y vías predichas entre casos y controles, así como coeficientes de correlación con EC. Se hizo especial hincapié en las enzimas involucradas en la producción de ácidos grasos de cadena corta y las proteasas con actividad de endoprolil peptidasa. 

En este repositorio se encuentral los códigos utilizados para realizar el análisis.

1. Detección de variantes biológicas con Dada2
2. Asignación taxonómica con la base de datos Dada2 y SILVA
3. Análisis de diversidad con Phyloseq y Vegan
4. Abundancia diferencial con DESeq2 incluido en phyloseq

El primer paso en el análisis es filtrar las lecturas por calidad. Se deben buscar los parámetros óptimos explorando el filtrado de sus datos con diferentes parámetros. Hay dos opciones disponibles para este paso: cuando las lecturas están emparejadas, use primero cutadapt.h para eliminar los cebadores, luego flash.sh para fusionar los pares y finalmente filtering_SE.R para eliminar secuencias de mala calidad. Si las lecturas son de un solo extremo, haga lo mismo pero omita el paso flash.sh. 

Después de filtrar las lecturas, el segundo paso es inferir las variantes de secuencias y realizar la asignación taxonómica utilizando DADA2 y la base de datos SILVA. Para hacer esto, use el script Inf_Chim_tax.R. Sin embargo, si sus datos hacen que el análisis se quede sin memoria, utilice el script InfChim_tax_2.R en su lugar. El primer script hace la inferencia y asignacion taxonómica de manera simultánea para todas las muestras y el segundo la hace muestra por muestra, ahorrando memoria RAM pero haciendo el análisis mas lento.
