# Evaluación del perfil taxonómico y funcional del microbioma de pacientes con Enfermedad Celíaca.

Estudio de casos y controles del microbioma en la Enfermedad Celíaca (EC) mediante el análisis de datos de secuenciación de alto rendimiento del gen 16S disponibles en bases de datos públicas. Se analizaron datos de nueve estudios que representaron 465 muestras (casos=168, controles = 192) tomadas de 4 partes diferentes del cuerpo (duodeno =120, saliva =33, heces =195 y faringe =54) para evaluar las diferencias en la diversidad alfa y beta, la abundancia diferencial a nivel de género entre casos y controles y los coeficientes de correlación de los taxa con EC. Se infirió el potencial metabólico del microbioma y se obtuvieron abundancias diferenciales de todos los  genes y vías predichas entre casos y controles, así como coeficientes de correlación con EC. Se hizo especial hincapié en las enzimas involucradas en la producción de ácidos grasos de cadena corta y las proteasas con actividad de endoprolil peptidasa. 

En este repositorio se encuentral los códigos utilizados para realizar el análisis.

##PASOS DEL ANÁLISIS##

***1. Detección de variantes biológicas con Dada2 y Asignación taxonómica con SILVA***

El primer paso en el análisis es filtrar las lecturas por calidad. Se deben buscar los parámetros óptimos explorando el filtrado de sus datos con diferentes parámetros. Hay dos opciones disponibles para este paso: cuando las lecturas están emparejadas, use primero cutadapt.h para eliminar los cebadores, luego flash.sh para fusionar los pares y finalmente filtering_SE.R para eliminar secuencias de mala calidad. Si las lecturas son de un solo extremo, haga lo mismo pero omita el paso flash.sh. 

Después de filtrar las lecturas, el segundo paso es inferir las variantes de secuencias y realizar la asignación taxonómica utilizando DADA2 y la base de datos SILVA. Para hacer esto, use el script Inf_Chim_tax.R. Sin embargo, si sus datos hacen que el análisis se quede sin memoria, utilice el script InfChim_tax_2.R en su lugar. El primer script hace la inferencia y asignacion taxonómica de manera simultánea para todas las muestras y el segundo la hace muestra por muestra, ahorrando memoria RAM pero haciendo el análisis mas lento.

***2. Análisis de diversidad con Phyloseq y Vegan***

Con los resultados anteriores es posible hacer un análisis de diversidad. El archivo Diversity.R contiene las instrucciones para fusionar los resultados de diferentes estudios analizados con las instrucciones anteriores, filtrar las muestras y los taxones de acuerdo con la información de las curvas de rarefacción y realizar la inferencia de diversidad alfa y beta y el análisis estadístico para comparar estos valores entre las muestras usando Shannon, Simpson, Chao1 y PCoA junto con Distancia ponderada de Unifrac. Este no es un script que se ejecute automáticamente porque se debe tener en cuenta el resultado del análisis en cada paso para decidir los parámetros del siguiente. Esto es solo una guía de los pasos en el análisis y debe ejecutarlo línea por línea. También contiene los códigos para hacer gráficos de buena calidad de los resultados usando ggplot2.

***3. Análisis de Abundancia diferencial con DESeq2 incluido en phyloseq***

Lo siguiente que se puede hacer después de analizar los cambios en la diversidad global de las muestras es enfocarse en taxones específicos que pueden estar asociados con la enfermedad o el tratamiento que se está estudiando. Eso podría evaluarse con un análisis de abundancia diferencial. Este archivo "DifferentialAbundance.R" contiene un conjunto de instrucciones para realizar el análisis de abundancia diferencial utilizando DESeq2. Al igual que el archivo anterior, este no es un script que se ejecute automáticamente, debe ejecutarlo línea por línea y revisar sus resultados en cada paso.

***4. Inferencia del potencial metabólico de las comunidadad microbianas con PICRUST2***

Finalmente, es posible utilizar la información taxonómica y de abundancia de taxones para inferir el potencial metabólico de las comunidades microbianas en las muestras. Este script "picrust2Qimme2.R" adapta los formatos de los archivos de salida obtenidos en el análisis anterior para importarlos a Qiime2 y ejecutar el análisis de inferencia metabólica usando el plugin PICRUST2 de Qiime2, este script también cambia el formatode  los archivos salida  obtenidos de PICRUST2 formato separado por tabulaciones que se puede leer en Excel o R.



# Assessing microbiome profiles of celiac disease patients

Case-control study of the microbiome in celiac disease through the analysis of 16S high-throughput sequencing data available in public databases. Data of nine studies accounting for 465 samples taken from 4 different parts of the body (duodenum, saliva, stool, and pharynx) was analyzed to assess the differences in alpha and beta diversity, differential abundance at genus level and correlation coefficients. The metabolic potential of the microbiome was inferred and differential abundance of genes and pathways, as well as correlation coefficients, were obtained. Special emphasis was stated on enzymes involved in Short-Chain Fatty Acids production and proteases with endoprolyl peptidase activity.

here you can find the code used to perform the analysis.

Biological variants detection with Dada2
Taxonomic Assignment with Dada2 and SILVA database
Diversity Analysis with Phyloseq and Vegan
Differential Abundance with DESeq2 included in phyloseq

The first step in the analysis is to filter the reads by quality. You should find the optimal parameters through exploring filtering your data with different parameters. Two options are available for this step: When the reads are paired-end use first cutadapt.h to remove the primers, then flash.sh to merge the pairs and finally filtering_SE.R. If the reads are single-end do the same but skipping the flash.sh step. After having the reads filtered the second step is to infer sequences variants and perform taxonomic assignment using DADA2 and SILVA database. To do this use the Inf_Chim_tax.R script. However, if your data make the analysis run out of memory use the InfChim_tax_2.R script instead.

Diversity.R file contains the instructions to merge the results from the different studies analyzed, filter the samples and taxa according to rarefaction curves information, and perform the alpha and beta diversity inference and statistical analysis of the samples using Shannon, Simpson, Chao1 and  PCoA coupled with Weighted Unifrac distance. This is not a script to run automatically because you have to be aware of the output of your analysis in each step to decide the parameters of the next one. This is just a guide of the steps in the analysis and you should run it line by line. It also contains the codes to make beautiful graphics of the results using ggplot2.

Nex thing you can do after analyzing changes in the global diversity of your samples is to focus in specific taxons that may be associated to the disease or treatment you are studying. That could be assessed with a differential abundance analysis. This file "DifferentialAbundance.R" contains a set of instructions to perform the differential abundance analysis using DESeq2. Like the previous file, this is not a script to run automatically, you should run it line by line and review your results in each step.

Finally, it is possible to use the taxonomic and abundance of taxa information to infer the metabolic potential of the microbial communities in the samples. This script "picrust2Qimme2" adapts the format files of the outputs obtained in the previous analysis to import them to Qiime2 and runt the metabolic inference analysis using the Qiime2 plugging PICRUST2, this script also changes the PICRUST2 output format of the obtained files to a tab-separated format that can be read in Excel or R.
