source /home/lmarcos/miniconda3/etc/profile.d/conda.sh 
#In QIIME 2 (qiime2-2018.4)
conda activate qiime2-2019.10 

biom head -i otu_duodeno_biom.biom 
sed -i 's/OTU ID/OTU_ID/g' otu_duodeno.biom 
biom head -i otu_duodeno.biom 

#Import feature table from exported biom

qiime tools import \
--input-path otu_duodeno_biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path duo_feature-table.qza

qiime feature-table summarize --i-table duo_feature-table.qza    --o-visualization duo_feature-table-summary

#Import the taxonomy table:
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path tax_duo.txt \
--output-path duo_taxonomy_featuredatatax.qza

#Imported duo_rep-seqs.fna as DNASequencesDirectoryFormat to duo_rep-seqs.qza
qiime tools import \
--input-path duo_rep-seqs.fna \
--type 'FeatureData[Sequence]' \
--output-path duo_rep-seqs.qza


#Run PICRUST2
qiime picrust2 full-pipeline --i-table duo_feature-table.qza\
 --i-seq duo_rep-seqs.qza \
 --p-threads 6 \
 --p-highly-verbose \
 --o-ko-metagenome ko-metagenome_duo \
 --o-ec-metagenome ec-metagenome_duo \
 --o-pathway-abundance pathway-abundance_duo 



qiime feature-table summarize \
   --i-table pathway-abundance_sal.qza \
   --o-visualization pathway_abundance_sal.qzv
   

# FeatureTable[Frequency] artifact will be exported as a BIOM v2.1.0 formatted file.

qiime tools export \
     --input-path ec-metagenome_duo.qza \
     --output-path ec-metagenome_duo.biom

qiime tools export \
     --input-path ko-metagenome_duo.qza \
     --output-path ko-metagenome_duo.biom

qiime tools export \
     --input-path pathway-abundance_duo.qza \
     --output-path pathway-abundance_duo.biom
	 
	 
#Change files format from biom format to tsv
biom convert -i ec-metagenome_duo.biom -o ec-metagenome_duo.tsv --to-tsv
biom convert -i ko-metagenome_duo.biom -o ko-metagenome_duo.tsv --to-tsv
biom convert -i pathway-abundance_duo.biom -o pathway-abundance_duo.tsv --to-tsv

#Add genes names to ko terms and pathways names to METACYC codes 
#Make sure the file do not have a header line  # Constructed from biom file
add_descriptions.py -i ec-metagenome_duo.tsv -m EC  -o ec-metagenome_duo_descrip.tsv;
add_descriptions.py -i ko-metagenome_duo.tsv -m KO  -o ko-metagenome_duo_descrip.tsv;
add_descriptions.py -i pathway-abundance_duo.tsv -m METACYC -o pathway-abundance_duo_descrip.tsv;
