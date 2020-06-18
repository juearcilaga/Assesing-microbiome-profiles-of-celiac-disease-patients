#!/bin/bash
source conda
source /home/lmarcos/miniconda3/etc/profile.d/conda.sh #CHANGE to your conda source path
conda activate

ls *_1.fastq.gz|while read line;
do
I=$(basename "${line}" _1.fastq.gz);
echo $I
results=$(flash ${I}_1.fastq.gz ${I}_2.fastq.gz -o ${I}_corrida.flash -M 300 --allow-outies -t 8 -z);
done
