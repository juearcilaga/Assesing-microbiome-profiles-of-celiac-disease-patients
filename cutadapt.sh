#!/bin/bashsource conda
source /home/lmarcos/miniconda3/etc/profile.d/conda.sh
conda activate

ls *.fastq.gz|while read line;
do
I=$(basename "${line}" .fastq.gz);
echo $I
results=$(cutadapt -m 100 -O 15 -q 20,20 -j 8 -g "GTGCCAGCMGCCGCGG"  -a "ACTYAAAKGAATTGACGG" -o ${I}_primers.fastq.gz  ${I}.fastq.gz >>cutadapt);
done

