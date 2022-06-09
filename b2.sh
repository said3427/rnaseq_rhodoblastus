#!/bin/sh
#Run by chmod u+x sh/b2.sh
#sh/b2.sh

#Runed bowtie2-build Annotation/newAnnotation/final/updated.fasta B2/index

echo "Analysis Start"

cd ~/Desktop/FYP

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {1..6};
do echo "${dataDir[i-1]} ###############################"
cd ~/Desktop/FYP
bowtie2 -x B2/index -1 RNAseq/rcorr/${dataDir[i-1]}/unfixrm_tri${dataDir[i-1]}_1.cor.fq -2 RNAseq/rcorr/${dataDir[i-1]}/unfixrm_tri${dataDir[i-1]}_2.cor.fq -S B2/${dataDir[i-1]}.sam -p 15 --rf
cd ./B2
samtools view -bS ${dataDir[i-1]}.sam > ${dataDir[i-1]}.bam
rm ${dataDir[i-1]}.sam
samtools sort ${dataDir[i-1]}.bam -o ${dataDir[i-1]}.bam
samtools index ${dataDir[i-1]}.bam
htseq-count -r pos -t CDS -i ID -n 7 -s reverse ${dataDir[i-1]}.bam updated.gff -c ${dataDir[i-1]}count.csv
done




