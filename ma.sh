#!/bin/sh
#Run by chmod u+x sh/ma.sh
#sh/ma.sh

echo "Analysis Start##################################"

#Map each pair
cd ~/Desktop/FYP

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

#for i in {1..6};
#do mkdir Mapped/${dataDir[i-1]};
#STAR --runThreadN 5 --genomeDir Indexed --runMode alignReads --outSAMtype BAM Unsorted --outFileNamePrefix Mapped/${dataDir[i-1]}/ --readFilesIn RNAseq/trimmed/tri${dataDir[i-1]}_1.fq RNAseq/trimmed/tri${dataDir[i-1]}_2.fq;
#done;

for i in {4..6};
do mkdir Mapped/${dataDir[i-1]}
#gunzip RNAseq/rcorr/${dataDir[i-1]}/unfixrm_tri${dataDir[i-1]}_1.cor.fq.gz -k
#gunzip RNAseq/rcorr/${dataDir[i-1]}/unfixrm_tri${dataDir[i-1]}_2.cor.fq.gz -k
STAR --runThreadN 5 --genomeDir Indexed10 --runMode alignReads --outSAMtype BAM Unsorted --outFileNamePrefix Mapped/${dataDir[i-1]}/ --readFilesIn RNAseq/rcorr/${dataDir[i-1]}/unfixrm_tri${dataDir[i-1]}_1.cor.fq RNAseq/rcorr/${dataDir[i-1]}/unfixrm_tri${dataDir[i-1]}_2.cor.fq
#rm RNAseq/rcorr/${dataDir[i-1]}/unfixrm_tri${dataDir[i-1]}_1.cor.fq
#rm RNAseq/rcorr/${dataDir[i-1]}/unfixrm_tri${dataDir[i-1]}_2.cor.fq
done;

# STAR --runThreadN 5 --genomeDir Indexed --runMode alignReads --outSAMtype BAM Unsorted --outFileNamePrefix Mapped/Aero_3/ --readFilesIn RNAseq/rcorr/Aero_3/unfixrm_triAero_3_1.cor.fq RNAseq/rcorr/Aero_3/unfixrm_triAero_3_2.cor.fq
# STAR --runThreadN 5 --genomeDir Indexed --runMode alignReads --outSAMtype BAM Unsorted --outFileNamePrefix Mapped/Ana_1/ --readFilesIn RNAseq/rcorr/Ana_1/unfixrm_triAna_1_1.cor.fq RNAseq/rcorr/Ana_1/unfixrm_triAna_1_2.cor.fq
