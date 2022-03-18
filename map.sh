#!/bin/sh
#Run by chmod u+x sh/map.sh
#sh/map.sh

echo "Analysis Start##################################"

#Map each pair
cd ~/Desktop/FYP

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {4..6};
do mkdir Mapped/${dataDir[i-1]};
STAR --runThreadN 5 --genomeDir Indexed --runMode alignReads --outSAMtype BAM Unsorted --outFileNamePrefix Mapped/${dataDir[i-1]}/ --readFilesIn RNAseq/trimmed/tri${dataDir[i-1]}_1.fq RNAseq/trimmed/tri${dataDir[i-1]}_2.fq;
done;


