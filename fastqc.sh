#!/bin/sh
#Run by chmod u+x sh/fastqc.sh
#sh/fastqc.sh
echo "Analysis Start"

#trim raw data
cd ~/Desktop/FYP/RNAseq/trimmed

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {1..6};
do fastqc tri${dataDir[i-1]}_1.fq;
fastqc tri${dataDir[i-1]}_2.fq;
done