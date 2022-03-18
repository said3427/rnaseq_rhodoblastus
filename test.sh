#!/bin/sh
#Run by chmod u+x sh/test.sh
#sh/test.sh
echo "Analysis Start"

#trim raw data
cd ~/Desktop/FYP/RNAseq

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {1..6};
do fastp -i raw_data/${dataDir[i-1]}/${dataDir[i-1]}_1.fq.gz -I raw_data/${dataDir[i-1]}/${dataDir[i-1]}_2.fq.gz -o trimmed/tri${dataDir[i-1]}_1.fq.gz -O trimmed/tri${dataDir[i-1]}_2.fq.gz -f 10 -F 10;
done