#!/bin/sh
#Run by chmod u+x sh/sort.sh
#sh/sort.sh
echo "Analysis Start"

#trim raw data
cd ~/Desktop/FYP/Mapped

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {1..6};
do samtools sort ${dataDir[i-1]}/Aligned.out.bam -o ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam;
samtools index ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam;
done