#!/bin/sh
#Run by chmod u+x sh/cn.sh
#sh/cn.sh
echo "Analysis Start"

#trim raw data
cd ~/Desktop/FYP/Mapped

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {1..6};
do mv ${dataDir[i-1]}/sorted.bam ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam;
done;