#!/bin/sh
#Run by chmod u+x sh/bamqc.sh
#sh/bamqc.sh
echo "Analysis Start"

#trim raw data
cd ~/Desktop/FYP/Mapped

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {1..6};
do mkdir ${dataDir[i-1]}/qc;
mkdir ${dataDir[i-1]}/qc/clipping;
echo ${dataDir[i-1]} "Clipping #################################";
clipping_profile.py -i ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam -o ${dataDir[i-1]}/qc/clipping/${dataDir[i-1]} -s "PE";
#mkdir ${dataDir[i-1]}/qc/coverage;
#echo ${dataDir[i-1]} "Coverage #################################";
#geneBody_coverage.py -i ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam -o ${dataDir[i-1]}/qc/coverage/${dataDir[i-1]} -r ag.bed;
mkdir ${dataDir[i-1]}/qc/insertion;
echo ${dataDir[i-1]} "Insertion #################################";
insertion_profile.py -i ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam -o ${dataDir[i-1]}/qc/insertion/${dataDir[i-1]} -s "PE";
mkdir ${dataDir[i-1]}/qc/gc;
echo ${dataDir[i-1]} "GC #################################";
read_GC.py -i ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam -o ${dataDir[i-1]}/qc/gc/${dataDir[i-1]};
echo ${dataDir[i-1]} "Stat #################################";
bam_stat.py -i ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam > ${dataDir[i-1]}/qc/${dataDir[i-1]}bam_stat.txt;
echo ${dataDir[i-1]} "Infer #################################";
infer_experiment.py -i ${dataDir[i-1]}/${dataDir[i-1]}sorted.bam -r ag.bed > ${dataDir[i-1]}/qc/${dataDir[i-1]}infer_experiment.txt;
done;

: '

bam_stat.py
clipping_profile.py
geneBody_coverage.py ** Runtime too long
infer_experiment.py
inner_distance.py**
insertion_profile.py
read_distribution.py**
read_GC.py
RNA_fragment_size.py**

do fastp -i raw_data/${dataDir[i-1]}/${dataDir[i-1]}_1.fq.gz -I raw_data/${dataDir[i-1]}/${dataDir[i-1]}_2.fq.gz -o trimmed/tri${dataDir[i-1]}_1.fq.gz -O trimmed/tri${dataDir[i-1]}_2.fq.gz -f 10 -F 10;
done;'