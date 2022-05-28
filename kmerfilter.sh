#!/bin/sh
#Run by chmod u+x sh/kmerfilter.sh
#sh/kmerfilter.sh
echo "Analysis Start"

#trim raw data
cd ~/Desktop/FYP/RNAseq

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {2..6};
do mkdir rcorr/${dataDir[i-1]}
echo ${dataDir[i-1]} "Rcorrector############################"
run_rcorrector.pl -1 trimmed/tri${dataDir[i-1]}_1.fq -2 trimmed/tri${dataDir[i-1]}_2.fq -od rcorr/${dataDir[i-1]}/ -t 14
cd ./rcorr/${dataDir[i-1]}
echo ${dataDir[i-1]} "Filering############################"
python ~/Desktop/FYP/RNAseq/FilterUncorrectabledPEfastq.py -1 tri${dataDir[i-1]}_1.cor.fq -2 tri${dataDir[i-1]}_2.cor.fq -s Aero_1
rm tri${dataDir[i-1]}_1.cor.fq
rm tri${dataDir[i-1]}_2.cor.fq
echo ${dataDir[i-1]} "Compress 1############################"
gzip unfixrm_tri${dataDir[i-1]}_1.cor.fq
echo ${dataDir[i-1]} "COmpress 2############################"
gzip unfixrm_tri${dataDir[i-1]}_2.cor.fq
cd ~/Desktop/FYP/RNAseq
done

#run_rcorrector.pl -1 trimmed/triAero_1_1.fq -2 trimmed/triAero_1_2.fq -od rcorr/
#python FilterUncorrectabledPEfastq.py -1 rcorr/triAero_1_1.cor.fq -2 rcorr/triAero_1_2.cor.fq -s Aero_1