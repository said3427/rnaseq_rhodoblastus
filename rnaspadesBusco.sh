#!/bin/sh
#Run by chmod u+x sh/test.sh
#sh/test.sh
echo "Analysis Start"

#trim raw data
cd ~/Desktop/FYP/RNAseq

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {1..6};
do fastp -i raw_data/${dataDir[i-1]}/${dataDir[i-1]}_1.fq.gz -I raw_data/${dataDir[i-1]}/${dataDir[i-1]}_2.fq.gz -o trimmed/tri${dataDir[i-1]}_1.fq -O trimmed/tri${dataDir[i-1]}_2.fq -f 10 -F 10;
done

rnaspades.py -t 10 -m 10 --only-assembler --ss rf \
--pe1-1 rcorr/Aero_1/unfixrm_triAero_1_1.cor.fq.gz \
--pe1-2 rcorr/Aero_1/unfixrm_triAero_1_2.cor.fq.gz \
--pe2-1 rcorr/Aero_2/unfixrm_triAero_2_1.cor.fq.gz \
--pe2-2 rcorr/Aero_2/unfixrm_triAero_2_2.cor.fq.gz \
--pe3-1 rcorr/Aero_3/unfixrm_triAero_3_1.cor.fq.gz \
--pe3-2 rcorr/Aero_3/unfixrm_triAero_3_2.cor.fq.gz \
--pe4-1 rcorr/Ana_1/unfixrm_triAna_1_1.cor.fq.gz \
--pe4-2 rcorr/Ana_1/unfixrm_triAna_1_2.cor.fq.gz \
--pe5-1 rcorr/Ana_2/unfixrm_triAna_2_1.cor.fq.gz \
--pe5-2 rcorr/Ana_2/unfixrm_triAna_2_2.cor.fq.gz \
--pe6-1 rcorr/Ana_3/unfixrm_triAna_3_1.cor.fq.gz \
--pe6-2 rcorr/Ana_3/unfixrm_triAna_3_2.cor.fq.gz \
-o Assembly/

busco -i RNAseq/Assembly/transcripts.fasta -o RNAseq/Assembly/qc/ -m transcriptome --auto-lineage-prok
busco -i trinity/trinity_output.Trinity.fasta -o trinity/qc/ -m transcriptome --auto-lineage-prok


rnaspades.py -t 10 -m 10 --only-assembler --ss rf --dataset Assembly2/samples.yaml -o Assembly2/