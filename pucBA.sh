#!/bin/sh
#Run by chmod u+x sh/test.sh
#sh/test.sh
echo "Analysis Start"

cd ~/Desktop/FYP/Annotation/newAnnotation/pucFasta

#Extract unique transcripts
cd-hit-est-2d -i Eve.fasta -i2 NCBI.fasta -o novel1.fasta -c 0.95 -n 8
cat Eve.fasta novel1.fasta > tem.fasta

cd-hit-est-2d -i tem.fasta -i2 NRRJ00000000.1.fasta -o novel2.fasta -c 0.95 -n 8
cat tem.fasta novel2.fasta > pucBA.fasta



grep -c ">" old.fa
grep -c ">" new.fa
grep -c ">" novelTest.fasta
