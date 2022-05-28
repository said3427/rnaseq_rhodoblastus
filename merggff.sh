#!/bin/sh
#Run by chmod u+x sh/test.sh
#sh/test.sh
echo "Analysis Start"

cd ~/Desktop/FYP/Annotation

#Extract transcript from old annotation
gffread -w newAnnotation/oldTranscripts.fa -g Annotation/JJZ00M.fna Annotation/out.gtf

#Extract unique transcripts
cd-hit-est-2d -i newAnnotation/oldTranscripts.fa -i2 newAnnotation/trinity_output.Trinity.fasta -o newAnnotation/novel.fasta -c 0.95 -n 8

experinment on differnt -c and -n
original: 2806
-c (-n 8)
0.99: 2630
0.95: 2599
0.90: 2561
0.80: 2501
-n (-c 0.95)
10: 2599
8: 2599
6: 2599

#Run Prokka
prokka --outdir newAnnotation/prokka --force --cpus 7 newAnnotation/novel.fasta

#Combine the novel transcripts with original transcripts
cat newAnnotation/oldTranscripts.fa newAnnotation/novel.fasta > newAnnotation/trinotateNew/updated.fasta

#Annotation and creation of xls in trinotate.sh

#Merge .gffs
#NOT USED gffcompare -o newAnnotation/mergegff/merged newAnnotation/prokka/PROKKA_05242022.gff Annotation/JJZ00.gff3
#run sh/gffmerge.r

#Testing real transcript numbers
gffread -w newAnnotation/transNum/old.fa -g Annotation/JJZ00M.fna Annotation/out.gtf
gffread -w newAnnotation/transNum/new.fa -g newAnnotation/prokka/PROKKA_05242022.fna newAnnotation/prokka/PROKKA_05242022.gff
cd ~/Desktop/FYP/Annotation/newAnnotation/transNum
cd-hit-est-2d -i old.fa -i2 new.fa -o novelTest.fasta -c 0.95 -n 8

grep -c ">" old.fa
grep -c ">" new.fa
grep -c ">" novelTest.fasta


#Send trinotate and merged gff DONE
#Extact CDS from newggff and CT-HIT-EST-2D, send the numbers DONE
#Create faste for pucBA and run prokkaa and add in gff

#Update the gff and fasta to get the final version
#Add in the pucBA into final gff and fasta
#Try to run to STAR
send gff, fasta, STAR commond
#Combine gene onotology terms





convert2bed -i gff < Annotation/JJZ00.gff > newAnnotation/new.bed

chmod +x gff3ToGenePred Annotation/JJZ00.gff temp.genePred
chmod +x genePredToBed temp.genePred newAnnotation/out.bed
rm temp.genePred

geneBody_coverage.py -r newAnnotation/out.bed -i newAnnotation/cov/Aero_1sorted.bam  -o cov