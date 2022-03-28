#!/bin/sh
#Run by chmod u+x sh/count.sh
#sh/count.sh
echo "Analysis Start"

#trim raw data
cd ~/Desktop/FYP

dataDir=(Aero_1 Aero_2 Aero_3 Ana_1 Ana_2 Ana_3)

for i in {1..6};
do htseq-count -r pos -t CDS -i ID Mapped/${dataDir[i-1]}/${dataDir[i-1]}sorted.bam Mapped/JJZ00M.gff -c Count/${dataDir[i-1]}count.csv;
done

#The -i ID used to change the deflut attribute name from 'gene ID' to 'ID'
#The JZZ00M.gff is modified JZZ00.gff, deleted the FASTA sequences at the end using sed