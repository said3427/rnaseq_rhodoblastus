#!/bin/sh
#Run by chmod u+x sh/trinoAuto.sh
#sh/trinoAuto.sh

echo "Analysis Start"

cd ~/Desktop/FYP/Annotation/newAnnotation/final

#Run trinotate on updated transcript fasta
Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

echo "blastx"

blastx -query updatedTra.fasta -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

echo "blastp"

blastp -query longest_orfs.pep -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

echo "pfam"

hmmscan --cpu 15 --domtblout PFAM.out Pfam-A.hmm longest_orfs.pep > pfam.log

echo "Complete"