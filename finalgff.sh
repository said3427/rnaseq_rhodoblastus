#!/bin/sh
#Run by chmod u+x sh/test.sh
#sh/test.sh

#Update the fasta and gff to get the final version
#Add in the pucBA into final gff and fasta

echo "Analysis Start"

cd ~/Desktop/FYP/Annotation

#Extract transcript from old annotation
gffread -w newAnnotation/final/old.fa -g Annotation/JJZ00M.fna Annotation/out.gtf
#Extrac transcript from new annotation
gffread -w newAnnotation/final/new.fa -g newAnnotation/prokka/PROKKA_05242022.fna newAnnotation/prokka/PROKKA_05242022.gff

#Extract novel transcripts
cd ~/Desktop/FYP/Annotation/newAnnotation/final
cd-hit-est-2d -i old.fa -i2 new.fa -o novel.fasta -c 0.95 -n 8

#Create Prokka annotation for novel transcripts
prokka --outdir prokka --prefix novel --force --cpus 7 novel.fasta
#Create Prokka annotaiton for pucBA genes
prokka --outdir prokka --prefix puc --force --cpus 7 pucBA.fasta

#Append novel transcripts to the end of genome
cat JJZ00.fna novel.fasta > tem.fasta
#Append pucBA genes to the end of genome
cat tem.fasta pucBA.fasta > updated.fasta

#Create transcript fasta
gffread -w pucTra.fasta -g pucBA.fasta prokka/puc.gff
cat old.fa new.fa > tem.fa
cat tem.fa pucTra.fasta > updatedTra.fasta

#Get gene_trans_map for updatedTra.fasta using sh/trans_map.R

#Create updated gff using sh/newgffmerge.R

#Run trinotate on updated transcript fasta
Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
TransDecoder.LongOrfs -t updatedTra.fasta
blastx -query updatedTra.fasta -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query longest_orfs.pep -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
hmmscan --cpu 15 --domtblout PFAM.out Pfam-A.hmm longest_orfs.pep > pfam.log
#Load the data into Trinotate SQLite Database
Trinotate Trinotate.sqlite init --gene_trans_map updatedTra.fasta.gene_trans_map --transcript_fasta updatedTra.fasta --transdecoder_pep longest_orfs.pep
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam PFAM.out
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
#Output annotation file
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls

#Run STAR (Said)

#Conver gff to gtf
agat_convert_sp_gff2gtf.pl --gff Annotation/newAnnotation/final/updated.gff -o Annotation/newAnnotation/final/updated.gtf

#Creat genome index
STAR --runThreadN 5 --runMode genomeGenerate --genomeSAindexNbases 10 --genomeDir finalmap/Indexed --genomeFastaFiles Annotation/newAnnotation/final/updated.fasta --sjdbGTFfile Annotation/newAnnotation/final/updated.gtf
