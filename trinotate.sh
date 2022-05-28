#!/bin/sh
#Run by chmod u+x sh/test.sh
#sh/test.sh
echo "Analysis Start"

cd ~/Desktop/FYP/trinotate

#Extracte the long open reading frames
TransDecoder.LongOrfs -t trinity_output.Trinity.fasta
#For updated fasta, runed the second line for each step
TransDecoder.LongOrfs -t Annotation/newAnnotation/trinotateNew/updated.fasta

#Search Trinity transcripts
blastx -query trinity_output.Trinity.fasta -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blast/blastx.outfmt6
blastx -query Annotation/newAnnotation/trinotateNew/updated.fasta -db trinotate/uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > Annotation/newAnnotation/trinotateNew/blastx.outfmt6

#Search Transdecoder-predicted proteins
blastp -query transdecoder/longest_orfs.pep -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blast/blastp.outfmt6
blastp -query Annotation/newAnnotation/trinotateNew/longest_orfs.pep -db trinotate/uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > Annotation/newAnnotation/trinotateNew/blastp.outfmt6

#Running HMMER to identify protein domains
hmmscan --cpu 15 --domtblout pfam/TrinotatePFAM.out Pfam-A.hmm transdecoder/longest_orfs.pep > pfam/pfam.log
hmmscan --cpu 15 --domtblout Annotation/newAnnotation/trinotateNew/PFAM.out trinotate/Pfam-A.hmm  Annotation/newAnnotation/trinotateNew/longest_orfs.pep >  Annotation/newAnnotation/trinotateNew/pfam.log

#Load the data into Trinotate SQLite Database
Trinotate Trinotate.sqlite init --gene_trans_map trinity_output.Trinity.fasta.gene_trans_map --transcript_fasta trinity_output.Trinity.fasta --transdecoder_pep transdecoder/longest_orfs.pep
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blast/blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam pfam/TrinotatePFAM.out
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blast/blastx.outfmt6

#run rscript trans_map.r
cat Annotation/newAnnotation/oldTranscripts.fa.gene_trans_map \
Annotation/newAnnotation/trinity_output.Trinity.fasta.gene_trans_map \
> Annotation/newAnnotation/trinotateNew/updated.fasta.gene_trans_map

cd ~/Desktop/FYP/Annotation/newAnnotation/trinotateNew
Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
Trinotate Trinotate.sqlite init --gene_trans_map updated.fasta.gene_trans_map --transcript_fasta updated.fasta --transdecoder_pep longest_orfs.pep
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam PFAM.out
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6

#Output annotation file
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls

#Transcrip quantification (NOT USED)
#Prepare rsem reference
#rsem-prepare-reference [options] reference_fasta_file(s) reference_name
rsem-prepare-reference --bowtie2 --transcript-to-gene-map trinotate/trinity_output.Trinity.fasta.gene_trans_map trinotate/trinity_output.Trinity.fasta rsem/rsem


#Calculate expression values
#rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name 
rsem-calculate-expression --strandedness reverse \
-p 15 \
--bowtie2 \
--append-names \
--time \
--paired-end \
RNAseq/rcorr/Aero_3/unfixrm_triAero_3_1.cor.fq \
RNAseq/rcorr/Aero_3/unfixrm_triAero_3_2.cor.fq \
rsem/rsem \
rsem/Aero_3/Aero_3

#Combine new transcript annotation with old genome annotation
PASA=/Users/maxwellhou/opt/anaconda3/envs/rnaseq/opt/PASA
base_dir=/Users/maxwellhou/Desktop/FYP/Annotation
cd ~/Desktop/FYP/Annotation

#Generate .sqlite for old annoation(NOT USED)
python3
db=gffutils.create_db("Annotation/JJZ00.gff3", "newAnnotation/old.sqlite")
exit()


docker run --rm -it -v /tmp:/tmp -v $base_dir:$base_dir \
pasapipeline/pasapipeline:latest \
bash -c 'cd /$base_dir && \
cd ./Users/maxwellhou/Desktop/FYP/Annotation && \
/usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
-c newAnnotation/alignAssembly.config -g Annotation/JJZ00M.fna \
-P Annotation/JJZ00.gff3'

#Test
docker run --rm -it -v /tmp:/tmp -v $base_dir:$base_dir \
pasapipeline/pasapipeline:latest \
bash -c 'cd /$base_dir && \
cd ./Users/maxwellhou/Desktop/FYP/Annotation/test && \
/usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config -g genome_sample.fasta \
-P orig_annotations_sample.gff3'


$PASA/scripts/Load_Current_Gene_Annotations.dbi \
-c newAnnotation/alignAssembly.config -g Annotation/JJZ00M.fna \
-P Annotation/JJZ00.gff3

