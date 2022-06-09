# Globle Transcriptomic Analysis on *Rhodoblastus acidophilus*
Repository for detailed codes and commands used in the project, aimed for a reproducible workflow.

Each section below is in the same sequence as in the methodology section. All command are to be run a command lines. For each sample command, sample name Aero_1 is used as an example. All custom scripts for this project are linked in the text below, like so [demo.code](demo.code).

## 1. Prerequisite packages and setup
All packages, except `R` packages, used in this project are stored in [condaEnvironment.yml](condaEnvironment.yml). This file can be used to directly create `conda` environment, as shown below, for all analysis illustrated.
```
conda env create -f condaEnvironment.yml
```
All shell scripts `.sh` should be configured for running using the command below.
```
chmod u+x sh/shellScriptName.sh
```

## 2. Raw data processing
The raw data is first quality checked using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
```
fastqc Aero_1_1.fq
```
Custom shell script [fastqc.sh](fastqc.sh) is used to automatically apply this step for all samples.
```
sh/fastqc.sh
```
All `FastQC` reports are combined together using [MultiQC](https://multiqc.info).
```
multiqc .
```
The `FastQC` reports show there are adaptor sequences remained from each sequence, judging from `per base sequence content`. Thus, each sequence are trimmed by length using [fastp](https://github.com/OpenGene/fastp).
```
fastp -i raw_data/Aero_1/Aero_1_1.fq.gz -I raw_data/Aero_1/Aero_1_2.fq.gz -o trimmed/triAero_1.fq -O trimmed/triAero_2.fq -f 10 -F 10
```
Custom shell script [trim.sh](trim.sh) is used to automatically apply this step for all samples.
```
sh/trim.sh
```
The trimmed sequences are checked again using `FastQC` for potential problems.

## 3. *De Novo* Transcriptomic Assembly
The sequences after adapter trimming where further marked for random sequencing errors using [Rcorrector](https://github.com/mourisl/Rcorrector). 
```
run_rcorrector.pl -1 trimmed/triAero_1_1.fq -2 trimmed/triAero_1_2.fq -od rcorr/Aero_1/ -t 14
```
 Then the marked errors were removed using [FilterUncorrectabledPEfastq.py](https://github.com/harvardinformatics/TranscriptomeAssemblyTools/blob/master/FilterUncorrectabledPEfastq.py) script.
 ```
python ~/Desktop/FYP/RNAseq/FilterUncorrectabledPEfastq.py -1 triAero_1_1.cor.fq -2 triAero_1_2.cor.fq -s Aero_1
 ```
Custom shell script [fmerfilter.sh](kmerfilter.sh) are used for automatically applying the above 2 steps for all samples, including compression of result files.
```
sh/fmerfilter.sh
```
 The processed sequences from different growth condition are then used together for de novo transcriptomic assembly using [Trinity](https://github.com/alexdobin/STAR).
 ```
 Trinity --seqType fq --max_memory 15G --CPU 15 --samples_file samples.txt --SS_lib_type RF --output Trinity
 ```
The --samples_file [samples.txt](samples.txt) specify all RNA sequences being used.

## 4. Identification and adidition of novel transcripts and pucBA complexes found in NCBI database to genome sequences and genome annotation

The assembled transcriptome and the pucBA operon genes found in NCBI database are annotated using [Prokka](https://github.com/tseemann/prokka).
```
prokka --outdir newAnnotation/prokka --force --cpus 7 newAnnotation/trinit.Trinity.out.fasta

prokka --outdir newAnnotation/prokka --force --cpus 7 newAnnotation/pucBA.fasta
```
- --force used here to write in pre-existing folder.

The transcript sequences are extracted from original annotated genome sequences, annotated transcriptome sequence, and identified pucBA operon genes using [GffRead](http://ccb.jhu.edu/software/stringtie/gff.shtml).
```
gffread -w newAnnotation/final/old.fa -g Annotation/JJZ00M.fna Annotation/out.gtf

gffread -w newAnnotation/final/new.fa -g newAnnotation/prokka/PROKKA_05242022.fna newAnnotation/prokka/PROKKA_05242022.gff
```
The novel transcripts are identified by comparing between the original transcripts and assembled transcriptome transcripts using [CD-HIT-EST-2D](http://www.bioinformatics.org/cd-hit/cd-hit-user-guide). The identified pucBA operon genes are also compared with original transcripts to remove duplicated transcripts.
```
cd-hit-est-2d -i old.fa -i2 new.fa -o novel.fasta -c 0.95 -n 8

cd-hit-est-2d -i old.fa -i2 pucBA.fa -o pucBA.fasta -c 0.95 -n 8
```

The novel transcripts and pucBA operon genes are added into the original genome sequences using `cat`.
```
cat old.fa new.fa > tem.fa

cat tem.fa pucBA.fasta > updatedTra.fasta
```
The novel transcripts are annotated again using `Prokka`. This is faster compare to writing a custom script to match and extract annotation from whole transcriptome annotation file.
```
prokka --outdir newAnnotation/prokka --force --cpus 7 newAnnotation/trinit.Trinity.out.fasta
```
The identified pucBA operon genes from NCBI databases are first compared between each other to remove duplicated genes with different names using `CD-HIT-EST-2d` and `cat`. Then, the genes are annotated using `Prokka`.
```
cd-hit-est-2d -i Eve.fasta -i2 NCBI.fasta -o novel1.fasta -c 0.95 -n 8

cat Eve.fasta novel1.fasta > tem.fasta

cd-hit-est-2d -i tem.fasta -i2 NRRJ00000000.1.fasta -o novel2.fasta -c 0.95 -n 8

cat tem.fasta novel2.fasta > pucBA.fasta

prokka --outdir newAnnotation/prokka --force --cpus 7 newAnnotation/pucBA.fasta
```
Now, the components of updated transcripts sequences, genome sequences, and genome annotation, are ready and assembled.

- The transcripts sequences are put together using `cat`.
```
cat old.fa novel.fa > tem.fa

cat tem.fa pucBA.fasta > updatedTra.fasta
```
- The genome sequences are put together using `cat`.
```
cat JJZ00.fna novel.fasta > tem.fasta

cat tem.fasta pucBA.fasta > updated.fasta
```
- The genome annotations are put together using custom r script [newgffmerge.R](newgffmerge.R), utilizing [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html).
```
Rscript newgffmerge.R
```

## 5. Alignment, count table generation, and DE analysis
The alignment of each reads to the updated annotated genome was done using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
```
bowtie2 -x B2/index -1 RNAseq/rcorr/Aero_1/unfixrm_triAero_1_1.cor.fq -2 RNAseq/rcorr/${dataDir[i-1]}/unfixrm_triAero_1_2.cor.fq -S B2/Aero_1.sam -p 15 --rf
```
The file conversion from `.sam` to sorted and indexed `.bam` is carried out using [samtools](https://www.htslib.org).
```
samtools view -bS Aero_1.sam > Aero_1.bam

rm Aero_1.sam

samtools sort Aero_1.bam -o Aero_1.bam

samtools index Aero_1.bam
```
The count table was then generated using the alignment file by using [HTseq](https://htseq.readthedocs.io/en/master/).
```
htseq-count -r pos -t CDS -i ID -n 7 -s reverse Aero_1.bam updated.gff -c $Aero_1count.csv
```
- Non-unique alignment reads excluded from counting, as HTseq default behaviour.

For the alignment, file conversion, and generation of count table, Custom shell script [b2.sh](b2.sh) is used to automatically apply this step for all samples.
```
sh/b2.sh
```
The alignment were quality checked using [RSeQC](http://rseqc.sourceforge.net). The `bam_stat` and `infer_expriment` QC reports are generated. 
```
bam_stat.py -i Aero_1.bam > qc/stat/Aero_1bam_stat.txt

infer_experiment.py -i Aero_1.bam -r ~/Desktop/FYP/Annotation/Annotation/ag.bed > qc/infer_experiment/Aero_1infer_experiment.txt
```
The counting matrix (read assignment) is quality checked using [MultiQC](https://multiqc.info).
```
MultiQC *
```

Finally, the statistical analysis including raw count normalization and differential analysis was done in R using DEseq2 following standard protocol using custom R script [Deseqanalysis.R](Deseqanalysis.R)
```
Rscript Deseq analysis.R
```
 
Here, the same procedure was run using the original annotated genome for comparison to identify any introduced error during the update of genome annotation and genome sequences.



## 6. Gene ontology annotation and pathway enrichment analysis

The `.gene_trans.map` is required for gene ontology annotation and file is created using custom r script [trans_map.R](trans_map.R).
```
Rscript trans_map.R
```
The updated transcripts sequences are annotated for gene ontology terms using [Trinotate](https://github.com/Trinotate/Trinotate.github.io/blob/master/index.asciidoc) and associated packages. The Trinotate is set to search through SwissProt for gene ontology terms and Pfam for protein domains identification (Bairoch, 2000; Punta et al., 2012).

- Initiate new Trinotate sqlite database.
```
Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
```
- Prepare protein database.
```
makeblastdb -in uniprot_sprot.pep -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```
- Generate more likely longest-ORF peptide candidates using [TransDecoder](https://github.com/TransDecoder/TransDecoder/releases).
```
TransDecoder.LongOrfs -t updatedTra.fasta
```
- Run database search.
```
blastx -query updatedTra.fasta -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query longest_orfs.pep -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
hmmscan --cpu 15 --domtblout PFAM.out Pfam-A.hmm longest_orfs.pep > pfam.log
```
- Load the data into Trinotate SQLite Database.
```
Trinotate Trinotate.sqlite init --gene_trans_map updatedTra.fasta.gene_trans_map --transcript_fasta updatedTra.fasta --transdecoder_pep longest_orfs.pep
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam PFAM.out
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
```
- Output annotation file.
```
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
```

The pathway enrichment analysis is carried out using [GOseq](https://bioconductor.org/packages/release/bioc/html/goseq.html) in r using custom R script [goseq.R](goseq.R).
```
Rscript goseq.R
```

## 7. Dendrogram construction
The R package [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [DECIPHER](http://www2.decipher.codes), and [dendextend](https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html) was used with agglomeration method set as Neighbor-Joining method for the construction. The code for its construction is stored in R script [Deseqanalysis.R](Deseqanalysis.R).
