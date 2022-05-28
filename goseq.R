library('rtracklayer')
library('stringr')

setwd('~/Desktop/FYP/Annotation/newAnnotation/final')

annotation = read.csv('trinotate_annotation_report.csv')
gff =readGFF('updated.gff')

#Extract all gene ID
assayed.genes = gff$ID

#Extract all gene length
gene.length = gff$end - gff$start

#Extract all GO terms
go = as.data.frame(matrix(NA, nrow = 0, ncol=2))



for(i in 1:dim(annotation)[1])
{
  geneName = annotation[i, 1]
  goTerms = unlist(sapply(str_split(annotation[i, 13:15], '`'), function(x) substr(x, 1, 10)))
  goTerms = goTerms[!duplicated(goTerms)]
  go = rbind(go, data.frame(geneName, goTerms))
}

