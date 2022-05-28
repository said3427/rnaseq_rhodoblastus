setwd('~/Desktop/FYP')

fullPuc = read.csv('fullPUC.csv')

fullPuc = fullPuc[-c(47:68),]



for(i in levels(as.factor(fullPuc$sourse)))
{
  fasta = as.data.frame(matrix(NA, nrow = 0, ncol=1))
  dftem = subset(fullPuc, fullPuc$sourse == i)
  for(j in 1:dim(dftem)[1])
  {
    fasta[nrow(fasta)+1, ] =paste('>', dftem[j, 2], sep = '')
    fasta[nrow(fasta)+1, ] = dftem[j, 3]
  }
  write.table(fasta, file = paste('pucFasta/', i, '.fasta', sep = ''), sep = '\t', quote = F, row.names = F, col.names = F, na = '')
}levels(as.factor(fullPuc$sourse))
