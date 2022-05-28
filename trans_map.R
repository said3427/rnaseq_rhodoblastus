setwd('~/Desktop/FYP/Annotation/newAnnotation/final')

tra = read.delim('updatedTra.fasta', header = F)
df = data.frame(matrix(NA, nrow = 0, ncol=2))

for(i in tra[, 1])
{
  if(substr(i, 1, 1) == '>')
  {
    tem = gregexpr(pattern = ' ', i)[[1]][1] #Identify the index of space
    
    if(tem < 0)
    {
      df[nrow(df) + 1, ] = rep(substr(i, 2, nchar(i)))
    }
    else
    {
      df[nrow(df) + 1, ] = rep(substr(i, 2, tem - 1))
    }
  }
}

write.table(df, file = 'updatedTra.fasta.gene_trans_map', sep = '\t', quote = F, row.names = F, col.names = F)

while(F)
{
  for(i in oldtra[, 1])
  {
    if(substr(i, 1, 2) == '>O')
      if(substr(i, 16, 20) == '_gene')
        df[nrow(df) + 1, ] = rep(substr(i, 2, 20))
    else
      df[nrow(df) + 1, ] = rep(substr(i, 2, 15))
  }
  
  write.table(df, file = 'newAnnotation/oldTranscripts.fa.gene_trans_map', sep = '\t', quote = F, row.names = F, col.names = F)
  
  ######################
  
  oldgff = read.delim('Annotation/JJZ00.gff3', skip = 5, header = F, sep = '\t', nrows = 9028)
  
  paste(substr(oldgff[1, 9], 4, 17), '_gene', sep = '')
  
  df = data.frame(matrix(NA, nrow = 4514, ncol=2))
  
  i = 1
  j = 1
  while(i < 9029)
  {
    df[j, ] =  paste(substr(oldgff[i, 9], 4, 17), '_gene', sep = '')
    j = j + 1
    i = i + 2
  }
  
  write.table(df, file = 'newAnnotation/oldTranscripts.fa.gene_trans_map', sep = '\t', quote = F, row.names = F, col.names = F)
}




