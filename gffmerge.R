setwd('~/Desktop/FYP/Annotation')

oldgff1 = read.delim('Annotation/JJZ00.gff', header = F, nrows = 5)
oldgff2 = read.delim('Annotation/JJZ00.gff', header = F, nrows = 9033-5, skip = 5)
oldgff3 = read.delim('Annotation/JJZ00.gff', header = F, skip = 9033)

newgff1 = read.delim('newAnnotation/prokka/PROKKA_05242022.gff', header = F, nrows = 2600-1, skip = 1)
newgff2 = read.delim('newAnnotation/prokka/PROKKA_05242022.gff', header = F, nrows = 15327-2600, skip = 2600)
newgff3 = read.delim('newAnnotation/prokka/PROKKA_05242022.gff', header = F, skip = 15327)

df1 = rbind(oldgff1, newgff1)
df2 = rbind(oldgff2, newgff2)
df3 = rbind(oldgff3, newgff3)

df1[nrow(df1), 2:9] = matrix(NA, nrow = nrow(df1), ncol=8)
df3[nrow(df3), 2:9] = matrix(NA, nrow = nrow(df3), ncol=8)

df = rbind(df1, df2, df3)

write.table(df, file = 'newAnnotation/merged.gff', sep = '\t', quote = F, row.names = F, col.names = F, na = '')
