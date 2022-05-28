library('rtracklayer')

setwd('~/Desktop/FYP/Annotation/newAnnotation/final')

oldgff =readGFF('JJZ00.gff')
newgff = readGFF('prokka/novel.gff')
pucgff = readGFF('prokka/puc.gff')

dim(oldgff)
dim(newgff)
dim(pucgff)

#Match the colname
colnames(oldgff)
colnames(newgff)
colnames(pucgff)

#Fix newgff by removing db_xref values to NA and rename as Parent
colOrder = match(colnames(newgff), colnames(oldgff))
colOrder
colnames(newgff)[which(is.na(colOrder))]
newgff[, which(is.na(colOrder))] = rep(NA, nrow(newgff))
colnames(newgff)[which(is.na(colOrder))] = 'Parent'
#Fix newgff by adding in a column called protein_id with NA values
newgff[, 'protein_id'] = rep(NA, nrow(newgff))

#Fix pucgff by adding colunm with NA and correct colnames
colOrder = match(colnames(pucgff), colnames(oldgff))
colOrder #lacking col 10 13 15 16
for(i in 1:4)
{
  pucgff[, ncol(pucgff) + 1] = rep(NA, nrow(pucgff))
}
colnames(pucgff)[15:18] = colnames(oldgff[c(10, 13, 15, 16)])

#Reorder columns
newgff = newgff[, match(colnames(newgff), colnames(oldgff))]
pucgff = pucgff[, match(colnames(pucgff), colnames(oldgff))]

#Combine gffs
merged = rbind(oldgff, newgff, pucgff)

#Export combined gff
export(merged, 'updated.gff', format = 'GFF3')

