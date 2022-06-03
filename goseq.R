library('DESeq2')
library('rtracklayer')
library('stringr')
library('goseq')
library('pheatmap')

#Directory where the htseq-cout output is
directory = "~/Desktop/FYP/B2/count"

#Extract general statistic from Htseq2

#Creating the dds
sampleFiles = grep("",list.files(directory),value=TRUE)[1:6]
sampleCondition = c(sub("(.*Aero).*","\\1",sampleFiles[1:3]), sub("(.*Ana).*","\\1",sampleFiles[4:6]))
sampleTable = data.frame(sampleName = c(substr(sampleFiles[1:3], 1, 6), substr(sampleFiles[4:6], 1, 5)),
                         fileName = sampleFiles,
                         condition = sampleCondition)
sampleTable$condition = factor(sampleTable$condition)

dds = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition)

#Set Aero as reference level
dds$condition <- relevel(dds$condition, ref = "Aero")

#DE analysis
dds = DESeq(dds) #!!!!!!!!!!!!!!!!!!!! Important
ntd = normTransform(dds) #normalized count in log2(n + 1)
res = results(dds) #DE analysis

summary(res)
sum(res$padj < 0.01, na.rm = T)

#Sort the dds based on P-value
#resOrdered = res[order(res$pvalue),]



#Plot the data
par(mfrow = c(1, 2))
plotMA(res, ylim = c(-5, 5), main = 'MA plot before shrinkage')
#invesitgate linear pattern ?gene-dependent dependencies ?solved by shrkinage (below)?

#Log fold change (LFC) shrinkage
coef = resultsNames(dds)[2]
resLFC = lfcShrink(dds, coef = coef, type = 'apeglm')
#resLFC = lfcShrink(dds, coef = coef, type = 'ashr')
plotMA(resLFC, ylim = c(-5, 5), main = 'MA plot after shrinkage')

#Plot dispersion
plotDispEsts(dds, main = 'dispersion estimates plot')

#QC for normolization and shrikage before and after normailzaiton
par(mfrow = c(2, 1))
log2dds = log2(assay(dds) + 1) #Needs to change this
boxplot(log2dds, main = 'Before normalization')
abline(h = median(log2dds), lty = 2)
boxplot(assay(ntd), main = 'After normalization')
abline(h = median(assay(ntd)), lty = 2)

par(mfrow = c(2, 6))
for(i in 1:6)
{
  tem = data.frame(before = log2dds[, i], after = assay(ntd)[, i])
  boxplot(tem, main = colnames(dds)[i])
  plot(density(tem$before), lty = 2, main = colnames(dds)[i])
  lines(density(tem$after))
}

#Volcano Plot
dim(assay(dds))

sig = subset(resLFC, resLFC$padj <= 0.05)
dim(sig)

upReg = subset(sig, sig$log2FoldChange > 1)
dim(upReg)

downReg = subset(sig, sig$log2FoldChange < -1)
dim(downReg)

nonSig1 = subset(resLFC, resLFC$padj > 0.05)
nonSig2 = subset(sig, sig$log2FoldChange >= -1 & sig$log2FoldChange <= 1)

plot(resLFC$log2FoldChange, -log10(resLFC$padj)) #used to set xlim and ylim

plot(nonSig1$log2FoldChange,
     -log10(nonSig1$padj),
     xlab = expression('Log'[2] * '(Fold Change)'),
     ylab = expression('-Log'[10] * '(Adjusted P-value)'),
     main = 'Volcano plot showing differently regulated genes in Aero vs Ana (|LFC| > 1)',
     xlim = c(-10, 10),
     ylim = c(0, 100),
     pch = 16, col = 'black')
points(nonSig2$log2FoldChange, -log10(nonSig2$padj), pch = 16, col = 'black')
points(upReg$log2FoldChange, -log10(upReg$padj), pch = 24, col = 'red', bg = 'red')
points(downReg$log2FoldChange, -log10(downReg$padj), pch = 25, col = 'blue', bg = 'blue')
abline(h = -log10(0.05), lty = 2)
abline(v = -1, lty = 2)
abline(v = 1, lty = 2)
text(8, 50, paste('UpReg: ', round((dim(upReg)[1]/dim(resLFC)[1]) * 100, 2), '%', sep = ''))
text(-8, 50, paste('DownReg: ', round((dim(downReg)[1]/dim(resLFC)[1]) * 100, 2), '%', sep = ''))
legend('topright',
       legend = c('NonReg', 'UpReg', 'DownReg'),
       pch = c(16, 24, 25),
       col = c('black', 'red', 'blue'),
       pt.bg = c('black', 'red', 'blue'))

#PucBA Heatmap and barchat
namecode = read.csv('~/Desktop/FYP/namecode.csv')
complexes = assay(ntd)[namecode[1:34, 1], ]
rownames(complexes) = namecode[1:34, 2]

complexes = complexes[order(rowSums(complexes), decreasing = T), ]

df = as.data.frame(colData(dds)[,"condition"])
colnames(df) = 'Condition'
rownames(df) = colnames(dds)
pheatmap(complexes,
         cluster_rows = FALSE,
         show_rownames = T,
         cluster_cols = FALSE,
         annotation_col = df,
         main = 'pucBA')



complexes = as.data.frame(resLFC[namecode[1:16, 1], ])
rownames(complexes) = namecode[1:16, 2]
complexes = complexes[order(complexes$log2FoldChange, decreasing = T, na.last = F), ]
complexes = subset(complexes, complexes$padj < 0.05)

#par(mar = c(2, 5, 15, 1)) #down, left, up, right
barplot(complexes$log2FoldChange,
        ylim = c(-1, 10),
        names.arg = rownames(complexes),
        horiz = F,
        las = 2,
        ylab = expression('Log'[2] * ' Fold Change'))

##########################Run Deseq analysis.R before this (use old alignemnt data)
puc.stat = as.data.frame(matrix(NA, nrow = 16, ncol = 5))
colnames(puc.stat) = c('Gene', 'Aeroic Count', 'Anaerobic Count', 'Aerobic Contibution', 'Anaerobic Contribution')
puc.stat$Gene = namecode[1:16, 2]
puc.stat$`Aeroic Count` = rowSums(assay(ntd)[namecode[1:16, 1], 1:3])
puc.stat$`Anaerobic Count` = rowSums(assay(ntd)[namecode[1:16, 1], 4:6])

puc.stat$`Aerobic Contibution`[seq(1, 16, 2)] = puc.stat$`Aeroic Count`[seq(1, 16, 2)]/sum(puc.stat$`Aeroic Count`[seq(1, 16, 2)])
puc.stat$`Aerobic Contibution`[seq(2, 16, 2)] = puc.stat$`Aeroic Count`[seq(2, 16, 2)]/sum(puc.stat$`Aeroic Count`[seq(2, 16, 2)])

puc.stat$`Anaerobic Contribution`[seq(1, 16, 2)] = puc.stat$`Anaerobic Count`[seq(1, 16, 2)]/sum(puc.stat$`Anaerobic Count`[seq(1, 16, 2)])
puc.stat$`Anaerobic Contribution`[seq(2, 16, 2)] = puc.stat$`Anaerobic Count`[seq(1, 16, 2)]/sum(puc.stat$`Anaerobic Count`[seq(2, 16, 2)])

puc.ratio = as.data.frame(matrix(NA, nrow = 9, ncol = 3))
colnames(puc.ratio) = c('Gene', 'Aerobic Ratio', 'Anaerobic Ratio')
puc.ratio$Gene = c('pucBA.7', 'pucBA.9', 'pucBA.10.1', 'pucBA.10.2', 'pucBA.23',
                   'pucBA90.1', 'pucBA90.2', 'pucBA90.3', 'Tolal pucBA')

puc.ratio$`Aerobic Ratio`[1:8] = puc.stat$`Aeroic Count`[seq(1, 16, 2)]/puc.stat$`Aeroic Count`[seq(2, 16, 2)]
puc.ratio$`Aerobic Ratio`[9] = sum(puc.stat$`Aeroic Count`[seq(1, 16, 2)])/sum(puc.stat$`Aeroic Count`[seq(2, 16, 2)])

puc.ratio$`Anaerobic Ratio`[1:8] = puc.stat$`Anaerobic Count`[seq(1, 16, 2)]/puc.stat$`Anaerobic Count`[seq(2, 16, 2)]
puc.ratio$`Anaerobic Ratio`[9] = sum(puc.stat$`Anaerobic Count`[seq(1, 16, 2)])/sum(puc.stat$`Anaerobic Count`[seq(2, 16, 2)])

#Scatter plot
max(puc.ratio$`Aerobic Ratio`)
max(puc.ratio$`Anaerobic Ratio`)

plot(puc.ratio$`Aerobic Ratio`,
     xaxt = 'n',
     xlab = 'pucBA Gene Pair Name',
     ylab = 'Ratio of pucB to pucA (pucB/pucA)')
axis(1, at = 1:9, labels = puc.ratio$Gene)
points(puc.ratio$`Anaerobic Ratio`, pch = 0)
abline(h = 1, lty = 2)
legend('topright',
       legend = c('Aerobic Condition', 'Anaerobic Condition'),
       pch = c(1, 2))

#Pie chart
par(mar = c(1, 1, 1, 1), mfrow = c(2, 2))
#Aerobic
pie(puc.stat$`Aerobic Contibution`[seq(1, 16, 2)],
    labels = paste(puc.stat$Gene[seq(1, 16, 2)],
                   ' (',
                   round(puc.stat$`Aerobic Contibution`[seq(1, 16, 2)], 3) * 100,
                   '%)',
                   sep = ''),
    init.angle = 90,
    main = 'Total pucB Composition in Aerobic Condition')

pie(puc.stat$`Aerobic Contibution`[seq(2, 16, 2)],
    labels = paste(puc.stat$Gene[seq(2, 16, 2)],
                   ' (',
                   round(puc.stat$`Aerobic Contibution`[seq(2, 16, 2)], 3) * 100,
                   '%)',
                   sep = ''),
    init.angle = 90,
    main = 'Total pucA Composition in Aerobic Condition')

#Anaerobic
pie(puc.stat$`Anaerobic Contribution`[seq(1, 16, 2)],
    labels = paste(puc.stat$Gene[seq(1, 16, 2)],
                   ' (',
                   round(puc.stat$`Anaerobic Contribution`[seq(1, 16, 2)], 3) * 100,
                   '%)',
                   sep = ''),
    init.angle = 90,
    main = 'Total pucB Composition in Anaerobic Condition')

pie(puc.stat$`Anaerobic Contribution`[seq(2, 16, 2)],
    labels = paste(puc.stat$Gene[seq(2, 16, 2)],
                   ' (',
                   round(puc.stat$`Anaerobic Contribution`[seq(2, 16, 2)], 3) * 100,
                   '%)',
                   sep = ''),
    init.angle = 90,
    main = 'Total pucA Composition in Anaerobic Condition')


write.csv(puc.stat, '~/Desktop/FYP/Report/puc.stat.csv', quote = F, row.names = T, col.names = T, na = '')
##########################End of using old alignment data





##############################
#GOseq analysis prepare
setwd('~/Desktop/FYP/Annotation/newAnnotation/final')
go = read.csv('goTerms.csv') #Extract pre-paried go term list
gff = readGFF('updated.gff') #Acquire gff for name and ID relationship
NCBIgff = readGFF('NCBI.gff')
transcriptomegff = readGFF('transcriptome.gff')

#Extract general information of gff
types = levels(as.factor(gff$type))
length(subset(gff$seqid, gff$type == 'gene')) #gene
length(subset(gff$seqid, gff$type == 'tRNA')) #tRNA
length(subset(gff$seqid, gff$type == 'rRNA')) #rRNA

types = levels(as.factor(NCBIgff$type))
length(subset(gff$seqid, NCBIgff$type == 'gene')) #gene
length(subset(gff$seqid, NCBIgff$type == 'tRNA')) #tRNA
length(subset(gff$seqid, NCBIgff$type == 'rRNA')) #rRNA

types = levels(as.factor(transcriptomegff$type))
length(subset(gff$seqid, transcriptomegff$type == 'CDS'))

#Creat vector for all genes and up de and down de genes
setwd("~/Desktop/FYP/B2")
assayed.genes = read.delim('Aero_1count.csv')[, 1] #All genes

up.de.genes = subset(resLFC, resLFC$padj < 0.05 & resLFC$log2FoldChange > 1)
up.de.genes = row.names(up.de.genes)
down.de.genes = subset(resLFC, resLFC$padj < 0.05 & resLFC$log2FoldChange < 1)
down.de.genes = row.names(down.de.genes)

up.gene.vector = as.integer(assayed.genes%in%up.de.genes)
names(up.gene.vector) = assayed.genes
down.gene.vector = as.integer(assayed.genes%in%down.de.genes)
names(down.gene.vector) = assayed.genes

#Extract gene length
all.length = data.frame(row.names = gff$ID, length = gff$end - gff$start)
assayed.length = as.vector(sapply(assayed.genes, function(x) all.length[x, ]))

#Extract GO terms
go$geneName = substr(go$geneName, 1, 14) #Correct ID with assed.genes
go[go == '.'] = NA #Correct no go term from . to NA

#Begin of GO analysis
#Fitting PWF
up.pwf = nullp(up.gene.vector, bias.data = assayed.length)
down.pwf = nullp(down.gene.vector, bias.data = assayed.length)

#Wallenius approximation
up.go.wall = goseq(pwf = up.pwf, gene2cat = go, use_genes_without_cat = F)
down.go.wall = goseq(pwf = down.pwf, gene2cat = go, use_genes_without_cat = F)

#Run revigo use abovig 2
revigo.up = up.go.wall[, 1:2]
revigo.down = down.go.wall[, 1:2]

setwd('~/Desktop/FYP/B2/revigo')
write.table(revigo.up, 'revigo.up.csv', sep = ' ', quote = F, row.names = F, col.names = F, na = '')
write.table(revigo.down, 'revigo.down.csv', sep = ' ', quote = F, row.names = F, col.names = F, na = '')

#Remove outlier (count = 1, p-value = 1)
up.go.wall = subset(up.go.wall, up.go.wall$over_represented_pvalue < 0.999 & up.go.wall$numInCat > 1)
down.go.wall = subset(down.go.wall, down.go.wall$over_represented_pvalue < 0.999 & down.go.wall$numInCat > 1)
# < 1 not used due to rounding issues

#Seperate by GO general catagoary
up.bp.go.wall = subset(up.go.wall, up.go.wall$ontology == 'BP')
up.cc.go.wall = subset(up.go.wall, up.go.wall$ontology == 'CC')
up.mf.go.wall = subset(up.go.wall, up.go.wall$ontology == 'MF')

down.bp.go.wall = subset(down.go.wall, down.go.wall$ontology == 'BP')
down.cc.go.wall = subset(down.go.wall, down.go.wall$ontology == 'CC')
down.mf.go.wall = subset(down.go.wall, down.go.wall$ontology == 'MF')

#Perform p-value adjustment for multitesting
up.bp.go.wall$padj = p.adjust(up.bp.go.wall$over_represented_pvalue, method = 'BH')
up.cc.go.wall$padj = p.adjust(up.cc.go.wall$over_represented_pvalue, method = 'BH')
up.mf.go.wall$padj = p.adjust(up.mf.go.wall$over_represented_pvalue, method = 'BH')

down.bp.go.wall$padj = p.adjust(down.bp.go.wall$over_represented_pvalue, method = 'BH')
down.cc.go.wall$padj = p.adjust(down.cc.go.wall$over_represented_pvalue, method = 'BH')
down.mf.go.wall$padj = p.adjust(down.mf.go.wall$over_represented_pvalue, method = 'BH')

#Calculate general statistic
go.stat = as.data.frame(matrix(NA, nrow = 0, ncol = 4))
colnames(go.stat) = c('Type', 'Total', 'P-value < 0.05', 'Padj < 0.05')
go.stat[nrow(go.stat) + 1, ] = c('UP.BP',
                                 length(up.bp.go.wall$category),
                                 length(subset(up.bp.go.wall$category, up.bp.go.wall$over_represented_pvalue < 0.05)),
                                 length(subset(up.bp.go.wall$category, up.bp.go.wall$padj < 0.05)))
go.stat[nrow(go.stat) + 1, ] = c('UP.CC',
                                 length(up.cc.go.wall$category),
                                 length(subset(up.cc.go.wall$category, up.cc.go.wall$over_represented_pvalue < 0.05)),
                                 length(subset(up.cc.go.wall$category, up.cc.go.wall$padj < 0.05)))
go.stat[nrow(go.stat) + 1, ] = c('UP.MF',
                                 length(up.mf.go.wall$category),
                                 length(subset(up.mf.go.wall$category, up.mf.go.wall$over_represented_pvalue < 0.05)),
                                 length(subset(up.mf.go.wall$category, up.mf.go.wall$padj < 0.05)))

go.stat[nrow(go.stat) + 1, ] = c('DOWN.BP',
                                 length(down.bp.go.wall$category),
                                 length(subset(down.bp.go.wall$category, down.bp.go.wall$over_represented_pvalue < 0.05)),
                                 length(subset(down.bp.go.wall$category, down.bp.go.wall$padj < 0.05)))
go.stat[nrow(go.stat) + 1, ] = c('DOWN.CC',
                                 length(down.cc.go.wall$category),
                                 length(subset(down.cc.go.wall$category, down.cc.go.wall$over_represented_pvalue < 0.05)),
                                 length(subset(down.cc.go.wall$category, down.cc.go.wall$padj < 0.05)))
go.stat[nrow(go.stat) + 1, ] = c('DOWN.MF',
                                 length(down.mf.go.wall$category),
                                 length(subset(down.mf.go.wall$category, down.mf.go.wall$over_represented_pvalue < 0.05)),
                                 length(subset(down.mf.go.wall$category, down.mf.go.wall$padj < 0.05)))

write.csv(go.stat, '~/Desktop/FYP/Report/go.stat.csv', quote = F, row.names = T, col.names = T, na = '')

#Select for significantly enriched and unenriched, using coorerction for multitsting
enriched.go = go.wall$category[p.adjust(go.wall$over_represented_pvalue, method = 'BH') < 0.05]
unenriched.go = go.wall$category[p.adjust(go.wall$under_represented_pvalue, method = 'BH') < 0.05]
#Shows no significance, containue with undajust p-value

#Plot
#Extract signiifcant enriched
go.wall$perChange = go.wall$numDEInCat/go.wall$numInCat
enriched.go = subset(go.wall, go.wall$over_represented_pvalue < 0.05)
#Sort based on percentatge enriched
enriched.go = enriched.go[order(enriched.go$perChange, decreasing = F), ]
#Plot
par(mar = c(1, 20, 1, 1))
barplot(enriched.go$perChange, names.arg = enriched.go$term, horiz = T, las = 2)

#Add in up or down, while consider overal sampling size
#Add up number
up.enriched.go = subset(up.go.wall, up.go.wall$padj < 0.05)
enriched.go$up = 0






#Plot up
up.go.wall$perChange = up.go.wall$numDEInCat/up.go.wall$numInCat #Calculate percentate of enrichment
up.go.wall$padj = p.adjust(up.go.wall$over_represented_pvalue, method = 'BH') #Add padj

up.enriched.go = subset(up.go.wall, up.go.wall$padj < 0.05)
up.enriched.go = up.enriched.go[order(up.enriched.go$perChange, decreasing = F), ]

par(mar = c(1, 20, 1, 1))
barplot(up.enriched.go$perChange, names.arg = up.enriched.go$term, horiz = T, las = 2)

#Using undajusted pvalue
up.enriched.go = subset(up.go.wall, up.go.wall$over_represented_pvalue < 0.05)
up.enriched.go = up.enriched.go[order(up.enriched.go$perChange, decreasing = F), ]

par(mar = c(1, 20, 1, 1))
barplot(up.enriched.go$perChange, names.arg = up.enriched.go$term, horiz = T, las = 2)

#Plot down
down.go.wall$perChange = down.go.wall$numDEInCat/down.go.wall$numInCat #Calculate percentate of enrichment
down.go.wall$padj = p.adjust(down.go.wall$over_represented_pvalue, method = 'BH') #Add padj

down.enriched.go = subset(down.go.wall, down.go.wall$over_represented_pvalue < 0.05)
down.enriched.go = down.enriched.go[order(down.enriched.go$perChange, decreasing = F), ]

par(mar = c(1, 20, 1, 1))
barplot(down.enriched.go$perChange, names.arg = down.enriched.go$term, horiz = T, las = 2)



#TO DO
remove cat with p-value = 1 & NumInCat = 1
seperate ontology catagorys into 3, run padjust after sepreation, because each catafoty has different p-value distribution
run revigo (use orignal p-values and all actagory)

for genreal summize of DE, use shrinked LFC > 1

seperate up and down, shiked LFC > 1

general statistic for each step

QC
Raw FastQC
Alignment Gernal
Assembly BUSCO
Annotation: nuber of gene annotated, number of rRNA...
