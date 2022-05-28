library('DESeq2')
library("RColorBrewer")
library('pheatmap')
library('apeglm')
library('DECIPHER')
library('dendextend')

setwd('~/Desktop/FYP')

#Directory where the htseq-cout output is
directory = "~/Desktop/FYP/CountR"

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



investigate = subset(res, res[1][, 1] < 1)

plotMA(investigate)

investigate = subset(investigate, (investigate[1][, 1] > 0.2) &  (investigate[1][, 1] < 0.5))
investigate = subset(investigate, investigate[2][, 1] < -1)
plotMA(investigate)
investigate = investigate[order(investigate$baseMean),]
plotMA(investigate)

#Plot dispersion
plotDispEsts(dds, main = 'dispersion estimates plot')





#QC for normolization and shrikage before and after normailzaiton
par(mfrow = c(2, 1))
log2dds = log2(assay(dds) + 1)
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








#Plot count changes for a single gene, here used the one with smalles padj
plotCounts(dds, gene = which.min(res$padj), intgroup = 'condition')

#subset pucBA genes
#pucBA.7
par(mfrow = c(1, 2))
plotCounts(dds, gene = 'OMJNHDJD_05480', intgroup = 'condition')
plotCounts(dds, gene = 'OMJNHDJD_05490', intgroup = 'condition')
tem = plotCounts(dds, gene = 'OMJNHDJD_05480', intgroup = 'condition', returnData = T)
boxplot(tem$count~tem$condition)
tem = plotCounts(dds, gene = 'OMJNHDJD_05490', intgroup = 'condition', returnData = T)
boxplot(tem$count~tem$condition)
#pucBA.9
plotCounts(dds, gene = 'OMJNHDJD_24380', intgroup = 'condition')
plotCounts(dds, gene = 'OMJNHDJD_24390', intgroup = 'condition')
tem = plotCounts(dds, gene = 'OMJNHDJD_24380', intgroup = 'condition', returnData = T)
boxplot(tem$count~tem$condition)
tem = plotCounts(dds, gene = 'OMJNHDJD_24390', intgroup = 'condition', returnData = T)
boxplot(tem$count~tem$condition)
#pucBA.10.1
plotCounts(dds, gene = 'OMJNHDJD_05480', intgroup = 'condition')
plotCounts(dds, gene = 'OMJNHDJD_05490', intgroup = 'condition')
#pucBA.10.2
plotCounts(dds, gene = 'OMJNHDJD_23260', intgroup = 'condition')
plotCounts(dds, gene = 'OMJNHDJD_23270', intgroup = 'condition')
#pucBA.23
plotCounts(dds, gene = 'OMJNHDJD_44670', intgroup = 'condition')
plotCounts(dds, gene = 'OMJNHDJD_44680', intgroup = 'condition')
#pucBA.90.1
plotCounts(dds, gene = 'OMJNHDJD_42040', intgroup = 'condition')
plotCounts(dds, gene = 'OMJNHDJD_42050', intgroup = 'condition')
#pucBA.90.2
plotCounts(dds, gene = 'OMJNHDJD_43310', intgroup = 'condition')
plotCounts(dds, gene = 'OMJNHDJD_43320', intgroup = 'condition')
#pucBA.90.3
plotCounts(dds, gene = 'OMJNHDJD_43210', intgroup = 'condition')
plotCounts(dds, gene = 'OMJNHDJD_43220', intgroup = 'condition')



#estimate dispersion trend and apply a variance stabilizing transformation
vsd = vst(dds, blind = F)
#For normalisez count only assay(vsd)

#transform to log2 count
rld = rlog(dds, blind = F)



#Heat map of sample-to-sample distance
sampleDists = dist(t(assay(vsd))) #Calculate distances
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = colnames(vsd)
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main = 'sample-to-sample distance',
         cex = 1.5)

#PCA plot
plotPCA(vsd, intgroup= "condition")

#Volcano plot
sig = subset(resLFC, resLFC$padj <= 0.05)
dim(sig) #1071

upReg = subset(sig, sig$log2FoldChange > 0)
dim(upReg) #572

downReg = subset(sig, sig$log2FoldChange < 0)
dim(downReg) #499

nonSig = subset(resLFC, resLFC$padj > 0.05)

plot(resLFC$log2FoldChange, -log10(resLFC$padj)) #used to set xlim and ylim

plot(nonSig$log2FoldChange,
     -log10(nonSig$padj),
     xlab = expression('Log'[2] * '(Fold Change)'),
     ylab = expression('-Log'[10] * '(Adjusted P-value)'),
     main = 'Volcano plot showing differently regulated genes in Aero vs Ana',
     xlim = c(-10, 10),
     ylim = c(0, 70),
     pch = 16, col = 'black')
points(upReg$log2FoldChange, -log10(upReg$padj), pch = 24, col = 'red', bg = 'red')
points(downReg$log2FoldChange, -log10(downReg$padj), pch = 25, col = 'blue', bg = 'blue')
abline(h = -log10(0.05), lty = 2)
abline(v = 0, lty = 2)
text(3, 50, paste('UpReg: ', round((dim(upReg)[1]/dim(resLFC)[1]) * 100, 2), '%', sep = ''))
text(-5, 50, paste('DownReg: ', round((dim(downReg)[1]/dim(resLFC)[1]) * 100, 2), '0%', sep = ''))
legend('topright',
       legend = c('NonReg', 'UpReg', 'DownReg'),
       pch = c(16, 24, 25),
       col = c('black', 'red', 'blue'),
       pt.bg = c('black', 'red', 'blue'))

#Volcano plot > 1
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
     ylim = c(0, 70),
     pch = 16, col = 'black')
points(nonSig2$log2FoldChange, -log10(nonSig2$padj), pch = 16, col = 'black')
points(upReg$log2FoldChange, -log10(upReg$padj), pch = 24, col = 'red', bg = 'red')
points(downReg$log2FoldChange, -log10(downReg$padj), pch = 25, col = 'blue', bg = 'blue')
abline(h = -log10(0.05), lty = 2)
abline(v = -1, lty = 2)
abline(v = 1, lty = 2)
text(3, 50, paste('UpReg: ', round((dim(upReg)[1]/dim(resLFC)[1]) * 100, 2), '%', sep = ''))
text(-5, 50, paste('DownReg: ', round((dim(downReg)[1]/dim(resLFC)[1]) * 100, 2), '0%', sep = ''))
legend('topright',
       legend = c('NonReg', 'UpReg', 'DownReg'),
       pch = c(16, 24, 25),
       col = c('black', 'red', 'blue'),
       pt.bg = c('black', 'red', 'blue'))

#Volcano plot > 2
sig = subset(resLFC, resLFC$padj <= 0.05)
dim(sig) #1071

upReg = subset(sig, sig$log2FoldChange > 2)
dim(upReg) #572

downReg = subset(sig, sig$log2FoldChange < -2)
dim(downReg) #499

nonSig1 = subset(resLFC, resLFC$padj > 0.05)
nonSig2 = subset(sig, sig$log2FoldChange >= -2 & sig$log2FoldChange <= 2)

plot(resLFC$log2FoldChange, -log10(resLFC$padj)) #used to set xlim and ylim

plot(nonSig1$log2FoldChange,
     -log10(nonSig1$padj),
     xlab = expression('Log'[2] * '(Fold Change)'),
     ylab = expression('-Log'[10] * '(Adjusted P-value)'),
     main = 'Volcano plot showing differently regulated genes in Aero vs Ana (|LFC| > 2)',
     xlim = c(-10, 10),
     ylim = c(0, 70),
     pch = 16, col = 'black')
points(nonSig2$log2FoldChange, -log10(nonSig2$padj), pch = 16, col = 'black')
points(upReg$log2FoldChange, -log10(upReg$padj), pch = 24, col = 'red', bg = 'red')
points(downReg$log2FoldChange, -log10(downReg$padj), pch = 25, col = 'blue', bg = 'blue')
abline(h = -log10(0.05), lty = 2)
abline(v = -2, lty = 2)
abline(v = 2, lty = 2)
text(3, 50, paste('UpReg: ', round((dim(upReg)[1]/dim(resLFC)[1]) * 100, 2), '%', sep = ''))
text(-5, 50, paste('DownReg: ', round((dim(downReg)[1]/dim(resLFC)[1]) * 100, 2), '0%', sep = ''))
legend('topright',
       legend = c('NonReg', 'UpReg', 'DownReg'),
       pch = c(16, 24, 25),
       col = c('black', 'red', 'blue'),
       pt.bg = c('black', 'red', 'blue'))

#Heatmap
namecode = read.csv('namecode.csv')[1:16,]
namecode$type = c('B', 'A')

ntd = normTransform(dds) #this gives log2(n + 1)

#pucB
tem = subset(namecode, namecode$type == 'B')
complexes = assay(ntd)[tem[, 1],]
rownames(complexes) = tem[, 2]
df = as.data.frame(colData(dds)[,"condition"])
colnames(df) = 'Condition'
rownames(df) = colnames(ntd)
pheatmap(complexes,
         cluster_rows = T,
         show_rownames = T,
         cluster_cols = T,
         annotation_col = df,
         main = 'pucB')


#pucA
tem = subset(namecode, namecode$type == 'A')
complexes = assay(ntd)[tem[, 1],]
rownames(complexes) = tem[, 2]
df = as.data.frame(colData(dds)[,"condition"])
colnames(df) = 'Condition'
rownames(df) = colnames(ntd)
pheatmap(complexes,
         cluster_rows = T,
         show_rownames = T,
         cluster_cols = T,
         annotation_col = df,
         main = 'pucA')


#pucBA
complexes = assay(ntd)[namecode[, 1],]
rownames(complexes) = namecode[, 2]
df = as.data.frame(colData(dds)[,"condition"])
colnames(df) = 'Condition'
rownames(df) = colnames(dds)
pheatmap(complexes,
         cluster_rows = FALSE,
         show_rownames = T,
         cluster_cols = FALSE,
         annotation_col = df,
         main = 'pucBA unclustered')

pheatmap(complexes,
         cluster_rows = T,
         show_rownames = T,
         cluster_cols = T,
         annotation_col = df,
         main = 'pucBA clustered')



select = order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df = as.data.frame(colData(dds)[,"condition"])
rownames(df) = colnames(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)







#Sequence distance of pucBA genes
seq = DNAStringSet(namecode$seq)
seq = OrientNucleotides(seq) #Orient sequences
seq = AlignSeqs(seq) #Align sequences
#BrowseSeqs(seq) #Browse sequence alighment in browser
dmatrix = DistanceMatrix(seq, type = 'dist')
tree = IdClusters(dmatrix, cutoff = 10, method = "NJ", showPlot = F, type = "dendrogram")
tree = as.dendrogram(tree)
labels(tree) = namecode[labels(tree), 2]
plot(tree, horiz = F, main = 'sequence dendrogram')

#Heat map in order of dendrogram with dendrogram
complexes = assay(ntd)[namecode[, 1],]
rownames(complexes) = namecode[, 2]
complexes = complexes[rev(match(labels(tree), rownames(complexes))),]#reorder according to dendrogram seequences
df = as.data.frame(colData(dds)[,"condition"])
colnames(df) = 'Condition'
rownames(df) = colnames(dds)
pheatmap(complexes,
         cluster_rows = FALSE,
         show_rownames = T,
         cluster_cols = FALSE,
         annotation_col = df,
         main = 'pucBA unclustered')

plot(tree, horiz = T, main = 'sequence dendrogram')

plot_horiz.dendrogram(tree)

#For poster
complexes = assay(ntd)[namecode[, 1], 4:6]
colnames(complexes) = c('Sample_1', 'Sample_2', 'Sample_3')
complexes = as.data.frame(complexes)
complexes['Globle Median',] = colMedians(assay(ntd))[1:3]
rownames(complexes)[1:16] = namecode[, 2]
tem = rowSums(complexes)
complexes = complexes[order(tem, decreasing = T),]
pheatmap(complexes,
         cluster_rows = FALSE,
         show_rownames = T,
         cluster_cols = FALSE,
         main = 'pucBA unclustered',
         angle_col = 0)

par(mar = c(5, 12, 5, 12))
plot(tree, horiz = T, main = 'sequence dendrogram')

plot_horiz.dendrogram(tree)

#full puc dendrogram
fullPuc = read.csv('fullPUC.csv')
aaseq  = AAStringSet(fullPuc$pro)
aaseq = AlignSeqs(aaseq)
dmatrix = DistanceMatrix(aaseq, type = 'dist')
tree = IdClusters(dmatrix, cutoff = 10, method = "NJ", showPlot = F, type = "dendrogram")
tree = as.dendrogram(tree)
labels(tree) = fullPuc[labels(tree), 2]
par(cex = 0.5)
plot(tree, horiz = T)

