#RNAseq
setwd("/Users/eveli/Desktop/rnaseq/")

#BiocManager::install("DESeq2")
#BiocManager::install("EnhancedVolcano")

library('DESeq2')
library("RColorBrewer")
library('pheatmap')
library('apeglm')
library("EnhancedVolcano")
library("ape")
library("Biostrings")

ff <- list.files( path = "./", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
# Column 4 is equivalent to htseq-count option -s reverse dUTP library prep
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) )
ff <- gsub( "ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "./", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1



sampleTable=data.frame(samples=colnames(counts))

sampleTable$Condition=c(rep("Aero",3),rep("Ana",3),rep("puBA10_35",3),rep("Positive",2),rep("puBA7_10",3),rep("puBA10_10",3))

# Anaerobic - Aerobic

countsAeroAna=counts[,1:6]
sampleTableAeroAna=sampleTable[1:6,]

sampleTableAeroAna$Condition = factor(sampleTableAeroAna$Condition)


# Distribution Raw
boxplot(log2(1+countsAeroAna),las=2)


ddsAeroAna = DESeqDataSetFromMatrix(colData = sampleTableAeroAna, countData = countsAeroAna, design= ~ Condition)

ddsAeroAna$Condition <- relevel(ddsAeroAna$Condition, ref = "Aero")

ddsAeroAna = DESeq(ddsAeroAna) #Normalization

countsNormAeroAna=DESeq2::counts(ddsAeroAna,normalize=T)

boxplot(log2(countsNormAeroAna+1),las=2)


resAeroAna = results(ddsAeroAna) #DE analysis

coefAeroAna = resultsNames(ddsAeroAna)[2]

resLFCAeroAna = lfcShrink(ddsAeroAna, coef = coefAeroAna, type = 'apeglm')

sampleDists = dist(t(countsNormAeroAna)) #Calculate distances
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = colnames(countsNormAeroAna)
colnames(sampleDistMatrix) = NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main = 'sample-to-sample distance',
         show_colnames = T,
         cex = 1.5)



plotPCA(normTransform(ddsAeroAna), intgroup= "Condition")

EnhancedVolcano(resAeroAna,lab = rownames(resAeroAna),  
                x = 'log2FoldChange',
                y = 'padj')


#gff=read_gff("JJZ00.gff")


# Crispr

countsCrispr=counts[,-1:-6]
sampleTableCrispr=sampleTable[-1:-6,]

sampleTableCrispr$Condition = factor(sampleTableCrispr$Condition)


# Distribution Raw
boxplot(log2(1+countsCrispr),las=2)


ddsCrispr = DESeqDataSetFromMatrix(colData = sampleTableCrispr, countData = countsCrispr, design= ~ Condition)

ddsCrispr$Condition <- relevel(ddsCrispr$Condition, ref = "Positive")

ddsCrispr = DESeq(ddsCrispr) #Normalization

countsNormCrispr=DESeq2::counts(ddsCrispr,normalize=T)

boxplot(log2(countsNormCrispr+1),las=2)


resCrispr = results(ddsCrispr) #DE analysis


sampleDists = dist(t(countsNormCrispr)) #Calculate distances
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = colnames(countsNormCrispr)
colnames(sampleDistMatrix) = NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main = 'sample-to-sample distance',
         show_colnames = T,
         cex = 1.5)



plotPCA(normTransform(ddsCrispr), intgroup= "Condition")

dataResultsFinal=data.frame()
# puBA10_10_vs_Positive
coefCrispr = resultsNames(ddsCrispr)[2]

resLFCCrispr = lfcShrink(ddsCrispr, coef = coefCrispr, type = 'apeglm')

EnhancedVolcano(resCrispr,lab = rownames(resCrispr),  
                x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1,pCutoff = 0.05)

index=which(abs(resLFCCrispr$log2FoldChange)>1 & resLFCCrispr$padj<0.05)
dataResults=data.frame(resLFCCrispr[index,c("log2FoldChange")])

dataResults$ID=rownames(dataResults)
dataResults$Condition="puBA10_10"
dataResultsFinal=rbind(dataResultsFinal,dataResults)
# puBA10_35_vs_Positive

coefCrispr = resultsNames(ddsCrispr)[3]

resLFCCrispr = lfcShrink(ddsCrispr, coef = coefCrispr, type = 'apeglm')

index=which(abs(resLFCCrispr$log2FoldChange)>1 & resLFCCrispr$padj<0.05)
dataResults=data.frame(resLFCCrispr[index,c("log2FoldChange")])

dataResults$ID=rownames(dataResults)
dataResults$Condition="puBA10_35"

dataResultsFinal=rbind(dataResultsFinal,dataResults)

EnhancedVolcano(resCrispr,lab = rownames(resCrispr),  
                x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1,pCutoff = 0.05)

# puBA7_10_vs_Positive
coefCrispr = resultsNames(ddsCrispr)[4]

resLFCCrispr = lfcShrink(ddsCrispr, coef = coefCrispr, type = 'apeglm')
index=which(abs(resLFCCrispr$log2FoldChange)>1 & resLFCCrispr$padj<0.05)
dataResults=data.frame(resLFCCrispr[index,c("log2FoldChange")])
dataResults$ID=rownames(dataResults)
dataResults$Condition="puBA7_10"

dataResultsFinal=rbind(dataResultsFinal,dataResults)


EnhancedVolcano(resCrispr,lab = rownames(resCrispr),  
                x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1,pCutoff = 0.05)



