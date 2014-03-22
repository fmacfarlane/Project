install.packages("DESeq2")
install.packages("RCurl")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("biomaRt")
install.packages("stringr")
require(DESeq2)
require(RCurl)
require(biomaRt)
require(stringr)
require(limma)


geneCounts<-read.delim("~/GitHub/project/RNAseqCounts.txt",head=T,sep="\t",skip=1, row.names=1)

nonZeroCounts<-geneCounts[rowSums(geneCounts[,6:28])>0,6:28]

treatments <- as.factor(str_sub(colnames(nonZeroCounts),-7,-7))

dds <- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
                              as.data.frame(treatments),
                              design=~treatments)

dds$treatments <- relevel(dds$treatments, "l" )

dds <- DESeq(dds)

DEresults <- results(dds)
#get results for padj
res <- DEresults[order(DEresults$padj),]
head(res)
# get annotations
require(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                host="www.ensembl.org", 
                path="/biomart/martservice")

gg4 <- useDataset("ggallus_gene_ensembl",mart=mart)

annot <- getBM(attributes=c("ensembl_gene_id","external_gene_id","affy_chicken"),filter="ensembl_gene_id",values=rownames(DEresults),mart=gg4)

annotResults <- merge(annot,DEresults,by.x="ensembl_gene_id",by.y="row.names")
colData(dds)


#heat maps + plots
#library("RColorBrewer")
#library("gplots")
#install.packages("gplots")
#rld <- rlogTransformation(dds)
#vsd <- varianceStabilizingTransformation(dds)
#select <-order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:30]
#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
#          Rowv = FALSE, Colv = FALSE, scale="none",
#          dendrogram="none", trace="none", margin=c(10,6))
#heatmap.2(assay(rld)[select,], col = hmcol,
#          Rowv = FALSE, Colv = FALSE, scale="none",
#          dendrogram="none", trace="none", margin=c(10, 6))
#heatmap.2(assay(vsd)[select,], col = hmcol,
#          Rowv = FALSE, Colv = FALSE, scale="none",
#          dendrogram="none", trace="none", margin=c(10, 6))
