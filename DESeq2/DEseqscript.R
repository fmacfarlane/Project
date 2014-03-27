#Install packages required to run the analysis
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


#Read in the data from the sequence count file
geneCounts<-read.delim("~/GitHub/project/DESeq2/RNAseqCounts.txt",head=T,sep="\t",skip=1, row.names=1)

nonZeroCounts<-geneCounts[rowSums(geneCounts[,6:28])>0,6:28]
#run background correction and normalisation on the data
treatments <- as.factor(str_sub(colnames(nonZeroCounts),-9,-5))

dds <- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
                              as.data.frame(treatments),
                              design=~treatments)
#relevel the data
dds$treatments <- relevel(dds$treatments, "minus" )
# run DESequencing on the data
dds <- DESeq(dds)

DEresults <- results(dds)
# get annotations and extra gene info from ensembl
require(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                host="www.ensembl.org", 
                path="/biomart/martservice")
#choose which database information has to come from
gg4 <- useDataset("ggallus_gene_ensembl",mart=mart)
#Download the relevant information from ensembl
annot <- getBM(attributes=c("ensembl_gene_id","external_gene_id","affy_chicken"),filter="ensembl_gene_id",values=rownames(DEresults),mart=gg4)
#merge the results dataframes
annotResults <- merge(annot,DEresults,by.x="ensembl_gene_id",by.y="row.names")
#sort the results by adjusted p value, lowest first
DEsorted <- head(DEresults[order(DEresults$padj),],desc=TRUE)
annotsorted <- head(annotResults[order(annotResults$padj),],desc=TRUE)#sorted for lowest padj
#view the results to find the lowest adjusted p value
View(annotsorted)