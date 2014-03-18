require(DESeq2)
require(RCurl)
require(biomaRt)

#remoteFile <- getURL(paste0("http://ts-ug-dev.lifesci.dundee.ac.uk/usr/local/share/BS32010-data/expression/data/RNAseqCounts.txt"))

geneCounts<-read.delim("/homes/fmacfarlane/bioinformatics/RNAseqCounts.txt",head=T,sep="\t",skip=1)

nonZeroCounts<-geneCounts[rowSums(geneCounts[,6:28])>0,6:28]

treatments <- as.factor(substr(colnames(nonZeroCounts),1,1))

dds <- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
                              as.data.frame(treatments),
                              design=~treatments)

dds$treatments <- relevel(dds$treatments,"L" )

dds <- DESeq(dds)

DEresults <- results(dds)

# get annotations
require(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                host="www.ensembl.org", 
                path="/biomart/martservice")

gg4 <- useDataset("ggallus_gene_ensembl",mart=mart)

annot <- getBM(attributes=c("ensembl_gene_id","external_gene_id","affy_chicken"),filter="ensembl_gene_id",values=rownames(DEresults),mart=gg4)

annotResults <- merge(annot,DEresults,by.x="ensembl_gene_id",by.y="row.names")

