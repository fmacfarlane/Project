#Limma script for project
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")
biocLite("chicken.db")
biocLite("affyPLM")
biocLite("affy")
biocLite("arrayQualityMetrics")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("GEOquery")
install.packages("affy")
install.packages("affyPLM")
install.packages("arrayQualityMetrics")
install.packages("annotate")
install.packages("chicken.db")
install.packages("chickenprobe")
install.packages("genefilter")
install.packages("limma")
#
library(ggplot2)
library(reshape2)
library(GEOquery)
library(affy)
library(affyPLM)
library(arrayQualityMetrics)
library(annotate)
library(chicken.db)
library(chickenprobe)
library(genefilter)
library(limma)

.ggboxIt <- function(dfin,pdfname)
{	# draw a ggplot boxplot of data frame produced by .summarise.it()
	gp <- ggplot(dfin, aes(x = name, ymin = lwisk, lower = lbox,
													 middle = mid, upper = ubox,
													 ymax = uwisk, fill = name))
  gp <- gp + geom_boxplot(stat="identity",show_guide=F) + coord_flip() 
	gp <- gp + theme(panel.background = element_blank())
	gp <- gp + labs(y = "Log2(Intensity)", x = "")
  ggsave(pdfname) 
}

.summariseIt <- function(dfraw,phenoData)
{ # create quantile data from a dataframe of depths
  # useful for munging into ggplot, however probably a plyr function for this too
	dfin <- melt(dfraw)
	colnames(dfin) <- c("probeId","sample","value")
	dfin$group<-factor(dfin$sample)
  stats <- as.data.frame(t(sapply(levels(dfin$sample),
            function(x){
              quantile(log2(as.numeric(dfin$value))[which(dfin$group==x)],
											 prob=c(0.05,0.25,0.5,0.75,0.95),na.rm=T)
            })))
	colnames(stats)<-c("lwisk","lbox","mid","ubox","uwisk")
  stats$name <- rownames(stats)
  rownames(stats) <- NULL
 	return(stats)
}

.getData <- function()
{
	baseDir <- "~/GitHub/project"
	workingDir=paste0(baseDir)


	filenames <- c("ROS1-_9.CEL", "ROS1+_10.CEL", "ROS2-_11.CEL", "ROS2+_12.CEL", "ROS3-_13.CEL", "ROS3+_14.CEL", "ROS4-_15.CEL", "ROS4+_16.CEL")
	samplenames <- c("ROS1-_9", "ROS1+_10", "ROS2-_11", "ROS2+_12", "ROS3-_13", "ROS3+_14", "ROS4-_15", "ROS4+_16")
	targets <- c("minus", "plus", "minus", "plus", "minus", "plus", "minus", "plus")

	phenodata<-as.data.frame(cbind(filenames,samplenames,targets))
	write.table(phenodata,paste(workingDir,"phenodata.txt",sep="/")
							,quote=F,row.name=F)
	celRAW <- ReadAffy(celfile.path=workingDir,compress=T,
										 phenoData=phenodata)
	#arrayQualityMetrics(expressionset=celRAW,
	#										outdir=paste0(baseDir,"celRAW_AQM"),
	#										force=T,do.logtransform=T)
}

.plotDensity <- function(exps,filename)
{
	pdf(filename)
	# Plot a density vs log intensity histogram for the unnormalised data
	d <- apply(exps,2,
				function(x){
					density(x)
				})
	xmax <- max(sapply(d,function(x)max(x$x)))
	xmin <- min(sapply(d,function(x)min(x$x)))
	ymax <- max(sapply(d,function(x)max(x$y)))
	plot(0,pch='',ylab='',xlab='',
			 xlim=c(xmin,round(xmax+1)),ylim=c(0,ymax))
	lapply(1:length(d),function(x) lines(d[[x]],col=x))
	dev.off()
}

# Perform probe-level metric calculations on the CEL files:
.doPLM <- function(celRAW)
{
	pdf("celRAWqc.pdf")
	#affyPLM is required to interrogate celRMA
	celRAWqc <- fitPLM(celRAW)

	# Create an image of GSM24662.CEL:
	image(celRAWqc, which=1, add.legend=TRUE)


	image(celRAWqc, which=4, add.legend=TRUE)

	# affyPLM also provides more informative boxplots
	# RLE (Relative Log Expression) plots should have

	RLE(celRAWqc, main="RLE")

	# We can also use NUSE (Normalised Unscaled Standard Errors).
	# The median standard error should be 1 for most genes.
	
	NUSE(celRAWqc, main="NUSE")
	dev.off()
}

.doCluster <- function(celRMA)
{
#Quick and dirty clustering much better use pvclust
#Clustering like normalization is a bit of a black art
	eset <- exprs(celRMA)
	distance <- dist(t(eset),method="maximum")
	clusters <- hclust(distance)
	plot(clusters)
}

.doFilter <- function(celRMA)
{
	celfiles.filtered <- nsFilter(celRMA, 
															require.entrez=FALSE, 
															remove.dupEntrez=FALSE, )
}


.doDE <- function(eset)
{
	samples <- eset$targets
# check the results of this
# convert into factors
	samples <- as.factor(samples)
# set up the experimental desi
	design <- model.matrix(~0 + samples)
	colnames(design) <- c("plus","minus")

# fit the linear model to the filtered expression set
	fit <- lmFit(exprs(eset), design)


# Now the contrast matrix is combined with the per-probeset linear model fit.

	plus_ebFit <- eBayes(fit)
  

#	ttab <- topTable(plus_ebFit, number=30000, coef=1, p.value=0.05)

#	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=5, p.value=0.05))
#	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=4, p.value=0.05))
#	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=3, p.value=0.05))
#	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=2, p.value=0.05))
# Get a list for probesets with a four fold change or more
	tTable <- topTable(plus_ebFit, number=30000, p.value=0.05)
 # return(tTable)

	annotation <- as.data.frame(select(chicken.db,
																		 rownames(tTable), 
																		 c("ENSEMBL","SYMBOL")))
	colnames(annotation) <- c("probeId","ensemblId","geneSymbol")
	#results <- merge(annotation, tTable,by.x="probeId",by.y="row.names")
	results <- merge(annotation[1:100,], tTable,by.x="probeId",by.y="row.names")
  head(results)
	write.table(results, "results.txt", sep="\t", quote=FALSE)
	return(results)
}

if(!exists("celResults"))
{
	celRAW <- .getData()
#	eset<-exprs(celRAW)
	celRMA <- rma(celRAW)
#	.ggboxIt(.summariseIt(log2(exprs(celRAW))),"sumRAW.pdf")
#	.ggboxIt(.summariseIt(exprs(celRMA)),"sumRMA.pdf")
#	.plotDensity(log2(exprs(celRAW)),"densityRAW.pdf")
#	.plotDensity(log2(exprs(celRMA)),"densityRMA.pdf")
  celFilt <- .doFilter(celRMA)
	celResults <- .doDE(celFilt$eset)

}


