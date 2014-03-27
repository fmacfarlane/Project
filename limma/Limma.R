#Limma script for project
#install nd require neccessary packages
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
biocLite("chickenprobe")
biocLite("genefilter")
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
library(BiocGenerics)
library(genefilter)
library(limma)
#function to plot the summarised boxplots
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
#function to summarise the raw and normalised data, which can then be plotted using the ggplot function
.summariseIt <- function(dfraw,phenoData)
{
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
#Input data from ROS.Cel files in project/limma directory
.getData <- function()
{#set current directory as working directory
	baseDir <- "~/GitHub/project/limma"
	workingDir=paste0(baseDir)

#input filenames and samplenames
	filenames <- c("ROS1-_9.CEL", "ROS1+_10.CEL", "ROS2-_11.CEL", "ROS2+_12.CEL", "ROS3-_13.CEL", "ROS3+_14.CEL", "ROS4-_15.CEL", "ROS4+_16.CEL")
	samplenames <- c("ROS1-_9", "ROS1+_10", "ROS2-_11", "ROS2+_12", "ROS3-_13", "ROS3+_14", "ROS4-_15", "ROS4+_16")
	#set targete to correspond to treated(plus) or untreated(minus) samples
  targets <- c("minus", "plus", "minus", "plus", "minus", "plus", "minus", "plus")
# combine the above filenames,samplenames and targets into a dataframe and export table as a text file
	phenodata<-as.data.frame(cbind(filenames,samplenames,targets))
	write.table(phenodata,paste(workingDir,"phenodata.txt",sep="/")
							,quote=F,row.name=F)
	celRAW <- ReadAffy(celfile.path=workingDir,compress=T,
										 phenoData=phenodata)
}
#Function to plot the distribution of signal intensities for the raw and normalised data
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

	# Create an image of the ROS data
	image(celRAWqc, which=1, add.legend=TRUE)


	image(celRAWqc, which=4, add.legend=TRUE)

	# affyPLM also provides more informative boxplots


	RLE(celRAWqc, main="RLE")

	# We can also use NUSE (Normalised Unscaled Standard Error)
	
	NUSE(celRAWqc, main="NUSE")
	dev.off()
}

.doCluster <- function(celRMA)
{
#performs clustering on the normalised data
	eset <- exprs(celRMA)
	distance <- dist(t(eset),method="maximum")
	clusters <- hclust(distance)
	plot(clusters)#these are plotted as boxplots
}
#We then filter the data
.doFilter <- function(celRMA)
{
	celfiles.filtered <- nsFilter(celRMA, 
															require.entrez=FALSE, 
															remove.dupEntrez=FALSE, )
}


.doDE <- function(eset)
{
	samples <- eset$targets
# convert into factors
	samples <- as.factor(samples)
# set up the experimental design
	design <- model.matrix(~0 + samples)
	colnames(design) <- c("minus","plus")

# fit the linear model to the filtered expression set
	fit <- lmFit(exprs(eset), design)


# Now the contrast matrix is combined with the per-probeset linear model fit.

	plus_ebFit <- eBayes(fit)
  

	ttab <- topTable(plus_ebFit, number=30000, coef=1)

	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=5))
	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=4))
	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=3))
	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=2))
 #Get a list for probesets with a four fold change or more
	tTable <- topTable(plus_ebFit, number=30000, p.value=0.05)
 # sets the p-value limit to be 0.05
#then create table of results with relevant gene information taken from the chicken database
	annotation <- as.data.frame(select(chicken.db,
																		 rownames(tTable), 
																		 c("ENSEMBL","SYMBOL")))
	colnames(annotation) <- c("probeId","ensemblId","geneSymbol")
	results <- merge(annotation, tTable,by.x="probeId",by.y="row.names")
	#results <- merge(annotation[1:100,], tTable,by.x="probeId",by.y="row.names")
  head(results)
	write.table(results, "results.txt", sep="\t", quote=FALSE)
	return(results)
}
#set the if loop so that if celResults exist the limma analysis will take place
if(!exists("celResults"))
{
	celRAW <- .getData()
	eset<-exprs(celRAW)
	celRMA <- rma(celRAW)#this is the rma normalisation step
	.ggboxIt(.summariseIt(log2(exprs(celRAW))),"sumRAW.pdf")
	.ggboxIt(.summariseIt(exprs(celRMA)),"sumRMA.pdf")#plot the cluster summarised data boxplots
	.plotDensity(log2(exprs(celRAW)),"densityRAW.pdf")
	.plotDensity(log2(exprs(celRMA)),"densityRMA.pdf")#plot the distribution of intensitites
  celFilt <- .doFilter(celRMA)#filter the data
	celResults <- .doDE(celFilt$eset)

}
# the results can then be ordered by adjusted p value to find the most significant genes
sorted <- head(celResults[order(celResults$adj.P.Val),])#find lowest adjusted p values
