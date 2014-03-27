#Install and require the neccessary packages
install.packages("ape")
install.packages("phangorn")
install.packages("picante")
require(ape)
require(phangorn)
require(picante)
#Read in the data from the fasta file.
x<-read.dna("homologuesclustalalign.fna.mfa",format="fasta")
#Compare the number of substitutions seperating any pair of sequences(distance)
d<-dist.dna(x, model="raw")
#Write this into a viewable table
write.table(as.matrix(d),"distances1.csv")
#Plotting trees from distance matrices
#
#choose method of construction(bionj in this case)
tr.bionj<-bionj(d)
#plot the tree
plot(tr.bionj)
rownames(x)
#root the tree at the chicken node
tr.bionjr<-root(tr.bionj,outgroup="Chicken", resolve.root=TRUE )
#plot the rooted tree
plot(tr.bionjr);add.scale.bar(length=0.001)
#
#
#BOOTSTRAPPING
#
#Fit the tree to the data
fit<-pml(tr.bionj,as.phyDat(x))
#optimise plot and set random seed for bootstrap process
fit=optim.pml(fit,T)
plot(fit)
set.seed(8)
#bootstrap the data and plot the results
bs<-bootstrap.pml(fit,bs=100,optNni=T)
treeBS<-plotBS(fit$tree, type="p", bs)
#calculate the evolutionary distinctiveness scores
orig<-evol.distinct(tr.bionjr,type="fair.proportion")
#view ed score results table
orig
#plot a histogram of ed score frequencies
hist(orig$w, main='ED scores', xlab='ED Score')
