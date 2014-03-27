README.md file for the Comparative genomics folder.
* The .fna files are the fasta sequences downloaded from ENSEMBL to be used in the phylogenetic analysis of the gene. Those that also read_region, were the genome regions compared with the comparative genomic techniques.
* The distance1.csv displays the distance matrix calculated within the tree_plot.R file.
* THe file homologuesclustalalign.fna.mfa contains the jalview alignment of the homologues genes.
* THe R script tree_plot.R was used to plot various trees using the information from the jalview alignment of the homolgues, the evolutionary distinctiveness scores were also calculated.
The folder wga_output hold the output files from running ACT,MUMMER and BlAST.
  - The act.png file shows a visualisation obtained once act had been run using the blastn data.
  - The blastn.tab holds the output from running BLASTN on the two genome region sequences.
  - The mummer.png is the plot obtained from running mummer on the to genome regions.
  -all other files with the prefix mummer are outoutresults from the mummer system.