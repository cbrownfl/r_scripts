# get from input 
data <- commandArgs(trailingOnly=TRUE)[1] # get the file name with the counts
mode <- commandArgs(trailingOnly=TRUE)[2] # do you want to use the raw counts or the percent of total reads?
num_reads <- commandArgs(trailingOnly=TRUE)[3] # file with the number of reads

# load the files into R
library(lme4)
library(qvalue)
data <- read.delim(data, header=T, sep="\t", row.names=1)
names <- row.names(data)
attach(data)

# define some variables 
row <- length(data[,1]) # number of functions or OTUs
col <- length(data[1,]) # number of samples
treatment <- gl(2,4,col) # which samples are cases and controls ? 
outcome <- gl(4,1,col) # which samples are paired ?

# make matricies 
counts <- matrix(0, nrow=row, ncol=col)	
pvalue <- array(0, dim=row)
pertotal <- matrix(0, nrow=row, ncol=col)
qvalue <- matrix(0, nrow=row, ncol=1)

# calcualte the percent of total reads
if (mode=="pertotal")
{
	num_reads <- read.delim(num_reads, header=F, sep=" ", row.names=1) # load the file that has the number of reads for each sample
	for (i in 1:row)
	{
		for (j in 1:col)
		{	
			pertotal[i:i, j:j] <- (data[i:i, j:j]/num_reads[j:j, 1])*100
		}
	} 
    data <- pertotal
}

# poisson!
for (i in 1: row)
{	
	counts[i,] <- t(data[i,])
	model <- glmer(counts[i,] ~ treatment +(0+treatment|outcome),family=poisson())
	pvalue[i] <- summary(model)@coefs[2,4]
}

# get the qvalues from the pvalues
qvalue <- qvalue(pvalue)$qvalues

#print out the data and pvalues
print <- data.frame(names, data, pvalue, qvalue)
write.table(print, file="out", row.name=FALSE, quote=FALSE, sep="\t")