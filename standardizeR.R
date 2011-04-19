#input to script
file = commandArgs(trailingOnly=TRUE)[1]
num_reads = commandArgs(trailingOnly=TRUE)[2]
num_cases = commandArgs(trailingOnly=TRUE)[3]
method = commandArgs(trailingOnly=TRUE)[4]
ttest = commandArgs(trailingOnly=TRUE)[5]
library(qvalue)

#read data
data <- read.delim(file=file, sep="\t", header=T, row.names=NULL)
header <- read.delim(file=file, sep="\t", header=F, row.names=NULL)
reads <- read.delim(file=num_reads, sep=" ", header=F, row.names=1)
num_cases <- as.numeric(num_cases)
first_cont <- num_cases+1

#format the data if there are annotations and a header
col <- length(data[1,])
row <- length(data[,1])
top <- matrix(0, nrow=1, ncol=col+1)
top[1:1, 1:col-1] <- header[1:1,2:col]
annotations <- data[,1:1]
data <- data[1:row,2:col]
col <- length(data[1,])
row <- length(data[,1])

#make matricies
pertotal <- matrix(0, nrow=row, ncol=col)
pertotal_adj <- matrix(0, nrow=row, ncol=col)
trans <- matrix(0, nrow=row, ncol=col)
stand <- matrix(0, nrow=row, ncol=col)
sdev <- matrix(0, nrow=row, ncol=col)
sdev_case <- matrix(0, nrow=row, ncol=col)
average <- matrix(0, nrow=row, ncol=col)
average_case <- matrix(0, nrow=row, ncol=col)
sdev_cont <- matrix(0, nrow=row, ncol=col)
average_cont <- matrix(0, nrow=row, ncol=col)

#log transform

## method 0: Log transform the percent of total reads.  Calculate standard deviation and average for each sample, separetely, across columns.
if (method==0)
{
    for (i in 1:row)
    {
	    for (j in 1:col)
	    {	
		    pertotal[i:i, j:j] <- ((data[i:i, j:j])/(reads[j,1]))*100
		    pertotal_adj[i:i, j:j] <- pertotal[i:i, j:j]*1000000
			trans[i:i, j:j] <- log(pertotal_adj[i:i, j:j]+1, 2)
		    sdev[1, j:j] <- sd(trans[1:row, j:j])
		    average[1, j:j] <- mean(trans[1:row, j:j])
	    }
    }
}

## method 1: Log transform the percent of total reads. Calculate standard deviation and average for cases and controls separately, acrosss rows.
if (method==1)
{
    for (i in 1:row)
    {
	    for (j in 1:col)
	    {	
		    pertotal[i:i, j:j] <- ((data[i:i, j:j])/(reads[j,1]))*100
		    trans[i:i, j:j] <- log(pertotal[i:i, j:j]+1, 10)
		    sdev_case[i,1] <- sd(trans[i:i, 1:num_cases])
		    average_case[i,1] <- mean(trans[i:i, 1:num_cases])
			sdev_cont[i,1] <- sd(trans[i:i, first_cont:col])
		    average_cont[i,1] <- mean(trans[i:i, first_cont:col])
	    }
    }
}

## method 2: Log transform the percent of total reads.  Calculate standard deviation and average for cases and controls together, across rows.
if (method==2)
{
    for (i in 1:row)
    {
	    for (j in 1:col)
	    {	
		    pertotal[i:i, j:j] <- ((data[i:i, j:j])/(reads[j,1]))*100
		    trans[i:i, j:j] <- log(pertotal[i:i, j:j]+1, 10)
		    sdev[i,1] <- sd(trans[i:i, 1:col])
		    average[i,1] <- mean(trans[i:i, 1:col])
	    }
    }
}

## method 3: Log transform the raw counts.  Calculate standard deviation and average for cases and controls, separately, across columns.  
if (method==3)
{
    for (i in 1:row)
    {
	    for (j in 1:col)
	    {	
		    trans[i:i, j:j] <- log(data[i:i, j:j]+1, 2)
		    sdev[1, j:j] <- sd(trans[1:row, j:j])
		    average[1, j:j] <- mean(trans[1:row, j:j])
	    }
    }
}

## method 4: Log transform the raw counts.  Calculate standard deviation and average for cases and controls, together, across rows.  
if (method==4)
{
    for (i in 1:row)
    {
	    for (j in 1:col)
	    {	
		    trans[i:i, j:j] <- log(data[i:i, j:j]+1, 10)
		    sdev[i,1] <- sd(trans[i:i, 1:col])
		    average[i,1] <- mean(trans[i:i, 1:col])
	    }
    }
}

##method 5: Don't bother doing the log transformation or standardization! Just do the t-test on percent of total reads!
if (method==5)
{
    for (i in 1:row)
    {
	    for (j in 1:col)
	    {	
		    pertotal[i:i, j:j] <- ((data[i:i, j:j])/(reads[j,1]))*100
		}
	}
}

#standardize

##method 0
if (method==0)
{
	for (i in 1:row)
	{
		for (j in 1:col)
		{
			stand[i:i, j:j] <- (trans[i:i, j:j] - average[1, j:j])/(sdev[1, j:j])
		}
	}
}

##method 1
if (method==1)
{
	for (i in 1:row)
	{
		for (j in 1:num_cases)
		{
			stand[i:i,j:j] <- (trans[i:i,j:j] - average_case[i:i, 1])/(sdev_case[i:i, 1])
		}
	}
	for (i in 1:row)
	{
		for (j in first_cont:col)
		{
			stand[i:i,j:j] <- (trans[i:i,j:j] - average_cont[i:i, 1])/(sdev_cont[i:i, 1])
		}
	}
}

##method 2
if (method==2)
{
	for (i in 1:row)
	{
		for (j in 1:col)
		{
			stand[i:i,j:j] <- (trans[i:i,j:j] - average[i:i, 1])/(sdev[i:i, 1])
		}
	}
}

##method 3
if (method==3)
{
	for (i in 1:row)
	{
		for (j in 1:col)
		{
			stand[i:i, j:j] <- (trans[i:i, j:j] - average[1, j:j])/(sdev[1, j:j])
		}
	}
}

##method 4
if (method==4)
{
	for (i in 1:row)
	{
		for (j in 1:col)
		{
			stand[i:i,j:j] <- (trans[i:i,j:j] - average[i:i, 1])/(sdev[i:i, 1])
		}
	}
}

##method 5
if (method==5)
{
			stand <- pertotal
}

#t-test
if (ttest=="paired")
{
	ttest <- matrix(0, nrow=row, ncol=1)
	for (i in 1:row)
	{
		ttest[i:i, 1] <- t.test(stand[i:i, 1:num_cases], stand[i:i, first_cont:col], paired=TRUE)$p.value
	}
}
if (ttest=="notpaired")
{
	ttest <- matrix(0, nrow=row, ncol=1)
	for (i in 1:row)
	{
		ttest[i:i, 1] <- t.test(stand[i:i, 1:num_cases], stand[i:i, first_cont:col], paired=FALSE)$p.value
	}
}

#determine if p and q values need to be plotted and saved

# if doing a ttest
if (ttest!="none")
{
	#make more matricies 
	qvalues <- matrix(0, nrow=row, ncol=1)
	qval <- qvalue(ttest[,1])$qvalues
	qvalues <- as.matrix(qval)

	#make the table
	ttestc <- col+1
	qvalc <- col+2
	table <- matrix(0, nrow=row, ncol=qvalc)
	table[,1:col] <- as.matrix(stand)
	table[, ttestc] <- ttest
	table[, qvalc] <- qvalues
	print <- data.frame(table, row.names=annotations)
	print <- print[order(print[,ttestc], print[,qvalc]),]
	
	#plot some data
	boxplot(data, main="raw counts")
	boxplot(pertotal, main="raw percent of total reads")
	boxplot(trans, main="log transformed data")
	boxplot(stand[,1:col], main="standardized data")
	hist(ttest, main="t-test: p values", breaks=50)
	hist(qvalues, main="q values", breaks=50)
	for (i in 1:row)
	{
		plot(stand[i:i,1:num_cases], col="red", type="l"); lines(stand[i:i, first_cont:col], col="blue", type="l")
	}
}

#if not doing a ttest
if (ttest=="none")
{	
	#make the table
	table <- matrix(0, nrow=row, ncol=col)
	table[,1:col] <- as.matrix(stand)
	print <- data.frame(table, row.names=annotations)
	
	#plot some data
	boxplot(data, main="raw counts")
	boxplot(pertotal, main="raw percent of total reads")
	boxplot(trans, main="log transformed data")
	boxplot(stand[,1:col], main="standardized data")
	for (i in 1:row)
	{
		plot(stand[i:i, 1:col], type="n"); lines(stand[i:i, 1:num_cases], col="red", type="l"); lines(stand[i:i, first_cont:col], col="blue", type="l")
	}
}

#print the table
write.table(print, file="out", row.name=TRUE, quote=FALSE, sep="\t")
#write.table(data.frame(pertotal), file="pertotal", row.names=TRUE, quote=FALSE, sep="\t")