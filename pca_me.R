#get this stuff from input
file = commandArgs(trailingOnly=TRUE)[1] #data file with the counts
samples = commandArgs(trailingOnly=TRUE)[2] #tab delim file with header and the number of samples for each test. Or "none."
num_cases = commandArgs(trailingOnly=TRUE)[3] #if there are cases and controls, how many cases are there? If not, 0.
ttest = commandArgs(trailingOnly=TRUE)[4] #are there pvalues and qvalues in the input table? yes or no.
trans = commandArgs(trailingOnly=TRUE)[5] #do you want to log transform the data before doing pca? yes or no.
threeD = commandArgs(trailingOnly=TRUE)[6] #do you want to plot 3D graphs?  yes or no.
print_table = commandArgs(trailingOnly=TRUE)[7] #do you want to print a table with the rotation coordinates? yes or no

#load the library and data
library(scatterplot3d)
data <- read.delim(file=file, sep="\t", header=T, row.names=1)
if (num_cases==0)
{
	samples <- read.delim(file=samples, sep="\t", header=T, row.names=1)
	samples <- as.matrix(samples)
}
row <- length(data[,1])
if (ttest=="yes")
{
	col <- length(data[1,]) - 2
}
if (ttest=="no")
{
	col <- length(data[1,])
}

if (num_cases!=0)
{
	num_cases <- as.numeric(num_cases)
	num_controls <- col - num_cases
}


#log transform the data
if (trans=="yes")
{
	logtrans <- matrix(0, ncol=col, nrow=row)
	for (i in 1:row)
	{
		for (j in 1:col)
		{
			logtrans[i:i, j:j] <- log(data[i:i, j:j] + 1, 10)
		}
	}
	data <- logtrans
}

#do pca
pca <- prcomp(data[, 1:col], scale=TRUE, center=TRUE, retx=TRUE)
plot <- as.matrix(pca$rotation)

#setup colors
if (num_cases!=0)
{
	color <- matrix(0, ncol=1, nrow=col)
	cont_start <- num_cases+1
	for (i in 1:num_cases)
	{
		color[i:i, 1] <- "red"
	}
	for (i in cont_start:col)
	{
		color[i:i, 1] <- "blue"
	}
}
if (num_cases==0)
{
	color <- matrix(0, ncol=1, nrow=col)
	num_tests <- length(samples[,1])
	start <- 1
	count <- 0
	for (i in 1:num_tests)
	{
		samples_in_test <- as.numeric(samples[i:i, 1])
		end <- samples_in_test + count
		for (j in start:end)
		{
			color[j:j, 1] <- as.character(samples[i:i, 2])
		}
		start <- end + 1
		count <- count + samples_in_test
	}
}

#save the pca rotation data
if (print_table=="yes")
{
	write.table(pca$rotation, file="pca", sep="\t")
}


#plot pca

#2D
par(xpd=T, mar=par()$mar+c(0,0,0,7))
plot(plot[,1], plot[,2], col=color, xlab="PC1", ylab="PC2")
legend(.22, .25, row.names(samples), fill=samples[,2])

#3D
if (threeD=="yes")
{
	for (i in 1:col)
	{
		for (j in 1:col)
		{
			for (k in 1:col)
			{
				if (i!=j)
				{
					if (i!=k)
					{
						if (j!=k)
							{
								scatterplot3d(plot[,i], plot[,j], plot[,k], color=color, xlab=i, ylab=j, zlab=k)
							}
					}
				}
			}
		}
	}
}