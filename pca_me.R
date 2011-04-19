#get this stuff from input
file = commandArgs(trailingOnly=TRUE)[1]
num_cases = commandArgs(trailingOnly=TRUE)[2]

#load the library and data
library(scatterplot3d)
data <- read.delim(file=file, sep="\t", header=T)
col <- length(data[1,]) - 2
num_cases <- as.numeric(num_cases)
num_controls <- col - num_cases

#do pca
pca <- prcomp(data[, 1:col], scale=TRUE, center=TRUE, retx=TRUE)
plot <- as.matrix(pca$rotation)

#setup colors
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

#plot pca
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
