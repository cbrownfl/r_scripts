#from input
counts = commandArgs(trailingOnly=TRUE)[1]
num_reads = commandArgs(trailingOnly=TRUE)[2]
cases = commandArgs(trailingOnly=TRUE)[3] # the number of cases in the dataset

#read in data
counts <- read.delim(file=counts, sep="\t", header=T, row.names=1)
num_reads <- read.delim(file=num_reads, sep=" ", header=F, row.names=1)
num_reads <- as.matrix(num_reads)
counts <- as.matrix(counts)
annotation <- rownames(counts)

#count samples and tests
num_tests <- length(counts[,1])
num_samples <- length(counts[1,])
cases <- as.numeric(cases)
controls <- num_samples - cases
first_cont <- cases + 1

#make matricies
case_average_count <- matrix(0, ncol=1, nrow=num_tests)
cont_average_count <- matrix(0, ncol=1, nrow=num_tests)
case_average_pertotal <- matrix(0, ncol=1, nrow=num_tests)
cont_average_pertotal <- matrix(0, ncol=1, nrow=num_tests)

#get percent of total reads
pertotal <- matrix(0, nrow=num_tests, ncol=num_samples)
for (i in 1:num_tests)
{
	for (j in 1:num_samples)
	{
		pertotal[i:i, j:j] <- (counts[i:i, j:j]/num_reads[j:j, 1])*100
	}
}

#get the average number of reads for cases and controls
avg_tot_cases <- mean(num_reads[1:cases, 1])
avg_tot_cont <- mean(num_reads[first_cont:num_samples, 1])

#plot the total number of reads for cases and controls and their averages
plot(c(1:num_samples), num_reads[1:num_samples, 1], type="l", ylab="number of reads", xlab="sample", main="total number of reads for cases and controls")
	abline(h=avg_tot_cases, col="red")
	abline(h=avg_tot_cont, col="blue")

#plot all of the samples
for (i in 1:num_tests)
{
	#get the average for cases and controls for both the counts and percent of total reads
	case_average_count[i:i, 1] <- mean(counts[i:i, 1:cases])
	case_average_pertotal[i:i, 1] <- mean(pertotal[i:i, 1:cases])
	cont_average_count[i:i, 1] <- mean(counts[i:i, first_cont:num_samples])
	cont_average_pertotal[i:i, 1] <- mean(pertotal[i:i, first_cont:num_samples])

	#plotting magic
	#plot total number of reads for cases and controls with a black line
	#title("main="counts=dashed, percent_total=solid", sub="total_reads= black, cases=red, controls=blue")
	#plot(c(1:num_samples), num_reads[1:num_samples, 1], type="l", yaxt="n", ylab="relative number of reads or percent of total reads", xlab="sample")
		#abline(h=avg_tot_cases, col="black", lty=20)
		#abline(h=avg_tot_cont, col="black")
		#par(new=TRUE)
	#plot counts with a dotted line
	plot(counts[i:i, 1:num_samples], type="n", xaxt="n", yaxt="n", ylab="", xlab="")
		lines(counts[i:i, 1:cases], col="red", type="l", lty=20, xaxt="n", yaxt="n", ylab="", xlab="")
		lines(c(first_cont:num_samples), counts[i:i, first_cont:num_samples], col="blue", type="l", lty=20, xaxt="n", yaxt="n", ylab="", xlab="")
		#abline(h=case_average_count[i:i, 1], col="red", lty=20)
		#abline(h=cont_average_count[i:i, 1], col="blue", lty=20)
		par(new=TRUE)
	#plot percent of total reads with a solid line
	plot(pertotal[i:i, 1:num_samples], main=annotation[i:i], cex.main=.5, type="n", xaxt="n", yaxt="n", ylab="", xlab="")
		lines(pertotal[i:i, 1:cases], col="red", type="l", xaxt="n", yaxt="n", ylab="", xlab="")
		lines(c(first_cont:num_samples), pertotal[i:i, first_cont:num_samples], col="blue", type="l", xaxt="n", yaxt="n", ylab="", xlab="")
		abline(h=case_average_pertotal[i:i, 1], col="red")
		abline(h=cont_average_pertotal[i:i, 1], col="blue")
}