library(dupRadar)

args = commandArgs(trailingOnly=TRUE)

file = args[1]
id = args[2]
gtf = args[3]
stranded = as.integer(args[4])
paired = as.logical(args[5])
threads = as.integer(args[6])
outdir = args[7]

dm = analyzeDuprates(file, gtf, stranded, paired, threads)


if(!all(is.nan(dm$dupRate))){

	bitmap(paste(outdir, '/', id, '.densplot.png', sep=""), type='png16m', height=1000, width=1000, unit='px', taa=4, gaa=4)
	duprateExpDensPlot(DupMat=dm, main=id)
	dev.off()

	bitmap(paste(outdir, '/', id, '.boxplot.png', sep=""), type='png16m', height=1000, width=1000, unit='px', taa=4, gaa=4)
	par(mar=c(10, 4, 4, 2) + 0.1)
	duprateExpBoxplot(DupMat=dm, main=id)
	dev.off()

	bitmap(paste(outdir, '/', id, '.plot.png', sep=""), type='png16m', height=1000, width=1000, unit='px', taa=4, gaa=4)
	duprateExpPlot(DupMat=dm, main=id)
	dev.off()

	bitmap(paste(outdir, '/', id, '.hist.png', sep=""), type='png16m', height=1000, width=1000, unit='px', taa=4, gaa=4)
	expressionHist(dm)
	dev.off()

} else {

	print("dupRate is NaN")

}
