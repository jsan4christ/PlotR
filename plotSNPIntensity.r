##########################################################################################################################
# plotSNPIntensity.r - Matt Ravenhall											 #
# Provided an array of SNPs has been provided at the command line as shown below, produces SNP intensity plots for each. #
# e.g. 'plotSNPIntensity.r CHR:BP' will produce Case, Control & Combined SNP intensity plots.				 #
# Note that this script assumed out.XY.format exists within the working directory and contains XY data.			 #
##########################################################################################################################

# Warning: This script should now take CHR:BP instead of BP input. Testing now to see if that works.

# Optional
xlim.fixed=2
ylim.fixed=2

# Grab Position (in BP) of SNP
args <- commandArgs(trailingOnly=TRUE)

# Bring in all SNPs
print('Parsing XY Data')
dat <- read.table('out.XY.FORMAT',header=T)

# Loop SNPs, currently probably not the most efficient approach but it'll work
for (eachsnp in args) {
	# Use i for BP, j for CHR
	i <- strsplit(eachsnp, ':')[[1]][2]
	j <- strsplit(eachsnp, ':')[[1]][1]

	# Isolate SNP of interest
	print(paste('Isolating XY data for SNP on CHR ',j,' at BP ',i,sep=''))
#	snp <- data.frame(dat[dat$CHR == j,]))
#	snp <- data.frame(t(dat[dat$POS == i,]))
	snp <- data.frame(t(dat[dat$POS == i & dat$CHR == j,]))
	print('Removing CHROM and POS data from isolated subset')
	snp <- subset(snp, row.names(snp) != 'CHROM' & row.names(snp) != 'POS')

	# Split XY columns
	names(snp) <- c('XY')
	print('Splitting XY data into appropriate columns')
	snp$X <- lapply(strsplit(as.character(snp$XY), ','), '[', 1)
	snp$Y <- lapply(strsplit(as.character(snp$XY), ','), '[', 2)

	# Add colours
	print('Parsing GT data (for point colours)')
	cols <- read.table('out.GT.FORMAT',header=T)
	print(paste('Isolating GT data for SNP on CHR ',j' at BP ',i,sep=''))
	snpcol <- data.frame(t(cols[cols$POS == i,]))
	print('Removing CHROM and POS data from isolated subset')
	snpcol <- subset(snpcol, row.names(snpcol) != 'CHROM' & row.names(snpcol) != 'POS')

	print('Splitting GT data into appropriate columns')
	names(snpcol) <- c('XY')
	snpcol$X <- lapply(strsplit(as.character(snpcol$XY), '/'), '[', 1)
	snpcol$Y <- lapply(strsplit(as.character(snpcol$XY), '/'), '[', 2)
	print('Parsing GT to single numeric format')
	snpcol$colID <- as.numeric(snpcol$X) + as.numeric(snpcol$Y)
	print('Accounting for missing GT calls')
	snpcol$colID[is.na(snpcol$colID)] <- 4

	# Identify IDs as Case/Control
	print('Parsing phenotypic data')
	pheno <- read.table('../RawData/Tanzania_GWAS-2.5M_b37_915.sample',header=T)
	print('Identifying Case & Control individuals')
	CaseIDs <- pheno$scanID[pheno$caseorcontrol=='CASE']
	ControlIDs <- pheno$scanID[pheno$caseorcontrol=='CONTROL']

	# Plot to file
	print(paste('Plotting figure as SNPIntensity_BP',i,'.png',sep=''))
	png(paste('SNPIntensity_',j,':',i,'_Combined.png',sep=''))
	plot(snp$X, snp$Y, xlim=c(0,xlim.fixed), ylim=c(0,ylim.fixed), col=rainbow(4)[1 + as.numeric(snpcol$colID)], pch=16, main=paste(j,':',i,' Intensity (Case+Control)',sep=''), xlab='X', ylab='Y')
	dev.off()

	print(paste('Plotting figure as SNPIntensity_BP',i,'_Case.png',sep=''))
	png(paste('SNPIntensity_',j,':',i,'_Case.png',sep=''))
	plot(snp$X[CaseIDs], snp$Y[CaseIDs], xlim=c(0,xlim.fixed), ylim=c(0,ylim.fixed), col=rainbow(4)[1 + as.numeric(snpcol$colID[CaseIDs])], pch=16, main=paste(j,':',i,' Intensity (Case Only)',sep=''), xlab='X', ylab='Y')
	dev.off()

	print(paste('Plotting figure as SNPIntensity_BP',i,'_Control.png',sep=''))
	png(paste('SNPIntensity_BP',j,':',i,'_Control.png',sep=''))
	plot(snp$X[ControlIDs], snp$Y[ControlIDs], xlim=c(0,xlim.fixed), ylim=c(0,ylim.fixed), col=rainbow(4)[1 + as.numeric(snpcol$colID[ControlIDs])], pch=16, main=paste(j,':',i,' Intensity (Control Only)',sep=''), xlab='X', ylab='Y')
	dev.off()

	print('Done')
}
