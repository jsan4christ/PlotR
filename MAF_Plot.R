plotMAFs = function(xname, yname, xfile, yfile, removeBA=FALSE) {
	# Input is the .frq file produced by calling plink --freq
	x = paste(xfile, sep="")
	y = paste(yfile, sep="")
	Con1 = read.table(x, header=T)
	Con2 = read.table(y, header=T)

	# Remove NAs
	Con1.noNA = Con1[complete.cases(Con1), ]
	Con2.noNA = Con2[complete.cases(Con2), ]

	# Remove 0s
	Con1.noNA = Con1.noNA[Con1.noNA$A1 != 0, ]
	Con1.noNA = Con1.noNA[Con1.noNA$A2 != 0, ]
	Con2.noNA = Con2.noNA[Con2.noNA$A1 != 0, ]
	Con2.noNA = Con2.noNA[Con2.noNA$A2 != 0, ]

	# Subset to SNP name and MAF only
	cordsCon1 = data.frame(Con1.noNA$SNP, Con1.noNA$MAF)
	cordsCon2 = data.frame(Con2.noNA$SNP, Con2.noNA$MAF)
	names(cordsCon1) = c("SNP", "Con1")
	names(cordsCon2) = c("SNP", "Con2")

	# Sort by SNP name (ensures easy merge later)
	cordsCon1.sorted = cordsCon1[order(cordsCon1$SNP),]
	cordsCon2.sorted = cordsCon2[order(cordsCon2$SNP),]
	names(cordsCon1.sorted) = c("SNP", "Con1")
	names(cordsCon2.sorted) = c("SNP", "Con2")

	# Remove non-shared SNPs (again, eases merge)
	cordsCon1.shared = cordsCon1.sorted[cordsCon1.sorted$SNP %in% cordsCon2.sorted$SNP, ]
	cordsCon2.shared = cordsCon2.sorted[cordsCon2.sorted$SNP %in% cordsCon1.sorted$SNP, ]

	# Merge MAFs for shared SNPs ready for plotting
	cordsCleaned = data.frame(cordsCon1.shared$SNP, cordsCon1.shared$Con1, cordsCon2.shared$Con2)
	names(cordsCleaned) = c("SNP", "Con1", "Con2")

	correl <- round(cor(cordsCleaned$Con1, cordsCleaned$Con2, method = "spearman"), digits=5)

	# Plot as solid transparent circles, allows measure of intensity
	plot(cordsCleaned$Con1, cordsCleaned$Con2, 
		xlab=xname, ylab=yname,
		xlim=range(0.0,0.5), ylim=range(0.0,0.5),
		pch=16, col=rgb(0,0,0,0.01),
		main=paste("rho =", correl))
	# add 1:1 dotted line
	abline(0, 1, lty=2, col="red")
}

png("MAF_Plots.png", width = 1200, height = 450)

# subplot 3 col, 1 row
par(mfrow=c(1,3), oma=c(0,0,2,0))
plotMAFs("Case", "Control 1", "CaseFreqs.frq", "Con1Freqs.frq")
plotMAFs("Case", "Control 2", "CaseFreqs.frq", "Con2Freqs.frq")
plotMAFs("Control 1", "Control 2", "Con1Freqs.frq", "Con2Freqs.frq")
title(main="Minor Allele Frequencies per SNP per dataset", outer=T)

dev.off()