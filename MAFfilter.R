prepfile <- function(infile) {
	# convert Plink --freq output file, return SNPs and their MAFs
	# take in raw file location, return read.table version
	infile1 = read.table(infile, header=T)

	# remove NAs
	infile2 = infile1[complete.cases(infile1), ]

	# remove 0s in A1 and A2
	infile3 = infile2[infile2$A1 != 0, ]
	infile3 = infile2[infile2$A2 != 0, ]
  
	# cut to SNP and MAF only
	infile4 = data.frame(infile3$SNP, infile3$MAF)

	# check names are in place for later reference
	names(infile4) = c("SNP", "MAF")

	# sort by SNP column
	infile5 = infile4[order(infile4$SNP), ]
  
  # remove row column that order adds
  infile6 = data.frame(infile5$SNP, infile5$MAF)
  
  # reset names
  names(infile6) = c("SNP", "MAF")

	return(infile6)
}

cleanMAFs <- function(file1, file2, difference=0.1) {
	# take .freq files from Plink output, plus parameter for max difference in MAFs
	# return list of SNPs with MAF max difference or more different

	Con1 <- prepfile(file1)
	Con2 <- prepfile(file2)

	names(Con1) = c("SNP", "Con1")
	names(Con2) = c("SNP", "Con2")

	# Remove non-shared SNPs (eases merge)
	Con1.shared = Con1[Con1$SNP %in% Con2$SNP, ]
	Con2.shared = Con2[Con2$SNP %in% Con1$SNP, ]

	# Merge MAFs for shared SNPs ready for plotting
	merged = data.frame(Con1.shared$SNP, Con1.shared$Con1, Con2.shared$Con2)
	names(merged) = c("SNP", "Con1", "Con2")

	# Initiate vector for SNPs to exclude
	SNPlist <- c()

	# for each row in merged, find the SNPs where the abs(Con1MAF-Con2MAF >= maxDifference)
	for (i in 1:nrow(merged)){
		if (abs(merged[i, 2] - merged[i, 3]) > difference) { # the co-ords here are probably wrong
			# add SNP to list to exclude later
			SNPlist <- c(SNPlist, as.character(merged[i,1]))
		}
	}
	return(SNPlist)
}

file1 = #file location
file2 = #file location
file3 = #file location

write.table(cleanMAFs(file1, file2), file ="MAFexclusion_1_p1", sep='\n', quote=F, row.names=F, col.names=F)
write.table(cleanMAFs(file1, file3), file ="MAFexclusion_2_p1", sep='\n', quote=F, row.names=F, col.names=F)
