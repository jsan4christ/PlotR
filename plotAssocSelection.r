# plotAssocSelection.r
# Matt Ravenhall
# Given an association results files and an iHS or XP-EHH results files, plus an option list of haplotypes,
# Plots assoc & selection signals for a genomic region (default full range of loci in assoc file)

# TODO: Add some form of annotation overlay to identify where genes etc. are located.

# Format of arguments
# assoc: variable containing a data frame
# assocpos: BP column name for assoc data frame (default: 'BP')
# assocval: P value column name for assoc data frame (default: 'P')
# ihs: variable containing an ihs data frame, with at least position and ihs columns
# ihspos: BP column name for iHS data frame (default: 'V2')
# ihsval: iHS value column name for iHS data frame (default: 'V7')
# xpehh: variable containing an xp-ehh data frame, with at least position and ihs columns
# xpehhpos: BP column name for XP-EHH data frame (default: 'pos')
# xpehhval: XP-EHH value column name for XP-EHH data frame (default: 'normxpehh')
# haplops: vector containing haplotype start and end positions in pairs (ie. hap1start, hap1end, hap2start, hap2end)
# x1: start of xlim
# x2: end of xlim
# title: Figure title

plotAssocSelection <- function(assoc, assocpos='BP', assocval='P', ihs=NA, ihspos='V2', ihsval='V7', xpehh=NA, xpehhpos='pos', xpehhval='normxpehh', haplops=NA, x1=min(assoc[,assocpos]), x2=max(assoc[,assocpos]), title=NA, plot=FALSE, plotname='defaultOut.png') {

	# Optionally plot the figure as a .png, note that book-ending this function with png() and dev.off() will do the same.
	if (plot) { png(plotname, height=750, width=1000, res=100) }

	# TODO: check that ihs OR xpehh have been included, currently we default at iHS

	options(scipen=5)

	# Plot base associations
	print('Assoc')
	plot(assoc[,assocpos], -log10(assoc[,assocval]), pch=20, cex=0.5, xlim=c(x1,x2), xlab='Position (bp)', ylab=expression(-log[10](italic(p))), main=title)

	# Add second plot layer only if ihs or xpehh file was provided

	if (is.data.frame(ihs)) {
		print('Notice: Assuming input is an iHS results file.')
		m2 <- '|iHS|'
		par(new=T)
		plot(ihs[,ihspos], abs(ihs[,ihsval]), pch=20, cex=0.5, xlim=c(x1,x2), xlab=NA, ylab=NA, axes=F, col='red')
	} else if (is.data.frame(xpehh)) {
		print('Notice: Assuming input is an XP-EHH results file.')
		m2 <- '|XP-EHH|'
		par(new=T)
		plot(xpehh[,xpehhpos], abs(xpehh[,xpehhval]), pch=20, cex=0.5, xlim=c(x1,x2), xlab=NA, ylab=NA, axes=F, col='red')
	}

	# Selection plot axis
	axis(side=4, col='red')
	mtext(side=4, line=3, m2, col='red')

	# HaploPS Segments
	if(is.vector(haplops)) {
		# check haplops is a vector of pairs
		if ( (length(haplops)/2)%%1 == 0 ) {
			print(paste(length(haplops)/2,' HaploPS pairs detected.',sep=''))
			for (i in 1:(length(haplops)/2)) {
				print(paste(haplops[(2*i)-1], haplops[2*i], sep=' '))
				segments(haplops[(2*i)-1], 3, x1=haplops[2*i], y1=3, col=rgb(0,0,1,0.4), lwd=1000, lend=1)
			}
			axis(side=1, col='blue')
		} else {
			print('Warning: Variable haplops is not a vector of paired integers.')
		}
	}

	# Legend
	legend('topleft', legend=c(expression(-log[10](italic(p))), m2, 'HaploPS'), pch=c(20,20,15), col=c('black','red','blue'))

	if (plot) { dev.off() }
}
