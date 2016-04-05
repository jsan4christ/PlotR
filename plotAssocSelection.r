# plotAssocSelection.r
# Matt Ravenhall
# Given an association results files and an iHS or XP-EHH results files, plus an option list of haplotypes,
# Plots assoc & selection signals for a genomic region (default full range of loci in assoc file)

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
# plot: boolean, if TRUE plots to default name png
# plotname: optional name for the default plot
# y1bot: optional, customise y1 bottom axis
# y1top: optional, customise y1 top axis
# y2bot: optional, customise y2 bottom axis
# y2top: optional, customise y2 top axis
# exons: optional vector of numeric pairs, for plotting exons
# y1line: add line for y1 axis (assoc)
# y2line: add line for y2 axis (ihs/xpehh)
# gene: optional vector of length 2, for plotting custom start and end of gene, if absent uses start and end of xons
# genename: optional string for gene name in legend
# allcex: all cex sizes
# allpch: all pch plot points

plotAssocSelection <- function(assoc, assocpos='BP', assocval='P', ihs=NA, ihspos='V2', ihsval='V7', xpehh=NA, xpehhpos='pos', xpehhval='normxpehh', haplops=NA, x1=min(assoc[,assocpos]), x2=max(assoc[,assocpos]), title=NA, plot=FALSE, plotname='defaultOut.png', y1bot=-log10(max(assoc$P)), y1top=-log10(min(assoc$P)), y2bot=NA, y2top=NA, exons=NA, y1line=NA, y2line=NA, gene=NA, genename='Gene',allcex=0.5,allpch=20) {

	options(scipen=5)
	ishaplops <- (is.vector(haplops) & !is.na(haplops))
	isexons <- (is.vector(exons) & !is.na(exons))
	isihs <- is.data.frame(ihs)
	isxpehh <- is.data.frame(xpehh)
	if (isihs & isxpehh) { print('Warning: Both iHS and XP-EHH files passed, defaulting to iHS only.') }
	
	if (plot) { png(plotname, height=750, width=1000, res=100) }

	# Plot base associations
	print('Assoc')
	plot(assoc[,assocpos], -log10(assoc[,assocval]), pch=allpch, cex=allcex, xlim=c(x1,x2), ylim=c(y1bot, y1top), xlab='Position (bp)', ylab=expression(-log[10](italic(p))), main=title)
	if (is.numeric(y1line)) { abline(h=y1line,col='black',lty=2) }

	# Add second plot layer only if ihs or xpehh file was provided
	if (isihs) {
		print('Notice: Assuming input is an iHS results file.')
		m2 <- '|iHS|'
		par(new=T)
		if (is.na(y2bot)) { y2bot <- min(abs(ihs[ihsval])) }
		if (is.na(y2top)) { y2top <- max(abs(ihs[ihsval])) }
		plot(ihs[,ihspos], abs(ihs[,ihsval]), pch=allpch, cex=allcex, xlim=c(x1,x2), ylim=c(y2bot, y2top), xlab=NA, ylab=NA, axes=F, col='red')
		if (is.numeric(y2line)) { abline(h=y2line,col='red',lty=2) }
	} else if (isxpehh) {
		print('Notice: Assuming input is an XP-EHH results file.')
		m2 <- '|XP-EHH|'
		par(new=T)
		if (is.na(y2bot)) { y2bot <- min(abs(xpehh[xpehhval])) }
		if (is.na(y2top)) { y2top <- max(abs(xpehh[xpehhval])) }
		plot(xpehh[,xpehhpos], abs(xpehh[,xpehhval]), pch=allpch, cex=allcex, xlim=c(x1,x2), ylim=c(y2bot, y2top), xlab=NA, ylab=NA, axes=F, col='red')
		if (is.numeric(y2line)) { abline(h=y2line,col='red',lty=2) }
	} else {
		print('Notice: No iHS or XP-EHH results file provided')
	}

	# Selection plot axis
	axis(side=4, col='red')
	mtext(side=4, line=3, m2, col='red')

	# HaploPS Segments
	if(ishaplops) {
		# check haplops is a vector of pairs
		if ( (length(haplops)/2)%%1 == 0 ) {
			print(paste(length(haplops)/2,' HaploPS pairs detected.',sep=''))
			for (i in 1:(length(haplops)/2)) { #exons
				print(paste(haplops[(2*i)-1], haplops[2*i], sep=' '))
				segments(haplops[(2*i)-1], 4, x1=haplops[2*i], y1=4, col=rgb(0,0,1,0.4), lwd=1000, lend=1)
			}
			axis(side=1, col='blue')
		} else {
			print('Warning: Variable haplops is not a vector of paired integers.')
		}
	}

	exonh=-0.2
	# Custom gene annotation
	if (isexons) {
		# Check gene locations are in pairs
		if ( (length(exons)/2)%%1 == 0 ) {
			print(paste(length(exons)/2,' exons detected.',sep=''))
			for (i in 1:(length(exons)/2)) {
				print(paste(exons[(2*i)-1], exons[2*i], sep=' '))
				segments(exons[(2*i)-1], exonh, x1=exons[2*i], y1=exonh, col='black', lwd=30, lend=1)
			}
			# whole gene
			if (length(gene) != 2) {
				print('Assuming gene length as start of exon 1 to end of exon n.')
				segments(exons[1],exonh,x1=exons[length(exons)],y1=exonh, col=rgb(0.1,0.5,0.1,0.6), lwd=30, lend=1)
			} else {
				print(paste('Plotting gene length as ',gene[1],' to ',gene[2],'.',sep=''))
				segments(gene[1],exonh,x1=gene[2],y1=exonh, col=rgb(0.1,0.5,0.1,0.6), lwd=30, lend=1)
			}
		}
	}

	# Legend
	# check for haplops and genes
	if (ishaplops) {
		if (isexons) {
			# Yes HaploPS & Yes Gene
			legend('topleft', legend=c(expression(-log[10](italic(p))), m2, 'HaploPS', genename), pch=c(20,20,15,15), col=c('black','red','blue',rgb(0.1,0.5,0.1,1)))
		} else {
			# Yes HaploPS, No gene
			legend('topleft', legend=c(expression(-log[10](italic(p))), m2, 'HaploPS'), pch=c(20,20,15), col=c('black','red','blue'))
		}
	} else {
		if (isexons) {
			# No HaploPS, Yes gene
			legend('topleft', legend=c(expression(-log[10](italic(p))), m2, genename), pch=c(20,20,15), col=c('black','red',rgb(0.1,0.5,0.1,1)))
		} else {
			# No HaploPS, No gene
			legend('topleft', legend=c(expression(-log[10](italic(p))), m2), pch=c(20,20), col=c('black','red'))
		}
	}
	if (plot) { dev.off() }
}
