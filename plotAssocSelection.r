# Currently under construction
# Aim to clean this up and produce a script which will overlay association and selection data at a given locus
# Take in genomic region & chromosome (xr1:start and xr2:end), selection (association, ihs, xpehh, haplos) as bools & paths,
# Output plot

plotTemplate <- function(assoc, ihs, xpehh, haplo, xrange1, xrange2, title, haplox1, haplox2) {

}


plotAsoc_XPEHH <- function(ASSOC, XP, xrange1, xrange2, title) {
  png('outFile.png',height=750,width=1000,res=150)
  plot(get(ASSOC)['BP'], -log10(get(ASSOC)['P']), pch=20, cex=0.5, xlim=c(xrange1,xrange2),xlab='Position (bp)',ylab=expression(-log[10](italic(p))), main=title)
  par(new=T); plot(get(XP)['pos'], abs(get(XP)['normxpehh']), pch=20, cex=0.5, xlim=c(xrange1,xrange2),xlab=NA,ylab=NA,axes=F,col='red')
  axis(side=4, col='red'); mtext(side=4,line=3,'|XP-EHH|', col='red')
  legend('topleft', legend=c(expression(-log[10](italic(p))), '|XP-EHH|'), pch=c(20, 20), col=c('black','red')); dev.off()
}

plotAssoc_iHS_HaploPS <- function() {
  png('outFile.png',height=750,width=1000,res=150)
  par(mar=c(5,5,2,5))
  plot(assoc.12$BP, -log10(assoc.12$P), pch=20, xlim=c(xr1,xr2), main=title, xlab='Position (bp)', ylab=expression(-log[10](italic(p))))
  par(new=T)
  plot(ihs$V2, abs(ihs$V7), xlim=c(xr1,xr2), col='red', axes=F,xlab=NA,ylab=NA,pch=20)
  segments(xxr1.2,3,x1=xxr2.2,y1=3, col=rgb(0,0,1,0.4), lwd=1000, lend=1)
  axis(side=4, col='red'); mtext(side=4,line=3,'|iHS|',col='red')
  legend('topleft', legend=c(expression(-log[10](italic(p))), '|iHS|', 'HaploPS'), pch=c(20, 20, 15), col=c('black','red','blue'))
  axis(side=1, col='blue')
  dev.off()
}
