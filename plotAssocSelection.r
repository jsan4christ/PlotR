# Currently under construction
# Take in genomic region & chromosome (xr1:start and xr2:end),& selection (association, ihs, xpehh, haplos)
# Output plot

# TODO: Include HaploPS, perhaps as a list of start and ends

plotAssocSelection <- function(assoc,assocpos='BP',assocval='P',ihs=NA,ihspos='V2',ihsval='V7',xpehh=NA,xpehhpos='pos',xpehhval='norm$
        print(max(assoc[,assocpos]))

        # TODO: check that ihs OR xpehh have been included

        # TODO: check non NA files exist, return error that these need to be parsed if they do not

        options(scipen=5)

        # Plot base associations
        print('Assoc')
        plot(assoc[,assocpos], -log10(assoc[,assocval]), pch=20, cex=0.5, xlim=c(x1,x2), xlab='Position (bp)', ylab=expression(-log[10](itali$

        # Add second plot layer only if ihs or xpehh file was provided

        if (is.data.frame(ihs)) {
                print('iHS')
                m2 <- '|iHS|'
                par(new=T)
                plot(ihs[,ihspos], abs(ihs[,ihsval]), pch=20, cex=0.5, xlim=c(x1,x2), xlab=NA, ylab=NA, axes=F, col='red')
        } else if (is.data.frame(xpehh)) {
                print('XP-EHH')
                m2 <- '|XP-EHH|'
                par(new=T)
                plot(xpehh[,xpehhpos], abs(xpehh[,xpehhval]), pch=20, cex=0.5, xlim=c(x1,x2), xlab=NA, ylab=NA, axes=F, col='red')
        }

        # Legend
        axis(side=4, col='red')
        mtext(side=4, line=3, m2, col='red')

        # HaploPS Segments

        ## TO DO ##

        # Legend
        legend('topleft', legend=c(expression(-log[10](italic(p))), m2, 'HaploPS'), pch=c(20,20,15), col=c('black','red','blue'))
}

# Precursor functions
plotAsoc_XPEHH <- function(ASSOC, XP, xrange1, xrange2, title) {
  png('outFile.png',height=750,width=1000,res=150)
  plot(ASSOC['BP'], -log10(ASSOC['P']), pch=20, cex=0.5, xlim=c(xrange1,xrange2),xlab='Position (bp)',ylab=expression(-log[10](italic(p))), main=title)
  par(new=T); plot(XP['pos'], abs(XP['normxpehh']), pch=20, cex=0.5, xlim=c(xrange1,xrange2),xlab=NA,ylab=NA,axes=F,col='red')
  axis(side=4, col='red'); mtext(side=4,line=3,'XP-EHH|', col='red')
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
