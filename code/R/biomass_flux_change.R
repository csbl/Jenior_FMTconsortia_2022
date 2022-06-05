
path <- '~/Desktop/fmt_16S/data'
tax.seqtab.lachno <- read.delim(paste0(path, "/lachno.abund.tsv"), sep='\t', header=TRUE, row.names=1)
test <- apply(tax.seqtab.lachno, 1, FUN='median')
test <- (test/sum(test)) * 100
test <- sort(test)



complete_max <- 159.665
rich_max <- 125.621
minimal_max <- 41.92
competitive_max <- 0
cooperative_max <- 110.119

complete_percent <- (complete_max / complete_max) * 100
rich_percent <- (rich_max / complete_max) * 100
minimal_percent <- (minimal_max / complete_max) * 100
competitive_percent <- (competitive_max / complete_max) * 100
cooperative_percent <- (cooperative_max / complete_max) * 100


png(filename='~/Desktop/Jenior_Consortia_2022/results/biomass_change.png', width=3, height=4, units='in', res=300)
par(mar=c(3,3,1.5,1), mgp=c(1.6, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(0, type='n', ylab='Optimal Biomass Flux (%)', xaxt='n', yaxt='n', xlab='', ylim=c(0,105), xlim=c(0.8,2.2),
     main='iCdR703 biomass in consortia-derived media', cex.main=0.7)
legend('bottomleft', legend=c('Complete media','Rich medium','Minimal medium','Cooperative consortia','Competitive consortia'), 
       col=c('black','forestgreen','gray40','firebrick3','dodgerblue4'), pch=16, pt.cex=1.2, cex=0.75, bty='n')
axis(side=1, at=c(1,2), cex.axis=0.1, lwd=2)
axis(side=2, at=seq(0,100,20), cex.axis=0.75, lwd=2)

lines(c(100,rich_percent), lwd=2.5, col='forestgreen', type='b', pch=19) # rich media
lines(c(100,competitive_percent), lwd=2.5, col='dodgerblue4', type='b', pch=19) # competitive
lines(c(100,cooperative_percent), lwd=2.5, col='firebrick3', type='b', pch=19) # cooperative
lines(c(100,minimal_percent), lwd=2.5, col='gray40', type='b', pch=19) # minimal media
points(x=1, y=100, pch=19) # complete media

par(xpd=TRUE)
text(x=c(1,2), y=-10, cex=0.8, labels=c('All metabolites\navailable','Specified media\ncondition'))
box()

dev.off()


