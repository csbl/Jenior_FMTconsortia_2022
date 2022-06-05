
# Format data
bproducta_comp <- c(159.67, 149.49, 147.87, 0.0, 0.0, 0.0)
bproducta_comp <- (bproducta_comp / bproducta_comp[1]) * 100.
bproducta_coop <- c(69.45, 80.04, 83.97, 83.97, 83.97, 83.97)
bproducta_coop <- (bproducta_coop / bproducta_coop[1]) * 100.
blongum_comp <- c(149.49, 145.63, 145.63, 121.0, 80, 80)
blongum_comp <- (blongum_comp / blongum_comp[1]) * 100.
blongum_coop <- c(69.45, 82.62, 92.25, 100.05, 107.59, 107.59)
blongum_coop <- (blongum_coop / blongum_coop[1]) * 100.
timepoints <- c(0,1,2,3,4,5)

# Specific axes of interaction
bproducta_index <- read.delim('~/Desktop/Jenior_Consortia_2022/data/Bproducta_CdR20291_competition.tsv', sep='\t', header=TRUE)
bproducta_index <- subset(bproducta_index, metabolite != 'H+')
bproducta_index$genre2_comp_index <- bproducta_index$genre2_comp_index * 100
bproducta_index$genre2_comp_index <- bproducta_index$genre2_comp_index * -1.0
bproducta_index <- bproducta_index[order(-bproducta_index$genre1_comp_index),] 
bproducta_index <- bproducta_index[,c('metabolite','genre1_comp_index','genre2_comp_index')]
blongum_index <- read.delim('~/Desktop/Jenior_Consortia_2022/data/Blongum_CdR20291_mutualism.tsv', sep='\t', header=TRUE)
blongum_index <- subset(blongum_index, metabolite != 'H+')
blongum_index[is.na(blongum_index)] <- 0
blongum_index$colsums <- blongum_index$genre1_coop_index + blongum_index$genre2_coop_index
blongum_index <- subset(blongum_index, colsums != 0)
blongum_index$colsums <- NULL
blongum_index$genre2_coop_index <- blongum_index$genre2_coop_index * -1.0
blongum_index <- blongum_index[,c('metabolite','genre1_coop_index','genre2_coop_index')]
blongum_index <- blongum_index[order(-blongum_index$genre1_coop_index),] 

#---------------------------------------------------------------------------------------------#

# Generate figure
bproducta_col <- 'chocolate2'
blongum_col <- 'cyan3'
cdiff_col <- 'darkorchid3'
bproducta_lab <- as.expression(bquote(paste('with ',italic('B. producta'),' GENRE')))
blongum_lab <- as.expression(bquote(paste('with ',italic('B. longum'),' GENRE')))

pdf(file='~/Desktop/Jenior_Consortia_2022/results/insilico_interaction1.pdf', width=4, height=4)
par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2, xaxs='i')
plot(timepoints, bproducta_comp, type='b', pch=16, xlim=c(0,5.2), ylim=c(0,180), cex=1.2, xaxt='n', yaxt='n',
     col='white', xlab='Simulation time point', ylab='iCdR703 Biomass Flux (%)')
axis(side=1, at=c(1:5), cex.axis=0.7, lwd=2)
axis(side=2, at=seq(0,180,20), cex.axis=0.7, lwd=2)
abline(h=100, col='gray40', lty=3)
lines(timepoints, blongum_comp, pch=20, col=blongum_col, type='b', cex=1.2)
lines(timepoints, blongum_coop, pch=20, col=blongum_col, type='b', cex=1.2, lty=2)
lines(timepoints, bproducta_comp, pch=18, col=bproducta_col, type='b', cex=1.2)
lines(timepoints, bproducta_coop, pch=18, col=bproducta_col, type='b', cex=1.2, lty=2)
legend('topleft', legend=c(blongum_lab,bproducta_lab), pt.cex=c(1.3,1.5), 
       col=c(blongum_col,bproducta_col), cex=1, text.font=3, pch=c(16,18))
legend('bottomleft', legend=c('Syntrophy','Competition'), pt.cex=0, lwd=2, cex=1, lty=c(2,1), pch=1)
box()
dev.off()

#cex.lab=1.1, cex.main=1, main=substitute(paste(bolditalic('C. difficile'), bold(' GENRE interactions')))

#--------------------------------------------------#

bproducta_lab <- c(as.expression(bquote(paste(italic('B. producta'),' contested'))), as.expression(bquote(paste(italic('C. difficile'),' contested'))))
blongum_lab <- c(as.expression(bquote(paste(italic('B. longum'), ' crossfed'))), as.expression(bquote(paste(italic('C. difficile'), ' crossfed'))))
bproducta_col <- 'chocolate2'
blongum_col <- 'cyan3'
cdiff_col <- 'darkorchid3'

pdf(file='~/Desktop/Jenior_Consortia_2022/results/insilico_interaction2.pdf', width=5, height=5)
layout(matrix(c(1,
                2), nrow=2, ncol=1, byrow=TRUE))

# B. longum cooperation
par(mar=c(0.5,3,0.5,1.5), las=1, mgp=c(1.7,0.7,0), lwd=2, xpd=FALSE)
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,nrow(blongum_index)), ylim=c(-4,22), xlab='', ylab='Cooperative Index', cex.lab=1.1)
abline(h=0, lty=2, col='gray40')
box()
axis(2, at=seq(-4,20,4), cex.axis=0.9, lwd=2)
barplot(blongum_index$genre1_coop_index, col=blongum_col, width=1, space=0, xaxt='n', yaxt='n', add=TRUE)
barplot(blongum_index$genre2_coop_index, col=cdiff_col, width=1, space=0, xaxt='n', yaxt='n', add=TRUE)
text(x=seq(0.5, nrow(blongum_index)-0.5, 1), y=blongum_index$genre1_coop_index + 6.5, 
     labels=blongum_index$metabolite, cex=0.8, srt=90, font=2)
legend('topright', legend=blongum_lab, pt.bg=c(blongum_col,cdiff_col), pch=22, pt.cex=1.8, pt.lwd=1.5, cex=0.8)
par(xpd=TRUE)
text(x=nrow(blongum_index)+1.5, y=9, 'Syntrophy', srt=270, cex=1.2, font=2)

# B. producta competition
par(mar=c(0.5,3,0.5,1.5), las=1, mgp=c(1.7,0.7,0), lwd=2, xpd=FALSE)
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,nrow(bproducta_index)), ylim=c(-4,16), 
     xlab='', ylab='Competitive Index', cex.lab=1.1)
abline(h=0, lty=2, col='gray40')
box()
axis(2, at=seq(-4,16,4), cex.axis=0.9, lwd=2)
barplot(bproducta_index$genre1_comp_index, col=bproducta_col, width=1, space=0, xaxt='n', yaxt='n', add=TRUE)
barplot(bproducta_index$genre2_comp_index, col=cdiff_col, width=1, space=0, xaxt='n', yaxt='n', add=TRUE)
text(x=seq(0.5, nrow(bproducta_index)-0.5, 1), y=bproducta_index$genre1_comp_index + 5, 
     labels=bproducta_index$metabolite, cex=0.8, srt=90, font=2)
legend('topright', legend=bproducta_lab, pt.bg=c(bproducta_col,cdiff_col), pch=22, pt.cex=1.8, pt.lwd=1.5, cex=0.8)
par(xpd=TRUE)
text(x=nrow(bproducta_index)+1, y=6, 'Competition', srt=270, cex=1.2, font=2)

dev.off()


