




# iCdR703 base growth demands
solo_demands <- read.delim(file='/home/mjenior/Desktop/solo_demands.tsv', sep='\t', header=TRUE)
#solo_demands$median <- log10(solo_demands$median + 1)
#solo_demands$q25 <- log10(solo_demands$q25 + 1)
#solo_demands$q75 <- log10(solo_demands$q75 + 1)
solo_demands$median <- round(solo_demands$median)
solo_demands <- subset(solo_demands, median > 1)

bproducta_demands <- read.delim(file='/home/mjenior/Desktop/bproducta_demands.tsv', sep='\t', header=TRUE)
#bproducta_demands$median <- log10(bproducta_demands$median + 1)
#bproducta_demands$q25 <- log10(bproducta_demands$q25 + 1)
#bproducta_demands$q75 <- log10(bproducta_demands$q75 + 1)
#solo_demands$q75 <- log10(solo_demands$q75 + 1)
bproducta_demands$median <- round(bproducta_demands$median)
bproducta_demands <- subset(bproducta_demands, median > 1)

png(filename='~/Desktop/solo_demands.png', units='in', width=5, height=3, res=300)
par(mar=c(7,3,0.5,0.5), las=1, mgp=c(1.9,0.7,0), xpd=FALSE, lwd=2)
barplot(as.vector(solo_demands$median), ylab='Substrate Uptake Rate', col=as.vector(solo_demands$color), 
        yaxt='n', ylim=c(0,1000), cex.lab=0.9,
        names.arg=as.vector(solo_demands$substrate), las=2, cex.names=0.7)
axis(side=2, at=seq(0,1000,200), cex.axis=0.8, lwd=2)
legend('topright', legend=c('Nucleotide','Amino acid','Carbohydrate','Ion','Vitamin','Other'), pt.bg=c('blue3','firebrick2','chartreuse2','darkorchid1','darkorange1','ivory'), 
       pch=22, pt.cex=1.5, cex=0.9, box.lwd=0)
#text(x=12, y=3.75, font=2, cex=0.8, labels='iCdR703 substrate usage in solo-culture')
dev.off()

png(filename='~/Desktop/bproducta_demands.png', units='in', width=5, height=3, res=300)
par(mar=c(7,3,0.5,0.5), las=1, mgp=c(1.9,0.7,0), xpd=FALSE, lwd=2)
barplot(as.vector(bproducta_demands$median), ylab='Substrate Uptake Rate', col=as.vector(bproducta_demands$color), 
        yaxt='n', ylim=c(0,1000), cex.lab=0.9,
        names.arg=as.vector(bproducta_demands$substrate), las=2, cex.names=0.7)
axis(side=2, at=seq(0,1000,200), cex.axis=0.8, lwd=2)
legend('topright', legend=c('Nucleotide','Amino acid','Carbohydrate','Ion','Vitamin','Other'), pt.bg=c('blue3','firebrick2','chartreuse2','darkorchid1','darkorange1','ivory'), 
       pch=22, pt.cex=1.5, cex=0.9, box.lwd=0)
#text(x=12, y=3.75, font=2, cex=0.8, labels='iCdR703 substrate usage in co-culture with B. producta')
dev.off()


#----------------------------------------------------#

# in silico 
pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/insilico_interaction.pdf', width=6, height=3)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

# B. producta and R20291 competition 
bproducta_biomass <- c(446.25, 430.17, 421.39, 418.93)
r20291_biomass <- c(159.67, 149.49, 147.87, 0.0)
bproducta_biomass <- (bproducta_biomass / bproducta_biomass[1]) * 100.
r20291_biomass <- (r20291_biomass / r20291_biomass[1]) * 100.
timepoints <- c(0,1,2,3)
#pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/bproducta_competition.pdf', width=4, height=3.5)
par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.9,0.7,0), xpd=FALSE, lwd=2.5)
plot(timepoints, bproducta_biomass, type='b', pch=16, xlim=c(0,3), ylim=c(0,100), cex=1.2,
     col='chocolate2', xlab='Time point', ylab='Growth Rate (%)', cex.axis=0.8)
lines(timepoints, r20291_biomass, pch=18, col='chartreuse3', type='b', cex=1.2)
legend('bottomleft', legend=c('B. producta','C. difficile'), pt.cex=1.2, lwd=2,
       col=c('chocolate2', 'chartreuse3'), cex=0.9, text.font=3, lty=1, pch=c(16,18))
#dev.off()

# S. thermophillus and R20291 cooperation
sthermophillus_biomass <- c(778.06, 778.06, 827.47, 832.77, 842.41, 842.98)
r20291_biomass <- c(69.45, 90.03, 97.14, 97.17, 99.96, 99.96)
sthermophillus_biomass <- (sthermophillus_biomass / sthermophillus_biomass[1]) * 100.
r20291_biomass <- (r20291_biomass / r20291_biomass[1]) * 100.
timepoints <- c(0,1,2,3,4,5)
#pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/sthermophillus_cooperation.pdf', width=4, height=3.5)
par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.9,0.7,0), xpd=FALSE, lwd=2.5)
plot(timepoints, sthermophillus_biomass, type='b', pch=16, xlim=c(0,5), ylim=c(100,160), cex=1.2,
     col='firebrick2', xlab='Time point', ylab='Growth Rate (%)', cex.axis=0.8)
lines(timepoints, r20291_biomass, pch=18, col='chartreuse3', type='b', cex=1.2)
legend('topleft', legend=c('S. thermophillus','C. difficile'), pt.cex=1.2, lwd=2,
       col=c('firebrick2', 'chartreuse3'), cex=0.9, text.font=3, lty=1, pch=c(16,18))

dev.off()

#-----------------------------------------------------------------------------------------------------#

# in vitro 

collateNorm <- function(col_num, data) {
  blank <- paste0('A', as.character(col_num))
  wells <- c(paste0('B', as.character(col_num)), paste0('C', as.character(col_num)), 
             paste0('D', as.character(col_num)), paste0('E', as.character(col_num)), 
             paste0('F', as.character(col_num)), paste0('G', as.character(col_num)))
  blank <- as.matrix(data[,blank])
  growth <- as.matrix(data[,wells])
  blank <- apply(blank, 1, quantile, probs=c(0.25,0.5,0.75))
  growth <- apply(growth, 1, quantile, probs=c(0.25,0.5,0.75))
  growth <- abs(growth - blank)
  growth <- t(growth)
  growth <- growth[1:169,] # 18 hours
  return(growth)}
growth <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/spent_media.2.R20291.tsv', sep='\t', header=TRUE)
dsm2950 <- collateNorm(3, growth) # Blautia producta
growth <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/spent_media.3.R20291.tsv', sep='\t', header=TRUE)
fresh <- collateNorm(1, growth) # no depletion
lmd9 <- collateNorm(9, growth) # Streptococcus thermophillus
rm(growth)

library(scales)
draw_growth <- function(data, clr) {
  polygon(c(1:nrow(data), rev(1:nrow(data))), c(data[,3], rev(data[,1])), col=alpha(clr, 0.25), border=NA)
  lines(data[,2], lwd=3, col=clr)}

pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/invitro_interaction.pdf', width=4.5, height=4.5)
par(mar=c(3,3,0.5,0.5), mgp=c(2, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(fresh[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', ylim=c(0,0.8), xlim=c(1,145))
axis(side=1, at=seq(1,176,16), labels=seq(0,20,2), cex.axis=0.8, lwd=2)
axis(side=2, at=seq(0,0.8,0.2), cex.axis=0.8, lwd=2)
draw_growth(fresh, 'black')
draw_growth(dsm2950, 'chocolate2')
draw_growth(lmd9, 'firebrick2')
legend('topleft', legend=c('Fresh media', 'B. producta conditioned', 'S. thermophillus conditioned'), pt.cex=1.2,
       col=c('black', 'chocolate2', 'firebrick2'), cex=0.9, lty=1, lwd=3)
dev.off()







plot_diff <- function(metabolite, x) {
  fresh <- log10(as.vector(subset(B_producta_DSM2950, Depleting_species == 'Fresh_media')[,metabolite]))
  bproducta <- log10(as.vector(subset(B_producta_DSM2950, Depleting_species == 'B_producta_DSM2950')[,metabolite]))
  
  stripchart(fresh, vertical=T, pch=21, at=x, bg='gray40', 
             cex=1.1, method='jitter', jitter=0.15, add=TRUE)
  stripchart(bproducta, vertical=T, pch=21, at=x+1, bg='chocolate2', 
             cex=1.1, method='jitter', jitter=0.15, add=TRUE)
  
  pval <- round(wilcox.test(fresh, bproducta, exact=FALSE)$p.value, 4)
  print(pval)
  ymax <- max(c(max(fresh),max(bproducta))) * 1.1
  
  if (pval <= 0.05) {
    segments(x0=x, y0=ymax, x1=x+1, lwd=2)
    text(x=x+0.5, y=ymax*1.05, '*', font=2, cex=1.5)
  }
}

plot_diff('Isoleucylisoleucine', 0.75) 


test_diff('Leucyl.leucine')
test_diff('Leucylvaline')
test_diff('Cyclo.Leu.Pro.')


test <- as.matrix(t(competition[,c(5,6)]))
colnames(test) <- competition$metabolite



png(filename='~/Desktop/K22_resubmission/metabolomics.png', units='in', width=4, height=3.5, res=300)
par(mar=c(5,3,0.5,0.5), las=1, mgp=c(1.6,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,12), xlim=c(0,8), cex.axis=0.8,
     ylab='Scaled Intensity (Log10)', xlab='', xaxt='n', cex.lab=1.1)
legend('bottomright', legend=c('Fresh media', 'B. producta depleted'), pt.cex=1.2,
       pt.bg=c('black', 'chocolate2'), cex=0.9, pch=21)
par(xpd=TRUE)
text(x=c(0.75, 2.75, 4.75, 6.75), y=-2.5, labels=c('Isoleucine','Leucine','Valine','Proline'), srt=45)
par(xpd=FALSE)
plot_diff('Isoleucylisoleucine', 0.25)
plot_diff('Leucyl.leucine', 2.25)
plot_diff('Leucylvaline', 4.25)
plot_diff('Cyclo.Leu.Pro.', 6.25)
dev.off()






