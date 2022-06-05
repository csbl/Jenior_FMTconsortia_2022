
insilicoCompetition <- function(strain, r20291, strName, strCol, leg='bottomleft') {
  strain <- (strain / strain[1]) * 100.
  r20291 <- (r20291 / r20291[1]) * 100.
  timepoints <- c(1:length(strain)-1)
  par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.9,0.7,0), xpd=FALSE, lwd=2.5)
  plot(timepoints, strain, type='b', pch=16, xlim=c(0,length(strain)-1), ylim=c(0,100), cex=1.2,
       col=strCol, xlab='Time point', ylab='Max Growth Rate (%)', cex.axis=0.8)
  lines(timepoints, r20291, pch=18, col='chartreuse3', type='b', cex=1.2)
  legend(leg, legend=c(strName,'C. difficile R20291'), pt.cex=1.1, lwd=2,
         col=c(strCol, 'chartreuse3'), cex=0.8, text.font=3, lty=1, pch=c(16,18))}

png(file='~/Desktop/insilico_competition.png', units='in', width=5, height=4.5, res=300)
par(mar=c(3,3,1,1), xpd=TRUE, lwd=2, las=1, font=2)
layout(matrix(c(1,2,
                3,4), nrow=2, ncol=2, byrow=TRUE))
bproducta <- c(446.25, 430.17, 421.39, 418.93)
r20291_bproducta <- c(159.67, 149.49, 147.87, 0.0)
insilicoCompetition(bproducta, r20291_bproducta, 'B. producta', 'brown3')
erectale <- c(453.64, 429.53, 424.19, 421.24, 420.36, 420.18)
r20291_erectale <- c(159.67, 159.67, 159.67, 159.67, 125.0, 0.0)
insilicoCompetition(erectale, r20291_erectale, 'E. rectale', 'cyan3')
bvulgatus <- c(436.67, 381.64, 369.13, 365.72)
r20291_bvulgatus <- c(159.67, 149.49, 144.08, 0.0)
insilicoCompetition(bvulgatus, r20291_bvulgatus, 'P. vulgatus', 'darkorchid3')
bpsedo <- c(947.96, 945.51, 944.27, 942.63)
r20291_bpsedo <- c(159.67, 149.49, 144.25, 0.0)
insilicoCompetition(bpsedo, r20291_bpsedo, 'B. psedocatenulatum', 'goldenrod3')
dev.off()




insilicoCooperation <- function(strain, r20291, strName, strCol, leg='bottomright') {
  strain <- (strain / strain[1]) * 100.
  r20291 <- (r20291 / r20291[1]) * 100.
  timepoints <- c(1:length(strain)-1)
  ymax <- round(max(c(strain, r20291)))
  par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.9,0.7,0), xpd=FALSE, lwd=2.5)
  plot(timepoints, strain, type='b', pch=16, xlim=c(0,length(strain)-1), ylim=c(100,ymax), cex=1.2,
       col=strCol, xlab='Time point', ylab='Max Growth Rate (%)', cex.axis=0.8)
  lines(timepoints, r20291, pch=18, col='chartreuse3', type='b', cex=1.2)
  legend(leg, legend=c(strName,'C. difficile R20291'), pt.cex=1.1, lwd=2,
         col=c(strCol, 'chartreuse3'), cex=0.8, text.font=3, lty=1, pch=c(16,18))}

png(file='~/Desktop/insilico_cooperation.png', units='in', width=5, height=4.5, res=300)
par(mar=c(3,3,1,1), xpd=TRUE, lwd=2, las=1, font=2)
layout(matrix(c(1,2,
                3,4), nrow=2, ncol=2, byrow=TRUE))
rintestinalis <- c(698.27, 739.68, 863.06, 863.69, 864.75, 830.75)
r20291_rintestinalis <- c(69.45, 91.73, 96.19, 96.19, 96.37, 96.91)
insilicoCooperation(rintestinalis, r20291_rintestinalis, 'R. intestinalis', 'firebrick2')
blongum <- c(935.38, 935.38, 1000.0, 1000.0, 1000.0, 1000.0)
r20291_blongum <- c(69.45, 82.62, 90.25, 94.05, 94.59, 94.76)
insilicoCooperation(blongum, r20291_blongum, 'B. longum', 'darkorange1', leg='right')
ecoli <- c(357.57, 357.57, 531.88, 540.06, 540.06, 540.06)
r20291_ecoli <- c(69.45, 73.22, 77.03, 81.06, 81.09, 83.88)
insilicoCooperation(ecoli, r20291_ecoli, 'E. coli', 'deepskyblue2')
sthermophillus <- c(778.06, 778.06, 827.47, 832.77, 842.41, 842.98)
r20291_sthermophillus <- c(69.45, 90.03, 97.14, 97.17, 99.96, 99.96)
insilicoCooperation(sthermophillus, r20291_sthermophillus, 'S. thermophillus', 'darkseagreen4', leg='right')
dev.off()













Starting flux Bacteroides caccae: 1000.0 Clostridium difficile R20291: 159.67
Iteration 1 | 12 contested metabolites | Bacteroides caccae: 1000.0 Clostridium difficile R20291: 149.49
Iteration 2 | 6 contested metabolites | Bacteroides caccae: 1000.0 Clostridium difficile R20291: 148.81
Iteration 3 | 6 contested metabolites | Bacteroides caccae: 1000.0 Clostridium difficile R20291: 148.25
Iteration 4 | 6 contested metabolites | Bacteroides caccae: 1000.0 Clostridium difficile R20291: 148.13
Iteration 5 | 5 contested metabolites | Bacteroides caccae: 1000.0 Clostridium difficile R20291: 148.11
Ending search after 5 iterations

Initial growth: Bacteroides caccae = 965.89 | Clostridium difficile R20291 = 69.45
Current growth: Bacteroides caccae = 965.89 | Clostridium difficile R20291 = 87.98
Current growth: Bacteroides caccae = 1000.0 | Clostridium difficile R20291 = 94.41
Current growth: Bacteroides caccae = 1000.0 | Clostridium difficile R20291 = 96.01
Current growth: Bacteroides caccae = 1000.0 | Clostridium difficile R20291 = 100.36
Current growth: Bacteroides caccae = 1000.0 | Clostridium difficile R20291 = 101.0
A total of 36 new metabolites added to media
15 metabolites passed from Bacteroides caccae to Clostridium difficile R20291
13 metabolites passed from Clostridium difficile R20291 to Bacteroides caccae
Optimal Bacteroides caccae objective flux increased by 34.11
Optimal Clostridium difficile R20291 objective flux increased by 31.55





