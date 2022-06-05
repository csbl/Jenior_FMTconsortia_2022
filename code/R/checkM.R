
readCheckMTab <- function(sample_name, screen=TRUE) {
  path <- '/home/mjenior/Desktop/active_projects/CDI_FMT_project/fmt_metaG/checkM_results/'
  sample <- paste0(path, sample_name)
  sample <- read.delim(sample, header=TRUE, sep='\t')
  sample_name <- gsub('.checkM.all.tsv','',sample_name)
  sample$sample <- rep(sample_name, nrow(sample))
  sample$bin <- gsub('bin.','',sample$bin)
  sample$bin <- gsub('.format.pruned','',sample$bin)
  
  if (screen == TRUE) {
    sample <- subset(sample, completeness > 50)
    sample <- subset(sample, contamination < 33)
    sample <- subset(sample, genes >= 1354)
    sample <- subset(sample, genome_size >= 1300000)}
  return(sample)}

#----------------------------------#

# Before screening

checkM_01044A <- readCheckMTab('01044A.checkM.all.tsv', screen=FALSE)
checkM_01090A <- readCheckMTab('01090A.checkM.all.tsv', screen=FALSE)
checkM_01092A <- readCheckMTab('01092A.checkM.all.tsv', screen=FALSE)
checkM_01093A <- readCheckMTab('01093A.checkM.all.tsv', screen=FALSE)
completeness_01 <- c(checkM_01044A$completeness, checkM_01090A$completeness, 
                     checkM_01092A$completeness, checkM_01093A$completeness)
contamination_01 <- c(checkM_01044A$contamination, checkM_01090A$contamination, 
                      checkM_01092A$contamination, checkM_01093A$contamination)
#checkM_01_pre <- rbind(checkM_01044A, checkM_01090A, checkM_01092A, checkM_01093A)
rm(checkM_01044A, checkM_01090A, checkM_01092A, checkM_01093A)

checkM_05053B <- readCheckMTab('05053B.checkM.all.tsv', screen=FALSE)
checkM_05080M <- readCheckMTab('05080M.checkM.all.tsv', screen=FALSE)
checkM_05042G <- readCheckMTab('05042G.checkM.all.tsv', screen=FALSE)
checkM_05098A <- readCheckMTab('05098A.checkM.all.tsv', screen=FALSE)
completeness_05 <- c(checkM_05053B$completeness, checkM_05080M$completeness, 
                     checkM_05042G$completeness, checkM_05098A$completeness)
contamination_05 <- c(checkM_05053B$contamination, checkM_05080M$contamination, 
                      checkM_05042G$contamination, checkM_05098A$contamination)
#checkM_05_pre <- rbind(checkM_05053B, checkM_05080M, checkM_05042G, checkM_05098A)
rm(checkM_05053B, checkM_05080M, checkM_05042G, checkM_05098A)

checkM_28045D <- readCheckMTab('28045D.checkM.all.tsv', screen=FALSE)
checkM_28043C <- readCheckMTab('28043C.checkM.all.tsv', screen=FALSE)
checkM_28047D <- readCheckMTab('28047D.checkM.all.tsv', screen=FALSE)
checkM_28045A <- readCheckMTab('28045A.checkM.all.tsv', screen=FALSE)
completeness_28 <- c(checkM_28045D$completeness, checkM_28043C$completeness, 
                     checkM_28047D$completeness, checkM_28045A$completeness)
contamination_28 <- c(checkM_28045D$contamination, checkM_28043C$contamination, 
                      checkM_28047D$contamination, checkM_28045A$contamination)
#checkM_28_pre <- rbind(checkM_28045D, checkM_28043C, checkM_28047D, checkM_28045A)
rm(checkM_28045D, checkM_28043C, checkM_28047D, checkM_28045A)

pre_completeness <- c(completeness_01, completeness_05, completeness_28)
pre_contamination <- c(contamination_01, contamination_05, contamination_28)
rm(completeness_01, completeness_05, completeness_28)
rm(contamination_01, contamination_05, contamination_28)

#----------------------#

# After screening 

checkM_01044A <- readCheckMTab('01044A.checkM.tsv')

checkM_01090A <- readCheckMTab('01090A.checkM.tsv')
checkM_01092A <- readCheckMTab('01092A.checkM.tsv')
checkM_01093A <- readCheckMTab('01093A.checkM.tsv')
completeness_01 <- c(checkM_01044A$completeness, checkM_01090A$completeness, 
                     checkM_01092A$completeness, checkM_01093A$completeness)
contamination_01 <- c(checkM_01044A$contamination, checkM_01090A$contamination, 
                     checkM_01092A$contamination, checkM_01093A$contamination)
checkM_01_post <- rbind(checkM_01044A, checkM_01090A, checkM_01092A, checkM_01093A)
rm(checkM_01044A, checkM_01090A, checkM_01092A, checkM_01093A)

checkM_05053B <- readCheckMTab('05053B.checkM.tsv')
checkM_05080M <- readCheckMTab('05080M.checkM.tsv')
checkM_05042G <- readCheckMTab('05042G.checkM.tsv')
checkM_05098A <- readCheckMTab('05098A.checkM.tsv')
completeness_05 <- c(checkM_05053B$completeness, checkM_05080M$completeness, 
                     checkM_05042G$completeness, checkM_05098A$completeness)
contamination_05 <- c(checkM_05053B$contamination, checkM_05080M$contamination, 
                      checkM_05042G$contamination, checkM_05098A$contamination)
checkM_05_post <- rbind(checkM_05053B, checkM_05080M, checkM_05042G, checkM_05098A)
rm(checkM_05053B, checkM_05080M, checkM_05042G, checkM_05098A)

checkM_28045D <- readCheckMTab('28045D.checkM.tsv')
checkM_28043C <- readCheckMTab('28043C.checkM.tsv')
checkM_28047D <- readCheckMTab('28047D.checkM.tsv')
checkM_28045A <- readCheckMTab('28045A.checkM.tsv')
completeness_28 <- c(checkM_28045D$completeness, checkM_28043C$completeness, 
                     checkM_28047D$completeness, checkM_28045A$completeness)
contamination_28 <- c(checkM_28045D$contamination, checkM_28043C$contamination, 
                      checkM_28047D$contamination, checkM_28045A$contamination)
checkM_28_post <- rbind(checkM_28045D, checkM_28043C, checkM_28047D, checkM_28045A)
rm(checkM_28045D, checkM_28043C, checkM_28047D, checkM_28045A)
rm(readCheckMTab)

post_completeness <- c(completeness_01, completeness_05, completeness_28)
post_contamination <- c(contamination_01, contamination_05, contamination_28)
rm(completeness_01, completeness_05, completeness_28)
rm(contamination_01, contamination_05, contamination_28)

completeness_pval <- round(t.test(x=pre_completeness, y=post_completeness)$p.value, 4)
contamination_pval <- round(t.test(x=pre_contamination, y=post_contamination)$p.value, 4)

#---------------------------------#

# Generate figure
pdf(file='~/Desktop/Jenior_Consortia_2022/results/checkM.pdf', width=5, height=4)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

par(mar=c(3.5,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,120), xlim=c(0,2), 
     ylab='Completeness (%)', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(pre_completeness, cex=0, lwd=1.8, at=0.5, col='white', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
boxplot(post_completeness, cex=0, lwd=1.8, at=1.5, col='gray40', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(0,120,20), cex.axis=0.8)
segments(x0=0.5, y0=110, x1=1.5)
text(x=1, y=115, '***', cex=1.1, font=2)
mtext(c(paste0('Pre-screen\n(',length(pre_contamination),')'), 
        paste0('Post-screen\n(',length(post_contamination),')')), at=c(0.5,1.5), side=1, cex=0.8, padj=1)

par(mar=c(3.5,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,20), xlim=c(0,2), 
     ylab='Contamination (%)', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(pre_contamination, cex=0, lwd=1.8, at=0.5, col='white', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
boxplot(post_contamination, cex=0, lwd=1.8, at=1.5, col='gray40', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(0,20,5), cex.axis=0.8)
mtext(c(paste0('Pre-screen\n(',length(pre_contamination),')'), 
        paste0('Post-screen\n(',length(post_contamination),')')), at=c(0.5,1.5), side=1, cex=0.8, padj=1)
segments(x0=0.5, y0=5, x1=1.5)
text(x=1, y=6, 'n.s.', cex=0.7)

dev.off()

