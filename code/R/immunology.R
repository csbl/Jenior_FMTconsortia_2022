
# Read in 
immune <- read.delim(file='~/Desktop/Jenior_Consortia_2022/data/JL66_regating.tsv', sep='\t', header=TRUE)
immune$mouse <- NULL
immune$well <- NULL

# subset data
control <- subset(immune, treatment_group == 'control')
control$treatment_group <- NULL
competitive <- subset(immune, treatment_group == 'competitive')
competitive$treatment_group <- NULL
cooperative <- subset(immune, treatment_group == 'cooperative')
cooperative$treatment_group <- NULL
rm(immune)

# Normalize by live cells
control_live <- control$Live
control$Live <- NULL
competitive_live <- competitive$Live
competitive$Live <- NULL
cooperative_live <- cooperative$Live
cooperative$Live <- NULL
for (x in colnames(control)) {
  control[,x] <- (control[,x] / control_live) * 100
  competitive[,x] <- (competitive[,x] / competitive_live) * 100
  cooperative[,x] <- (cooperative[,x] / cooperative_live) * 100}
rm(x, control_live, competitive_live, cooperative_live)

# Test differences
control_competitive <- c()
control_cooperative <- c()
competitive_cooperative <- c()
for (x in colnames(control)) {
  pval <- round(wilcox.test(control[,x], competitive[,x], exact=FALSE)$p.value, 4)
  control_competitive <- c(control_competitive, pval)
  pval <- round(wilcox.test(control[,x], cooperative[,x], exact=FALSE)$p.value, 4)
  control_cooperative <- c(control_cooperative, pval)
  pval <- round(wilcox.test(competitive[,x], cooperative[,x], exact=FALSE)$p.value, 4)
  competitive_cooperative <- c(competitive_cooperative, pval)}
pvalues <- as.data.frame(cbind(control_competitive, control_cooperative, competitive_cooperative))
colnames(pvalues) <- c('control_competitive','control_cooperative','competitive_cooperative')
rownames(pvalues) <- colnames(control)
rm(x, pval, control_competitive, control_cooperative, competitive_cooperative)

# Generate figures
immune_plot <- function(celltype, formatted_name='default', sig_lab='default') {
  
  pval1 <- pvalues[celltype,'control_competitive']
  pval2 <- pvalues[celltype,'control_cooperative']
  pval3 <- pvalues[celltype,'competitive_cooperative']
  
  control_curr <- control[,celltype]
  competitive_curr <- competitive[,celltype]
  cooperative_curr <- cooperative[,celltype]
  
  ymax <- max(max(c(control_curr, competitive_curr, cooperative_curr)) * 1.2)
  axis_labs <- round(seq(0, ymax, length.out=5), 1)
  
  if (formatted_name == 'default') {
    formatted_name <- paste(celltype, '(% of live)')
    } else {formatted_name <- paste(formatted_name, '(% of live)')}
  
  par(mar=c(4.5,3,0.5,0.5), las=1, mgp=c(1.8,0.7,0), xpd=FALSE, lwd=2)
  plot(0, type='n', ylim=c(0, ymax), xlim=c(0,3), 
       ylab=formatted_name, xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
  boxplot(control_curr, cex=0, at=0.5, col='gray40', ylab='', staplewex=0.6,
          boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
  boxplot(competitive_curr, cex=0, at=1.5, col='dodgerblue4', ylab='', staplewex=0.6, 
          boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
  boxplot(cooperative_curr, cex=0, at=2.5, col='firebrick3', ylab='', staplewex=0.6, 
          boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
  
  axis(side=2, at=axis_labs, cex.axis=0.7, lwd=2)
  
  segments(x0=0.5, y0=ymax*0.87, x1=1.5)
  segments(x0=0.5, y0=ymax*0.92, x1=2.5)
  segments(x0=1.5, y0=ymax*0.97, x1=2.5)
  par(xpd=TRUE)
  text(x=c(0.5,1.5,2.5), y=-(ymax*0.22), srt=40, cex=1, 
       labels=c('Control\ngavage','Competitive\nconsortia','Cooperative\nconsortia'))
  
  if (sig_lab == 'default') {sig_lab <- '*'}
  if (pval1 <= 0.05) {text(x=1, y=ymax*0.893, sig_lab, font=2, cex=1.2)
  } else {text(x=1, y=ymax*0.893, 'n.s.', cex=0.7)}
  
  if (pval2 <= 0.05) {text(x=1.5, y=ymax*0.943, sig_lab, font=2, cex=1.2)
  } else {text(x=1.5, y=ymax*0.943, 'n.s.', cex=0.7)}
  
  if (pval3 <= 0.05) {text(x=2, y=ymax*0.993, sig_lab, font=2, cex=1.2)
  } else {text(x=2, y=ymax*0.993, 'n.s.', cex=0.7)}
}

pdf(file='~/Desktop/Jenior_Consortia_2022/results/immunology.pdf', width=2, height=8)
layout(matrix(c(1,2,3), nrow=3, ncol=1, byrow=TRUE))
immune_plot('Eosinophils', sig_lab='**')
immune_plot('Neutrophils')
immune_plot('ILC3', 'ILC3 cells')
dev.off()

pdf(file='~/Desktop/Jenior_Consortia_2022/results/tregs.pdf', width=2, height=3)
immune_plot('Tregs')
dev.off()

