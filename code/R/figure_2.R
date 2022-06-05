
#------------------------------------------------------------------------------------------------------------#

# Old figure 1 - metaphlan results from metagenomics sequencing
# move to new figure 2

# Read in data
phylum <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/phylum/phylum.tsv', sep='\t', header=TRUE, row.names=1)
family <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/family/family.tsv', sep='\t', header=TRUE, row.names=1)
genus <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/genus/genus.tsv', sep='\t', header=TRUE, row.names=1)
species <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/species.tsv', sep='\t', header=TRUE, row.names=1)

# Calculate 'Other' category
calc_other <- function(shared) {
  shared[shared < 1] <- 0
  shared <- shared[, which(colSums(shared) > 0)]
  shared$Other <- 100 - rowSums(shared)
  return(shared)}
phylum <- calc_other(phylum)
family <- calc_other(family)
genus <- calc_other(genus)
species <- calc_other(species)
rm(calc_other)

# Group by metadata
metadata <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metadata.tsv', sep='\t', header=TRUE)
metadata$status <- NULL
metadata$sample_id <- NULL
metadata$outcome <- NULL
phylum <- merge(metadata, phylum, by.x='name', by.y='row.names')
rownames(phylum) <- phylum$name
phylum$name <- NULL
family <- merge(metadata, family, by.x='name', by.y='row.names')
rownames(family) <- family$name
family$name <- NULL
genus <- merge(metadata, genus, by.x='name', by.y='row.names')
rownames(genus) <- genus$name
genus$name <- NULL
species <- merge(metadata, species, by.x='name', by.y='row.names')
rownames(species) <- species$name
species$name <- NULL


# Subset data
phylum_donor1 <- subset(phylum, donor == 1)
phylum_donor1$donor <- NULL
phylum_donor5 <- subset(phylum, donor == 5)
phylum_donor5$donor <- NULL
phylum_donor28 <- subset(phylum, donor == 28)
phylum_donor28$donor <- NULL
phylum <- rbind(phylum_donor1, rep(0, ncol(phylum_donor1)), phylum_donor5, rep(0, ncol(phylum_donor1)), phylum_donor28)
rm(phylum_donor1, phylum_donor5, phylum_donor28)
family_donor1 <- subset(family, donor == 1)
family_donor1$donor <- NULL
family_donor5 <- subset(family, donor == 5)
family_donor5$donor <- NULL
family_donor28 <- subset(family, donor == 28)
family_donor28$donor <- NULL
family <- rbind(family_donor1, rep(0, ncol(family_donor1)), family_donor5, rep(0, ncol(family_donor1)), family_donor28)
rm(family_donor1, family_donor5, family_donor28)


species_donor1 <- subset(species, donor == 1)
species_donor1$donor <- NULL
species_donor5 <- subset(species, donor == 5)
species_donor5$donor <- NULL
species_donor28 <- subset(species, donor == 28)
species_donor28$donor <- NULL
species <- rbind(species_donor1, species_donor5, species_donor28)
genus_donor1 <- subset(genus, donor == 1)
genus_donor1$donor <- NULL
genus_donor5 <- subset(genus, donor == 5)
genus_donor5$donor <- NULL
genus_donor28 <- subset(genus, donor == 28)
genus_donor28$donor <- NULL
genus <- rbind(genus_donor1, genus_donor5, genus_donor28)

# Transform broader classifications to shared format
phylum <- as.data.frame(t(phylum))
family <- as.data.frame(t(family))
genus <- as.data.frame(t(genus))

# Test for differences at species-level
pvals <- c()
for (x in c(1:ncol(species_donor1))) {
  test <- c(species_donor5[,x], species_donor28[,x])
  pvals[x] <- round(wilcox.test(species_donor1[,x], test, exact=FALSE)$p.value, 3)}
sig_taxa <- c()
for (y in 1:length(pvals)) {if (pvals[y]<=0.01) {sig_taxa <- c(sig_taxa, colnames(species_donor1)[y])}}
super_donor <- c()
normal_donors <- c()
for (z in sig_taxa) {
  test1 <- median(species_donor1[,z])
  test2 <- median(c(species_donor5[,z], species_donor28[,z]))
  ifelse(test1 > test2, super_donor <- c(super_donor, z), normal_donors <- c(normal_donors, z))}
super_donor_abund <- species[,super_donor]
normal_donor_abund <- species[,normal_donors]
rm(test, test1, test2, sig_taxa, pvals, x, y, z, super_donor, normal_donors)

# Subset for plotting
super_donor_abund$type <- c('super','super','super','super','normal','normal','normal','normal','normal','normal','normal','normal')
normal_donor_abund$type <- c('super','super','super','super','normal','normal','normal','normal','normal','normal','normal','normal')
super_donor_hi <- subset(super_donor_abund, type == 'super')
super_donor_hi$type <- NULL
super_donor_lo <- subset(normal_donor_abund, type == 'super')
super_donor_lo$type <- NULL
normal_donor_hi <- subset(normal_donor_abund, type == 'normal')
normal_donor_hi$type <- NULL
normal_donor_lo <- subset(super_donor_abund, type == 'normal')
normal_donor_lo$type <- NULL


# Bar chart funcion
auto_barplot <- function(shared, col_pal, lgnd=TRUE) {
  layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
  par(mar=c(2,2,0.5,0.5), mgp=c(1, 0.25, 0), new=FALSE, xpd=FALSE)
  lty.o <- par('lty')
  par(lty=0)
  barplot(as.matrix(shared), col=col_pal, yaxt='n', xaxt='n', cex.lab=0.7,
          ylim=c(0,100), ylab='Relative Abundance (%)', space=0)
  par(lty=lty.o)
  box(lwd=1.5)
  axis(side=2, at=seq(0,100,20), labels=c('0','20','40','60','80','100'), tick=FALSE, las=1, cex.axis=0.7)
  abline(h=c(20,40,60,80), lty=2, lwd=0.5)
  par(xpd=TRUE)
  text(x=seq(0.5, ncol(shared)-0.5, 1), y=-5, labels=c('1','2','3','4','','1','2','3','4','','1','2','3','4'), 
       cex=0.7, font=2, col=c('#6AA84F','#6AA84F','#6AA84F','#6AA84F', 'black', '#B45F06','#B45F06','#6AA84F','#6AA84F', 'black', '#6AA84F','#B45F06','#6AA84F','#B45F06'))
  text(x=c(2,7,12), y=-11, cex=0.75, font=2, labels=c('Donor A','Donor B','Donor C'))
  par(xpd=FALSE)
  par(mar=c(0,0,0,0))
  plot(0, type='n', ylim=c(0,10), xlim=c(5,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
  if (lgnd == TRUE) {
    legend('topleft', legend=rev(rownames(shared)), pt.bg=col_pal, pch=22, pt.cex=1.5, cex=0.6, bty='n')
    par(font=2)
    legend('bottomleft', legend=c('Success','Failure'), pch=22, pt.cex=0, cex=0.6, bty='n', text.col=c('#6AA84F','#B45F06'))}}

# Generate figures
library(viridis)
taxa_order <- c("Verrucomicrobia", "Bacteroidetes", "Actinobacteria", "Firmicutes", "Other")
phylum <- phylum[taxa_order,]
png(filename='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/phylum.png', units='in', width=5, height=3, res=300)
auto_barplot(phylum, viridis(n=nrow(phylum)))
dev.off()
rm(phylum)

auto_barplot(genus, viridis(n=nrow(genus)))

#require('colortools')
#library(colortools)
#sequential('darkblue', percentage=11, verbose=FALSE, v=0.9)
taxa_order <- c("Verrucomicrobiaceae",
                "Ruminococcaceae","Lachnospiraceae","Eubacteriaceae","Clostridiaceae","Veillonellaceae","Erysipelotrichaceae","Oscillospiraceae","Streptococcaceae",
                "Bifidobacteriaceae","Coriobacteriaceae",
                "Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Prevotellaceae",
                "Other")
family <- family[taxa_order,]
reds <- c('firebrick4','firebrick','firebrick3','firebrick1')
greens <- c('forestgreen','chartreuse3')
blues <- c('#0202E6','#3535E6','#4E4EE6','#6767E6','#8181E6','#9A9AE6','#B3B3E6','#CCCCE6')
new_palette <- c('goldenrod1', rev(blues), rev(greens), rev(reds), 'gray30')
#png(filename='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/family.png', units='in', width=5, height=3, res=300)
pdf(file='~/Desktop/repos/Jenior_FMT_2021/results/family.pdf', width=5, height=3)
auto_barplot(family, new_palette, lgnd=FALSE)
dev.off()


pdf(file='~/Desktop/repos/Jenior_FMT_2021/results/family_legend.pdf', width=2, height=3)
#png(filename='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/family_legend.png', 
#    units='in', width=2, height=3, res=300)
par(mar=c(0,0,0,0))
plot(0, type='n', ylim=c(0,50), xlim=c(0,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend(x=1.6, y=49, legend="Other (<1% each)", pt.bg='gray30', pch=22, pt.cex=1.5, cex=0.6, bty='n')
legend(x=1.6, y=45.5, legend=rev(c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Prevotellaceae")), pt.bg=reds, pch=22, pt.cex=1.5, cex=0.6, bty='n')
legend(x=1.6, y=35, legend=rev(c("Bifidobacteriaceae","Coriobacteriaceae")), pt.bg=greens, pch=22, pt.cex=1.5, cex=0.6, bty='n')
legend(x=1.6, y=29, legend=rev(c("Ruminococcaceae","Lachnospiraceae","Eubacteriaceae","Clostridiaceae","Veillonellaceae","Erysipelotrichaceae","Oscillospiraceae","Streptococcaceae")), 
       pt.bg=blues, pch=22, pt.cex=1.5, cex=0.6, bty='n')
legend(x=1.6, y=10, legend="Verrucomicrobiaceae", pt.bg='goldenrod1', pch=22, pt.cex=1.5, cex=0.6, bty='n')
segments(x0=1.55, y0=c(44.5,34,28,8.75), x1=1.55, y1=c(35.5,29.5,10.5,6.75), lwd=2)
text(x=1.45, y=c(40.25,32,19.5,7.75), labels=c('Bacteroidetes','Actinobacteria','Firmicutes','Verrucomicrobia'), 
     adj=1, font=1, cex=0.5)
dev.off()
rm(family, taxa_order, new_palette)
rm(auto_barplot)


# Species-level ordination
require('vegan')
library(vegan)
species_dist <- vegdist(species, method='bray')
species_nmds <- as.data.frame(metaMDS(species_dist, k=2, trymax=25)$points)
species$donor <- c('super','super','super','super',
                   'normal','normal','normal','normal','normal','normal','normal','normal')
nmds_pval <- adonis(species_dist ~ donor, data=species, perm=999, method='bray')
nmds_pval <- round(nmds_pval$aov.tab[[6]][1], 3)
flux_x <- (abs(max(species_nmds$MDS1)) - abs(min(species_nmds$MDS1))) / 2
flux_y <- (abs(max(species_nmds$MDS2)) - abs(min(species_nmds$MDS2))) / 2
species_nmds$MDS1 <- species_nmds$MDS1 - flux_x
species_nmds$MDS2 <- species_nmds$MDS2 - flux_y
super_nmds_points <- subset(species_nmds, rownames(species_nmds) %in% c('fmt01044A','fmt01090A','fmt01092A','fmt01093A')) #
normal_nmds_points <- subset(species_nmds, rownames(species_nmds) %in% c('fmt05042G','fmt05053B','fmt05080M','fmt05098A','fmt28043C','fmt28045A','fmt28045D','fmt28047D')) #

require('scales')
library(scales)
#png(filename='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species_nmds.png', units='in', width=4.5, height=4, res=300)
#pdf(file='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species_nmds.pdf', width=4.5, height=4)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=species_nmds$MDS1, y=species_nmds$MDS2, xlim=c(-0.4,0.4), ylim=c(-0.4, 0.4),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=super_nmds_points$MDS1, y=super_nmds_points$MDS2, bg=alpha('chartreuse3',0.8), pch=21, cex=1.7)
points(x=normal_nmds_points$MDS1, y=normal_nmds_points$MDS2, bg=alpha('chocolate2',0.8), pch=21, cex=1.7)
legend('topright', legend=c('Super Donor', 'Normal Donor'), 
       pt.bg=c('chartreuse3', 'chocolate2'), pch=21, pt.cex=1.4, cex=0.7, box.lwd=2)
#dev.off()

#-------------------------------------------------------------------------------------------------------------------------------#

# Supervised machine learning
library(randomForest)
set.seed(9861)

# Super vs Normal donors

# Prep for machine learning
species <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/species.tsv', sep='\t', header=TRUE, row.names=1)
metadata <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metadata.tsv', sep='\t', header=TRUE)
metadata$status <- NULL
metadata$sample_id <- NULL
metadata$outcome <- NULL
species <- merge(metadata, species, by.x='name', by.y='row.names')
rownames(species) <- species$name
species$name <- NULL
super_donor_abund <- subset(species, donor == 1)
super_donor_abund$donor <- NULL
normal_donor_abund <- subset(species, donor %in% c(5,28))
normal_donor_abund$donor <- NULL
# Diversity
library(vegan)
super_div <- as.vector(diversity(super_donor_abund, index='invsimpson'))
normal_div <- as.vector(diversity(normal_donor_abund, index='invsimpson'))
pval <- wilcox.test(super_div, normal_div, exact=FALSE)$p.value
png(filename='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/super_inv_simpson.png', units='in', width=3, height=5, res=300)
ymax <- round(max(c(max(super_div), max(normal_div))) * 1.2)
par(mar=c(3, 3, 0.5, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
plot(0, type='n', xlab='', ylab='Inv. Simpson Diversity', xaxt='n', yaxt='n', 
     xlim=c(0,1), ylim=c(0,ymax))
stripchart(at=0.25, vertical=TRUE, jitter(super_div, amount=1e-5), lwd=2,
           pch=21, bg='blue3', method='jitter', jitter=0.12, cex=2, add=TRUE)
stripchart(at=0.75, vertical=TRUE, jitter(normal_div, amount=1e-5), lwd=2,
           pch=21, bg='green3', method='jitter', jitter=0.12, cex=2, add=TRUE)
axis(side=2, at=seq(0,ymax,ymax/5), cex.axis=0.8, lwd=2)
segments(x0=0.1, y0=median(super_div), x1=0.4, lwd=3)
segments(x0=0.6, y0=median(normal_div), x1=0.9, lwd=3)
segments(x0=0.25, y0=ymax*0.92, x1=0.75, lwd=2)
text(x=0.5, y=ymax*0.95, 'n.s.')
box(lwd=2)
par(xpd=TRUE)
text(x=0.25, y=-(ymax*0.07), labels='Super\nDonor', cex=1.2)
text(x=0.75, y=-(ymax*0.07), labels='Normal\nDonors', cex=1.2)
par(xpd=FALSE)
dev.off()


# Run random forest
species$donor <- gsub(1, 'super', species$donor)
species$donor <- gsub(5, 'normal', species$donor)
species$donor <- gsub(28, 'normal', species$donor)
condition <- as.factor(species$donor)
species$donor <- NULL
species <- droplevels(species)
rf_obj <- randomForest(condition ~ ., data=species, importance=TRUE, err.rate=TRUE, na.action=na.roughfix, 
                       ntree=10000, mtry=15)
# Number of trees: 1500
# No. of variables tried at each split: 15
# OOB estimate of  error rate: 0%
# Confusion matrix:
#             normal super       class.error
# normal      8      0           0
# super       0      4           0
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- subset(rf_obj, rf_obj > (abs(min(rf_obj)))) # significance
rf_mda <- as.data.frame(rf_mda)
rf_mda$feature <- rownames(rf_mda)

# Test for differences
sig <- c()
increased <- c()
abs_diff <- c()
for (x in rf_mda$feature) {
  abs_diff <- c(abs_diff, round(abs(median(super_donor_abund[,x]) - median(normal_donor_abund[,x])), 2))
  pval <- round(wilcox.test(super_donor_abund[,x], normal_donor_abund[,x], exact=FALSE)$p.value, 3)
  if (pval < 0.05) {
    sig <- c(sig, '*')
  } else {sig <- c(sig, 'non_sig')}
  
  super_med <- median(super_donor_abund[,x])
  normal_med <- median(normal_donor_abund[,x])
  if (super_med < normal_med) {
    increased <- c(increased, 'normal')
  } else {increased <- c(increased, 'super')}
}

# Reformat
colnames(rf_mda) <- c('mda', 'species')
rf_mda$mda <- round(rf_mda$mda, 3)
rf_mda$increased_in <- increased
rf_mda$wilcox_sig <- sig
rf_mda$med_abs_diff <- abs_diff
rf_mda_super <- subset(rf_mda, med_abs_diff > 0)


# Read in saved results
rf_mda_super <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/super_normal.sig_RF.tsv', sep='\t', header=TRUE)
#super_donor_abund <- super_donor_abund[, rf_mda_super$species]
#normal_donor_abund <- normal_donor_abund[, rf_mda_super$species]
#colnames(super_donor_abund) <- gsub('_', ' ', colnames(super_donor_abund))
#colnames(normal_donor_abund) <- gsub('_', ' ', colnames(normal_donor_abund))

# Format
rf_mda_super$species <- gsub('_', ' ', rf_mda_super$species)
rf_mda_super <- rf_mda_super[order(rf_mda_super$mda),]
point_col <- c()
for (x in rf_mda_super$increased_in) {
  
  if (x == 'normal') {
    point_col <- c(point_col, 'black')
  } else {point_col <- c(point_col, 'white')}
}
rf_mda_super$point_col <- point_col
rf_mda_super$mda <- rf_mda_super$mda * 2.0
rf_mda_super <- subset(rf_mda_super, mda > 5)

# Fig 1C
#pdf(file='~/Desktop/repos/Jenior_FMT_2021/results/figure_1c.pdf', width=3.5, height=6)
png(file='~/Desktop/repos/Jenior_FMT_2021/results/figure_1c.png', units='in', width=3.5, height=6, res=300)
par(mar=c(2.5, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(rf_mda_super$mda,  xlab='Mean Decrease Accuracy', xlim=c(0,max(rf_mda_super$mda)*1.3),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=0.1, cex=0.8)
points(x=rf_mda_super$mda, y=c(1:nrow(rf_mda_super)), pch=23, bg=rf_mda_super$point_col, cex=1.2)
legend('bottomright', legend=c('Sig. increased in:','Donor A', 'Donors B & C'), col=c('white', 'black', 'black'),
       pt.bg=c('white', 'white', 'black'), pch=23, pt.cex=1.2, cex=0.7, bg='white')
text(x=0, y=seq(1.3,nrow(rf_mda_super)+0.3,1), labels=rf_mda_super$species, cex=0.75, pos=4, font=4, col=as.character(rf_mda_super$tax_col))
#rect(xleft=3.55, xright=5, ytop=4.7, ybottom=3.5, col='white', border='white')
#text(x=3.5, y=4, labels='Out of Bag = 0%', pos=4, cex=0.7)
dev.off()

#---------------------------------#

# Success vs Failure - supplement

# Prep for machine learning
species <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/species.tsv', sep='\t', header=TRUE, row.names=1)
metadata <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metadata.tsv', sep='\t', header=TRUE)
metadata$sample_id <- NULL
metadata$status <- NULL
metadata$donor <- NULL
species <- merge(metadata, species, by.x='name', by.y='row.names')
rownames(species) <- species$name
species$name <- NULL
success_donor_abund <- subset(species, outcome == 'success')
success_donor_abund$outcome <- NULL
failure_donor_abund <- subset(species, outcome == 'failure')
failure_donor_abund$outcome <- NULL

# Diversity
library(vegan)
success_div <- as.vector(diversity(success_donor_abund, index='invsimpson'))
failure_div <- as.vector(diversity(failure_donor_abund, index='invsimpson'))
pval <- wilcox.test(success_div, normal_div, exact=FALSE)$p.value
png(filename='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/success_inv_simpson.png', units='in', width=3, height=5, res=300)
ymax <- round(max(c(max(success_div), max(failure_div))) * 1.2)
par(mar=c(3, 3, 0.5, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
plot(0, type='n', xlab='', ylab='Inv. Simpson Diversity', xaxt='n', yaxt='n', 
     xlim=c(0,1), ylim=c(0,ymax))
stripchart(at=0.25, vertical=TRUE, jitter(success_div, amount=1e-5), lwd=2,
           pch=21, bg='blue3', method='jitter', jitter=0.12, cex=2, add=TRUE)
stripchart(at=0.75, vertical=TRUE, jitter(failure_div, amount=1e-5), lwd=2,
           pch=21, bg='green3', method='jitter', jitter=0.12, cex=2, add=TRUE)
axis(side=2, at=seq(0,ymax,ymax/5), cex.axis=0.8, lwd=2)
segments(x0=0.1, y0=median(success_div), x1=0.4, lwd=3)
segments(x0=0.6, y0=median(failure_div), x1=0.9, lwd=3)
segments(x0=0.25, y0=ymax*0.92, x1=0.75, lwd=2)
text(x=0.5, y=ymax*0.95, 'n.s.')
box(lwd=2)
par(xpd=TRUE)
text(x=0.25, y=-(ymax*0.07), labels='Successful\nDonations', cex=0.9)
text(x=0.75, y=-(ymax*0.07), labels='Failed\nDonations', cex=0.9)
par(xpd=FALSE)
dev.off()


# Run random forest
condition <- as.factor(species$outcome)
species$outcome <- NULL
species <- droplevels(species)
rf_obj <- randomForest(condition ~ ., data=species, importance=TRUE, err.rate=TRUE, na.action=na.roughfix, ntree=5000, mtry=15)
# Number of trees: 1500
# No. of variables tried at each split: 15
# OOB estimate of  error rate: 66.67%
#Confusion matrix:
#               failure success   class.error
# failure       0       4         1.0
# success       4       4         0.5
rf_mda <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- subset(rf_mda, rf_mda > (abs(min(rf_mda)))) # significance
rf_mda <- as.data.frame(rf_mda)
rf_mda$feature <- rownames(rf_mda)

# Test for differences
sig <- c()
increased <- c()
abs_diff <- c()
for (x in rf_mda$feature) {
  abs_diff <- c(abs_diff, round(abs(median(success_donor_abund[,x]) - median(failure_donor_abund[,x])), 2))
  pval <- round(wilcox.test(success_donor_abund[,x], failure_donor_abund[,x], exact=FALSE)$p.value, 3)
  if (pval < 0.05) {
    sig <- c(sig, '*')
  } else {sig <- c(sig, 'non_sig')}
  
  success_med <- median(success_donor_abund[,x])
  failure_med <- median(failure_donor_abund[,x])
  if (success_med < failure_med) {
    increased <- c(increased, 'failure')
  } else {increased <- c(increased, 'success')}
}

# Reformat
colnames(rf_mda) <- c('mda', 'species')
rf_mda$mda <- round(rf_mda$mda, 3)
rf_mda$increased_in <- increased
rf_mda$wilcox_sig <- sig
rf_mda$med_abs_diff <- abs_diff
rf_mda_success <- subset(rf_mda, med_abs_diff > 0)
rm(x, species, metadata, condition, rf_mda, rf_obj, pval, sig, increased, 
   abs_diff, super_med, normal_med, success_med, failure_med)

# Save results
#write.table(rf_mda_super, file='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/super_normal.sig_RF.tsv', quote=FALSE, sep='\t', row.names=FALSE)
#write.table(rf_mda_success, file='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/success_failure.sig_RF.tsv', quote=FALSE, sep='\t', row.names=FALSE)
rf_mda_super <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/super_normal.sig_RF.tsv', sep='\t', header=TRUE)
#rf_mda_super$species <- gsub('_unclassified', '', rf_mda_super$species)
rf_mda_super$species <- gsub('_', ' ', rf_mda_super$species)
rf_mda_success <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/success_failure.sig_RF.tsv', sep='\t', header=TRUE)
#rf_mda_success$species <- gsub('_unclassified', '', rf_mda_success$species)
rf_mda_success$species <- gsub('_', ' ', rf_mda_success$species)


# Flux samples from informative features
png(filename='~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/figure_S__.png', 
    units='in', width=2.5, height=4, res=300)
par(mar=c(2.5, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(all_aucrf$importance,  xlab='Importance', xlim=c(0,0.6),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=0.1, cex=0.8)
text(x=-0.025, y=seq(1.3,10.3,1), labels=all_aucrf$species, cex=0.75, pos=4, font=4,
     col=c('#B3B3E6','#B3B3E6','chartreuse3','chartreuse3','firebrick1','#3535E6','firebrick1','#9A9AE6','#B3B3E6','firebrick'))
text(x=0.275, y=6.3, labels='unclassified', cex=0.75, pos=4, font=2, col='#3535E6')
points(x=all_aucrf$importance, y=c(1:10), pch=23, 
       bg=rev(c('white','white','white','white','black','white','black','white','black','white')), cex=1.2)
legend('bottomright', legend=c('Increased in:','Donor A', 'Donors B & C','p-value < 0.05'), col=c('white', 'black', 'black','white'),
       pt.bg=c('white', 'white', 'black', 'white'), pch=23, pt.cex=1.2, cex=0.7)
dev.off()


#-------------------------------------------------------------------------------------------------------------------#

# Inv. Simpson Diversity
library(vegan)

# Super vs Normal
species <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/species.tsv', sep='\t', header=TRUE, row.names=1)
metadata <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metadata.tsv', sep='\t', header=TRUE)
metadata$status <- NULL
metadata$sample_id <- NULL
metadata$outcome <- NULL
species <- merge(metadata, species, by.x='name', by.y='row.names')
rownames(species) <- species$name
species$name <- NULL
super_donor_abund <- subset(species, donor == 1)
super_donor_abund$donor <- NULL
normal_donor_abund <- subset(species, donor %in% c(5,28))
normal_donor_abund$donor <- NULL
super_div <- as.vector(diversity(super_donor_abund, index='invsimpson'))
normal_div <- as.vector(diversity(normal_donor_abund, index='invsimpson'))
div_pval <- round(wilcox.test(super_div, normal_div, exact=FALSE)$p.value, 4)

# Success vs Failure
species <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metaphlan/species/species.tsv', sep='\t', header=TRUE, row.names=1)
metadata <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metadata.tsv', sep='\t', header=TRUE)
metadata$status <- NULL
metadata$sample_id <- NULL
metadata$donor <- NULL
species <- merge(metadata, species, by.x='name', by.y='row.names')
rownames(species) <- species$name
species$name <- NULL
success_donor_abund <- subset(species, outcome == 'success')
success_donor_abund$outcome <- NULL
failure_donor_abund <- subset(species, outcome == 'failure')
failure_donor_abund$outcome <- NULL
success_div <- as.vector(diversity(success_donor_abund, index='invsimpson'))
failure_div <- as.vector(diversity(failure_donor_abund, index='invsimpson'))
success_div_pval <- round(wilcox.test(success_div, failure_div, exact=FALSE)$p.value, 4)

