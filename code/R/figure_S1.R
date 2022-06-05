
# Beta-diversity 

# Read in data
species <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metaphlan/species/species.tsv', sep='\t', header=TRUE, row.names=1)
metadata <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metadata.tsv', sep='\t', header=TRUE)
metadata$status <- NULL
metadata$sample_id <- NULL
species <- merge(metadata, species, by.x='name', by.y='row.names')
rownames(species) <- species$name
species$name <- NULL
rm(metadata)

# Subset data
super <- subset(species, donor == 1)
super$donor <- NULL
super$outcome <- NULL
normal <- subset(species, donor %in% c(5,28))
normal$donor <- NULL
normal$outcome <- NULL
success <- subset(species, outcome == 'success')
success$donor <- NULL
success$outcome <- NULL
failure <- subset(species, outcome == 'failure')
failure$donor <- NULL
failure$outcome <- NULL
rm(species)

# Calculate diversity
library(vegan)
super_div <- as.vector(diversity(super, index='invsimpson'))
normal_div <- as.vector(diversity(normal, index='invsimpson'))
success_div <- as.vector(diversity(success, index='invsimpson'))
failure_div <- as.vector(diversity(failure, index='invsimpson'))

# Calucate differences
pval1 <- round(wilcox.test(super_div, normal_div, exact=FALSE)$p.value, 4)
pval2 <- round(wilcox.test(success_div, failure_div, exact=FALSE)$p.value, 4)
pval3 <- round(wilcox.test(super_div, failure_div, exact=FALSE)$p.value, 4)

# Calucalte summary stats
super_q <- as.vector(c(quantile(super_div, 0.25), quantile(super_div, 0.5), quantile(super_div, 0.75)))
normal_q <- as.vector(c(quantile(normal_div, 0.25), quantile(normal_div, 0.5), quantile(normal_div, 0.75)))
success_q <- as.vector(c(quantile(success_div, 0.25), quantile(success_div, 0.5), quantile(success_div, 0.75)))
failure_q <- as.vector(c(quantile(failure_div, 0.25), quantile(failure_div, 0.5), quantile(failure_div, 0.75)))

# Plot results
plotDiversity <- function(group1, group2, groups, clrs, pval) {
  ymax <- round(max(c(max(group1), max(group2))) * 1.2)
  par(mar=c(2.5, 3, 0.5, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=3, las=1)
  barplot(c(group1[2], group2[2]), col=clrs, yaxt='n', ylim=c(0,ymax), 
          ylab='Inv. Simpson Diversity', names.arg=groups, cex.names=1.1)
  axis(side=2, at=seq(0,ymax,ymax/5), cex.axis=0.8, lwd=3)
  text(x=1.9, y=ymax*0.9, bquote(paste(italic(p), ' = ', .(pval))), cex=0.8)
  box(lwd=3)
  segments(x0=0.4, y0=group1[1], x1=1, lwd=3)
  segments(x0=0.4, y0=group1[3], x1=1, lwd=3)
  segments(x0=0.7, y0=group1[1], y1=group1[3], lwd=3)
  segments(x0=1.6, y0=group2[1], x1=2.2, lwd=3)
  segments(x0=1.6, y0=group2[3], x1=2.2, lwd=3)
  segments(x0=1.9, y0=group2[1], y1=group2[3], lwd=3)}

super_col <- '#8e7cc3'
normal_col <- '#e69138'
success_col <- '#6aa84f'
failure_col <- '#990000'

pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_S1AB.pdf', width=2.7, height=5)
layout(matrix(c(1, 2), nrow=2, ncol=1, byrow=TRUE))
plotDiversity(super_q, normal_q, groups=c('Super','Normal'), clrs=c(super_col,normal_col), pval1)
plotDiversity(success_q, failure_q, groups=c('Success','Failure'), clrs=c(success_col,failure_col), pval2)
dev.off()

#-------------------------------------------------------------------------------------------------------------------------------#

# Supervised machine learning - Super vs Normal

# Prep for machine learning
species <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metaphlan/species/species.tsv', sep='\t', header=TRUE, row.names=1)
metadata <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metadata.tsv', sep='\t', header=TRUE)
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

# Run random forest
species$donor <- gsub(1, 'super', species$donor)
species$donor <- gsub(5, 'normal', species$donor)
species$donor <- gsub(28, 'normal', species$donor)
condition <- as.factor(species$donor)
species$donor <- NULL
species <- droplevels(species)
library(randomForest)
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
#rf_mda <- subset(rf_obj, rf_obj > (abs(min(rf_obj)))) # significance
rf_mda <- subset(rf_obj, rf_obj > 0)
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

# Reformat and save
colnames(rf_mda) <- c('mda', 'species')
rf_mda$mda <- round(rf_mda$mda, 3)
rf_mda$increased_in <- increased
rf_mda$wilcox_sig <- sig
rf_mda$med_abs_diff <- abs_diff
rf_mda_super <- subset(rf_mda, med_abs_diff > 0)
write.table(rf_mda, file='~/Desktop/repos/Jenior_Consortia_2022/data/metaphlan/species/super_normal.RF.tsv', quote=FALSE, sep='\t', row.names=FALSE)






#------------#

# Supervised machine learning - Success vs Failure

# Prep for machine learning
species <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metaphlan/species/species.tsv', sep='\t', header=TRUE, row.names=1)
metadata <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metadata.tsv', sep='\t', header=TRUE)
metadata$sample_id <- NULL
metadata$donor <- NULL
species <- merge(metadata, species, by.x='name', by.y='row.names')
rownames(species) <- species$name
species$name <- NULL
success_donor_abund <- subset(species, outcome == 'success')
success_donor_abund$outcome <- NULL
success_donor_abund$status <- NULL
failure_donor_abund <- subset(species, outcome == 'failure')
failure_donor_abund$outcome <- NULL
failure_donor_abund$status <- NULL
super_donor_abund <- subset(species, status == 'super')
super_donor_abund$outcome <- NULL
super_donor_abund$status <- NULL
normal_donor_abund <- subset(species, status == 'normal')
normal_donor_abund$outcome <- NULL
normal_donor_abund$status <- NULL
species$status <- NULL

# Reformat
condition <- as.factor(species$outcome)
species$outcome <- NULL
species <- droplevels(species)
species <- as.data.frame(t(species))
for (x in colnames(species)) {species[,x] <- as.numeric(as.character(species[,x]))}
species <- subset(species, rowSums(species) > 0)

# Create 3 quantiles
for (x in colnames(species)) {
  grouping <- as.integer(cut(rank(species[,x], ties.method='first'), 
                             quantile(rank(species[,x], ties.method='first'), 
                                      probs=0:3/3, na.rm=TRUE), include.lowest=TRUE))
  species[,x] <- as.factor(grouping)}; rm(grouping)
species <- as.data.frame(t(species))

# Subset by super vs normal results
rf_mda_super <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metaphlan/species/super_normal.RF.tsv', sep='\t', header=TRUE)
species <- species[,rf_mda_super$species]
rm(rf_mda_super)

# Run random forest
library(randomForest)
set.seed(9861)
rf_obj <- randomForest(condition ~ ., data=species, importance=TRUE, err.rate=TRUE, na.action=na.roughfix, ntree=5000, mtry=15)
print(rf_obj)
# OOB estimate of  error rate: 18.75%
# Confusion matrix:
#   failure success class.error
# failure       8       0       0.000
# success       3       5       0.375
rf_mda <- importance(rf_obj, type=1, scale=TRUE)
#rf_mda <- subset(rf_mda, rf_mda > abs(min(rf_mda)))
rf_mda <- subset(rf_mda, rf_mda > 0)
rf_mda <- as.data.frame(rf_mda)
rf_mda$feature <- rownames(rf_mda)
rf_mda <- rf_mda[order(rf_mda$MeanDecreaseAccuracy),] 

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
rf_mda_success <- subset(rf_mda, mda > 1)
rm(x, species, metadata, condition, rf_obj, rf_mda, pval, sig, increased, 
   abs_diff, success_med, failure_med)

# Save results
#write.table(rf_mda_success, file='~/Desktop/repos/Jenior_Consortia_2022/data/metaphlan/species/success_failure.sig_RF.tsv', quote=FALSE, sep='\t', row.names=FALSE)

# Read manually reformatted results
rf_mda_success <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metaphlan/species/success_failure.sig_RF.tsv', sep='\t', header=TRUE)
#rf_mda_success$species <- gsub('_unclassified', '', rf_mda_success$species)
rf_mda_success$species <- gsub('_', ' ', rf_mda_success$species)
rf_mda_success$increased_in <- gsub('success', 'black', rf_mda_success$increased_in)
rf_mda_success$increased_in <- gsub('failure', 'white', rf_mda_success$increased_in)

# MDA plot
pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_1c.pdf', width=3.75, height=5)
par(mar=c(2.5, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(rf_mda_success$mda,  xlab='Mean Decrease Accuracy (%)', xlim=c(0,ceiling(max(rf_mda_success$mda))),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=0.1, cex=0.8)
points(x=rf_mda_success$mda, y=c(1:nrow(rf_mda_success)), pch=23, bg=rf_mda_success$increased_in, cex=1.2)
text(x=-0.025, y=seq(1.3,nrow(rf_mda_success)+0.3,1), labels=rf_mda_success$species, cex=0.75, pos=4, font=4, col=as.character(rf_mda_success$tax_col))
legend('right', legend=c('Success','Failure'), col='black', bg='white',
       pt.bg=c('black','white'), pch=23, pt.cex=1.3, cex=0.9)
legend('bottomright', legend=c('Bacteroidetes','Firmicutes','Proteobacteria','Actinobacteria'), text.font=2, bg='white',
       text.col=c('firebrick','dodgerblue','darkorange','forestgreen'), pt.cex=0, cex=0.9)
dev.off()

#-------------------------------------------------------------------------------------------------------------------------------#

# Plot abundance of important pathways

# Read in data
# Pathway abundances
s01044A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.01044A.tsv', sep='\t', header=TRUE)
s01090A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.01090A.tsv', sep='\t', header=TRUE)
s01092A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.01092A.tsv', sep='\t', header=TRUE)
s01093A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.01093A.tsv', sep='\t', header=TRUE)
n05042G <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.05042G.tsv', sep='\t', header=TRUE)
n05053B <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.05053B.tsv', sep='\t', header=TRUE)
n05080M <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.05080M.tsv', sep='\t', header=TRUE)
n05098A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.05098A.tsv', sep='\t', header=TRUE)
n28043C <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.28043C.tsv', sep='\t', header=TRUE)
n28045A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.28045A.tsv', sep='\t', header=TRUE)
n28045D <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.28045D.tsv', sep='\t', header=TRUE)
n28047D <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/abundance/pathway_abundances.28047D.tsv', sep='\t', header=TRUE)

# Assemble tables
super <- merge(s01044A, s01090A, by='pathway')
super <- merge(super, s01092A, by='pathway')
super <- merge(super, s01093A, by='pathway')
rownames(super) <- super$pathway
super$pathway <- NULL
colnames(super) <- c('s01044A','s01090A','s01092A','s01093A')
super <- as.data.frame(t(super))
normal <- merge(n05042G, n05053B, by='pathway')
normal <- merge(normal, n05080M, by='pathway')
normal <- merge(normal, n05098A, by='pathway')
normal <- merge(normal, n28043C, by='pathway')
normal <- merge(normal, n28045A, by='pathway')
normal <- merge(normal, n28045D, by='pathway')
normal <- merge(normal, n28047D, by='pathway')
rownames(normal) <- normal$pathway
normal$pathway <- NULL
colnames(normal) <- c('n05042G','n05053B','n05080M','n05098A','n28043C','n28045A','n28045D','n28047D')
normal <- as.data.frame(t(normal))
success <- merge(s01044A, s01090A, by='pathway')
success <- merge(success, s01092A, by='pathway')
success <- merge(success, s01093A, by='pathway')
success <- merge(success, n05098A, by='pathway')
success <- merge(success, n05080M, by='pathway')
success <- merge(success, n28043C, by='pathway')
success <- merge(success, n28045D, by='pathway')
rownames(success) <- success$pathway
success$pathway <- NULL
colnames(success) <- c('s01044A','s01090A','s01092A','s01093A','n05098A','n05080M','n28043C','n28045D')
success <- as.data.frame(t(success))
failure <- merge(n05042G, n05053B, by='pathway')
failure <- merge(failure, n28045A, by='pathway')
failure <- merge(failure, n28047D, by='pathway')
rownames(failure) <- failure$pathway
failure$pathway <- NULL
colnames(failure) <- c('n05042G','n05053B','n28045A','n28047D')
failure <- as.data.frame(t(failure))
rm(s01044A,s01090A,s01092A,s01093A,n05098A,n05080M,
   n28043C,n28045D,n05042G,n05053B,n28045A,n28047D)

# Focus on non-metabolic pathways
metabolic_paths <- c("2-Oxocarboxylic_acid_metabolism","ABC_transporters","Alanine,_aspartate_and_glutamate_metabolism",
                     "Amino_sugar_and_nucleotide_sugar_metabolism","Arginine_and_proline_metabolism","Ascorbate_and_aldarate_metabolism",
                     "beta-Alanine_metabolism","Carbohydrate_digestion_and_absorption","Glyoxylate_and_dicarboxylate_metabolism",
                     "Carbon_metabolism","Citrate_cycle_(TCA_cycle)","Cyanoamino_acid_metabolism","Cysteine_and_methionine_metabolism",
                     "D-Alanine_metabolism","D-Arginine_and_D-ornithine_metabolism","D-Glutamine_and_D-glutamate_metabolism","Degradation_of_aromatic_compounds",
                     "Fatty_acid_degradation","Fatty_acid_metabolism","Fructose_and_mannose_metabolism","Galactose_metabolism",
                     "Glutathione_metabolism","Glycerolipid_metabolism","Glycerophospholipid_metabolism","Glycine,_serine_and_threonine_metabolism","Glycolysis_Gluconeogenesis",
                     "Glycosaminoglycan_degradation","Histidine_metabolism","Lysine_degradation","Glycosaminoglycan_biosynthesis",
                     "Microbial_metabolism_in_diverse_environments","Mineral_absorption","Lipoic_acid_metabolism",
                     "Nicotinate_and_nicotinamide_metabolism","Nitrogen_metabolism","Biosynthesis_of_unsaturated_fatty_acids",
                     "Pentose_and_glucuronate_interconversions","Pentose_phosphate_pathway","Phenylalanine_metabolism",
                     "Phosphonate_and_phosphinate_metabolism","Phosphotransferase_system_(PTS)",
                     "Polycyclic_aromatic_hydrocarbon_degradation","Propanoate_metabolism","Protein_digestion_and_absorption",
                     "Purine_metabolism","Pyrimidine_metabolism","Pyruvate_metabolism","Retinol_metabolism","Riboflavin_metabolism",
                     "Selenocompound_metabolism","Sphingolipid_metabolism","Starch_and_sucrose_metabolism",
                     "Sulfur_metabolism","Synthesis_and_degradation_of_ketone_bodies","Taurine_and_hypotaurine_metabolism",
                     "Thiamine_metabolism","Toluene_degradation","Tryptophan_metabolism","Tyrosine_metabolism",
                     "Valine,_leucine_and_isoleucine_degradation","Vitamin_B6_metabolism","Flavonoid_biosynthesis",
                     "Vitamin_digestion_and_absorption","Xylene_degradation","Arginine_biosynthesis")
super[,metabolic_paths] <- NULL
normal[,metabolic_paths] <- NULL
success[,metabolic_paths] <- NULL
failure[,metabolic_paths] <- NULL

# Subsample
library(vegan)
super <- as.data.frame(t(super))
normal <- as.data.frame(t(normal))
success <- as.data.frame(t(success))
failure <- as.data.frame(t(failure))
sub_donor <- ceiling(min(c(min(as.vector(colSums(super))),min(as.vector(colSums(normal))))) * 0.85)
for (x in 1:ncol(super)) {super[,x] <- as.vector(rrarefy(super[,x], sample=sub_donor))}
for (y in 1:ncol(normal)) {normal[,y] <- as.vector(rrarefy(normal[,y], sample=sub_donor))}
sub_donor <- ceiling(min(c(min(as.vector(colSums(success))),min(as.vector(colSums(failure))))) * 0.85)
for (x in 1:ncol(success)) {success[,x] <- as.vector(rrarefy(success[,x], sample=sub_donor))}
for (y in 1:ncol(failure)) {failure[,y] <- as.vector(rrarefy(failure[,y], sample=sub_donor))}

# Build final tables
donor <- merge(super, normal, by='row.names')
rownames(donor) <- donor$Row.names
donor$Row.names <- NULL
donor <- subset(donor, rowSums(donor) >= 3)
donor <- as.data.frame(t(donor))
colnames(donor) <- make.names(colnames(donor))
outcome <- merge(success, failure, by='row.names')
rownames(outcome) <- outcome$Row.names
outcome$Row.names <- NULL
outcome <- subset(outcome, rowSums(outcome) >= 3)
outcome <- as.data.frame(t(outcome))
colnames(outcome) <- make.names(colnames(outcome))
super <- as.data.frame(t(super))
normal <- as.data.frame(t(normal))
success <- as.data.frame(t(success))
failure <- as.data.frame(t(failure))

# Supervised machine learning
library(randomForest)
set.seed(906801)
condition <- as.factor(c('super','super','super','super','normal','normal','normal','normal','normal','normal','normal','normal'))
rf_donor <- randomForest(condition ~ ., data=donor, ntree=5000, mtry=15,
                         importance=TRUE, replace=FALSE, err.rate=TRUE)
# OOB estimate of  error rate: 16.67%
# Confusion matrix:
#   normal super class.error
# normal      8     0         0.0
# super       2     2         0.5
rf_donor <- importance(rf_donor, type=1, scale=FALSE)
rf_donor <- subset(rf_donor, rf_donor > abs(min(rf_donor)))
rf_donor <- as.data.frame(rf_donor)
rf_donor$feature <- rownames(rf_donor)
colnames(rf_donor) <- c('mda', 'pathway')
rf_donor <- rf_donor[order(-rf_donor$mda),]
if (nrow(rf_donor > 20)) {rf_donor <- rf_donor[1:20,]}
donor <- donor[,which(colnames(donor) %in% rf_donor$pathway)]

# Prep for radar plot
super <- donor[c('s01044A','s01090A','s01092A','s01093A'),]
super <- as.data.frame(apply(super, MARGIN=2, median))
normal <- donor[c('n05042G','n05053B','n28045A','n28047D','n05098A','n05080M','n28043C','n28045D'),]
normal <- as.data.frame(apply(normal, MARGIN=2, median))
donor <- merge(super, normal, by='row.names')
rownames(donor) <- donor$Row.names
donor$Row.names <- NULL
colnames(donor) <- c('super', 'normal')

# Calculate max of each group
library(matrixStats)
donor_max <- round(donor)
donor_max$max <- as.vector(rowMaxs(as.matrix(donor_max)))
donor_max$super <- NULL
donor_max$normal <- NULL

# Calculate percentages
donor$super_perc <- (donor$super / (donor$super + donor$normal)) * 100.0
donor$normal_perc <- (donor$normal / (donor$super + donor$normal)) * 100.0
donor$super <- NULL
donor$normal <- NULL
colnames(donor) <- c('super', 'normal')

# Radar plot tables
donor <- as.data.frame(t(donor))
donor <- rbind(rep(floor(min(donor))-5, ncol(donor)) , rep(ceiling(max(donor))+5, ncol(donor)) , donor)
donor <- rbind(rep(0, ncol(donor)) , rep(100, ncol(donor)) , donor)
#write.table(donor, file='~/Desktop/repos/Jenior_Consortia_2022/data/supp_donor_radar_abund.tsv', append = FALSE, quote=FALSE, sep='\t')

# Read in previous results for better figure reproduceability
donor <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/supp_donor_radar_abund.tsv', sep='\t', header=TRUE, row.names=1)

# Rank by degree of difference
donor <- as.data.frame(t(donor))
donor$diff <- donor$super - donor$normal
donor <- donor[order(-donor$diff),] 
donor$diff <- NULL
donor <- as.data.frame(t(donor))
colnames(donor) <- gsub('_', ' ', colnames(donor))
colnames(donor) <- gsub('\\.', '\n', colnames(donor))

# Generate figures
library(fmsb)
library(scales)
super_col <- '#8e7cc3'
normal_col <- '#e69138'
success_col <- '#6aa84f'
failure_col <- '#990000'

# Generate figure
pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_S1C.pdf', width=5, height=5)
par(mar=c(1,1,1,1), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(donor, pcol=c(normal_col, super_col), pfcol=c(alpha(normal_col,0.7),alpha(super_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey60', cglty=5, cglwd=1)
text(x=-0.2, y=0.99,'100%', cex=0.5, font=1, srt=10)
text(x=-0.04, y=0.22, '0%', cex=0.4, font=1, srt=17)
legend('center', legend=c('Super Donor','Normal Donors'), col=c(super_col,normal_col),
       pt.bg=c(alpha(super_col,0.7),alpha(normal_col,0.7)), pch=22, pt.cex=1.5, cex=0.7, bg='white')
text(x=0, y=-0.2, 'OOB-error = 16.67%', cex=0.7, font=2)
dev.off()

#-----------#

# Supervised machine learning
condition <- as.factor(c('success','success','success','success','success','success','success','success','failure','failure','failure','failure'))
rf_donor <- randomForest(condition ~ ., data=outcome, ntree=5000, mtry=15,
                         importance=TRUE, replace=FALSE, err.rate=TRUE)
# OOB estimate of  error rate: 58.33%
# Confusion matrix:
#   failure success class.error
# failure       0       4       1.000
# success       3       5       0.375
rf_donor <- importance(rf_donor, type=1, scale=FALSE)
rf_donor <- subset(rf_donor, rf_donor > abs(min(rf_donor)))
rf_donor <- as.data.frame(rf_donor)
rf_donor$feature <- rownames(rf_donor)
colnames(rf_donor) <- c('mda', 'pathway')
rf_donor <- rf_donor[order(-rf_donor$mda),]
if (nrow(rf_donor > 20)) {rf_donor <- rf_donor[1:20,]}
outcome <- outcome[,which(colnames(outcome) %in% rf_donor$pathway)]

# Prep for radar plot
success <- outcome[c('s01044A','s01090A','s01092A','s01093A','n05098A','n05080M','n28043C','n28045D'),]
success <- as.data.frame(apply(success, MARGIN=2, median))
failure <- outcome[c('n05042G','n05053B','n28045A','n28047D'),]
failure <- as.data.frame(apply(failure, MARGIN=2, median))
outcome <- merge(success, failure, by='row.names')
rownames(outcome) <- outcome$Row.names
outcome$Row.names <- NULL
colnames(outcome) <- c('success', 'failure')

# Calculate max of each group
library(matrixStats)
outcome_max <- round(outcome)
outcome_max$max <- as.vector(rowMaxs(as.matrix(outcome_max)))

# Calculate percentages
outcome$success_perc <- (outcome$success / (outcome$success + outcome$failure)) * 100.0
outcome$failure_perc <- (outcome$failure / (outcome$success + outcome$failure)) * 100.0
outcome$success <- NULL
outcome$failure <- NULL
colnames(outcome) <- c('success', 'failure')

# Radar plot tables
outcome <- as.data.frame(t(outcome))
outcome <- rbind(rep(floor(min(outcome))-5, ncol(outcome)) , rep(ceiling(max(outcome))+5, ncol(outcome)) , outcome)
#write.table(outcome, file='~/Desktop/repos/Jenior_Consortia_2022/data/supp_outcome_radar_abund.tsv', append = FALSE, quote=FALSE, sep='\t')

# Read in previous results for better figure reproduceability
outcome <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/supp_outcome_radar_abund.tsv', sep='\t', header=TRUE, row.names=1)

# Rank by degree of difference
outcome <- as.data.frame(t(outcome))
outcome$diff <- outcome$success - outcome$failure
outcome <- outcome[order(-outcome$diff),] 
outcome$diff <- NULL
outcome <- as.data.frame(t(outcome))
colnames(outcome) <- gsub('_', ' ', colnames(outcome))
colnames(outcome) <- gsub('\\.', '\n', colnames(outcome))

# Generate figure
pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_S1D.pdf', width=4, height=4)
par(mar=c(1,1,1,1), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(outcome, pcol=c(failure_col, success_col), pfcol=c(alpha(failure_col,0.7),alpha(success_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.5, cglcol='grey60', cglty=5, cglwd=1)
text(x=-0.45, y=0.8,'87%', cex=0.5, font=1, srt=30)
text(x=-0.1, y=0.2, '13%', cex=0.5, font=1, srt=30)
legend('center', legend=c('Success','Failure'), col=c(success_col,failure_col),
       pt.bg=c(alpha(success_col,0.7),alpha(failure_col,0.7)), pch=22, pt.cex=1.4, cex=0.6, bg='white')
text(x=0, y=-0.2, 'OOB-error = 50.0%', cex=0.5, font=2)
dev.off()
