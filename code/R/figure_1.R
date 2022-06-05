
# Read in data
family <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metaphlan/family/family.tsv', sep='\t', header=TRUE, row.names=1)

# Calculate 'Other' category
family[family < 1] <- 0
family <- family[, which(colSums(family) > 0)]
family$Other <- 100 - rowSums(family)

# Group by metadata
metadata <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metadata.tsv', sep='\t', header=TRUE)
metadata$status <- NULL
metadata$sample_id <- NULL
metadata$outcome <- NULL
family <- merge(metadata, family, by.x='name', by.y='row.names')
rownames(family) <- family$name
family$name <- NULL

# Subset data
family_donor1 <- subset(family, donor == 1)
family_donor1$donor <- NULL
family_donor5 <- subset(family, donor == 5)
family_donor5$donor <- NULL
family_donor28 <- subset(family, donor == 28)
family_donor28$donor <- NULL
family <- rbind(family_donor1, rep(0, ncol(family_donor1)), family_donor5, rep(0, ncol(family_donor1)), family_donor28)
family <- as.data.frame(t(family))
rm(family_donor1, family_donor5, family_donor28)

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

pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_1B.pdf', width=5, height=3)
auto_barplot(family, new_palette, lgnd=FALSE)
dev.off()

pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_1B_legend.pdf', width=2, height=3)
par(mar=c(0,0,0,0))
plot(0, type='n', ylim=c(0,50), xlim=c(0,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend(x=2, y=49, legend="Other (<1% each)", pt.bg='gray30', pch=22, pt.cex=1.5, cex=0.6, bty='n')
legend(x=2, y=45.5, legend=rev(c("Bacteroidaceae","Porphyromonadaceae","Rikenellaceae","Prevotellaceae")), pt.bg=reds, pch=22, pt.cex=1.5, cex=0.6, bty='n')
legend(x=2, y=35, legend=rev(c("Bifidobacteriaceae","Coriobacteriaceae")), pt.bg=greens, pch=22, pt.cex=1.5, cex=0.6, bty='n')
legend(x=2, y=29, legend=rev(c("Ruminococcaceae","Lachnospiraceae","Eubacteriaceae","Clostridiaceae","Veillonellaceae","Erysipelotrichaceae","Oscillospiraceae","Streptococcaceae")), 
       pt.bg=blues, pch=22, pt.cex=1.5, cex=0.6, bty='n')
legend(x=2, y=10, legend="Verrucomicrobiaceae", pt.bg='goldenrod1', pch=22, pt.cex=1.5, cex=0.6, bty='n')
segments(x0=1.95, y0=c(44.5,34,28,8.75), y1=c(35.5,29.5,10.5,6.75), lwd=2)
text(x=1.82, y=c(40.25,32,19.5,7.75), labels=c('Bacteroidota\n(Bacteroidetes)','Actinomycetota\n(Actinobacteria)','Bacillota\n(Firmicutes)','Verrucomicrobiota'), 
     adj=1, font=2, cex=0.5)
dev.off()

#-------------------------------------------------------------------------------------------------------------------------------#

# FUNCTIONAL ANALYSIS

# Metabolic pathway abundance

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
# Super vs Normal
super <- merge(s01044A, s01090A, by='pathway')
super <- merge(super, s01092A, by='pathway')
super <- merge(super, s01093A, by='pathway')
rownames(super) <- super$pathway
super$pathway <- NULL
colnames(super) <- c('s01044A','s01090A','s01092A','s01093A')
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
# Success vs Failure
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
failure <- merge(n05042G, n05053B, by='pathway')
failure <- merge(failure, n28045A, by='pathway')
failure <- merge(failure, n28047D, by='pathway')
rownames(failure) <- failure$pathway
failure$pathway <- NULL
colnames(failure) <- c('n05042G','n05053B','n28045A','n28047D')
rm(s01044A,s01090A,s01092A,s01093A,n05098A,n05080M,n28043C,n28045D,n05042G,n05053B,n28045A,n28047D)

# Focus on metabolic pathways
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
super <- super[metabolic_paths,]
normal <- normal[metabolic_paths,]
success <- success[metabolic_paths,]
failure <- failure[metabolic_paths,]

# Subsample
library(vegan)
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

# Supervised machine learning
library(randomForest)
set.seed(906801)
condition <- as.factor(c('super','super','super','super','normal','normal','normal','normal','normal','normal','normal','normal'))
rf_donor <- randomForest(condition ~ ., data=donor, ntree=5000, mtry=15,
                         importance=TRUE, replace=FALSE, err.rate=TRUE)
print(rf_donor)
#OOB estimate of  error rate: 8.33%
#Confusion matrix:
#  normal super class.error
#normal      8     0        0.00
#super       1     3        0.25
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
#write.table(donor, file='~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar_abund.tsv', append = FALSE, quote=FALSE, sep='\t')

# Read in previous results for better figure reproduceability
donor <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar_abund.tsv', sep='\t', header=TRUE, row.names=1)

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
pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_1C.pdf', width=5, height=5)
par(mar=c(1,1,1,1), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(donor, pcol=c(normal_col, super_col), pfcol=c(alpha(normal_col,0.7),alpha(super_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey60', cglty=5, cglwd=1)
text(x=-0.2, y=0.99,'80%', cex=0.6, font=1, srt=10)
text(x=-0.04, y=0.22, '20%', cex=0.4, font=1, srt=17)
legend('center', legend=c('Super Donor','Normal Donors'), col=c(super_col,normal_col),
       pt.bg=c(alpha(super_col,0.7),alpha(normal_col,0.7)), pch=22, pt.cex=1.5, cex=0.7, bg='white')
text(x=0, y=-0.2, 'OOB-error = 8.33%', cex=0.7, font=2)
dev.off()

#-----------#

# Supervised machine learning
condition <- as.factor(c('success','success','success','success','success','success','success','success','failure','failure','failure','failure'))
rf_donor <- randomForest(condition ~ ., data=outcome, ntree=5000, mtry=15,
                         importance=TRUE, replace=FALSE, err.rate=TRUE)
print(rf_donor)
# OOB estimate of  error rate: 41.67%
# Confusion matrix:
#   failure success class.error
# failure       1       3        0.75
# success       2       6        0.25
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
#write.table(outcome, file='~/Desktop/repos/Jenior_Consortia_2022/data/outcome_radar_abund.tsv', append = FALSE, quote=FALSE, sep='\t')

# Read in previous results for better figure reproduceability
outcome <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/outcome_radar_abund.tsv', sep='\t', header=TRUE, row.names=1)

# Rank by degree of difference
outcome <- as.data.frame(t(outcome))
outcome$diff <- outcome$success - outcome$failure
outcome <- outcome[order(-outcome$diff),] 
outcome$diff <- NULL
outcome <- as.data.frame(t(outcome))
colnames(outcome) <- gsub('_', ' ', colnames(outcome))
colnames(outcome) <- gsub('\\.', '\n', colnames(outcome))

# Generate figure
pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_1D.pdf', width=4, height=4)
par(mar=c(1,1,1,1), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(outcome, pcol=c(failure_col, success_col), pfcol=c(alpha(failure_col,0.7),alpha(success_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey60', cglty=5, cglwd=1)
text(x=-0.45, y=0.8,'80%', cex=0.5, font=1, srt=30)
text(x=-0.1, y=0.2, '20%', cex=0.5, font=1, srt=30)
legend('center', legend=c('Success','Failure'), col=c(success_col,failure_col),
       pt.bg=c(alpha(success_col,0.7),alpha(failure_col,0.7)), pch=22, pt.cex=1.4, cex=0.6, bg='white')
text(x=0, y=-0.2, 'OOB-error = 41.67%', cex=0.5, font=2)
dev.off()
