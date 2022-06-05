
# Read in data
metag <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/compiled_mapping.pathways.tsv', sep='\t', header=TRUE)

# Subset to known pathways
metag <- as.data.frame(subset(metag, metag$pathway != 'unknown'))
metag$keggID <- NULL
metag$gene <- NULL

# Subset by number of genes in each pathway
keep <- as.data.frame(table(metag$pathway))
keep <- as.character(subset(keep, Freq > 10)$Var1)
metag <- subset(metag, pathway %in% keep)
rm(keep)

# Agregate and reformat
metag <- aggregate(. ~ pathway, metag, sum)
rownames(metag) <- make.names(metag$pathway)
metag$pathway <- NULL
metag <- droplevels(metag)

# Subsample groups evenly
library(vegan)
metag_super <- metag[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]
metag_normal <- metag[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]
sub_level <- ceiling(min(c(as.vector(colSums(metag_super)), as.vector(colSums(metag_normal)))) * 0.9)
for (x in 1:ncol(metag_super)) {metag_super[,x] <- as.vector(rrarefy(metag_super[,x], sample=sub_level))}
for (y in 1:ncol(metag_normal)) {metag_normal[,y] <- as.vector(rrarefy(metag_normal[,y], sample=sub_level))}
metag_super <- as.data.frame(t(metag_super))
metag_normal <- as.data.frame(t(metag_normal))

# Subset by colsums
test <- colSums(metag_super) + colSums(metag_normal)
metag_super <- metag_super[,as.vector(which(test > 10))]
metag_normal <- metag_normal[,as.vector(which(test > 10))]

# Subset by Wilcoxon test significance
pvals <- c()
for (x in 1:ncol(metag_super)) {pvals[x] <- wilcox.test(metag_super[,x], metag_normal[,x], exact=FALSE)$p.value}
metag_super_path <- metag_super[,which(pvals <= 0.05)]
metag_normal_path <- metag_normal[,which(pvals <= 0.05)]
metag_super_unchanged <- metag_super[,which(pvals > 0.05)]
metag_normal_unchanged <- metag_normal[,which(pvals > 0.05)]

# Calculate median of each group
metag_super_med <- apply(metag_super_path, 2, median)
metag_normal_med <- apply(metag_normal_path, 2, median)
metag_med <- cbind(as.numeric(as.character(metag_super_med)),as.numeric(as.character(metag_normal_med)))
rownames(metag_med) <- colnames(metag_super_path)
colnames(metag_med) <- c('super','normal')
metag_med <- as.data.frame(metag_med)
metag_med$abs_diff <- abs(metag_med$super - metag_med$normal)
metag_med <- metag_med[order(-metag_med$abs_diff),]
metag_med$abs_diff <- NULL

# Calculate percentages
metag_med$super_percent <- (metag_med$super / (metag_med$super + metag_med$normal)) * 100.0
metag_med$normal_percent <- (metag_med$normal / (metag_med$super + metag_med$normal)) * 100.0
metag_med$super <- NULL
metag_med$normal <- NULL

# Prep for radar plot
metag_med <- as.data.frame(t(metag_med))
metag_med <- rbind(rep(100, ncol(metag_med)) , rep(0, ncol(metag_med)) , metag_med)
metag_legend <- as.data.frame(cbind(rev(paste0(seq(1,ncol(metag_med),1),')')), colnames(metag_med)))
colnames(metag_legend) <- c('number', 'pathway')
colnames(metag_med) <- rev(seq(1,ncol(metag_med),1))
metag_legend$pathway  <- gsub('_', ' ', metag_legend$pathway)
metag_legend$pathway <- gsub(' .CAMP.', '', metag_legend$pathway)
metag_legend$pathway <- gsub('C5.', '', metag_legend$pathway)
metag_legend$pathway <- gsub('X2.', '', metag_legend$pathway)
metag_legend$pathway <- gsub(' . Vibrio cholerae', '', metag_legend$pathway)
metag_legend$pathway <- gsub('Cationic a', 'A', metag_legend$pathway)
metag_legend$pathway <- gsub('Biosynthesis of siderophore group nonribosomal peptides', 'Siderophore peptide biosynthesis', metag_legend$pathway)
metag_legend$pathway <- gsub('Phosphonate and phosphinate', 'Phosphoenolpyruvate', metag_legend$pathway)
metag_legend$pathway <- gsub('Limonene and pinene', 'Terpene', metag_legend$pathway)
metag_legend$pathway <- gsub('Biosynthesis of amino acids', 'Amino acid biosynthesis', metag_legend$pathway)
metag_legend$pathway <- gsub('Biosynthesis of ansamycins', 'Ansamycin biosynthesis', metag_legend$pathway)
metag_legend$pathway <- gsub('Penicillin and cephalosporin biosynthesis', 'Antiobiotic biosynthesis', metag_legend$pathway)

# Radar plot
# https://www.r-graph-gallery.com/142-basic-radar-chart.html
require('fmsb')
library(fmsb)
require('scales')
library(scales)
super_col <- 'deepskyblue3'
normal_col <- 'firebrick3'
path_labels <- c("Amino acid\nbiosynthesis\n(12889)","Carbon\nmetabolism\n(4570)","Antimicrobial peptide\nresistance\n(701)",
                 "Bacterial\nchemotaxis\n(729)","Sulfur\nrelay system\n(842)","Nucleotide\nexcision repair\n(707)",
                 "Siderophore peptide\nbiosynthesis\n(30)","Pantothenate & CoA\nbiosynthesis\n(199)",
                 "Sphingolipid\nmetabolism\n(50)","Propanoate\nmetabolism\n(17)","Oxocarboxylic acid\nmetabolism\n(53)")
colnames(metag_med) <- path_labels

pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/pathway_radar.pdf', width=7.5, height=5)
par(mar=c(0,1,0,1), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(metag_med, pcol=c(super_col,normal_col), pfcol=c(alpha(super_col,0.7),alpha(normal_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey30', cglty=5, cglwd=1)
#text(x=c(-0.91,0.91,0.09,0), y=c(-0.04,-0.04,0.93,-0.93), labels='100%', cex=0.6, font=1)
legend('bottomright', legend=c('Super Donor', 'Normal Donors'), col=c(super_col,normal_col),
       pt.bg=c(alpha(super_col,0.7),alpha(normal_col,0.7)), pch=22, pt.cex=1.5, cex=0.8, bty='n')
dev.off()

# Abundance bar plots
pathway_bars <- function(super_abund, normal_abund, path_name) {
  pval <- round(wilcox.test(super_abund, normal_abund, exact=FALSE)$p.value, 4)
  super_abund <- as.vector(quantile(super_abund, c(0.25, 0.5, 0.75)))
  normal_abund <- as.vector(quantile(normal_abund, c(0.25, 0.5, 0.75)))
  xmax <- max(c(super_abund[3], normal_abund[3])) * 1.3
  par(mar=c(2.6,4,2,1), xpd=FALSE, las=1, mgp=c(1.6,0.5,0), lwd=2.2)
  barplot(c(normal_abund[2], super_abund[2]), xlim=c(0,xmax), xaxs='i', xlab='Gene Abundance', cex.names=1.2,
          horiz=TRUE, col=c(normal_col, super_col), names.arg=c('Normal\nDonor','Super\nDonor'), xaxt='n',
          main=path_name, cex.main=1.1)
  axis(side=1, cex.axis=0.9, lwd=2)
  box()
  segments(x0=normal_abund[1], x1=normal_abund[3], y0=0.7)
  segments(x0=normal_abund[1], y0=0.55, y1=0.85)
  segments(x0=normal_abund[3], y0=0.55, y1=0.85)
  segments(x0=super_abund[1], x1=super_abund[3], y0=1.9)
  segments(x0=super_abund[1], y0=1.75, y1=2.05)
  segments(x0=super_abund[3], y0=1.75, y1=2.05)
  line_pos <- max(c(super_abund[3], normal_abund[3])) * 1.1
  segments(x0=line_pos, y0=0.7, y1=1.9, lwd=2.5)
  sig_pos <- max(c(super_abund[3], normal_abund[3])) * 1.2
  if (pval <= 0.001) {text(x=sig_pos, y=1.3, '***', font=2, cex=2, srt=270)}
  if (pval <= 0.01) {text(x=sig_pos, y=1.3, '**', font=2, cex=2, srt=270)}
  if (pval <= 0.05) {text(x=sig_pos, y=1.3, '*', font=2, cex=2, srt=270)}
  if (pval > 0.05) {text(x=sig_pos, y=1.3, 'n.s.', cex=1.2, srt=270)} 
  }

pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/pathway_abundances.pdf', width=9, height=6)
layout(matrix(c(1,2,3,
                4,5,6,
                7,8,9), nrow=3, ncol=3, byrow=TRUE))
pathway_bars(metag_super_path$Bacterial_chemotaxis, metag_normal_path$Bacterial_chemotaxis, 'Bacterial chemotaxis')
#pathway_bars(metag_super_path$Biosynthesis_of_amino_acids, metag_normal_path$Biosynthesis_of_amino_acids, 'Amino acid biosynthesis')
pathway_bars(metag_super_path$Biosynthesis_of_siderophore_group_nonribosomal_peptides, metag_normal_path$Biosynthesis_of_siderophore_group_nonribosomal_peptides, 'Siderophore peptide biosynthesis')
pathway_bars(metag_super_path$C5.Branched_dibasic_acid_metabolism, metag_normal_path$C5.Branched_dibasic_acid_metabolism, 'Branched dibasic acid metabolism')
#pathway_bars(metag_super_path$Carbon_metabolism, metag_normal_path$Carbon_metabolism, 'Carbon metabolism')
pathway_bars(metag_super_path$Cationic_antimicrobial_peptide_.CAMP._resistance, metag_normal_path$Cationic_antimicrobial_peptide_.CAMP._resistance, 'Cationic antimicrobial peptide resistance')
pathway_bars(metag_super_path$Ether_lipid_metabolism, metag_normal_path$Ether_lipid_metabolism, 'Ether lipid metabolism')
pathway_bars(metag_super_path$Limonene_and_pinene_degradation, metag_normal_path$Limonene_and_pinene_degradation, 'Terpene degradation')
#pathway_bars(metag_super_path$Nucleotide_excision_repair, metag_normal_path$Nucleotide_excision_repair, 'Nucleotide excision repair')
#pathway_bars(metag_super_path$Phosphonate_and_phosphinate_metabolism, metag_normal_path$Phosphonate_and_phosphinate_metabolism, 'Phosphoenolpyruvate metabolism')
pathway_bars(metag_super_path$Pantothenate_and_CoA_biosynthesis, metag_normal_path$Pantothenate_and_CoA_biosynthesis, 'Pantothenate and CoA biosynthesis')
pathway_bars(metag_super_path$Propanoate_metabolism, metag_normal_path$Propanoate_metabolism, 'Propanoate metabolism')
pathway_bars(metag_super_path$Sphingolipid_metabolism, metag_normal_path$Sphingolipid_metabolism, 'Sphingolipid metabolism')
dev.off()

#-------------------------------------------------------------------------------------------------------------#

# Break down important pathways

# Read in data
metag <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/compiled_mapping.pathways.tsv', sep='\t', header=TRUE)

subset_metag <- function(metaG, path) {
  metaG <- as.data.frame(subset(metaG, metaG$pathway != 'unknown'))
  sub_metag <- as.data.frame(subset(metaG, metaG$pathway == path))
  sub_metag$gene <- NULL
  sub_metag$pathway <- NULL
  #sub_metag <- aggregate(. ~ keggID, sub_metag, sum)
  rownames(sub_metag) <- sub_metag$keggID
  sub_metag$keggID <- NULL
  return(sub_metag)}

# Subset to each significant pathway group
metag_aa <- subset_metag(metag, 'Biosynthesis_of_amino_acids')
metag_sphingolipid <- subset_metag(metag, 'Sphingolipid_metabolism')
metag_pantothenate <- subset_metag(metag, 'Pantothenate_and_CoA_biosynthesis')
metag_limonene <- subset_metag(metag, 'Limonene_and_pinene_degradation')
metag_ether <- subset_metag(metag, 'Ether_lipid_metabolism')
metag_siderophore <- subset_metag(metag, 'Biosynthesis_of_siderophore_group_nonribosomal_peptides')
metag_propanoate <- subset_metag(metag, 'Propanoate_metabolism')
rm(metag)

# Break into the 2 groups
metag_aa_super <- as.data.frame(t(metag_aa[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]))
metag_aa_normal <- as.data.frame(t(metag_aa[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]))
rm(metag_aa)
metag_sphingolipid_super <- as.data.frame(t(metag_sphingolipid[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]))
metag_sphingolipid_normal <- as.data.frame(t(metag_sphingolipid[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]))
rm(metag_sphingolipid)
metag_pantothenate_super <- as.data.frame(t(metag_pantothenate[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]))
metag_pantothenate_normal <- as.data.frame(t(metag_pantothenate[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]))
rm(metag_pantothenate)
metag_limonene_super <- as.data.frame(t(metag_limonene[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]))
metag_limonene_normal <- as.data.frame(t(metag_limonene[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]))
rm(metag_limonene)
metag_ether_super <- as.data.frame(t(metag_ether[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]))
metag_ether_normal <- as.data.frame(t(metag_ether[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]))
rm(metag_ether)
metag_siderophore_super <- as.data.frame(t(metag_siderophore[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]))
metag_siderophore_normal <- as.data.frame(t(metag_siderophore[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]))
rm(metag_siderophore)
metag_propanoate_super <- as.data.frame(t(metag_propanoate[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]))
metag_propanoate_normal <- as.data.frame(t(metag_propanoate[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]))
rm(metag_propanoate)





#-------------------------------------------------------------------------------------------------------------#

# Read in data
metag <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/compiled_mapping.pathways.tsv', sep='\t', header=TRUE)

# Subset to known pathways
metag$gene <- make.names(metag$gene)
metag$keggID <- NULL
metabolic <- subset(metag, pathway == 'Metabolic_pathways')
metabolic$pathway <- NULL
carbon <- subset(metag, pathway == 'Carbon_metabolism')
carbon$pathway <- NULL
amino_acid <- subset(metag, pathway == 'Biosynthesis_of_amino_acids')
amino_acid$pathway <- NULL
other <- subset(metag, !pathway %in% c('Metabolic_pathways','Carbon_metabolism'))
other$pathway <- NULL
rm(metag)

# Subsample, aggregate, and reformat
library(vegan)
sub_sample <- function(metagenome) { 
  cols <- c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv',
            'fmt.metaG.01093A.reads2genes.norm.tsv','fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv',
            'fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv',
            'fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')
  sub_lvl <- round(min(colSums(metagenome[,cols])) * 0.9)
  for (x in cols) {metagenome[,x] <- as.vector(rrarefy(metagenome[,x], sample=sub_lvl))}
  metagenome <- aggregate(. ~ gene, metagenome, sum)
  rownames(metagenome) <- metagenome$gene
  metagenome$gene <- NULL
  metagenome <- droplevels(metagenome)
  metagenome <- as.data.frame(t(metagenome))
  return(metagenome)}

metabolic <- sub_sample(metabolic)
carbon <- sub_sample(carbon)
amino_acid <- sub_sample(amino_acid)
#other <- sub_sample(other)
#metag <- sub_sample(metag)

# Screen to avoid stack overflow
#keep <- c()
#y <- 1
#for (x in 1:ncol(metag)) {
#  if (var(metag[,x]) > 100) {
#    keep[y] <- colnames(metag)[x]
#    y <- y+1}}
#metag <- metag[,keep]
#rm(keep, x, y)

# Machine learning
library(randomForest)
condition <- as.factor(c('super','super','super','super',
                         'normal','normal','normal','normal',
                         'normal','normal','normal','normal'))
rf_obj <- randomForest(condition ~ ., data=amino_acid, importance=TRUE, err.rate=TRUE, ntree=5000, mtry=50)
# OOB estimate of  error rate: 0%
#   normal super class.error
# normal      8     0           0
# super       0     4           0
rf_mda <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- as.data.frame(rf_mda)
rf_mda$feature <- rownames(rf_mda)
rf_mda <- rf_mda[order(-rf_mda$MeanDecreaseAccuracy),] 
rf_mda <- subset(rf_mda, rf_mda$MeanDecreaseAccuracy > (abs(min(rf_mda$MeanDecreaseAccuracy)))) # significance
#write.table(rf_mda, file='~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/pathway_RF.tsv', quote=FALSE, sep='\t')
rm(rf_obj, condition)


#rf_mda <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/pathway_RF.tsv', sep='\t', header=TRUE)
#rf_mda <- subset(rf_mda, MeanDecreaseAccuracy > 6.0)
metag <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/compiled_mapping.pathways.tsv', sep='\t', header=TRUE)
#metag <- as.data.frame(subset(metag, metag$pathway != 'unknown'))
#metag <- as.data.frame(subset(metag, metag$gene != 'hypothetical_protein'))
metag$gene <- make.names(metag$gene)
gene_metadata <- metag[,c('gene','pathway')]
gene_metadata <- unique(gene_metadata)
metag$keggID <- NULL
metag$pathway <- NULL
metag <- aggregate(. ~ gene, metag, sum)
rownames(metag) <- metag$gene
metag$gene <- NULL
metag <- droplevels(metag)
metag <- as.data.frame(t(metag))
metag <- metag[, as.character(rf_mda$feature)]
rm(rf_mda)

# Subsample for plotting comparison
amino_acid <- sub_sample(as.data.frame(t(amino_acid)))


# Subset by donor type
metag$donor <- c('super','super','super','super','normal','normal','normal','normal','normal','normal','normal','normal')
metag_super <- subset(metag, donor == 'super')
metag_super$donor <- NULL
metag_normal <- subset(metag, donor == 'normal')
metag_normal$donor <- NULL
rm(metag)

# Test for differences at species-level
pvals <- c()
for (x in c(1:ncol(metag_super))) {
  pvals[x] <- round(wilcox.test(metag_super[,x], metag_normal[,x], exact=FALSE)$p.value, 3)}
sig_genes <- c()
for (y in 1:length(pvals)) {if (pvals[y]<=0.01) {sig_genes <- c(sig_genes, colnames(metag_super)[y])}}
metag_super <- metag_super[, sig_genes]
metag_normal <- metag_normal[, sig_genes]
gene_metadata <- subset(gene_metadata, gene %in% sig_genes)
rm(x, y, pvals, sig_genes)

# Abundance bar plots
auto_bars <- function(super_abund, normal_abund, path_name) {
  pval <- round(wilcox.test(super_abund, normal_abund, exact=FALSE)$p.value, 4)
  super_abund <- as.vector(quantile(super_abund, c(0.25, 0.5, 0.75)))
  normal_abund <- as.vector(quantile(normal_abund, c(0.25, 0.5, 0.75)))
  xmax <- max(c(super_abund[3], normal_abund[3])) * 1.3
  par(mar=c(2.6,4,2,1), xpd=FALSE, las=1, mgp=c(1.6,0.5,0), lwd=2.2)
  barplot(c(normal_abund[2], super_abund[2]), xlim=c(0,xmax), xaxs='i', xlab='Gene Abundance', cex.names=1.2,
          horiz=TRUE, col=c('green3', 'blue3'), names.arg=c('Normal\nDonor','Super\nDonor'), xaxt='n',
          main=path_name, cex.main=0.8)
  axis(side=1, cex.axis=0.9, lwd=2)
  box()
  segments(x0=normal_abund[1], x1=normal_abund[3], y0=0.7)
  segments(x0=normal_abund[1], y0=0.55, y1=0.85)
  segments(x0=normal_abund[3], y0=0.55, y1=0.85)
  segments(x0=super_abund[1], x1=super_abund[3], y0=1.9)
  segments(x0=super_abund[1], y0=1.75, y1=2.05)
  segments(x0=super_abund[3], y0=1.75, y1=2.05)
  line_pos <- max(c(super_abund[3], normal_abund[3])) * 1.1
  segments(x0=line_pos, y0=0.7, y1=1.9, lwd=2.5)
  sig_pos <- max(c(super_abund[3], normal_abund[3])) * 1.2
  if (pval <= 0.001) {text(x=sig_pos, y=1.3, '***', cex=2, srt=270)}
  if (pval <= 0.01) {text(x=sig_pos, y=1.3, '**', cex=2, srt=270)}
  if (pval <= 0.05) {text(x=sig_pos, y=1.3, '*', cex=2, srt=270)}
  if (pval > 0.05) {text(x=sig_pos, y=1.3, 'n.s.', cex=1.2, srt=270)} 
}

for (x in 1:ncol(metag_super)) {
  file_name <- paste0('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_analysis/amino_acid/', colnames(metag_super)[x], '.png')
  png(filename=file_name, units='in', width=4, height=3, res=300)
  auto_bars(metag_super[,x], metag_normal[,x], colnames(metag_super)[x])
  dev.off()}
rm(x)


#------------------------------------------------------------------------------------------------------#

# Completeness analysis
s01044A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/01044A/completeness_report.tsv', sep='\t', header=TRUE)
s01090A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/01090A/completeness_report.tsv', sep='\t', header=TRUE)
s01092A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/01092A/completeness_report.tsv', sep='\t', header=TRUE)
s01093A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/01093A/completeness_report.tsv', sep='\t', header=TRUE)
n05042G <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/05042G/completeness_report.tsv', sep='\t', header=TRUE)
n05053B <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/05053B/completeness_report.tsv', sep='\t', header=TRUE)
n05080M <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/05080M/completeness_report.tsv', sep='\t', header=TRUE)
n05098A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/05098A/completeness_report.tsv', sep='\t', header=TRUE)
n28043C <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/28043C/completeness_report.tsv', sep='\t', header=TRUE)
n28045A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/28045A/completeness_report.tsv', sep='\t', header=TRUE)
n28045D <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/28045D/completeness_report.tsv', sep='\t', header=TRUE)
n28047D <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/28047D/completeness_report.tsv', sep='\t', header=TRUE)

super <- merge(s01044A, s01090A, by='pathway')
super <- merge(super, s01092A, by='pathway')
super <- merge(super, s01093A, by='pathway')
normal <- merge(n05042G, n05053B, by='pathway')
normal <- merge(normal, n05080M, by='pathway')
normal <- merge(normal, n05098A, by='pathway')
normal <- merge(normal, n28043C, by='pathway')
normal <- merge(normal, n28045A, by='pathway')
normal <- merge(normal, n28045D, by='pathway')
normal <- merge(normal, n28047D, by='pathway')
success <- merge(s01044A, s01090A, by='pathway')
success <- merge(success, s01092A, by='pathway')
success <- merge(success, s01093A, by='pathway')
success <- merge(success, n05098A, by='pathway')
success <- merge(success, n05080M, by='pathway')
success <- merge(success, n28043C, by='pathway')
success <- merge(success, n28045D, by='pathway')
failure <- merge(n05042G, n05053B, by='pathway')
failure <- merge(failure, n28045A, by='pathway')
failure <- merge(failure, n28047D, by='pathway')
rm(s01044A, s01090A, s01092A, s01093A)
rm(n05042G, n05053B, n05080M, n05098A)
rm(n28043C, n28045A, n28045D, n28047D)

library(Matrix)
completeness <- merge(super, normal, by='pathway')
rownames(completeness) <- completeness$pathway
completeness$pathway <- NULL
keep_paths <- c("2-Oxocarboxylic_acid_metabolism","ABC_transporters","Acarbose_and_validamycin_biosynthesis","Acridone_alkaloid_biosynthesis","Aflatoxin_biosynthesis",
                "Alanine,_aspartate_and_glutamate_metabolism","Aldosterone_synthesis_and_secretion","Aldosterone-regulated_sodium_reabsorption","alpha-Linolenic_acid_metabolism",
                "Amino_sugar_and_nucleotide_sugar_metabolism","Aminoacyl-tRNA_biosynthesis","Aminobenzoate_degradation","AMPK_signaling_pathway","Amyotrophic_lateral_sclerosis",
                "Anthocyanin_biosynthesis","Antifolate_resistance","Antigen_processing_and_presentation","Apelin_signaling_pathway","Apoptosis_-_multiple_species",
                "Arabinogalactan_biosynthesis_-_Mycobacterium","Arachidonic_acid_metabolism","Arginine_and_proline_metabolism","Arginine_biosynthesis",
                "Ascorbate_and_aldarate_metabolism","Atrazine_degradation","Bacterial_chemotaxis","Bacterial_invasion_of_epithelial_cells","Bacterial_secretion_system",
                "Basal_cell_carcinoma","Basal_transcription_factors","Base_excision_repair","Benzoate_degradation","Benzoxazinoid_biosynthesis","beta-Alanine_metabolism",
                "beta-Lactam_resistance","Betalain_biosynthesis","Bile_secretion","Biofilm_formation_-_Escherichia_coli","Biofilm_formation_-_Pseudomonas_aeruginosa",
                "Biofilm_formation_-_Vibrio_cholerae","Biosynthesis_of_12-,_14-_and_16-membered_macrolides","Biosynthesis_of_amino_acids","Biosynthesis_of_ansamycins",
                "Biosynthesis_of_enediyne_antibiotics","Biosynthesis_of_secondary_metabolites","Biosynthesis_of_siderophore_group_nonribosomal_peptides",
                "Biosynthesis_of_type_II_polyketide_backbone","Biosynthesis_of_type_II_polyketide_products","Biosynthesis_of_unsaturated_fatty_acids",
                "Biosynthesis_of_vancomycin_group_antibiotics","Biosynthesis_of_various_secondary_metabolites_-_part_1","Biosynthesis_of_various_secondary_metabolites_-_part_2",
                "Biosynthesis_of_various_secondary_metabolites_-_part_3","Biotin_metabolism","Bisphenol_degradation","Brassinosteroid_biosynthesis","Butanoate_metabolism",
                "C-type_lectin_receptor_signaling_pathway","C5-Branched_dibasic_acid_metabolism","Calcium_signaling_pathway","cAMP_signaling_pathway","Caprolactam_degradation",
                "Carbapenem_biosynthesis","Carbohydrate_digestion_and_absorption","Carbon_fixation_pathways_in_prokaryotes","Carbon_metabolism","Cardiac_muscle_contraction",
                "Carotenoid_biosynthesis","Cationic_antimicrobial_peptide_(CAMP)_resistance","Cell_adhesion_molecules","Cell_cycle","Cellular_senescence",
                "Central_carbon_metabolism_in_cancer","cGMP-PKG_signaling_pathway","Chemokine_signaling_pathway","Chloroalkane_and_chloroalkene_degradation",
                "Chlorocyclohexane_and_chlorobenzene_degradation","Cholesterol_metabolism","Cholinergic_synapse","Circadian_entrainment","Circadian_rhythm",
                "Citrate_cycle_(TCA_cycle)","Clavulanic_acid_biosynthesis","Complement_and_coagulation_cascades","Cortisol_synthesis_and_secretion",
                "Cutin,_suberine_and_wax_biosynthesis","Cyanoamino_acid_metabolism","Cysteine_and_methionine_metabolism","Cytokine-cytokine_receptor_interaction",
                "Cytosolic_DNA-sensing_pathway","D-Alanine_metabolism","D-Arginine_and_D-ornithine_metabolism","D-Glutamine_and_D-glutamate_metabolism",
                "Degradation_of_aromatic_compounds","Dilated_cardiomyopathy","Dioxin_degradation","Diterpenoid_biosynthesis","DNA_replication","Dopaminergic_synapse",
                "Dorso-ventral_axis_formation","ECM-receptor_interaction","EGFR_tyrosine_kinase_inhibitor_resistance","Endocrine_and_other_factor-regulated_calcium_reabsorption",
                "Endocytosis","Epithelial_cell_signaling_in_Helicobacter_pylori_infection","ErbB_signaling_pathway","Estrogen_signaling_pathway","Ether_lipid_metabolism",
                "Ethylbenzene_degradation","Fanconi_anemia_pathway","Fatty_acid_biosynthesis","Fatty_acid_degradation","Fatty_acid_elongation","Fatty_acid_metabolism",
                "Fc_epsilon_RI_signaling_pathway","Fc_gamma_R-mediated_phagocytosis","Flagellar_assembly","Flavone_and_flavonol_biosynthesis","Flavonoid_biosynthesis",
                "Fluorobenzoate_degradation","Folate_biosynthesis","FoxO_signaling_pathway","Fructose_and_mannose_metabolism","Furfural_degradation","GABAergic_synapse",
                "Galactose_metabolism","Geraniol_degradation","Glucagon_signaling_pathway","Glucosinolate_biosynthesis","Glutamatergic_synapse","Glutathione_metabolism",
                "Glycerolipid_metabolism","Glycerophospholipid_metabolism","Glycine,_serine_and_threonine_metabolism","Glycolysis_Gluconeogenesis",
                "Glycosaminoglycan_biosynthesis_-_chondroitin_sulfate_dermatan_sulfate","Glycosaminoglycan_biosynthesis_-_heparan_sulfate_heparin",
                "Glycosaminoglycan_biosynthesis_-_keratan_sulfate","Glycosaminoglycan_degradation","Glycosphingolipid_biosynthesis_-_ganglio_series",
                "Glycosphingolipid_biosynthesis_-_globo_and_isoglobo_series","Glycosphingolipid_biosynthesis_-_lacto_and_neolacto_series",
                "Glycosylphosphatidylinositol_(GPI)-anchor_biosynthesis","Glyoxylate_and_dicarboxylate_metabolism","GnRH_secretion","GnRH_signaling_pathway",
                "Growth_hormone_synthesis,_secretion_and_action","Hedgehog_signaling_pathway","HIF-1_signaling_pathway","Hippo_signaling_pathway","Histidine_metabolism",
                "Homologous_recombination","IL-17_signaling_pathway","Indole_alkaloid_biosynthesis","Indole_diterpene_alkaloid_biosynthesis",
                "Inflammatory_mediator_regulation_of_TRP_channels","Inositol_phosphate_metabolism","Insulin_secretion","Insulin_signaling_pathway",
                "Intestinal_immune_network_for_IgA_production","Isoflavonoid_biosynthesis","Isoquinoline_alkaloid_biosynthesis","JAK-STAT_signaling_pathway",
                "Legionellosis","Limonene_and_pinene_degradation","Linoleic_acid_metabolism","Lipoarabinomannan_(LAM)_biosynthesis","Lipoic_acid_metabolism",
                "Lipopolysaccharide_biosynthesis","Lysine_biosynthesis","Lysine_degradation","Lysosome","Mannose_type_O-glycan_biosynthesis","MAPK_signaling_pathway",
                "Melanogenesis","Metabolic_pathways","Methane_metabolism","Microbial_metabolism_in_diverse_environments","Mineral_absorption","Mismatch_repair",
                "Monobactam_biosynthesis","Monoterpenoid_biosynthesis","mRNA_surveillance_pathway","Mucin_type_O-glycan_biosynthesis","N-Glycan_biosynthesis",
                "Naphthalene_degradation","Neomycin,_kanamycin_and_gentamicin_biosynthesis","Neuroactive_ligand-receptor_interaction","Neurotrophin_signaling_pathway",
                "NF-kappa_B_signaling_pathway","Nicotinate_and_nicotinamide_metabolism","Nitrogen_metabolism","Nitrotoluene_degradation","NOD-like_receptor_signaling_pathway",
                "Non-homologous_end-joining","Nonribosomal_peptide_structures","Notch_signaling_pathway","Novobiocin_biosynthesis","Nucleotide_excision_repair",
                "O-Antigen_nucleotide_sugar_biosynthesis","One_carbon_pool_by_folate","Other_glycan_degradation","Other_types_of_O-glycan_biosynthesis",
                "Oxidative_phosphorylation","Pantothenate_and_CoA_biosynthesis","Parathyroid_hormone_synthesis,_secretion_and_action","Pathogenic_Escherichia_coli_infection",
                "PD-L1_expression_and_PD-1_checkpoint_pathway_in_cancer","Penicillin_and_cephalosporin_biosynthesis","Pentose_and_glucuronate_interconversions",
                "Pentose_phosphate_pathway","Peptidoglycan_biosynthesis","Peroxisome","Pertussis","Phagosome","Phenazine_biosynthesis","Phenylalanine_metabolism",
                "Phenylalanine,_tyrosine_and_tryptophan_biosynthesis","Phenylpropanoid_biosynthesis","Phosphatidylinositol_signaling_system","Phospholipase_D_signaling_pathway",
                "Phosphonate_and_phosphinate_metabolism","Phosphotransferase_system_(PTS)","PI3K-Akt_signaling_pathway","Plant-pathogen_interaction",
                "Polycyclic_aromatic_hydrocarbon_degradation","Polyketide_sugar_unit_biosynthesis","Porphyrin_and_chlorophyll_metabolism","Primary_immunodeficiency",
                "Prodigiosin_biosynthesis","Progesterone-mediated_oocyte_maturation","Prolactin_signaling_pathway","Propanoate_metabolism","Proteasome",
                "Protein_digestion_and_absorption","Protein_export","Proximal_tubule_bicarbonate_reclamation","Purine_metabolism","Pyrimidine_metabolism",
                "Pyruvate_metabolism","Quorum_sensing","Rap1_signaling_pathway","Ras_signaling_pathway","Regulation_of_actin_cytoskeleton","Renin_secretion",
                "Renin-angiotensin_system","Retinol_metabolism","Retrograde_endocannabinoid_signaling","Riboflavin_metabolism","Ribosome","RNA_degradation","RNA_polymerase",
                "RNA_transport","Salmonella_infection","Secondary_bile_acid_biosynthesis","Selenocompound_metabolism","Sesquiterpenoid_and_triterpenoid_biosynthesis",
                "Shigellosis","Signaling_pathways_regulating_pluripotency_of_stem_cells","Sphingolipid_metabolism","Sphingolipid_signaling_pathway",
                "Staphylococcus_aureus_infection","Starch_and_sucrose_metabolism","Staurosporine_biosynthesis","Steroid_biosynthesis","Steroid_degradation",
                "Stilbenoid,_diarylheptanoid_and_gingerol_biosynthesis","Streptomycin_biosynthesis","Styrene_degradation","Sulfur_metabolism","Sulfur_relay_system",
                "Synthesis_and_degradation_of_ketone_bodies","Taurine_and_hypotaurine_metabolism","Terpenoid_backbone_biosynthesis","Tetracycline_biosynthesis",
                "Thiamine_metabolism","Toll_and_Imd_signaling_pathway","Toll-like_receptor_signaling_pathway","Toluene_degradation",
                "Tropane,_piperidine_and_pyridine_alkaloid_biosynthesis","Tryptophan_metabolism","Tuberculosis","Two-component_system","Tyrosine_metabolism",
                "Ubiquinone_and_other_terpenoid-quinone_biosynthesis","Ubiquitin_mediated_proteolysis","Valine,_leucine_and_isoleucine_biosynthesis",
                "Valine,_leucine_and_isoleucine_degradation","Vancomycin_resistance","Various_types_of_N-glycan_biosynthesis","VEGF_signaling_pathway",
                "Vibrio_cholerae_infection","Vitamin_B6_metabolism","Vitamin_digestion_and_absorption","Wnt_signaling_pathway","Xylene_degradation",
                "Yersinia_infection","Zeatin_biosynthesis")
metabolic_paths <- c("2-Oxocarboxylic_acid_metabolism","ABC_transporters","Alanine,_aspartate_and_glutamate_metabolism",
                     "Amino_sugar_and_nucleotide_sugar_metabolism","Arginine_and_proline_metabolism","Ascorbate_and_aldarate_metabolism",
                     "beta-Alanine_metabolism","Carbohydrate_digestion_and_absorption",
                     "Carbon_metabolism","Citrate_cycle_(TCA_cycle)","Cyanoamino_acid_metabolism","Cysteine_and_methionine_metabolism",
                     "D-Alanine_metabolism","D-Arginine_and_D-ornithine_metabolism","D-Glutamine_and_D-glutamate_metabolism","Degradation_of_aromatic_compounds",
                     "Fatty_acid_degradation","Fatty_acid_metabolism","Fructose_and_mannose_metabolism","Galactose_metabolism",
                     "Glutathione_metabolism","Glycerolipid_metabolism","Glycerophospholipid_metabolism","Glycine,_serine_and_threonine_metabolism","Glycolysis_Gluconeogenesis",
                     "Glycosaminoglycan_degradation","Histidine_metabolism","Lysine_degradation",
                     "Microbial_metabolism_in_diverse_environments","Mineral_absorption",
                     "Nicotinate_and_nicotinamide_metabolism","Nitrogen_metabolism",
                     "Pentose_and_glucuronate_interconversions","Pentose_phosphate_pathway","Phenylalanine_metabolism",
                     "Phosphonate_and_phosphinate_metabolism","Phosphotransferase_system_(PTS)",
                     "Polycyclic_aromatic_hydrocarbon_degradation","Propanoate_metabolism","Protein_digestion_and_absorption",
                     "Purine_metabolism","Pyrimidine_metabolism","Pyruvate_metabolism","Retinol_metabolism","Riboflavin_metabolism",
                     "Selenocompound_metabolism","Sphingolipid_metabolism","Starch_and_sucrose_metabolism",
                     "Sulfur_metabolism","Synthesis_and_degradation_of_ketone_bodies","Taurine_and_hypotaurine_metabolism",
                     "Thiamine_metabolism","Toluene_degradation","Tryptophan_metabolism","Tyrosine_metabolism",
                     "Valine,_leucine_and_isoleucine_degradation","Vitamin_B6_metabolism",
                     "Vitamin_digestion_and_absorption","Xylene_degradation")
completeness <- completeness[metabolic_paths,]
rownames(super) <- super$pathway
super$pathway <- NULL
super <- super[metabolic_paths,]
rownames(normal) <- normal$pathway
normal$pathway <- NULL
normal <- normal[metabolic_paths,]
rownames(success) <- success$pathway
success$pathway <- NULL
success <- success[metabolic_paths,]
rownames(failure) <- failure$pathway
failure$pathway <- NULL
failure <- failure[metabolic_paths,]
completeness <- as.data.frame(t(completeness))
colnames(completeness) <- make.names(colnames(completeness))
keep <- c()
for (x in colnames(completeness)) {
  if (nnzero(completeness[,x]) >= 4) {
    keep <- c(keep, x)}}
rm(x)
completeness <- completeness[,keep]
super <- as.data.frame(t(super))
colnames(super) <- make.names(colnames(super))
super <- super[,keep]
super <- droplevels(super)
normal <- as.data.frame(t(normal))
colnames(normal) <- make.names(colnames(normal))
normal <- normal[,keep]
normal <- droplevels(normal)
success <- as.data.frame(t(success))
colnames(success) <- make.names(colnames(success))
success <- success[,keep]
success <- droplevels(success)
failure <- as.data.frame(t(failure))
colnames(failure) <- make.names(colnames(failure))
failure <- failure[,keep]
failure <- droplevels(failure)
rm(keep)


library(randomForest)
metadata <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metadata.tsv', sep='\t', header=TRUE)
donor <- as.factor(metadata$status)
outcome <- as.factor(metadata$outcome)
rm(metadata)

rf_obj <- randomForest(donor ~ ., data=completeness, ntree=10000, mtry=50,
                       importance=TRUE, replace=FALSE, err.rate=TRUE)
# OOB estimate of  error rate: 0%
# Confusion matrix:
#           normal super class.error
# normal      8     0           0
# super       0     4           0
rf_feat <- importance(rf_obj, type=1, scale=FALSE)
super_rf_importance <- subset(rf_feat, rf_feat > abs(min(rf_feat)))
super_rf_importance <- as.data.frame(super_rf_importance)
super_rf_importance$feature <- rownames(super_rf_importance)
colnames(super_rf_importance) <- c('mda', 'pathway')
super_rf_importance <- super_rf_importance[order(-super_rf_importance$mda),]
rm(rf_obj, rf_feat)

rf_obj <- randomForest(outcome ~ ., data=completeness, ntree=10000, mtry=50,
                       importance=TRUE, replace=FALSE, err.rate=TRUE)
#OOB estimate of  error rate: 50%
#Confusion matrix:
#  failure success class.error
#failure       0       4        1.00
#success       2       6        0.25
rf_feat <- importance(rf_obj, type=1, scale=FALSE)
outcome_rf_importance <- subset(rf_feat, rf_feat > abs(min(rf_feat)))
outcome_rf_importance <- as.data.frame(outcome_rf_importance)
outcome_rf_importance$feature <- rownames(outcome_rf_importance)
colnames(outcome_rf_importance) <- c('mda', 'pathway')
outcome_rf_importance <- outcome_rf_importance[order(-outcome_rf_importance$mda),]
rm(rf_obj, rf_feat)


super_completeness <- completeness[,super_rf_importance$pathway]
super_completeness <- droplevels(super_completeness)
super_completeness <- as.data.frame(t(super_completeness))
super <- super_completeness[,c("X01044A","X01090A","X01092A","X01093A")]
super <- as.data.frame(t(super))
normal <- super_completeness[,c("X05042G","X05053B","X05080M","X05098A","X28043C","X28045A","X28045D","X28047D")]
normal <- as.data.frame(t(normal))
rm(super_completeness)

outcome_completeness <- completeness[,outcome_rf_importance$pathway]
outcome_completeness <- droplevels(outcome_completeness)
outcome_completeness <- as.data.frame(t(outcome_completeness))
success <- outcome_completeness[,c("X01044A","X01090A","X01092A","X01093A","X05080M","X05098A","X28043C","X28045D")]
success <- as.data.frame(t(success))
failure <- outcome_completeness[,c("X05042G","X05053B","X28045A","X28047D")]
failure <- as.data.frame(t(failure))
rm(outcome_completeness)

super_pvals <- c()
for (x in 1:ncol(super)) {super_pvals[x] <- wilcox.test(super[,x], normal[,x], exact=FALSE)$p.value}
super_pvals <- p.adjust(super_pvals, method='BH')

outcome_pvals <- c()
for (x in 1:ncol(success)) {outcome_pvals[x] <- wilcox.test(success[,x], failure[,x], exact=FALSE)$p.value}
outcome_pvals <- p.adjust(outcome_pvals, method='BH')

super_med <- c()
normal_med <- c()
for (x in 1:ncol(super)) {
  super_med[x] <- median(super[,x])
  normal_med[x] <- median(normal[,x])
  
  ymax <- round(max(c(max(super[,x]), max(normal[,x]))) * 1.2)
  formattedName <- gsub('_', ' ', colnames(super)[x])
  formattedName <- gsub('\\.', '-', formattedName)
  mda <- paste('MDA:', super_rf_importance$mda[x])

  fileName <- paste0('~/Desktop/repos/Jenior_Consortia_2022/results/pathway_completeness/', colnames(super)[x], '.png')
  png(filename=fileName, units='in', width=3, height=5, res=300)
  
  par(mar=c(3, 3, 2, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
  plot(0, type='n', xlab='', ylab='KEGG Pathway Completeness (%)', xaxt='n', yaxt='n', 
       xlim=c(0,1), ylim=c(0,ymax), main=formattedName, cex.main=0.8)
  stripchart(at=0.25, vertical=TRUE, jitter(super[,x], amount=1e-5), lwd=2,
             pch=21, bg='deepskyblue3', method='jitter', jitter=0.12, cex=2, add=TRUE)
  stripchart(at=0.75, vertical=TRUE, jitter(normal[,x], amount=1e-5), lwd=2,
             pch=21, bg='firebrick3', method='jitter', jitter=0.12, cex=2, add=TRUE)
  axis(side=2, at=seq(0,ymax,ymax/5), cex.axis=0.8, lwd=2)
  segments(x0=0.1, y0=median(super[,x]), x1=0.4, lwd=3)
  segments(x0=0.6, y0=median(normal[,x]), x1=0.9, lwd=3)
  segments(x0=0.25, y0=ymax*0.92, x1=0.75, lwd=2)
  box(lwd=2)
  
  par(xpd=TRUE)
  text(x=0.25, y=-(ymax*0.07), labels='Super\nDonor', cex=1.2)
  text(x=0.75, y=-(ymax*0.07), labels='Normal\nDonors', cex=1.2)
  par(xpd=FALSE)

  if (super_pvals[x] <= 0.05) {
    text(x=0.5, y=ymax*0.95, '*', cex=1.3, font=2)
  } else {text(x=0.5, y=ymax*0.95, 'n.s.')}
  
  dev.off()
}

library(vioplot)
super_all <- as.vector(t(super))
normal_all <- as.vector(t(normal))
pval <- wilcox.test(super_all, normal_all, exact=FALSE)$p.value
ymax <- 100

fileName <- '~/Desktop/repos/Jenior_Consortia_2022/results/pathway_completeness/all_pathways.png'
png(filename=fileName, units='in', width=3, height=4.5, res=300)
par(mar=c(3, 3, 2, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
vioplot(super_all, normal_all, col=c('blue3', 'green3'), main='All Pathways',
        ylim=c(0, ymax), ylab='KEGG Pathway Completeness (%)', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')


png(filename='~/Desktop/repos/Jenior_Consortia_2022/results/pathway_completeness/catabolic_pathways.png', 
    units='in', width=3, height=4.5, res=300)
par(mar=c(3,3,2,0.5), las=1, mgp=c(1.5,0.7,0), xpd=FALSE, lwd=2.5)
plot(0, type='n', ylim=c(0,70), xlim=c(0,1.5), yaxt='n', yaxs='i', main='Catabolic Pathways',
     ylab='KEGG Pathway Completeness (%)', xlab='', xaxt='n', cex.lab=1.1)
axis(side=2, cex.axis=0.8, lwd=2.5)
par(xpd=TRUE)
text(x=c(0.4,1.15), y=-5, labels=c('Super\ndonor','Normal\ndonors'), cex=1.2)
par(xpd=FALSE)
barplot(c(median(super_all),median(normal_all)), width=0.6, lwd=2.5, col=c('deepskyblue3','firebrick3'), yaxt='n', add=TRUE)
segments(x0=0.2, y0=quantile(super_all,0.25), x1=0.6, lwd=2.5)
segments(x0=0.2, y0=quantile(super_all,0.75), x1=0.6, lwd=2.5)
segments(x0=0.4, y0=quantile(super_all,0.25), y1=quantile(super_all,0.75), lwd=2.5)
segments(x0=0.95, y0=quantile(normal_all,0.25), x1=1.35, lwd=2.5)
segments(x0=0.95, y0=quantile(normal_all,0.75), x1=1.35, lwd=2.5)
segments(x0=1.15, y0=quantile(normal_all,0.25), y1=quantile(normal_all,0.75), lwd=2.5)
dev.off()





pval <- wilcox.test(super_med, normal_med, exact=FALSE)$p.value
ymax <- 100


fileName <- '~/Desktop/repos/Jenior_Consortia_2022/results/pathway_completeness/median_pathways.png'
png(filename=fileName, units='in', width=3, height=4.5, res=300)
par(mar=c(3, 3, 2, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
vioplot(super_med, normal_med, col=c('blue3', 'green3'), main='All Pathways',
        ylim=c(0, ymax), ylab='KEGG Pathway Completeness (%)', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0,ymax,ymax/5), cex.axis=0.8, lwd=2.5)
segments(x0=1, y0=ymax*0.92, x1=2, lwd=2)
if (pval <= 0.05) {
  text(x=1.5, y=ymax*0.95, '*', cex=1.3, font=2)
} else {text(x=1.5, y=ymax*0.95, 'n.s.')}
par(xpd=TRUE)
text(x=1, y=-(ymax*0.08), labels='Super\nDonor', cex=1.2)
text(x=2, y=-(ymax*0.08), labels='Normal\nDonors', cex=1.2)
par(xpd=FALSE)
dev.off()










for (x in 1:ncol(super)) {
  ymax <- round(max(c(max(super[,x]), max(normal[,x]))) * 1.2)
  formattedName <- gsub('_', ' ', colnames(super)[x])
  formattedName <- gsub('\\.', '-', formattedName)
  mda <- paste('MDA:', super_rf_importance$mda[x])
  
  fileName <- paste0('~/Desktop/repos/Jenior_Consortia_2022/results/pathway_completeness/', colnames(super)[x], '.png')
  png(filename=fileName, units='in', width=3, height=5, res=300)
  
  par(mar=c(3, 3, 2, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
  plot(0, type='n', xlab='', ylab='KEGG Pathway Completeness (%)', xaxt='n', yaxt='n', 
       xlim=c(0,1), ylim=c(0,ymax), main=formattedName, cex.main=0.8)
  stripchart(at=0.25, vertical=TRUE, jitter(super[,x], amount=1e-5), lwd=2,
             pch=21, bg='blue3', method='jitter', jitter=0.12, cex=2, add=TRUE)
  stripchart(at=0.75, vertical=TRUE, jitter(normal[,x], amount=1e-5), lwd=2,
             pch=21, bg='green3', method='jitter', jitter=0.12, cex=2, add=TRUE)
  axis(side=2, at=seq(0,ymax,ymax/5), cex.axis=0.8, lwd=2)
  segments(x0=0.1, y0=median(super[,x]), x1=0.4, lwd=3)
  segments(x0=0.6, y0=median(normal[,x]), x1=0.9, lwd=3)
  segments(x0=0.25, y0=ymax*0.92, x1=0.75, lwd=2)
  legend('bottomright', mda, bty='n', cex=1.2, pt.cex=0)
  box(lwd=2)
  
  par(xpd=TRUE)
  text(x=0.25, y=-(ymax*0.07), labels='Super\nDonor', cex=1.2)
  text(x=0.75, y=-(ymax*0.07), labels='Normal\nDonors', cex=1.2)
  par(xpd=FALSE)
  
  if (super_pvals[x] <= 0.05) {
    text(x=0.5, y=ymax*0.95, '*', cex=1.3, font=2)
  } else {text(x=0.5, y=ymax*0.95, 'n.s.')}
  
  dev.off()
}

super_all <- as.vector(t(super))
normal_all <- as.vector(t(normal))
pval <- wilcox.test(super_all, normal_all, exact=FALSE)$p.value
ymax <- round(max(c(super_all, normal_all)) * 1.2)

library(vioplot)
fileName <- '~/Desktop/repos/Jenior_Consortia_2022/results/pathway_completeness/all_pathways.png'
png(filename=fileName, units='in', width=3, height=4.5, res=300)
par(mar=c(3, 3, 2, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
vioplot(super_all, normal_all, col=c('blue3', 'green3'), main='All Pathways',
        ylim=c(0, ymax), ylab='KEGG Pathway Completeness (%)', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0,ymax,ymax/5), cex.axis=0.8, lwd=2.5)
segments(x0=1, y0=ymax*0.92, x1=2, lwd=2)
if (pval <= 0.05) {
  text(x=1.5, y=ymax*0.95, '*', cex=1.3, font=2)
} else {text(x=1.5, y=ymax*0.95, 'n.s.')}
par(xpd=TRUE)
text(x=1, y=-(ymax*0.08), labels='Super\nDonor', cex=1.2)
text(x=2, y=-(ymax*0.08), labels='Normal\nDonors', cex=1.2)
par(xpd=FALSE)
dev.off()

#------------------------------------------------------------------------------------------------------#

# Pathway abundances
s01044A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/01044A/pathway_abundances.tsv', sep='\t', header=TRUE)
s01090A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/01090A/pathway_abundances.tsv', sep='\t', header=TRUE)
s01092A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/01092A/pathway_abundances.tsv', sep='\t', header=TRUE)
s01093A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/01093A/pathway_abundances.tsv', sep='\t', header=TRUE)
n05042G <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/05042G/pathway_abundances.tsv', sep='\t', header=TRUE)
n05053B <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/05053B/pathway_abundances.tsv', sep='\t', header=TRUE)
n05080M <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/05080M/pathway_abundances.tsv', sep='\t', header=TRUE)
n05098A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/05098A/pathway_abundances.tsv', sep='\t', header=TRUE)
n28043C <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/28043C/pathway_abundances.tsv', sep='\t', header=TRUE)
n28045A <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/28045A/pathway_abundances.tsv', sep='\t', header=TRUE)
n28045D <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/28045D/pathway_abundances.tsv', sep='\t', header=TRUE)
n28047D <- read.delim('~/Desktop/active_projects/CDI_FMT_project/fmt_metaG/pathway_completeness/28047D/pathway_abundances.tsv', sep='\t', header=TRUE)

# Assemble tables
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

# Focus pathways
success <- success[keep_paths,]
failure <- failure[keep_paths,]
super <- super[keep_paths,]
normal <- normal[keep_paths,]

# Subsample
library(vegan)
sub_outcome <- ceiling(min(c(min(as.vector(colSums(success))), min(as.vector(colSums(failure))))) * 0.9)
for (x in 1:ncol(success)) {success[,x] <- as.vector(rrarefy(success[,x], sample=sub_outcome))}
for (y in 1:ncol(failure)) {failure[,y] <- as.vector(rrarefy(failure[,y], sample=sub_outcome))}
sub_donor <- ceiling(min(c(min(as.vector(colSums(super))),min(as.vector(colSums(normal))))) * 0.9)
for (x in 1:ncol(super)) {super[,x] <- as.vector(rrarefy(super[,x], sample=sub_donor))}
for (y in 1:ncol(normal)) {normal[,y] <- as.vector(rrarefy(normal[,y], sample=sub_donor))}

# Build final tables
outcome <- merge(success, failure, by='row.names')
rownames(outcome) <- outcome$Row.names
outcome$Row.names <- NULL
outcome <- subset(outcome, rowSums(outcome) >= 12)
outcome <- as.data.frame(t(outcome))
colnames(outcome) <- make.names(colnames(outcome))
donor <- merge(super, normal, by='row.names')
rownames(donor) <- donor$Row.names
donor$Row.names <- NULL
donor <- subset(donor, rowSums(donor) >= 12)
donor <- as.data.frame(t(donor))
colnames(donor) <- make.names(colnames(donor))

# Supervised macine learning
library(randomForest)
set.seed(906801)
condition <- as.factor(c('success','success','success','success','success','success','success','success','failure','failure','failure','failure'))
rf_outcome <- randomForest(condition ~ ., data=outcome, ntree=5000, mtry=15,
                         importance=TRUE, replace=FALSE, err.rate=TRUE)
# OOB estimate of  error rate: 58.33%
# Confusion matrix:
#           failure success class.error
# failure       0       4       1.000
# success       3       5       0.375
rf_outcome <- importance(rf_outcome, type=1, scale=FALSE)
rf_outcome <- subset(rf_outcome, rf_outcome > abs(min(rf_outcome)))
rf_outcome <- as.data.frame(rf_outcome)
rf_outcome$feature <- rownames(rf_outcome)
colnames(rf_outcome) <- c('mda', 'pathway')
rf_outcome <- rf_outcome[order(-rf_outcome$mda),]
if (nrow(rf_outcome >= 18)) {rf_outcome <- rf_outcome[1:17,]}
outcome <- outcome[,rf_outcome$pathway]
condition <- as.factor(c('super','super','super','super','normal','normal','normal','normal','normal','normal','normal','normal'))
rf_donor <- randomForest(condition ~ ., data=donor, ntree=5000, mtry=15,
                       importance=TRUE, replace=FALSE, err.rate=TRUE)
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
if (nrow(rf_donor >= 18)) {rf_donor <- rf_donor[1:17,]}
donor <- donor[,rf_donor$pathway]

# Prep for radar plot
success <- outcome[c('s01044A','s01090A','s01092A','s01093A','n05098A','n05080M','n28043C','n28045D'),]
success <- as.data.frame(apply(success, MARGIN=2, median))
failure <- outcome[c('n05042G','n05053B','n28045A','n28047D'),]
failure <- as.data.frame(apply(failure, MARGIN=2, median))
outcome <- merge(success, failure, by='row.names')
rownames(outcome) <- outcome$Row.names
outcome$Row.names <- NULL
colnames(outcome) <- c('success', 'failure')
super <- donor[c('s01044A','s01090A','s01092A','s01093A'),]
super <- as.data.frame(apply(super, MARGIN=2, median))
normal <- donor[c('n05042G','n05053B','n28045A','n28047D','n05098A','n05080M','n28043C','n28045D'),]
normal <- as.data.frame(apply(normal, MARGIN=2, median))
donor <- merge(super, normal, by='row.names')
rownames(donor) <- donor$Row.names
donor$Row.names <- NULL
colnames(donor) <- c('super', 'normal')

# Calculate percentages
outcome$success_perc <- (outcome$success / (outcome$success + outcome$failure)) * 100.0
outcome$failure_perc <- (outcome$failure / (outcome$success + outcome$failure)) * 100.0
outcome$success <- NULL
outcome$failure <- NULL
colnames(outcome) <- c('success', 'failure')
donor$super_perc <- (donor$super / (donor$super + donor$normal)) * 100.0
donor$normal_perc <- (donor$normal / (donor$super + donor$normal)) * 100.0
donor$super <- NULL
donor$normal <- NULL
colnames(donor) <- c('super', 'normal')

# Eliminate host-specific columns
outcome <- as.data.frame(t(outcome))
outcome$Salivary_secretion <- NULL
outcome$Ferroptosis <- NULL
outcome$Pancreatic_secretion <- NULL
outcome$Pathways_in_cancer <- NULL
outcome$Primary_bile_acid_biosynthesis <- NULL
outcome$Kaposi_sarcoma.associated_herpesvirus_infection <- NULL
outcome$Alzheimer_disease <- NULL
outcome$Apoptosis <- NULL
outcome$Choline_metabolism_in_cancer <- NULL
outcome$Estrogen_signaling_pathway <- NULL
outcome$Systemic_lupus_erythematosus <- NULL
outcome$Peroxisome <- NULL
colnames(outcome) <- gsub('_', ' ', colnames(outcome))
colnames(outcome) <- gsub('\\.', ', ', colnames(outcome))
donor <- as.data.frame(t(donor))
donor$Salivary_secretion <- NULL
donor$Ferroptosis <- NULL
donor$Pancreatic_secretion <- NULL
donor$Pathways_in_cancer <- NULL
donor$Primary_bile_acid_biosynthesis <- NULL
donor$Kaposi_sarcoma.associated_herpesvirus_infection <- NULL
donor$Alzheimer_disease <- NULL
donor$Apoptosis <- NULL
donor$Choline_metabolism_in_cancer <- NULL
donor$Estrogen_signaling_pathway <- NULL
donor$Systemic_lupus_erythematosus <- NULL
donor$Peroxisome <- NULL
colnames(donor) <- gsub('_', ' ', colnames(donor))
colnames(donor) <- gsub('\\.', ', ', colnames(donor))

# Radar plot tables
outcome <- rbind(rep(100, ncol(outcome)) , rep(0, ncol(outcome)) , outcome)
write.table(outcome, file='~/Desktop/repos/Jenior_Consortia_2022/data/outcome_radar.tsv', append = FALSE, quote=FALSE, sep='\t')
donor <- rbind(rep(100, ncol(donor)) , rep(0, ncol(donor)) , donor)
write.table(donor, file='~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar.tsv', append = FALSE, quote=FALSE, sep='\t')


#-----------------------------#


# Read in previous results for better figure reporduceability
outcome <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/outcome_radar.tsv', sep='\t', header=TRUE, row.names=1)
colnames(outcome) <- gsub('_', ' ', colnames(outcome))
colnames(outcome) <- gsub('\\.', '\n', colnames(outcome))
donor <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar.tsv', sep='\t', header=TRUE, row.names=1)
colnames(donor) <- gsub('_', ' ', colnames(donor))
colnames(donor) <- gsub('\\.', '\n', colnames(donor))

# Generate figures
library(fmsb)
library(scales)
success_col <- 'chartreuse3'
failure_col <- 'turquoise3'
super_col <- 'deepskyblue3'
normal_col <- 'firebrick3'

png(file='~/Desktop/repos/Jenior_Consortia_2022/results/donor_pathway_radar.png', units='in', width=5, height=5, res=300)
par(mar=c(0,2,0,2), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(donor, pcol=c(super_col,normal_col), pfcol=c(alpha(super_col,0.7),alpha(normal_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey30', cglty=5, cglwd=1)
text(x=0, y=0, 'OOB = 8.33%', font=2, cex=0.8)
legend('bottomright', legend=c('Super Donor', 'Normal Donors'), col=c(super_col,normal_col),
       pt.bg=c(alpha(super_col,0.7),alpha(normal_col,0.7)), pch=22, pt.cex=1.5, cex=0.8, bty='n')
dev.off()

png(file='~/Desktop/repos/Jenior_Consortia_2022/results/outcome_pathway_radar.png', units='in', width=5, height=5, res=300)
par(mar=c(0,2,0,2), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(outcome, pcol=c(success_col,failure_col), pfcol=c(alpha(success_col,0.7),alpha(failure_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey30', cglty=5, cglwd=1)
text(x=0, y=0, 'OOB = 58.33%', font=2, cex=0.8)
legend('bottomright', legend=c('Successful FMT', 'Failed FMT'), col=c(success_col,failure_col),
       pt.bg=c(alpha(success_col,0.7),alpha(failure_col,0.7)), pch=22, pt.cex=1.5, cex=0.8, bty='n')
dev.off()








