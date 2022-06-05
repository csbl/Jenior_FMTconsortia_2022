
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
require('vegan')
library(vegan)
metag_super <- metag[,c('fmt.metaG.01044A.reads2genes.norm.tsv','fmt.metaG.01090A.reads2genes.norm.tsv','fmt.metaG.01092A.reads2genes.norm.tsv','fmt.metaG.01093A.reads2genes.norm.tsv')]
metag_normal <- metag[,c('fmt.metaG.05042G.reads2genes.norm.tsv','fmt.metaG.05053B.reads2genes.norm.tsv','fmt.metaG.05080M.reads2genes.norm.tsv','fmt.metaG.05098A.reads2genes.norm.tsv','fmt.metaG.28043C.reads2genes.norm.tsv','fmt.metaG.28045A.reads2genes.norm.tsv','fmt.metaG.28045D.reads2genes.norm.tsv','fmt.metaG.28047D.reads2genes.norm.tsv')]
sub_level <- ceiling(min(c(as.vector(colSums(metag_super)), as.vector(colSums(metag_normal)))) * 0.85)
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
library(fmsb)
library(scales)
super_col <- 'deepskyblue3'
normal_col <- 'firebrick3'
path_labels <- c("Amino acid\nbiosynthesis\n(12889)","Carbon\nmetabolism\n(4570)","Antimicrobial peptide\nresistance\n(701)",
                 "Bacterial\nchemotaxis\n(729)","Sulfur\nrelay system\n(842)","Nucleotide\nexcision repair\n(707)",
                 "Siderophore peptide\nbiosynthesis\n(30)","Pantothenate & CoA\nbiosynthesis\n(199)",
                 "Sphingolipid\nmetabolism\n(50)","Propanoate\nmetabolism\n(17)","Oxocarboxylic acid\nmetabolism\n(53)")
colnames(metag_med) <- path_labels

#pdf(file='~/Desktop/repos/Jenior_FMT_2021/results/pathway_radar.pdf', width=7.5, height=5)
png(file='~/Desktop/repos/Jenior_FMT_2021/results/pathway_radar.png', units='in', width=5, height=5, res=300)
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

pdf(file='~/Desktop/repos/Jenior_FMT_2021/results/pathway_abundances.pdf', width=9, height=6)
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
rm(s01044A, s01090A, s01092A, s01093A)
rm(n05042G, n05053B, n05080M, n05098A)
rm(n28043C, n28045A, n28045D, n28047D)

library(Matrix)
completeness <- merge(super, normal, by='pathway')
rownames(completeness) <- completeness$pathway
completeness$pathway <- NULL

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
rm(keep)


library(randomForest)
metadata <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/metadata.tsv', sep='\t', header=TRUE)
donor <- as.factor(metadata$status)
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

super_completeness <- completeness[,super_rf_importance$pathway]
super_completeness <- droplevels(super_completeness)
super_completeness <- as.data.frame(t(super_completeness))
super <- super_completeness[,c("X01044A","X01090A","X01092A","X01093A")]
super <- as.data.frame(t(super))
normal <- super_completeness[,c("X05042G","X05053B","X05080M","X05098A","X28043C","X28045A","X28045D","X28047D")]
normal <- as.data.frame(t(normal))
rm(super_completeness)



library(vioplot)
super_all <- as.vector(t(super))
normal_all <- as.vector(t(normal))
pval <- wilcox.test(super_all, normal_all, exact=FALSE)$p.value
ymax <- 100

fileName <- '~/Desktop/repos/Jenior_Consortia_2022/all_pathways.png'
png(filename=fileName, units='in', width=3, height=4.5, res=300)
par(mar=c(3, 3, 2, 0.5), mgp=c(2.1, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
vioplot(super_all, normal_all, col=c('blue3', 'green3'), main='All Pathways',
        ylim=c(0, ymax), ylab='KEGG Pathway Completeness (%)', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')

fileName <- '~/Desktop/repos/Jenior_FMT_2021/results/pathway_completeness/all_pathways.png'
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
rm(s01044A,s01090A,s01092A,s01093A,n05098A,n05080M,n28043C,n28045D,n05042G,n05053B,n28045A,n28047D)

# Focus on metabolic pathways
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
super <- super[metabolic_paths,]
normal <- normal[metabolic_paths,]


# Subsample
library(vegan)
sub_donor <- ceiling(min(c(min(as.vector(colSums(super))),min(as.vector(colSums(normal))))) * 0.85)
for (x in 1:ncol(super)) {super[,x] <- as.vector(rrarefy(super[,x], sample=sub_donor))}
for (y in 1:ncol(normal)) {normal[,y] <- as.vector(rrarefy(normal[,y], sample=sub_donor))}

# Build final tables
donor <- merge(super, normal, by='row.names')
rownames(donor) <- donor$Row.names
donor$Row.names <- NULL
donor <- subset(donor, rowSums(donor) >= 12)
donor <- as.data.frame(t(donor))
colnames(donor) <- make.names(colnames(donor))

# Supervised macine learning
library(randomForest)
set.seed(906801)
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

# Read in previous results for better figure reporduceability
donor <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar_abund.tsv', sep='\t', header=TRUE, row.names=1)
#donor <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar_abund_all.tsv', sep='\t', header=TRUE, row.names=1)
#donor$Apoptosis <- NULL
#donor$Choline_metabolism.in_cancer <- NULL
#donor$Systemic_lupus.erythematosus <- NULL

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
super_col <- 'deepskyblue3'
normal_col <- 'firebrick3'
png(file='~/Desktop/repos/Jenior_Consortia_2022/results/figure_2A.png', units='in', width=5, height=5, res=300)
par(mar=c(0,2,0,2), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(donor, pcol=c(normal_col, super_col), pfcol=c(alpha(normal_col,0.7),alpha(super_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey60', cglty=5, cglwd=1)
text(x=-0.2, y=0.99,'80%', cex=0.5, font=1, srt=10)
text(x=-0.04, y=0.22, '20%', cex=0.4, font=1, srt=17)
legend('bottomright', legend=c('Normal Donors','Super Donor'), col=c(normal_col,super_col),
       pt.bg=c(alpha(normal_col,0.7), alpha(super_col,0.7)), pch=22, pt.cex=1.5, cex=0.8, bty='n')
#text(x=0, y=1.4, 'Gene Abundance by Pathway', font=2, cex=1)
dev.off()

#--------------#

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
rm(s01044A,s01090A,s01092A,s01093A,n05098A,n05080M,n28043C,n28045D,n05042G,n05053B,n28045A,n28047D)

# Break into the 2 groups and plot
library(vegan)
library(vioplot)
pathwayPlot <- function(pathway) {
  super_path <- as.vector(super[,pathway])
  normal_path <- as.vector(normal[,pathway])
  absDiff <- abs(median(super_path) - median(normal_path))
  print(absDiff)

  pval <- wilcox.test(super_path, normal_path, exact=FALSE)$p.value
  plot_title <- gsub('_', ' ', pathway)
  plot_title <- gsub('and', '\\&', plot_title)
  ymax <- max(c(max(super_path), max(normal_path))) * 1.2
  
  par(mar=c(3, 3.2, 2, 0.5), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i', lwd=2, las=1)
  vioplot(normal_path, super_path, col=c('firebrick3','deepskyblue3'), main=plot_title, cex.main=0.9,
          ylim=c(0, ymax), ylab='Gene Abundance', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  axis(side=2, at=round(seq(0,ymax,ymax/5)), cex.axis=0.7, lwd=2.5)
  #segments(x0=1, y0=ymax*0.92, x1=2, lwd=2)
  #if (pval <= 0.05) {
  #  text(x=1.5, y=ymax*0.96, '*', cex=1.3, font=2)
  #} else {text(x=1.5, y=ymax*0.96, 'n.s.', cex=0.8)}
  par(xpd=TRUE)
  text(x=1, y=-(ymax*0.12), labels='Normal\nDonors', font=2)
  text(x=2, y=-(ymax*0.12), labels='Super\nDonor', font=2)
  par(xpd=FALSE)}

# Overrepresented in Normal
png(file='~/Desktop/repos/Jenior_Consortia_2022/results/normal_associated.png', units='in', width=6, height=4, res=300)
layout(matrix(c(1,2,3,
                1,2,3,
                4,5,3,
                4,5,3), nrow=4, ncol=3, byrow=TRUE))
pathwayPlot('Glutathione_metabolism')
pathwayPlot('Valine,_leucine_and_isoleucine_degradation')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), xlab='', ylab='')
box(col='white', lwd=5)
pathwayPlot('Amino_sugar_and_nucleotide_sugar_metabolism')
pathwayPlot('Starch_and_sucrose_metabolism')
dev.off()

# Overrepresented in Super
png(file='~/Desktop/repos/Jenior_Consortia_2022/results/super_associated.png', units='in', width=6, height=4, res=300)
layout(matrix(c(1,2,3,
                1,2,4,
                5,6,4,
                5,6,7), nrow=4, ncol=3, byrow=TRUE))
pathwayPlot('Purine_metabolism')
pathwayPlot('Arginine_and_proline_metabolism')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), xlab='', ylab='')
box(col='white', lwd=5)
pathwayPlot('Cysteine_and_methionine_metabolism')
pathwayPlot('Fructose_and_mannose_metabolism')
pathwayPlot('Glycosaminoglycan_degradation')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), xlab='', ylab='')
box(col='white', lwd=5)
dev.off()

#------------------------------------------------------------------------------------------------------#

# Metabolic pathway completeness

# Pathway completeness
s01044A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_01044A.tsv', sep='\t', header=TRUE)
s01090A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_01090A.tsv', sep='\t', header=TRUE)
s01092A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_01092A.tsv', sep='\t', header=TRUE)
s01093A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_01093A.tsv', sep='\t', header=TRUE)
n05042G <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_05042G.tsv', sep='\t', header=TRUE)
n05053B <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_05053B.tsv', sep='\t', header=TRUE)
n05080M <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_05080M.tsv', sep='\t', header=TRUE)
n05098A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_05098A.tsv', sep='\t', header=TRUE)
n28043C <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_28043C.tsv', sep='\t', header=TRUE)
n28045A <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_28045A.tsv', sep='\t', header=TRUE)
n28045D <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_28045D.tsv', sep='\t', header=TRUE)
n28047D <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_28047D.tsv', sep='\t', header=TRUE)

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
rm(s01044A,s01090A,s01092A,s01093A,n05098A,n05080M,n28043C,n28045D,n05042G,n05053B,n28045A,n28047D)

# Focus on metabolic pathways
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
super <- super[metabolic_paths,]
normal <- normal[metabolic_paths,]

# Build final tables
donor <- merge(super, normal, by='row.names')
rownames(donor) <- donor$Row.names
donor$Row.names <- NULL
donor <- subset(donor, rowSums(donor) > 0)
donor <- as.data.frame(t(donor))
colnames(donor) <- make.names(colnames(donor))

# Supervised mahcine learning
library(randomForest)
set.seed(906801)
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
rf_donor <- rf_donor[1:20,]
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

# Radar plot tables
donor <- as.data.frame(t(donor))
donor <- rbind(rep(floor(min(donor))-5, ncol(donor)) , rep(ceiling(max(donor))+5, ncol(donor)) , donor)
#donor <- rbind(rep(0, ncol(donor)) , rep(ceiling(max(donor))+5, ncol(donor)) , donor)
#donor$Lysosome <- NULL
#donor$Systemic_lupus_erythematosus <- NULL
#write.table(donor, file='~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar_complete.tsv', append = FALSE, quote=FALSE, sep='\t')
#write.table(donor, file='~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar_complete_all.tsv', append = FALSE, quote=FALSE, sep='\t')

# Read in previous results for better figure reporduceability
donor <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar_complete.tsv', sep='\t', header=TRUE, row.names=1)
#donor <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar_complete_all.tsv', sep='\t', header=TRUE, row.names=1)

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
super_col <- 'deepskyblue3'
normal_col <- 'firebrick3'
png(file='~/Desktop/repos/Jenior_Consortia_2022/results/donor_pathway_completeness.png', units='in', width=5, height=5, res=300)
par(mar=c(0,2,0,2), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(donor, pcol=c(super_col,normal_col), pfcol=c(alpha(super_col,0.7),alpha(normal_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.5, cglcol='grey60', cglty=5, cglwd=1)
text(x=0, y=0, '91.67%\nAccuracy', font=2, cex=0.7)
text(x=-0.1, y=c(0.94,0.22), c('80%','20%'), cex=0.5, font=1)
legend('bottomright', legend=c('Super Donor', 'Normal Donors'), col=c(super_col,normal_col),
       pt.bg=c(alpha(super_col,0.7),alpha(normal_col,0.7)), pch=22, pt.cex=1.5, cex=0.8, bty='n')
text(x=0, y=1.4, 'Pathway Completeness', font=2, cex=1)
dev.off()


