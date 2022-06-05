
# Read in and format data
metabolome <- read.delim('~/Desktop/Jenior_Consortia_2022/data/combined_metabolome.tsv', sep='\t', header=TRUE)
exclude <- c('ZDXPYRJPNDTMRX-VKHMYHEASA-N','CQIUKKVOEOPUDV-IYSWYEEDSA-N',
             'DAXYUDFNWXHGBE-VKCGGMIFSA-N','LELBFTMXCIIKKX-UHFFFAOYNA-N',
             'SPLUHYDIQKURQW-UHFFFAOYSA-N','HJBWJAPEBGSQPR-GQCTYLIASA-N')
metabolome <- subset(metabolome, !INCHIKEY %in% exclude)
metabolome$MetaboliteName <- make.names(metabolome$MetaboliteName)
#inchikey <- metabolome[,c(1:2)]
metabolome$INCHIKEY <- NULL
rownames(metabolome) <- metabolome$MetaboliteName
metabolome$MetaboliteName <- NULL
metabolome$blank_1 <- NULL
metabolome$blank_2 <- NULL
metabolome$blank_3 <- NULL
metabolome$blank_4 <- NULL
metabolome <- as.data.frame(t(metabolome))
metadata <- read.delim('~/Desktop/Jenior_Consortia_2022/data/metabolome_metadata.tsv', sep='\t', header=TRUE)
metabolome <- as.data.frame(merge(metadata, metabolome, by.x='Label', by.y='row.names'))
rownames(metabolome) <- metabolome$Label
metabolome$Label <- NULL
rm(exclude, metadata)

#-------------------------------------------------------------#

# Define functions

# Create subset of decreased growth only
decreasedSubset <- function(data, species) {
  keep <- subset(data, R20291_secondary_growth == 'no')
  keep <- subset(keep, Depleting_species %in% c(species,'Fresh_media'))
  keep$R20291_secondary_growth <- NULL
  keep <- droplevels(keep)
  return(keep)}

# Test for significant differences - decreased
sigDecreased <- function(data) {
  group1 <- subset(data, Depleting_species == 'Fresh_media')
  group1$Depleting_species <- NULL
  group2 <- subset(data, Depleting_species != 'Fresh_media')
  group2$Depleting_species <- NULL
  keep <- c()
  for (x in colnames(group1)) {
    test1 <- as.vector(group1[,x])
    test2 <- as.vector(group2[,x])
    pval <- as.numeric(wilcox.test(test1, test2, exact=FALSE)$p.value)
    if (is.na(pval)) {pval = 1}
    if (pval <= 0.05) {if (median(test1) > median(test2)) {keep <- c(keep, x)}}
  }
  group2 <- group2[, keep]
  return(group2)}

findDecrease <- function(data, species='default') {
  decreased <- decreasedSubset(data, species)
  sig_decreased <- sigDecreased(decreased)
  return(sig_decreased)}


# Test for edges of competition
identifyCompetition <- function(data1, data2) {
  group1 <- colnames(data1)
  group2 <- colnames(data2)
  keep <- intersect(colnames(data1), colnames(data2))
  group1 <- data1[,keep]
  group2 <- data2[,keep]
  competition <- rbind(group1, group2)
  return(competition)}

#----------------------------------------------------------#

# Find consumed subtrate
Bproducta_decreased <- findDecrease(metabolome, 'B_producta_DSM2950')
Bpseudocatenulatum_decreased <- findDecrease(metabolome, 'B_pseudocatenulatum_ATCC27919')
Cdifficile_decreased <- findDecrease(metabolome, 'C_difficile_R20291')
rm(decreasedSubset, metabolome)

# Check edges of competition
Bproducta_competition <- identifyCompetition(Bproducta_decreased, Cdifficile_decreased)
Bpseudocatenulatum_competition <- identifyCompetition(Bpseudocatenulatum_decreased, Cdifficile_decreased)

# Identify those unique to each organism
Bproducta_competition_only <- setdiff(colnames(Bproducta_competition), colnames(Bpseudocatenulatum_competition))
Bproducta_competition_only <- Bproducta_competition[,Bproducta_competition_only]
Bpseudocatenulatum_competition_only <- setdiff(colnames(Bpseudocatenulatum_competition), colnames(Bproducta_competition))
Bpseudocatenulatum_competition_only <- Bpseudocatenulatum_competition[,Bpseudocatenulatum_competition_only]

# Find overlap
all_competition <- intersect(colnames(Bproducta_competition), colnames(Bpseudocatenulatum_competition))
all_competition <- Cdifficile_decreased[,all_competition]

#-------------------------------------------------------#

# Venn diagrams for competitive consortia with supporting strip charts for specific catagories
library(VennDiagram)
shared_12 <- intersect(colnames(Bproducta_decreased), colnames(Bpseudocatenulatum_decreased))
contested_1 <- intersect(colnames(Bproducta_decreased), colnames(Cdifficile_decreased))
contested_2 <- intersect(colnames(Bpseudocatenulatum_decreased), colnames(Cdifficile_decreased))
contested_12 <- intersect(shared_12, colnames(Cdifficile_decreased))

png(file='~/Desktop/repos/Jenior_Consortia_2022/results/invitro_competition.png', units='in', width=5, height=4, res=300)
draw.quad.venn(area1=length(contested_1), 
               area2=length(contested_2), 
               area3=length(contested_3), 
               area4=length(contested_4), 
               n12=length(contested_12), 
               n23=length(contested_23), 
               n13=length(contested_13), 
               n14=length(contested_14),
               n24=length(contested_24),
               n34=length(contested_34),
               n123=length(contested_123), 
               n124=length(contested_124), 
               n234=length(contested_234), 
               n134=length(contested_134), 
               n1234=length(contested_1234), 
               category=c('B. producta','B. pseudocatenulatum','P. distasonis','E. rectale'), 
               cat.fontfamily='sans', cat.fontface='italic', fontfamily='sans', cat.cex=1,
               fill=c('springgreen2','tan1','thistle2','paleturquoise3'), lwd=0.0001)
dev.off()
