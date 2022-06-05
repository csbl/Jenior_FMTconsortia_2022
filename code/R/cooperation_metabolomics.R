
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

# Create subset of primary growth only
primarySubset <- function(data, species) {
  keep <- subset(data, R20291_secondary_growth == 'no')
  keep <- subset(keep, Depleting_species %in% c(species,'Fresh_media'))
  keep$R20291_secondary_growth <- NULL
  keep <- droplevels(keep)
  return(keep)}
# Create subset of secondary growth only
secondarySubset <- function(data, species) {
  keep <- subset(data, Depleting_species == species)
  keep$Depleting_species <- NULL
  keep <- droplevels(keep)
  return(keep)}

# Test for significant differences - Primary
sigIncreased <- function(data) {
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
    if (pval <= 0.05) {if (median(test1) < median(test2)) {keep <- c(keep, x)}}
  }
  group2 <- group2[, keep]
  return(group2)}

findIncrease <- function(data, species='default') {
  increased <- primarySubset(data, species)
  sig_increased <- sigIncreased(increased)
  return(sig_increased)}


# Test for significant differences - Secondary
testCrossfeeding <- function(data1, data2, species) {
  group1 <- subset(data1, R20291_secondary_growth == 'yes')
  group1$R20291_secondary_growth <- NULL
  group1 <- subset(group1, Depleting_species == species)
  group1$Depleting_species <- NULL
  keep <- c()
  for (x in colnames(data2)) {
    test1 <- as.vector(group1[,x])
    test2 <- as.vector(data2[,x])
    pval <- as.numeric(wilcox.test(test1, test2, exact=FALSE)$p.value)
    if (is.na(pval)) {pval = 1}
    if (pval <= 0.05) {if (median(test1) < median(test2)) {keep <- c(keep, x)}}
  }
  group1 <- group1[, keep]
  group2 <- data2[, keep]
  increased <- rbind(group1, group2)
  return(increased)}

#-------------------------------------------------------------#

# vs Fresh media
Blongum_increased <- findIncrease(metabolome, 'B_longum_ATCC55813')
Sthermophilus_increased <- findIncrease(metabolome, 'S_thermophilus_LMD9')

# check for potential cross-feeding
Blongum_cooperation <- testCrossfeeding(metabolome, Blongum_increased, 'B_longum_ATCC55813')
Sthermophilus_cooperation <- testCrossfeeding(metabolome, Sthermophilus_increased, 'S_thermophilus_LMD9')
rm(Blongum_increased, Sthermophilus_increased, metabolome)

# Identify those unique to each organism
Blongum_cooperation_only <- setdiff(colnames(Blongum_cooperation), colnames(Sthermophilus_cooperation))
Blongum_cooperation_only <- Blongum_cooperation[,Blongum_cooperation_only]
Sthermophilus_cooperation_only <- setdiff(colnames(Sthermophilus_cooperation), colnames(Blongum_cooperation))
Sthermophilus_cooperation_only <- Sthermophilus_cooperation[,Sthermophilus_cooperation_only]

# Find overlap
all_cooperation <- intersect(colnames(Blongum_cooperation), colnames(Sthermophilus_cooperation))
all_cooperation <- Cdifficile_decreased[,all_cooperation]

#-------------------------------------------------------#

# Venn diagrams for competitive consortia with supporting strip charts for specific catagories
library(VennDiagram)
shared_12 <- intersect(colnames(Blongum_increased), colnames(Sthermophilus_increased))
crossfed_1 <- intersect(colnames(Blongum_increased), colnames(Cdifficile_increased))
crossfed_2 <- intersect(colnames(Sthermophilus_increased), colnames(Cdifficile_increased))



png(file='~/Desktop/repos/Jenior_Consortia_2022/results/invitro_cooperation.png', units='in', width=5, height=4, res=300)
draw.quad.venn(area1=length(crossfed_1), 
               area2=length(crossfed_2), 
               area3=length(crossfed_3), 
               area4=length(crossfed_4), 
               n12=length(crossfed_12), 
               n23=length(crossfed_23), 
               n13=length(crossfed_13), 
               n14=length(crossfed_14),
               n24=length(crossfed_24),
               n34=length(crossfed_34),
               n123=length(crossfed_123), 
               n124=length(crossfed_124), 
               n234=length(crossfed_234), 
               n134=length(crossfed_134), 
               n1234=length(crossfed_1234), 
               category=c('B. longum','S. thermophilus','R. intestinalis','E. coli'), 
               cat.fontfamily='sans', cat.fontface='italic', fontfamily='sans', cat.cex=1,
               fill=c('springgreen2','tan1','thistle2','paleturquoise3'), lwd=0.0001)
dev.off()

