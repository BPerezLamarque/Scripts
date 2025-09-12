######  Before running this script ######  
#### The following packages need to be installed:

install.packages("ggplot2")
install.packages("vegan")
install.packages("nlme")
install.packages("ape")
install.packages("reshape2")
install.packages("igraph")


######  Step 1: Open and filter the data #####

rm(list=ls())

# Set your working directory (! You need to change the path)
setwd("/Users/bperez/Nextcloud/Recherche/Toulouse/Enseignement/Master_TULIP-GS/2025/Project/data/")

# Open the metadata
metadata <- read.table("metadata_project_La_Reunion.csv", sep=";", header = TRUE, dec=",")

# Open the OTU table of 18S metabarcoding 
OTU_table <- read.table("OTU_table_La_Reunion_18S_OTU97.csv", sep=";", header = TRUE, dec=",")

View(OTU_table) # OTUs are on rows and samples on columns 

OTU_infos <- colnames(OTU_table)[1:8] # the 8 first columns give information on the OTUs


# Select only the samples from Grand Brûlé 
list_samples <- metadata$SampleID[which(metadata$Site=="Grand brule")]


# Extract in the full OTU table
OTU_table_filtered <- OTU_table[,c(OTU_infos, list_samples)]


write.table(OTU_table_filtered, "OTU_table_Grand_brule_18S_OTU97.csv", sep=";", row.names=FALSE, quote=FALSE)


## You can imagine to only keep OTU belonging to Glomeromycotina (to study plant-Glomeromycotina interactions) 
OTU_table_filtered <- OTU_table[grep("Glomeromycotina", OTU_table$taxonomy),]

## If you want to only keep OTU belonging to Glomeromycotina + Mucoromycotina
OTU_table_filtered <- OTU_table[grep("Glomeromycotina|Mucoromycotina", OTU_table$taxonomy),]




######  Step 2: Rarefaction curves and conversion in relative abundances #####

rm(list=ls())

library(ggplot2)


setwd("/Users/bperez/Nextcloud/Recherche/Toulouse/Enseignement/Master_TULIP-GS/2025/Project/data/")


OTU_table <- read.table("OTU_table_Grand_brule_18S_OTU97.csv", sep=";", header=TRUE)


OTU_infos <- colnames(OTU_table)[1:8] # the 8 first columns give information on the OTUs

list_samples <- colnames(OTU_table)[grep("L", colnames(OTU_table))] # get the list of samples



### 1) Plot rarefaction curves at the sample level :

sort(colSums(OTU_table[,list_samples]))

# remove the samples with less than 1000 reads 
list_samples <- list_samples[colSums(OTU_table[,list_samples])>=1000]
OTU_table <- OTU_table[,c(OTU_infos, list_samples)]

max(colSums(OTU_table[,list_samples])) # maximum number of reads

table_raref <- c()    # mean OTU richness

list_OTUs <- as.character(OTU_table$OTU)

# Make a loop to perform the rarefaction every 2500 reads
for (sample in list_samples){
  
  for (nb_reads in seq(1000, 30000, 2500) ){
    
    if (sum(OTU_table[ ,sample])>nb_reads){ # only perform a rarefaction at nb_reads if the sample has more than nb_reads reads
      
      number <- 0
      
      # Subsample 10 times nb_reads reads in each sample 
      for (j in 1:10){
        list_sampled_OTUs <- sample(size= nb_reads, x = list_OTUs, prob = OTU_table[ ,sample], replace=TRUE)
        number <- number + length(unique(list_sampled_OTUs))
      }
      
      number <- number/10 # average number of OTU in the sample when subsampling at nb_reads reads
      
      table_raref <- rbind(table_raref, c(sample, nb_reads, number))
    }
    
  }
}

# Convert the table into a data frame
table_raref <- data.frame(table_raref)
colnames(table_raref) <- c("sample", "number_reads","OTU_richness")
table_raref$number_reads <- as.numeric(as.character(table_raref$number_reads))
table_raref$OTU_richness <- as.numeric(as.character(table_raref$OTU_richness))

pdf(paste0("Rarefaction_curve_fungi_Grand_Brule.pdf"), width=7, height= 4)
ggplot(table_raref, aes(x=number_reads, y=OTU_richness, color=sample)) + 
        geom_line(linewidth=1.2) + theme_classic() + ylab("OTU richness") + 
        xlab("Number of reads")
dev.off()




### 2) Convert the number of reads into relative abundances and filter possible contaminants

# Consider that on OTU present in less than 3 reads is a cross-contamination


for (sample in list_samples){
  
  OTU_table[,sample][OTU_table[,sample]<3] <- 0
  
  # Compute relative abundance
  OTU_table[,sample] <- OTU_table[,sample]/sum(OTU_table[,sample])
  
  # Apply a minimum threshold to consider an OTU as truly present (e.g. 0.1%)
  OTU_table[which(OTU_table[,sample]<0.001), sample] <- 0
  
  # Compute relative abundance again
  OTU_table[,sample] <- OTU_table[,sample]/sum(OTU_table[,sample])
  
}

# Keep only OTUs that are still present

OTU_table <- OTU_table[rowSums(OTU_table[list_samples])>0,]


# Write the final OTU table for diversity analyses

write.table(OTU_table, "OTU_table_fungi_Grand_brule_final.csv", sep=";", row.names=FALSE, quote=FALSE)



######  Step 3: Taxonomic profiles #####

rm(list=ls())

library(ggplot2)
library(phyloseq)
library(scales)


setwd("/Users/bperez/Nextcloud/Recherche/Toulouse/Enseignement/Master_TULIP-GS/2025/Project/data/")

OTU_table <- read.table("OTU_table_fungi_Grand_brule_final.csv", sep=";", header=TRUE)

metadata <- read.table("metadata_project_La_Reunion.csv", sep=";", header = TRUE, dec=",")

# Order the samples per plant species

for (i in 1:ncol(OTU_table)){
  
  if (colnames(OTU_table)[i] %in% metadata$SampleID){
    colnames(OTU_table)[i] <- paste0(metadata$Plant.species[metadata$SampleID==colnames(OTU_table)[i]], "_",colnames(OTU_table)[i])
  }
  
}


list_samples <- colnames(OTU_table)[grep("L", colnames(OTU_table))]



# 1) Prepare the taxonomy

OTU_table$phylum <- OTU_table$taxonomy

OTU_table$phylum <- sapply(strsplit(split="|Fungi|", OTU_table$phylum, fixed=TRUE), "[[", 2)

# This step requires manual cleaning of the different order... 

OTU_table$phylum[grep("|Agaricomycotina|", OTU_table$phylum, fixed=TRUE)] <- "Agaricomycotina"
OTU_table$phylum[grep("|Mucoromycotina|", OTU_table$phylum, fixed=TRUE)] <- "Mucoromycotina"
OTU_table$phylum[grep("|Glomeromycotina|", OTU_table$phylum, fixed=TRUE)] <- "Glomeromycotina"
OTU_table$phylum[grep("|Taphrinomycotina|", OTU_table$phylum, fixed=TRUE)] <- "Taphrinomycotina"
OTU_table$phylum[grep("|Pucciniomycotina|", OTU_table$phylum, fixed=TRUE)] <- "Pucciniomycotina"
OTU_table$phylum[grep("Chytridiomycota|Incertae_Sedis|", OTU_table$phylum, fixed=TRUE)] <- "Chytridiomycotina"
OTU_table$phylum[grep("|Mortierellomycotina|", OTU_table$phylum, fixed=TRUE)] <- "Mortierellomycotina"
OTU_table$phylum[grep("|Saccharomycotina|", OTU_table$phylum, fixed=TRUE)] <- "Saccharomycotina"
OTU_table$phylum[grep("Cryptomycota|LKM11", OTU_table$phylum, fixed=TRUE)] <- "Cryptomycota_LKM11"
OTU_table$phylum[grep("|Pezizomycotina|", OTU_table$phylum, fixed=TRUE)] <- "Pezizomycotina"
OTU_table$phylum[grep("|Ustilaginomycotina|", OTU_table$phylum, fixed=TRUE)] <- "Ustilaginomycotina"
OTU_table$phylum[grep("|Entomophthoromycotina|", OTU_table$phylum, fixed=TRUE)] <- "Entomophthoromycotina"
OTU_table$phylum[grep("Neocallimastigomycota|", OTU_table$phylum, fixed=TRUE)] <- "Neocallimastigomycota"


table(OTU_table$phylum)


# Option: Do the same for the classes or orders... 


# 2) Plot the bar charts

taxonomy <- as.matrix(data.frame(taxa_OTU=OTU_table$OTU, Phylum=OTU_table$phylum))

otu_table_all = phyloseq(otu_table(OTU_table[,list_samples], taxa_are_rows = TRUE), tax_table(taxonomy))


pdf(paste0("Barplots_fungal_communities_Grand_brule.pdf"), width=15, height=7)

p <- plot_bar(otu_table_all,fill="Phylum") + 
  xlab("") + ylab("Relative abundance")+ scale_y_continuous(labels = percent_format())+
  theme_minimal() + theme(legend.title = element_blank()) +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme(panel.grid.major.y = element_blank())+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(p)

dev.off()





######  Step 4: Alpha diversity #####

rm(list=ls())

library(ggplot2)
library(vegan)
library(nlme)


setwd("/Users/bperez/Nextcloud/Recherche/Toulouse/Enseignement/Master_TULIP-GS/2025/Project/data/")

OTU_table <- read.table("OTU_table_fungi_Grand_brule_final.csv", sep=";", header=TRUE)

list_samples <- colnames(OTU_table)[grep("L", colnames(OTU_table))]

metadata <- read.table("metadata_project_La_Reunion.csv", sep=";", header = TRUE, dec=",")


# 1) Compute OTU richness and Shannon index for each sample

alpha_div <- c()

for (sample in list_samples){
  
  abundances <- OTU_table[,sample]
  alpha_div <- rbind(alpha_div, c(sample, 
                                  metadata$Plant.species[which(metadata$SampleID==sample)], # which plant species
                                  metadata$origin[which(metadata$SampleID==sample)],  # which origin (invasive, exotic or native)
                                  length(which(abundances>0)), # species richness
                                  diversity(abundances, index = "shannon") ) ) # Shannon index
  
}

alpha_div <- data.frame(alpha_div)
colnames(alpha_div) <- c("sample", "species", "origin", "OTU_richness", "Shannon")
alpha_div$OTU_richness <- as.numeric(as.character(alpha_div$OTU_richness))
alpha_div$Shannon <- as.numeric(as.character(alpha_div$Shannon))


# Plot the alpha diversities 

pdf(paste0("Alpha_div_fungal_communities_Grand_brule.pdf"), width=6.5, height=4)

ggplot(alpha_div, aes(x= species, y=OTU_richness, col = origin)) + geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Plant species") + ylab("OTU richness")
  
ggplot(alpha_div, aes(x= species, y=Shannon, col = origin)) + geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Plant species") + ylab("Shannon index")

dev.off()


# 2) Test whether invasive plants tend to associate with less fungal species than native ones


# test the effect of plant species on alpha diversity

kruskal.test(alpha_div$OTU_richness ~ alpha_div$species)
kruskal.test(alpha_div$Shannon ~ alpha_div$species)


# Use linear models and mixed linear models 

summary(lm(alpha_div$Shannon ~ alpha_div$origin ))

summary(lme(Shannon ~ origin, data = alpha_div, random = ~ 1|species  ))

summary(lm(alpha_div$OTU_richness ~ alpha_div$origin ))
summary(lme(OTU_richness ~ origin, data = alpha_div, random = ~ 1|species  ))




######  Step 5: Beta diversity #####

rm(list=ls())

library(ggplot2)
library(vegan)
library(nlme)
library(ape)



setwd("/Users/bperez/Nextcloud/Recherche/Toulouse/Enseignement/Master_TULIP-GS/2025/Project/data/")

OTU_table <- read.table("OTU_table_fungi_Grand_brule_final.csv", sep=";", header=TRUE)

list_samples <- colnames(OTU_table)[grep("L", colnames(OTU_table))]

metadata <- read.table("metadata_project_La_Reunion.csv", sep=";", header = TRUE, dec=",")


# 1) Compute beta diversities using Bray-Curtis distances and perform a principal coordinate analysis:

# Compute beta diversities using Bray-Curtis distances

dist_bray <- vegdist(t(OTU_table[,list_samples]), "bray") 

# Perform PCoA

res <- pcoa(dist_bray, correction="none", rn=NULL)

origin=metadata$origin[match(list_samples,metadata$SampleID)]
species=metadata$Plant.species[match(list_samples,metadata$SampleID)]

data_pcoa <- data.frame(res$vectors, origin, species)


# Plot PCoA

pdf(paste0("PCoA_fungal_communities_La_Reunion.pdf"), width=6.5, height=4)

# colors based on plant species
ggplot(data_pcoa,aes(x=Axis.1,y=Axis.2, color=species, fill=species))+xlab(paste0("PCoA1 (", round(res$values[,3][1]*100, 1),"%)"))+ylab(paste0("PCoA2 (", round(res$values[,3][2]*100, 1),"%)")) +labs(title=" ")+ theme_bw() + 
  geom_point(alpha=0.7,size=3)

# colors based on origin
ggplot(data_pcoa,aes(x=Axis.1,y=Axis.2, color=origin, fill=origin))+xlab(paste0("PCoA1 (", round(res$values[,3][1]*100, 1),"%)"))+ylab(paste0("PCoA2 (", round(res$values[,3][2]*100, 1),"%)")) +labs(title=" ")+ theme_bw() + 
  geom_point(alpha=0.7,size=3)

dev.off()


# 2) Test whether samples belonging to the same plant species tend to have similar composition 

# Test the dispersion of the beta diversity
mod2 <- betadisper(dist_bray, species)
print(anova(mod2))

mod2 <- betadisper(dist_bray, origin)
print(anova(mod2))


# PermANOVA

permanova <- adonis2(dist_bray ~ species, permutations=10000) 
permanova

permanova <- adonis2(dist_bray ~ origin, permutations=10000) 
permanova

permanova <- adonis2(dist_bray ~ species + origin, permutations=10000, by = "terms") 
permanova



######  Step 6: Network representation #####

rm(list=ls())

library(reshape2)
library(igraph)


setwd("/Users/bperez/Nextcloud/Recherche/Toulouse/Enseignement/Master_TULIP-GS/2025/Project/data/")

OTU_table <- read.table("OTU_table_fungi_Grand_brule_final.csv", sep=";", header=TRUE)

list_samples <- colnames(OTU_table)[grep("L", colnames(OTU_table))]

metadata <- read.table("metadata_project_La_Reunion.csv", sep=";", header = TRUE, dec=",")

# Format the OTU table for igraph

melt_network <- melt(as.matrix(OTU_table[,list_samples]))
melt_network <- melt_network[which(melt_network$value>0),]

g <- graph_from_edgelist(as.matrix(melt_network[,1:2]), directed = FALSE)
edge.attributes(g)$weight <- melt_network[,3]

V(g)$type <- rep("fungi", length(V(g)))
V(g)$type[which(names(V(g)) %in% metadata$SampleID[metadata$origin=="exotic"])] <-  "exotic"
V(g)$type[which(names(V(g)) %in% metadata$SampleID[metadata$origin=="native"])] <-  "native"
V(g)$type[which(names(V(g)) %in% metadata$SampleID[metadata$origin=="invasive"])] <-  "invasive"



# Define the colors, sizes and shapes of each category
shape <- c("circle","square","circle", "circle")
size <- c(3, 1, 3, 3)
col <- c( "#ba4a00", "#616a6b", "#3498db", "#2ecc71")
names(size) <- names(shape) <- names(col) <- c("exotic", "fungi", "invasive", "native")


# Plot the network

pdf("Network_fungal_communities_La_Réunion.pdf")
set.seed(3)
plot(g, vertex.color = col[V(g)$type], vertex.shape = shape[V(g)$type], vertex.label=NA, vertex.size=size[V(g)$type],
       edge.width=1+sqrt(100*E(g)$weight),  layout= layout_with_fr(g, weights = sqrt(100*E(g)$weight), niter = 5000 ) )
dev.off()

# You can try several layouts and different ways of computing weights










