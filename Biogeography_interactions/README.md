**Guidelines to measure cophylogenetic signal and phylogenetic congruence in a bipartite interaction network:**
====

<br>


It contains functions to analyze cophylogenetic systems, in particular to measure cophylogenetic signal, phylogenetic congruence, phylogenetic signal, or phylosymbiosis. 

It also allows to generate “null expectations” on the expected patterns based on the species' biogeography, i.e. whether the observed phylogenetic patterns can be explained by biogeographic conservatism alone.


<br>

**Citation:** Ibos Y., Sagnard A., Perez-Lamarque B., Assessing the influence of biogeography on (co)phylogenetic signals in species interactions, 2026, in prep.


**Contact:** benoit.perez-lamarque@utoulouse.fr

<br> 


# Example of R script: 


```r

library(ape)
library(phytools)
library(vegan)

load("functions_cophylo.R")


# Open the trees 
tree_B= read.tree("tree_A.tre")
tree_A = read.tree("tree_B.tre")
    
dist_B <- cophenetic(tree_B)
dist_A <- cophenetic(tree_A)
    

# open the interaction network
network = read.table(paste0("network_",name,"_",seed,".csv"),sep=";",header= TRUE)
# species on clade A are on columns, and species on clade B on rows
    

# Indicate the biogeography of clade A (columns) in the following table 
table_biogeo = data.frame(species= colnames(network),biogeo= bgeo)
    
nperm <- 10000

    
## Run PACo
D <-prepare_paco_data_test(dist_A, dist_B, network)
PACo_result <- PACo_test(D, nperm = nperm, seed = NA, nullmodel = 2, symmetric = TRUE)
PACo_result_biogeo <- PACo_test(D, nperm = nperm, seed = NA, nullmodel = 3, symmetric = TRUE, table_biogeo=table_biogeo) # when constraining the biogeography
    
# Parafit 
parafit_result <- parafit_test(dist_A, dist_B, network, nperm = nperm, nullmodel = 2)
parafit_result_biogeo <- parafit_test(dist_A, dist_B, network, nperm = nperm, nullmodel = 3, table_biogeo = table_biogeo) # when constraining the biogeography
    
# Mantel 
mantel_result <- phylosignal_network(network, tree_A, method="Jaccard_binary", nperm = nperm, correlation = "Pearson", only_A = TRUE, permutation ="shuffle")
mantel_result_biogeo <- phylosignal_network(network, tree_A, method="Jaccard_binary", nperm = nperm, correlation = "Pearson", only_A = TRUE, permutation ="biogeo", table_biogeo=table_biogeo) # when constraining the biogeography
    
# Summarize all the results
res <- c(name, seed, round(PACo_result$R2_obs,5), PACo_result$p, PACo_result_biogeo$p, 
             round(parafit_result$ParaFitGlobal,0), parafit_result$p, parafit_result_biogeo$p, 
             round(mantel_result$mantel_cor_A,4), mantel_result$p_upper_A, mantel_result_biogeo$p_upper_A)
names(res) <- c("name", "seed", "R2_paco", "p_paco", "p_paco_biogeo", "R_parafit", "p_parafit", "p_parafit_biogeo", "R_mantel", "p_mantel", "p_mantel_biogeo") 

print(res)
    
# Plot the null expectations in the patterns based on the biogeography
plot_H0(PACo_result, PACo_result_biogeo, "paco")
plot_H0(parafit_result, parafit_result_biogeo, "parafit")
plot_H0(mantel_result, mantel_result_biogeo, "mantel")

```