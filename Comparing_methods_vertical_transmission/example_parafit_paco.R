#### This script contains a quick tutorial on how to use the amended functions of ParaFit and PACo originally from the R-packages ape (Paradis et al., 2004) and paco (Hutchinson et al., 2017))
#### The amended functions now better handle the cases where the number of microbial strains is low (<4), which can generate bugs in the original functions during the permutations. 


#### Contact: Benoît Perez-Lamarque (benoit.perez.lamarque@gmail.com)

#### For more details about this script, see "Comparing different computational approaches for detecting long-term vertical transmission in host-associated microbiota", Perez-Lamarque B & Morlon H, in prep.


### Running ParaFit ####

library(ape)

source("functions_parafit_paco.R")

data(gopher.D)
data(lice.D)
data(HP.links)

# permutations with the null model 1
res <- parafit_test(gopher.D, lice.D, HP.links, nperm=10000, test.links=FALSE, nullmodel = 1)
res$ParaFitGlobal
res$p.global

# permutations with the null model 2
res <- parafit_test(gopher.D, lice.D, HP.links, nperm=10000, test.links=FALSE, nullmodel = 2)
res$ParaFitGlobal
res$p.global



### Running PACo ####

library(paco)

source("functions_parafit_paco.R")

data(gopherlice)
require(ape)

gdist <- cophenetic(gophertree)
ldist <- cophenetic(licetree)

D <- prepare_paco_data(gdist, ldist, gl_links)
D <- add_pcoord_test(D)


# permutations with the null model 1 (the same as method="r0" in the original function paco)
D <- PACo_test(D, nperm=10000, seed=42, nullmodel=1)
print(D$gof)

# permutations with the null model 2
D <- PACo_test(D, nperm=10000, seed=42, nullmodel=2)
print(D$gof)


