##### Step 1A: Make simulations trait matching #####

rm(list=ls())

library(ape)
library(phytools)
library(mvMORPH)


setwd("/path/simulated_data/")

name='simul_1'
nb_rep <- 250
seed=1

summary_stats <- c()


connectance <- c()

for (name in c('simul_1', 'simul_2', 'simul_3', 'simul_4',  'simul_5', 'simul_6', 'simul_7', 'simul_8')){
  
  print(name)
  
  if (name %in% c('simul_1', 'simul_2', 'simul_3', 'simul_4')){
    lambda_poisson=1.5
  }
  if (name %in% c('simul_5', 'simul_6', 'simul_7', 'simul_8')){
    lambda_poisson=1
  }
  
  if (name %in% c("simul_1","simul_5")){
    min_species <- 10
    max_species <- 50
  }
  if (name %in% c("simul_2","simul_6")){
    min_species <- 51
    max_species <- 100
  }
  if (name %in% c("simul_3","simul_7")){
    min_species <- 101
    max_species <- 150
  }
  if (name %in% c("simul_4","simul_8")){
    min_species <- 151
    max_species <- 200
  }
  
  for (seed in 1:nb_rep){
    
    set.seed(seed)
    
      nb_hosts <- floor(runif(1, min_species, max_species))
      nb_parasites <- floor(runif(1, min_species, max_species))
      
      tree_hosts <- NULL # in case the birth-death model leads to full extinctions
      tree_parasites <- NULL
      
      while (is.null(tree_hosts)) tree_hosts <- pbtree(n=nb_hosts, b = 0.1, d=0.003, extant.only=T)
      tree_hosts$tip.label <- gsub("t", "h", tree_hosts$tip.label)
      while (is.null(tree_parasites)) tree_parasites <- pbtree(n=nb_parasites, b = 0.3, d=0.003, extant.only=T)
      tree_parasites$tip.label <- gsub("t", "p", tree_parasites$tip.label)
      
      
      traits_hosts <- mvSIM(tree = tree_hosts, nsim=1, model="BM1", param = list(ntraits=1, sigma=1, theta=10))
      traits_parasites <- mvSIM(tree = tree_parasites, nsim=1, model="BM1", param = list(ntraits=1, sigma=1, theta=10))
      
      
      # Hosts are on rows and parasites on columns
      network <- matrix(0, nrow=Ntip(tree_hosts), ncol=Ntip(tree_parasites))
      rownames(network) <- tree_hosts$tip.label
      colnames(network) <- tree_parasites$tip.label
      
      
      for (i in 1:Ntip(tree_parasites)){
        
        nb_hosts <- rpois(n=1, lambda=lambda_poisson)
        
        hosts <- sample(1:Ntip(tree_hosts), size=nb_hosts, prob=1/abs(traits_hosts-traits_parasites[i]))
        
        network[hosts, i] <- 1
        
      }
      
      write.table(network, paste0("network_", name, "_", seed, ".csv"), sep=";")
      write.tree(tree_hosts, paste0("tree_hosts_", name, "_", seed, ".tre"))
      write.tree(tree_parasites, paste0("tree_parasites_", name, "_", seed, ".tre"))
      
      # Prepare for empress
      
      list_links <- reshape2::melt(network)
      list_links <- list_links[list_links$value>0,]
      write.table(paste0(list_links$Var2,":", list_links$Var1), paste0("links_", name, "_", seed, ".txt"), col.names=F, row.names=F, quote=F)
      
      network <- network[rowSums(network)>0,]
      network <- network[,colSums(network)>0]
      connectance <- c(connectance, length(which(network>0))/nrow(network)/ncol(network))
      
      # remove the tips with no interactions
      tree_hosts_interact <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[which(!tree_hosts$tip.label %in% list_links$Var1)])
      write.tree(tree_hosts_interact, paste0("tree_hosts_only_int_", name, "_", seed, ".tre"))
      tree_parasites_interact <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[which(!tree_parasites$tip.label %in% list_links$Var2)])
      write.tree(tree_parasites_interact, paste0("tree_parasites_only_int_", name, "_", seed, ".tre"))
      
      summary_stats <- rbind(summary_stats, c(Ntip(tree_hosts), nb_parasites, Ntip(tree_hosts_interact), Ntip(tree_parasites_interact), max(node.depth.edgelength(tree_hosts)), max(node.depth.edgelength(tree_parasites))))
      
      
      # save the trait values (for plotting - for simul_1 only)
      if (name=="simul_1"){
        write.table(traits_hosts, paste0("traits_hosts_", name, "_", seed, ".csv"), sep=";")
        write.table(traits_parasites, paste0("traits_parasites_", name, "_", seed, ".csv"), sep=";")
      }
    
  }
}

hist(connectance, breaks=100)

summary_stats <- data.frame(summary_stats)
colnames(summary_stats) <- c("nb_hosts", "nb_parasites", "nb_hosts_interact", "nb_parasites_interact","age_hosts", "age_parasites")

pdf(paste0("../summary_stats_simulations_lambda_",lambda_poisson,".pdf"), width=3.5, height=3)
hist(summary_stats$nb_hosts, ylab = "Count", xlab="Number of host species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))
hist(summary_stats$nb_parasites, ylab = "Count", xlab="Number of parasite species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))

hist(summary_stats$nb_hosts_interact, ylab = "Count", xlab="Number of interacting host species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))
hist(summary_stats$nb_parasites_interact, ylab = "Count", xlab="Number of interacting parasite species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))

hist(summary_stats$age_hosts, ylab = "Count", xlab="Age of the host clade (in Myr)", main="", col="antiquewhite", breaks=20, xlim=c(0,95))
hist(summary_stats$age_parasites, ylab = "Count", xlab="Age of the parasite clade (in Myr)", main="", col="antiquewhite", breaks=20, xlim=c(0,95))
dev.off()





##### Step 1B: Make simulations vicariance #####

rm(list=ls())

library(ape)
library(phytools)
library(TreeSim)
library(mvMORPH)


setwd("/path/simulated_data_geo/")


# Three areas 
# from 5 Myr vicariance to 20 Myr 

# Number of species to start with: 10 to 20


name='simul_geo_1'
nb_rep <- 250
seed=1

for (lambda_poisson in c(1.5, 1)){
  
  summary_stats <- c()
  
  for (name in c('simul_geo_1', 'simul_geo_2', 'simul_geo_3', 'simul_geo_4', 'simul_geo_5', 'simul_geo_6', 'simul_geo_7', 'simul_geo_8')){
  
  if (lambda_poisson==1.5){
    list_names=c('simul_geo_1', 'simul_geo_2', 'simul_geo_3', 'simul_geo_4')
  }
  if (lambda_poisson==1){
    list_names=c('simul_geo_5', 'simul_geo_6', 'simul_geo_7', 'simul_geo_8')
  }
  
  for (name in list_names){
    print(name)
    
    
    nb_min=10
    nb_max=20
    
    if (name %in% c("simul_geo_1","simul_geo_5")){
      age=5
    }
    if (name %in% c("simul_geo_2","simul_geo_6")){
      age=10
    }
    if (name %in% c("simul_geo_3","simul_geo_7")){
      age=15
    }
    if (name %in% c("simul_geo_4","simul_geo_8")){
      age=20
    }
    
    for (seed in 1:nb_rep){
      
      set.seed(seed)
      
      nb_hosts_before <- sample(seq(nb_min, nb_max), size=1)
      nb_parasites_before <- sample(seq(nb_min, nb_max), size=1)
      
      tree_hosts_before <- pbtree(n=nb_hosts_before, b = 0.1, d=0.003, extant.only=T)
      tree_parasites_before <- pbtree(n=nb_parasites_before, b = 0.1, d=0.003, extant.only=T)
      
      # Attribute randomly each species to a geographic area
      area_hosts <- sample(c("A","B","C"), size=nb_hosts_before, replace = TRUE)
      names(area_hosts) <- tree_hosts_before$tip.label
      area_parasites <- sample(c("A","B","C"), size=nb_parasites_before, replace = TRUE)
      names(area_parasites) <- tree_parasites_before$tip.label
      
      # simulate speciation at vicariance
      tree_hosts_vicariance <- tree_hosts_before
      for (i in 1:Ntip(tree_hosts_before)){
        if (runif(1)<0.5){ # 50% of chances of experiencing a speciation at the vicariance event
          ntip <- Ntip(tree_hosts_vicariance)+1
          tree_hosts_vicariance <- bind.tip(tree_hosts_vicariance, tip.label=paste0("t", ntip), edge.length=0, where=which(tree_hosts_vicariance$tip.label==tree_hosts_before$tip.label[i]), position=0)
          area_hosts <- c(area_hosts, sample(c("A","B","C")[which(c("A","B","C")!=area_hosts[tree_hosts_before$tip.label[i]])], size=1))
          names(area_hosts)[ntip] <- paste0("t", ntip)
        }
      }
      area_hosts <- area_hosts[tree_hosts_vicariance$tip.label]
      
      tree_parasites_vicariance <- tree_parasites_before
      for (i in 1:Ntip(tree_parasites_before)){
        if (runif(1)<0.5){ # 50% of chances of experiencing a speciation at the vicariance event
          ntip <- Ntip(tree_parasites_vicariance)+1
          tree_parasites_vicariance <- bind.tip(tree_parasites_vicariance, tip.label=paste0("t", ntip), edge.length=0, where=which(tree_parasites_vicariance$tip.label==tree_parasites_before$tip.label[i]), position=0)
          area_parasites <- c(area_parasites, sample(c("A","B","C")[which(c("A","B","C")!=area_parasites[tree_parasites_before$tip.label[i]])], size=1))
          names(area_parasites)[ntip] <- paste0("t", ntip)
        }
      }
      area_parasites <- area_parasites[tree_parasites_vicariance$tip.label]
      
      # simulate hosts lineages after vicariance
      tree_hosts <- tree_hosts_vicariance
      for (i in 1:length(area_hosts)){
        tip <- names(area_hosts)[i]
        area <- area_hosts[i]
        
        nb_species_area <- 0
        while (nb_species_area<1){
          tryCatch({tree_species_area <- sim.bd.age(age=age, numbsim=1, lambda = 0.1, mu = 0.003, mrca=FALSE,complete=FALSE)[[1]]}, error=function(e){tree_species_area=0})
          if (is.phylo(tree_species_area)){nb_species_area <- Ntip(tree_species_area)}
        }
        nb_species_total_area <- length(grep(area, tree_hosts$tip.label))
        tree_species_area$tip.label <- paste0("h", area, (nb_species_total_area+1):(nb_species_total_area+Ntip(tree_species_area)))
        
        tree_hosts <- bind.tree(tree_hosts, tree_species_area, where = which(tree_hosts$tip.label==tip), position = 0)
        
      }
      
      # simulate parasite lineages after vicariance
      tree_parasites <- tree_parasites_vicariance
      for (i in 1:length(area_parasites)){
        tip <- names(area_parasites)[i]
        area <- area_parasites[i]
        
        nb_species_area <- 0
        while (nb_species_area<1){
          tryCatch({tree_species_area <- sim.bd.age(age=age, numbsim=1, lambda = 0.1, mu = 0.003, mrca=FALSE,complete=FALSE)[[1]]}, error=function(e){tree_species_area=0})
          if (is.phylo(tree_species_area)){nb_species_area <- Ntip(tree_species_area)}
        }
        nb_species_total_area <- length(grep(area, tree_parasites$tip.label))
        tree_species_area$tip.label <- paste0("p", area, (nb_species_total_area+1):(nb_species_total_area+Ntip(tree_species_area)))
        
        tree_parasites <- bind.tree(tree_parasites, tree_species_area, where = which(tree_parasites$tip.label==tip), position = 0)
        
      }

      
      # Hosts are on rows and parasites on columns
      network <- matrix(0, nrow=Ntip(tree_hosts), ncol=Ntip(tree_parasites))
      rownames(network) <- tree_hosts$tip.label
      colnames(network) <- tree_parasites$tip.label
      
      # Interactions at random within each region
      for (i in 1:Ntip(tree_parasites)){
        
        if (length(grep("pA", tree_parasites$tip.label[i]))==1){
          nb_hosts <- rpois(n=1, lambda=lambda_poisson)
          list_hosts <- tree_hosts$tip.label[grep("hA", tree_hosts$tip.label)]
          if (length(list_hosts)>0){
            hosts <- unique(sample(list_hosts, size=nb_hosts, replace=TRUE))
            network[hosts, i] <- 1
          }
        }
        if (length(grep("pB", tree_parasites$tip.label[i]))==1){
          nb_hosts <- rpois(n=1, lambda=lambda_poisson)
          list_hosts <- tree_hosts$tip.label[grep("hB", tree_hosts$tip.label)]
          if (length(list_hosts)>0){
            hosts <- unique(sample(list_hosts, size=nb_hosts, replace=TRUE))
            network[hosts, i] <- 1
          }
        }
        if (length(grep("pC", tree_parasites$tip.label[i]))==1){
          nb_hosts <- rpois(n=1, lambda=lambda_poisson)
          list_hosts <- tree_hosts$tip.label[grep("hC", tree_hosts$tip.label)]
          if (length(list_hosts)>0){
            hosts <- unique(sample(list_hosts, size=nb_hosts, replace=TRUE))
            network[hosts, i] <- 1
          }
        }
      }
      
      write.table(network, paste0("network_", name, "_", seed, ".csv"), sep=";")
      write.tree(tree_hosts, paste0("tree_hosts_", name, "_", seed, ".tre"))
      write.tree(tree_parasites, paste0("tree_parasites_", name, "_", seed, ".tre"))
      
      # Prepare for empress
      
      list_links <- reshape2::melt(network)
      list_links <- list_links[list_links$value>0,]
      write.table(paste0(list_links$Var2,":", list_links$Var1), paste0("links_", name, "_", seed, ".txt"), col.names=F, row.names=F, quote=F)
      
      # remove the tips with no interactions
      tree_hosts_interact <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[which(!tree_hosts$tip.label %in% list_links$Var1)])
      write.tree(tree_hosts_interact, paste0("tree_hosts_only_int_", name, "_", seed, ".tre"))
      tree_parasites_interact <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[which(!tree_parasites$tip.label %in% list_links$Var2)])
      write.tree(tree_parasites_interact, paste0("tree_parasites_only_int_", name, "_", seed, ".tre"))
      
      summary_stats <- rbind(summary_stats, c(Ntip(tree_hosts), Ntip(tree_parasites), Ntip(tree_hosts_interact), Ntip(tree_parasites_interact), max(node.depth.edgelength(tree_hosts)), max(node.depth.edgelength(tree_parasites))))
      
    }
  }
  
  }
  
  summary_stats <- data.frame(summary_stats)
  colnames(summary_stats) <- c("nb_hosts", "nb_parasites", "nb_hosts_interact", "nb_parasites_interact","age_hosts", "age_parasites")
  
  max(summary_stats$nb_hosts)
  max(summary_stats$nb_parasites)
  
  max(summary_stats$nb_hosts_interact)
  max(summary_stats$nb_parasites_interact)
  
  min(summary_stats$nb_hosts_interact)
  min(summary_stats$nb_parasites_interact)
  
  pdf(paste0("../summary_stats_simulations_geo_lambda_",lambda_poisson,".pdf"), width=3.5, height=3)
  hist(summary_stats$nb_hosts, ylab = "Count", xlab="Number of host species", main="", col="antiquewhite", breaks=25, xlim=c(0,250))
  hist(summary_stats$nb_parasites, ylab = "Count", xlab="Number of parasite species", main="", col="antiquewhite", breaks=25, xlim=c(0,250))
  
  hist(summary_stats$nb_hosts_interact, ylab = "Count", xlab="Number of interacting host species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))
  hist(summary_stats$nb_parasites_interact, ylab = "Count", xlab="Number of interacting parasite species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))
  
  hist(summary_stats$age_hosts, ylab = "Count", xlab="Age of the host clade (in Myr)", main="", col="antiquewhite", breaks=20, xlim=c(0,95))
  hist(summary_stats$age_parasites, ylab = "Count", xlab="Age of the parasite clade (in Myr)", main="", col="antiquewhite", breaks=20, xlim=c(0,95))
  dev.off()
  
}


##### Step 1C: Make simulations co-diversification #####

rm(list=ls())

library(ape)
library(phytools)
library(HOME)


setwd("/path/simulated_data_codiv/")

name='simul_codiv_1'
nb_rep <- 250
seed=1

summary_stats <- c()


for (name in c('simul_codiv_1', 'simul_codiv_2', 'simul_codiv_3', 'simul_codiv_4')){
  
  print(name)
  
  if (name=="simul_codiv_1"){
    min_species <- 10
    max_species <- 50
  }
  if (name=="simul_codiv_2"){
    min_species <- 51
    max_species <- 100
  }
  if (name=="simul_codiv_3"){
    min_species <- 101
    max_species <- 150
  }
  if (name=="simul_codiv_4"){
    min_species <- 151
    max_species <- 200
  }
  
  for (seed in 1:nb_rep){
    figure=FALSE
  
    set.seed(seed)
    
    nb_hosts <- floor(runif(1, min_species, max_species))
 	  nb_transfers <- floor(runif(1, 0, nb_hosts/2))
    
    tree_hosts <- pbtree(n=nb_hosts, b = 0.1, d=0.003, extant.only=T)
    tree_hosts$tip.label <- gsub("t", "h", tree_hosts$tip.label)
  
    sim_microbiota(name, name_index=seed, provided_tree=tree_hosts, simul=nb_transfers, delta=0.001, figure=figure, scale=FALSE)
  
    
    tree_parasites <- read.tree(paste0("simulations/symbiont_tree_",name,"_",seed,".tre"))
    tree_parasites$tip.label <- gsub("h", "p", tree_parasites$tip.label)
    
    
    # rate of 10% of losses at present
	list_drop <- c()
	for (i in 1:Ntip(tree_parasites)){
	if (runif(1)<0.1) list_drop <- c(list_drop, tree_parasites$tip.label[i])
	}
	tree_parasites <- drop.tip(tree_parasites, tip=list_drop)
	
	    
    
    # Hosts are on rows and parasites on columns
    network <- matrix(0, nrow=Ntip(tree_hosts), ncol=Ntip(tree_parasites))
    tree_parasites$tip.label <- gsub("p", "h", tree_parasites$tip.label)
    rownames(network) <- tree_hosts$tip.label
    colnames(network) <- tree_parasites$tip.label
    
    
    for (i in 1:Ntip(tree_parasites)){
      network[strsplit(split="-", tree_parasites$tip.label[i])[[1]][1],i] <- 1
    }
    
    rownames(network) <- tree_hosts$tip.label <- paste0("h",1:Ntip(tree_hosts))
    colnames(network) <- tree_parasites$tip.label <- paste0("p",1:Ntip(tree_parasites))
    
    write.table(network, paste0("network_", name, "_", seed, ".csv"), sep=";")
    write.tree(tree_hosts, paste0("tree_hosts_", name, "_", seed, ".tre"))
    write.tree(tree_parasites, paste0("tree_parasites_", name, "_", seed, ".tre"))
    
    # Prepare for empress
    
    list_links <- reshape2::melt(network)
    list_links <- list_links[list_links$value>0,]
    write.table(paste0(list_links$Var2,":", list_links$Var1), paste0("links_", name, "_", seed, ".txt"), col.names=F, row.names=F, quote=F)
    
    
    # remove the tips with no interactions
    tree_hosts_interact <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[which(!tree_hosts$tip.label %in% list_links$Var1)])
    write.tree(tree_hosts_interact, paste0("tree_hosts_only_int_", name, "_", seed, ".tre"))
    tree_parasites_interact <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[which(!tree_parasites$tip.label %in% list_links$Var2)])
    write.tree(tree_parasites_interact, paste0("tree_parasites_only_int_", name, "_", seed, ".tre"))
    
    summary_stats <- rbind(summary_stats, c(Ntip(tree_hosts), Ntip(tree_parasites), Ntip(tree_hosts_interact), Ntip(tree_parasites_interact), max(node.depth.edgelength(tree_hosts)), max(node.depth.edgelength(tree_parasites))))

    
  }
}

summary_stats <- data.frame(summary_stats)
colnames(summary_stats) <- c("nb_hosts", "nb_parasites", "nb_hosts_interact", "nb_parasites_interact","age_hosts", "age_parasites")

pdf("../summary_stats_simulations_codiv.pdf", width=3.5, height=3)
hist(summary_stats$nb_hosts, ylab = "Count", xlab="Number of host species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))
hist(summary_stats$nb_parasites, ylab = "Count", xlab="Number of parasite species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))

hist(summary_stats$nb_hosts_interact, ylab = "Count", xlab="Number of interacting host species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))
hist(summary_stats$nb_parasites_interact, ylab = "Count", xlab="Number of interacting parasite species", main="", col="antiquewhite", breaks=20, xlim=c(0,200))

hist(summary_stats$age_hosts, ylab = "Count", xlab="Age of the host clade (in Myr)", main="", col="antiquewhite", breaks=20, xlim=c(0,95))
hist(summary_stats$age_parasites, ylab = "Count", xlab="Age of the parasite clade (in Myr)", main="", col="antiquewhite", breaks=20, xlim=c(0,95))
dev.off()




##### Step 1D: Prepare eMPRess with random bifurcations ####

rm(list=ls())

setwd("/path/simulated_data/")

nb_rep=250

name='simul_1'
seed=1

for (name in c('simul_1', 'simul_2', 'simul_3', 'simul_4', 'simul_5', 'simul_6', 'simul_7', 'simul_8')){
  print(name)
  for (seed in 1:nb_rep){
    
    set.seed(seed)
    
    network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";")
    network <- network[which(rowSums(network)>0),]
    network <- network[,which(colSums(network)>0)]
    
    tree_hosts <- read.tree(paste0("tree_hosts_only_int_", name, "_", seed, ".tre"))
    tree_parasites <- read.tree(paste0("tree_parasites_only_int_", name, "_", seed, ".tre"))
    
    position <- min(c(0.001, tree_parasites$edge.length))
    
    while (!all(colSums(network)==1)){
      ind_species <- which(colSums(network)>1)[sample(length(which(colSums(network)>1)), size=1)]
      species <- colnames(network)[ind_species]
      nb_parasite_names <- max(as.numeric(gsub("p","",colnames(network))))+1

      network <- cbind(network, rep(0, nrow(network)))
      colnames(network)[ncol(network)] <- paste0("p",nb_parasite_names)
      list_index <- which(network[,species]>0)
      index <- list_index[sample(length(list_index),size=1)]
      network[index, ncol(network)] <- 1
      network[index, species] <- 0
      
      tree_parasites <- bind.tip(tree_parasites, tip.label=paste0("p",nb_parasite_names), edge.length=0.001, where=which(tree_parasites$tip.label==species), position=position)
      
    }
    
    if (!is.ultrametric(tree_parasites)) tree_parasites <- force.ultrametric(tree_parasites, method="extend")
    
    list_links <- reshape2::melt(as.matrix(network))
    list_links <- list_links[list_links$value>0,]
    
    write.table(paste0(list_links$Var2,":", list_links$Var1), paste0("links_bifurcations_eMPRess_", name, "_", seed, ".txt"), col.names=F, row.names=F, quote=F)
    write.tree(tree_parasites, paste0("tree_parasites_bifurcations_eMPRess_", name, "_", seed, ".tre"))
    
  }
}


rm(list=ls())

setwd("/path/simulated_data_geo/")

nb_rep=250

name='simul_geo_1'
seed=1

for (name in c('simul_geo_1', 'simul_geo_2', 'simul_geo_3', 'simul_geo_4', 'simul_geo_5', 'simul_geo_6', 'simul_geo_7', 'simul_geo_8')){
  print(name)
  for (seed in 1:nb_rep){
    
    set.seed(seed)
    
    network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";")
    network <- network[which(rowSums(network)>0),]
    network <- network[,which(colSums(network)>0)]
    
    tree_hosts <- read.tree(paste0("tree_hosts_only_int_", name, "_", seed, ".tre"))
    tree_parasites <- read.tree(paste0("tree_parasites_only_int_", name, "_", seed, ".tre"))
    
    position <- min(c(0.001, tree_parasites$edge.length))
    
    while (!all(colSums(network)==1)){
      ind_species <- which(colSums(network)>1)[sample(length(which(colSums(network)>1)), size=1)]
      species <- colnames(network)[ind_species]
      area <- substr(species, 2, 2)
      nb_parasite_names <- max(as.numeric(gsub(paste0("p",area),"",colnames(network)[grep(paste0("p",area), colnames(network))])))+1
      
      network <- cbind(network, rep(0, nrow(network)))
      colnames(network)[ncol(network)] <- paste0("p",area,nb_parasite_names)
      list_index <- which(network[,species]>0)
      index <- list_index[sample(length(list_index),size=1)]
      network[index, ncol(network)] <- 1
      network[index, species] <- 0
      
      tree_parasites <- bind.tip(tree_parasites, tip.label=paste0("p",area,nb_parasite_names), edge.length=0.001, where=which(tree_parasites$tip.label==species), position=position)
      
    }
    
    if (!is.ultrametric(tree_parasites)) tree_parasites <- force.ultrametric(tree_parasites, method="extend")
    
    list_links <- reshape2::melt(as.matrix(network))
    list_links <- list_links[list_links$value>0,]
    
    write.table(paste0(list_links$Var2,":", list_links$Var1), paste0("links_bifurcations_eMPRess_", name, "_", seed, ".txt"), col.names=F, row.names=F, quote=F)
    write.tree(tree_parasites, paste0("tree_parasites_bifurcations_eMPRess_", name, "_", seed, ".tre"))
    
  }
}






##### Step 2A : Fit ParaFit and PACo for trait matching #####



rm(list=ls())

library(ape)
library(phytools)
library(mvMORPH)
library(paco)
library(parallel)


#setwd("/path/simulated_data/")
setwd("/path/simulated_data/")

source("../global_fit/functions_parafit.R")

name='simul_1'
seed=1

nb_cores <- 25
nb_rep <- 250

nullmodel=1

compute_global_fit <- function(seed, name){
  
  print(seed)
  
  res_global_fit <- c()
  
  network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";", header=TRUE)
  tree_hosts <- read.tree(paste0("tree_hosts_", name, "_", seed, ".tre"))
  tree_parasites <- read.tree(paste0("tree_parasites_", name, "_", seed, ".tre"))
  
  
  # Only keep interacting species
  
  network <- network[,which(colSums(network)>0)]
  network <- network[which(rowSums(network)>0),]
  tree_hosts <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[!tree_hosts$tip.label %in% rownames(network)])
  tree_parasites <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[!tree_parasites$tip.label %in% colnames(network)])
  
  
  ## Try PACo
  host_dist <- cophenetic(tree_hosts)
  parasite_dist <- cophenetic(tree_parasites)
  D <- prepare_paco_data(host_dist, parasite_dist, network)
  D <- add_pcoord(D)
  if (nullmodel==1) res_paco <- PACo_test(D, nperm = 10000, seed=42, nullmodel=1, symmetric=TRUE)   #res_paco <- PACo(D, nperm=10000, seed=42, method="r0") # problem with permutations
  if (nullmodel==2) res_paco <- PACo_test(D, nperm = 10000, seed=42, nullmodel=2, symmetric=TRUE)
  
  ## Try ParaFit
  if (nullmodel==1) res_parafit <- parafit_test(host_dist, parasite_dist, network, nperm = 10000, test.links = FALSE, correction = "none", nullmodel=1)
  if (nullmodel==2) res_parafit <- parafit_test(host_dist, parasite_dist, network, nperm = 10000, test.links = FALSE, correction = "none", nullmodel=2)
  
  
  res <- c(name, seed, Ntip(tree_hosts), Ntip(tree_parasites), round(res_paco$gof$ss,4), round(res_paco$gof$p,4), round(res_parafit$ParaFitGlobal,4), round(res_parafit$p.global,4)) 
  print(res)
  res_global_fit <- rbind(res_global_fit, res)
  
  res_global_fit <- data.frame(res_global_fit)
  colnames(res_global_fit) <- c("name", "seed", "N_hosts","N_parasites", "stat_paco", "p_paco", "stat_parafit", "p_parafit")
  
  write.table(res_global_fit, paste0("../global_fit/results_paco_parafit/results_paco_parafit_",name,"_", seed, "_nullmodel_",nullmodel,".csv"), sep=";", row.names=F)
  
}

mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_1', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_2', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_3', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_4', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_5', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_6', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_7', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_8', mc.preschedule = FALSE)


# summary all results
res_global_fit <- c()
for (name in c('simul_1', 'simul_2', 'simul_3', 'simul_4', 'simul_5', 'simul_6', 'simul_7', 'simul_8')){
  for (seed in 1:nb_rep){
    res_global_fit <- rbind(res_global_fit, read.table(paste0("../global_fit/results_paco_parafit/results_paco_parafit_",name,"_", seed, "_nullmodel_",nullmodel,".csv"), sep=";", header=TRUE))
  }
}

res_global_fit$N_hosts <- as.numeric(as.character(res_global_fit$N_hosts))
res_global_fit$N_parasites <- as.numeric(as.character(res_global_fit$N_parasites))
res_global_fit$stat_paco <- as.numeric(as.character(res_global_fit$stat_paco))
res_global_fit$p_paco <- as.numeric(as.character(res_global_fit$p_paco))
res_global_fit$stat_parafit <- as.numeric(as.character(res_global_fit$stat_parafit))
res_global_fit$p_parafit <- as.numeric(as.character(res_global_fit$p_parafit))


write.table(res_global_fit, paste0("../global_fit/results_all/results_paco_parafit_all_nullmodel_",nullmodel,".csv"), sep=";")

hist(res_global_fit$p_paco, breaks=20)
hist(res_global_fit$p_parafit, breaks=20)


table(res_global_fit$p_paco<0.05, res_global_fit$name)
table(res_global_fit$p_parafit<0.05, res_global_fit$name)


table(res_global_fit$p_paco<0.05, res_global_fit$name)/rbind(table(res_global_fit$name),table(res_global_fit$name))
table(res_global_fit$p_parafit<0.05, res_global_fit$name)/rbind(table(res_global_fit$name),table(res_global_fit$name))





##### Step 2B : Fit ParaFit and PACo for vicariance #####

rm(list=ls())

library(ape)
library(phytools)
library(mvMORPH)
library(paco)
library(parallel)


setwd("/path/simulated_data_geo/")

source("../global_fit/functions_parafit.R") # for null model 2

name='simul_geo_1'
seed=1

nb_cores <- 25
nb_rep <- 250

nullmodel=2

compute_global_fit <- function(seed, name){
  
  print(seed)
  
  res_global_fit <- c()
  
  network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";", header=TRUE)
  tree_hosts <- read.tree(paste0("tree_hosts_", name, "_", seed, ".tre"))
  tree_parasites <- read.tree(paste0("tree_parasites_", name, "_", seed, ".tre"))
  
  
  # Only keep interacting species
  
  network <- network[,which(colSums(network)>0)]
  network <- network[which(rowSums(network)>0),]
  tree_hosts <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[!tree_hosts$tip.label %in% rownames(network)])
  tree_parasites <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[!tree_parasites$tip.label %in% colnames(network)])
  
  
  ## Try PACo
  host_dist <- cophenetic(tree_hosts)
  parasite_dist <- cophenetic(tree_parasites)
  D <- prepare_paco_data(host_dist, parasite_dist, network)
  D <- add_pcoord(D)
  if (nullmodel==1) res_paco <- PACo(D, nperm=10000, seed=42, method="r0", symmetric=TRUE)
  if (nullmodel==2) res_paco <- PACo_test(D, nperm = 10000, seed=42, nullmodel=2, symmetric=TRUE)
  
  ## Try ParaFit
  if (nullmodel==1) res_parafit <- parafit(host_dist, parasite_dist, network, nperm = 10000, test.links = FALSE, correction = "none")
  if (nullmodel==2) res_parafit <- parafit_test(host_dist, parasite_dist, network, nperm = 10000, test.links = FALSE, correction = "none", nullmodel=2)
  
  
  res <- c(name, seed, Ntip(tree_hosts), Ntip(tree_parasites), round(res_paco$gof$ss,4), round(res_paco$gof$p,4), round(res_parafit$ParaFitGlobal,4), round(res_parafit$p.global,4)) 
  print(res)
  res_global_fit <- rbind(res_global_fit, res)
  
  res_global_fit <- data.frame(res_global_fit)
  colnames(res_global_fit) <- c("name", "seed", "N_hosts","N_parasites", "stat_paco", "p_paco", "stat_parafit", "p_parafit")
  
  write.table(res_global_fit, paste0("../global_fit/results_paco_parafit_geo/results_paco_parafit_",name,"_", seed, "_nullmodel_",nullmodel,".csv"), sep=";", row.names=F)
  
}

mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_geo_1', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_geo_2', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_geo_3', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_geo_4', mc.preschedule = FALSE)


# summary all results
res_global_fit <- c()
for (name in c('simul_geo_1', 'simul_geo_2', 'simul_geo_3', 'simul_geo_4')){
  for (seed in 1:nb_rep){
    res_global_fit <- rbind(res_global_fit, read.table(paste0("../global_fit/results_paco_parafit_geo/results_paco_parafit_",name,"_", seed, "_nullmodel_",nullmodel,".csv"), sep=";", header=TRUE))
  }
}

res_global_fit$N_hosts <- as.numeric(as.character(res_global_fit$N_hosts))
res_global_fit$N_parasites <- as.numeric(as.character(res_global_fit$N_parasites))
res_global_fit$stat_paco <- as.numeric(as.character(res_global_fit$stat_paco))
res_global_fit$p_paco <- as.numeric(as.character(res_global_fit$p_paco))
res_global_fit$stat_parafit <- as.numeric(as.character(res_global_fit$stat_parafit))
res_global_fit$p_parafit <- as.numeric(as.character(res_global_fit$p_parafit))


write.table(res_global_fit, paste0("../global_fit/results_all/results_paco_parafit_geo_all_nullmodel_",nullmodel,".csv"), sep=";")

hist(res_global_fit$p_paco, breaks=20)
hist(res_global_fit$p_parafit, breaks=20)


table(res_global_fit$p_paco<0.05, res_global_fit$name)
table(res_global_fit$p_parafit<0.05, res_global_fit$name)


table(res_global_fit$p_paco<0.05, res_global_fit$name)/nb_rep
table(res_global_fit$p_parafit<0.05, res_global_fit$name)/nb_rep




##### Step 2C : Fit ParaFit and PACo for codiversification #####

rm(list=ls())

library(ape)
library(phytools)
library(mvMORPH)
library(paco)
library(parallel)


setwd("/path/simulated_data_codiv/")

source("../global_fit/functions_parafit.R") # for null model 2

name='simul_codiv_1'
seed=1

nb_cores <- 25
nb_rep <- 250

nullmodel=1

compute_global_fit <- function(seed, name){
  
  print(seed)
  
  res_global_fit <- c()
  
  network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";", header=TRUE)
  tree_hosts <- read.tree(paste0("tree_hosts_", name, "_", seed, ".tre"))
  tree_parasites <- read.tree(paste0("tree_parasites_", name, "_", seed, ".tre"))
  
  
  # Only keep interacting species
  
  network <- network[,which(colSums(network)>0)]
  network <- network[which(rowSums(network)>0),]
  tree_hosts <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[!tree_hosts$tip.label %in% rownames(network)])
  tree_parasites <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[!tree_parasites$tip.label %in% colnames(network)])
  
  
  ## Try PACo
  host_dist <- cophenetic(tree_hosts)
  parasite_dist <- cophenetic(tree_parasites)
  D <- prepare_paco_data(host_dist, parasite_dist, network)
  D <- add_pcoord(D)
  if (nullmodel==1) res_paco <- PACo(D, nperm=10000, seed=42, method="r0", symmetric=TRUE)
  if (nullmodel==2) res_paco <- PACo_test(D, nperm = 10000, seed=42, nullmodel=2, symmetric=TRUE)
  
  ## Try ParaFit
  if (nullmodel==1) res_parafit <- parafit(host_dist, parasite_dist, network, nperm = 10000, test.links = FALSE, correction = "none")
  if (nullmodel==2) res_parafit <- parafit_test(host_dist, parasite_dist, network, nperm = 10000, test.links = FALSE, correction = "none", nullmodel=2)
  
  
  res <- c(name, seed, Ntip(tree_hosts), Ntip(tree_parasites), round(res_paco$gof$ss,4), round(res_paco$gof$p,4), round(res_parafit$ParaFitGlobal,4), round(res_parafit$p.global,4)) 
  print(res)
  res_global_fit <- rbind(res_global_fit, res)
  
  res_global_fit <- data.frame(res_global_fit)
  colnames(res_global_fit) <- c("name", "seed", "N_hosts","N_parasites", "stat_paco", "p_paco", "stat_parafit", "p_parafit")
  
  write.table(res_global_fit, paste0("../global_fit/results_paco_parafit_codiv/results_paco_parafit_",name,"_", seed, "_nullmodel_",nullmodel,".csv"), sep=";", row.names=F)
  
}

mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_codiv_1', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_codiv_2', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_codiv_3', mc.preschedule = FALSE)
mclapply(1:nb_rep, compute_global_fit, mc.cores=nb_cores, name='simul_codiv_4', mc.preschedule = FALSE)


# summary all results
res_global_fit <- c()
for (name in c('simul_codiv_1', 'simul_codiv_2', 'simul_codiv_3', 'simul_codiv_4')){
  for (seed in 1:nb_rep){
    res_global_fit <- rbind(res_global_fit, read.table(paste0("../global_fit/results_paco_parafit_codiv/results_paco_parafit_",name,"_", seed, "_nullmodel_",nullmodel,".csv"), sep=";", header=TRUE))
  }
}

res_global_fit$N_hosts <- as.numeric(as.character(res_global_fit$N_hosts))
res_global_fit$N_parasites <- as.numeric(as.character(res_global_fit$N_parasites))
res_global_fit$stat_paco <- as.numeric(as.character(res_global_fit$stat_paco))
res_global_fit$p_paco <- as.numeric(as.character(res_global_fit$p_paco))
res_global_fit$stat_parafit <- as.numeric(as.character(res_global_fit$stat_parafit))
res_global_fit$p_parafit <- as.numeric(as.character(res_global_fit$p_parafit))


write.table(res_global_fit, paste0("../global_fit/results_all/results_paco_parafit_codiv_all_nullmodel_",nullmodel,".csv"), sep=";")


hist(res_global_fit$p_paco, breaks=20)
hist(res_global_fit$p_parafit, breaks=20)


table(res_global_fit$p_paco<0.05, res_global_fit$name)
table(res_global_fit$p_parafit<0.05, res_global_fit$name)


table(res_global_fit$p_paco<0.05, res_global_fit$name)/nb_rep
table(res_global_fit$p_parafit<0.05, res_global_fit$name)/nb_rep






##### Step 2D : Plot results Global fit approaches ####



rm(list=ls())

setwd("/path/global_fit/")

library(ggplot2)

nullmodel=1
type="all"
type_2 <- "trait matching"

res_all <- c()
res_R2_paco <- c()

for (nullmodel in c(1,2)){
  for (type in c("all", "geo_new_all", "codiv_all")){
    
    if (type=="all") type_2 <- "trait matching"
    if (type=="geo_new_all") type_2 <- "vicariance"
    if (type=="codiv_all") type_2 <- "codiversification"
    
    
    print(nullmodel)
    print(type)
    print(type_2)
    
    res_global_fit <- read.table(paste0("results_paco_parafit_",type,"_symmetric_nullmodel_",nullmodel,".csv"), header=TRUE, sep=";")
    
    res_global_fit$type <- type_2
    
    res_connectance <- read.table(paste0("summary_stats_connectance_ratio.csv"), header=TRUE, sep=";")
    res_connectance_type <- res_connectance[grep(type_2, res_connectance$simul),]
    
    res_global_fit$connectance <- NA
    res_global_fit$Ratio_one_to_one <- NA
    for (i in 1:nrow(res_global_fit)){
      index <- intersect(which(res_connectance_type$name==res_global_fit$name[i]), which(res_connectance_type$seed==res_global_fit$seed[i]))
      res_global_fit$Ratio_one_to_one[i] <- res_connectance_type$Ratio_one_to_one[index]
      res_global_fit$connectance[i] <- res_connectance_type$Connectance[index]
    }
    
    res_global_fit$N_tot <- res_global_fit$N_hosts+res_global_fit$N_parasites
    res_global_fit$significance_paco <- res_global_fit$p_paco<0.05
    res_global_fit$significance_parafit <- res_global_fit$p_parafit<0.05
    
    res_global_fit$R2_paco <- 1-res_global_fit$stat_paco
    
    res_global_fit <- res_global_fit[which(res_global_fit$stat_paco<=1),]
    
    res_global_fit_paco_sig <- res_global_fit[which(res_global_fit$p_paco<0.05),]

    freq_paco <- c(length(which(res_global_fit_paco_sig$R2_paco<=0.25)),
                   length(intersect(which(res_global_fit_paco_sig$R2_paco>0.25), which(res_global_fit_paco_sig$R2_paco<=0.5))),
                   length(which(res_global_fit_paco_sig$R2_paco>0.50)))
    freq_paco <- freq_paco/sum(freq_paco)
    
    res_R2_paco <- rbind(res_R2_paco, c(nullmodel, type_2, round(freq_paco, 3)))
    
    #### stats 
    print("PACo")
    res_global_fit$significance_paco_numeric <- as.numeric(res_global_fit$significance_paco)
    model <- glm(formula = significance_paco_numeric ~ N_tot, family = binomial(link = "logit"), data=res_global_fit)
    print(summary(model))
    model <- glm(formula = significance_paco_numeric ~ connectance, family = binomial(link = "logit"), data=res_global_fit)
    print(summary(model))
    model <- glm(formula = significance_paco_numeric ~ Ratio_one_to_one, family = binomial(link = "logit"), data=res_global_fit)
    print(summary(model))
    
    print("ParaFit")
    res_global_fit$significance_parafit_numeric <- as.numeric(res_global_fit$significance_parafit)
    model <- glm(formula = significance_parafit_numeric ~ N_tot, family = binomial(link = "logit"), data=res_global_fit)
    print(summary(model))
    model <- glm(formula = significance_parafit_numeric ~ connectance, family = binomial(link = "logit"), data=res_global_fit)
    print(summary(model))
    model <- glm(formula = significance_parafit_numeric ~ Ratio_one_to_one, family = binomial(link = "logit"), data=res_global_fit)
    print(summary(model))
    
    ### PACo
    pdf(paste0("results_PACo_",type_2,"_stats_",nullmodel,".pdf"), height=3.5, width=5)
    print(ggplot(res_global_fit)+geom_point(aes(x=N_tot, y=R2_paco, col=significance_paco, shape=significance_paco))+
            theme_bw()+ylim(c(0,1))+xlab("Total number of species")+ ylab(expression(paste(R^2, " (PACo)")))+
            scale_color_manual(values=c("#1a5276", "#b03a2e"))+ guides(col=guide_legend(title="Significance\nPACO:"), shape="none"))
    
    print(ggplot(res_global_fit)+geom_point(aes(x=connectance, y=R2_paco, col=significance_paco, shape=significance_paco))+
            theme_bw()+ylim(c(0,1))+xlab("Connectance")+ ylab(expression(paste(R^2, " (PACo)")))+
            scale_color_manual(values=c("#1a5276", "#b03a2e"))+ guides(col=guide_legend(title="Significance\nPACO:"), shape="none"))
    
    print(ggplot(res_global_fit)+geom_point(aes(x=Ratio_one_to_one, y=R2_paco, col=significance_paco, shape=significance_paco))+
            theme_bw()+ylim(c(0,1))+xlab("Ratio of one-to-one interactions")+ ylab(expression(paste(R^2, " (PACo)")))+
            scale_color_manual(values=c("#1a5276", "#b03a2e"))+ guides(col=guide_legend(title="Significance\nPACO:"), shape="none"))
    dev.off()
    
    pdf(paste0("results_ratio_PACo_",type_2,"_stats_",nullmodel,".pdf"), height=3.5, width=3.5)
    
    print(ggplot(res_global_fit)+geom_boxplot(aes(x=significance_paco, y=Ratio_one_to_one, col=significance_paco))+
            theme_bw()+ylim(c(0,1))+xlab("Significance PACo")+ ylab(expression(paste("Ratio of one-to-one interactions")))+
            scale_color_manual(values=c("#1a5276", "#b03a2e"))+ theme(legend.position="none"))

    dev.off()
    
    ### ParaFit
    pdf(paste0("results_ParaFit_",type_2,"_stats_",nullmodel,".pdf"), height=3.5, width=5)
    print(ggplot(res_global_fit)+geom_point(aes(x=N_tot, y=stat_parafit, col=significance_parafit, shape=significance_parafit))+
            theme_bw()+xlab("Total number of species")+ ylab("ParaFit global statistic")+scale_y_continuous(trans='log10')+
            scale_color_manual(values=c("#1a5276", "#b03a2e"))+ guides(col=guide_legend(title="Significance\nParaFit:"), shape="none"))
    
    print(ggplot(res_global_fit)+geom_point(aes(x=connectance, y=stat_parafit, col=significance_parafit, shape=significance_parafit))+
            theme_bw()+xlab("Connectance")+ ylab("ParaFit global statistic")+scale_y_continuous(trans='log10')+
            scale_color_manual(values=c("#1a5276", "#b03a2e"))+ guides(col=guide_legend(title="Significance\nParaFit:"), shape="none"))
    
    print(ggplot(res_global_fit)+geom_point(aes(x=Ratio_one_to_one, y=stat_parafit, col=significance_parafit, shape=significance_parafit))+
            theme_bw()+xlab("Ratio of one-to-one interactions")+ ylab("ParaFit global statistic")+scale_y_continuous(trans='log10')+
            scale_color_manual(values=c("#1a5276", "#b03a2e"))+ guides(col=guide_legend(title="Significance\nParaFit:"), shape="none"))
    dev.off()
    
    pdf(paste0("results_ratio_ParaFit_",type_2,"_stats_",nullmodel,".pdf"), height=3.5, width=3.5)
    
    print(ggplot(res_global_fit)+geom_boxplot(aes(x=significance_parafit, y=Ratio_one_to_one, col=significance_parafit))+
            theme_bw()+ylim(c(0,1))+xlab("Significance ParaFit")+ ylab(expression(paste("Ratio of one-to-one interactions")))+
            scale_color_manual(values=c("#1a5276", "#b03a2e"))+ theme(legend.position="none"))
    
    dev.off()
    

    res_global_fit$null_model <- nullmodel
    res_all <- rbind(res_all, res_global_fit)
    
  }
}


res_R2_paco

res_all_nm <- res_all

res_all <- res_all[res_all$null_model==1,]
res_all$type_res_paco <- paste0(res_all$type,"_", res_all$significance_paco)
res_all$type_res_paco <- as.factor(res_all$type_res_paco)
res_all$type_res_paco <- factor(res_all$type_res_paco, levels = levels(res_all$type_res_paco)[c(4,3,6,5,2,1)])

pdf("summary_Paco_stats.pdf", width=3, height=2.25)
ggplot(res_all, aes(y=R2_paco, x=type_res_paco, fill=type_res_paco))+geom_violin(trim=TRUE, draw_quantiles=0.5)+
  theme_bw()+ylim(c(0,1))+xlab("")+ ylab(expression(paste(R^2, " (PACo)")))+
  scale_fill_manual(values=c("#b03a2e","#1a5276","#b03a2e","#1a5276","#b03a2e","#1a5276"))+ theme(legend.position = "none")
dev.off()


res_all$type_res_parafit <- paste0(res_all$type,"_", res_all$significance_parafit)
res_all$type_res_parafit <- as.factor(res_all$type_res_parafit)
res_all$type_res_parafit <- factor(res_all$type_res_parafit, levels = levels(res_all$type_res_parafit)[c(4,3,6,5,2,1)])

pdf("summary_ParaFit_stats.pdf", width=3, height=2.25)
ggplot(res_all, aes(y=stat_parafit, x=type_res_parafit, fill=type_res_parafit))+geom_violin(trim=TRUE, draw_quantiles=0.5)+
  theme_bw()+xlab("")+ ylab("ParaFit global statistic")+scale_y_continuous(trans='log10')+
  scale_fill_manual(values=c("#b03a2e","#1a5276","#b03a2e","#1a5276","#b03a2e","#1a5276"))+ theme(legend.position = "none")
dev.off()


nb_rep=250
table(res_all$p_paco<0.05, res_all$name)/nb_rep
table(res_all$p_parafit<0.05, res_all$name)/nb_rep


# nm2
res_all_nm2 <- res_all_nm[res_all_nm$null_model==2,]
table(res_all_nm2$p_parafit<0.05, res_all_nm2$name)/nb_rep
table(res_all_nm2$p_paco<0.05, res_all_nm2$name)/nb_rep





##### Step 3: Plot figure  - trait matching #####

rm(list=ls())

library(ape)
library(phytools)
library(mvMORPH)
library(paco)


setwd("/path/simulated_data/")


res_global_fit <- read.table(paste0("results_paco_parafit_all.csv"), sep=";", header=TRUE)


## Automatic selection

res_global_fit_signif <- res_global_fit[intersect(which(res_global_fit$p_parafit<0.05), which(res_global_fit$p_paco<0.05)),]

# remove the 31 (almost polytomy)
res_global_fit_signif <- res_global_fit_signif[-31,]

name=res_global_fit_signif$name[which.min(c(res_global_fit_signif$N_hosts+res_global_fit_signif$N_parasites))]
seed=res_global_fit_signif$seed[which.min(c(res_global_fit_signif$N_hosts+res_global_fit_signif$N_parasites))]

network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";", header=TRUE)
tree_hosts <- read.tree(paste0("tree_hosts_", name, "_", seed, ".tre"))
tree_parasites <- read.tree(paste0("tree_parasites_", name, "_", seed, ".tre"))
traits_hosts <- read.table(paste0("traits_hosts_", name, "_", seed, ".csv"), sep=";", header=TRUE)
traits_parasites <- read.table(paste0("traits_parasites_", name, "_", seed, ".csv"), sep=";", header=TRUE)



network <- network[,which(colSums(network)>0)]
network <- network[which(rowSums(network)>0),]
tree_hosts <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[!tree_hosts$tip.label %in% rownames(network)])
tree_parasites <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[!tree_parasites$tip.label %in% colnames(network)])

assoc <- reshape2::melt(as.matrix(network))
assoc <- assoc[assoc$value>0,]

pdf("figure/cophylo_traits.pdf", width=4, height=4)
cophyloplot(tree_hosts, tree_parasites, assoc = assoc, use.edge.length = TRUE, show.tip.label = F, space=50, 
            length.line=0, lwd=1.75, col="#626567")
dev.off()



# Plot traits

traits_hosts <- traits_hosts[tree_hosts$tip.label, 1]
names(traits_hosts) <- tree_hosts$tip.label

traits_parasites <- traits_parasites[tree_parasites$tip.label, 1]
names(traits_parasites) <- tree_parasites$tip.label

map_hosts <- contMap(tree_hosts, traits_hosts, res=200, fsize=c(0.6,0.8), outline=F, lims=NULL,lwd=5)
map_hosts <- setMap(map_hosts,c("#FFFFB2","#FECC5C","#FD8D3C","#E31A1C")) ## change color scheme

map_parasites <- contMap(tree_parasites, traits_parasites, res=200, fsize=c(0.6,0.8), outline=F, lims=NULL,lwd=5)
map_parasites <- setMap(map_parasites,c("#FFFFB2","#FECC5C","#FD8D3C","#E31A1C")) ## change color scheme


pdf("figure/trait_values.pdf", width=3, height=4)
plot(map_hosts,fsize=c(0.7,0.8), leg.txt="trait value", outline=F, lims=c(min(c(traits_hosts, traits_parasites)), max(c(traits_hosts, traits_parasites))), lwd=5)
add.scale.bar()
plot(map_parasites,fsize=c(0.7,0.8), leg.txt="trait value", outline=F, lims=c(min(c(traits_hosts, traits_parasites)), max(c(traits_hosts, traits_parasites))),lwd=5)
add.scale.bar()
dev.off()









##### Step 3: Plot figure  - vicariance #####

rm(list=ls())

library(ape)
library(phytools)
library(mvMORPH)
library(paco)


setwd("/path/simulated_data_geo/")


res_global_fit <- read.table(paste0("../global_fit/results_paco_parafit_geo_all_nullmodel_1.csv"), sep=";", header=TRUE)


## Automatic selection

res_global_fit_signif <- res_global_fit[intersect(which(res_global_fit$p_parafit<0.05), which(res_global_fit$p_paco<0.05)),]

res_global_fit_signif <- res_global_fit_signif[res_global_fit_signif$name!="simul_geo_1", ]

# remove the 31 (almost polytomy)
res_global_fit_signif <- res_global_fit_signif[-61,]

name=res_global_fit_signif$name[which.min(c(res_global_fit_signif$N_hosts+res_global_fit_signif$N_parasites))]
seed=res_global_fit_signif$seed[which.min(c(res_global_fit_signif$N_hosts+res_global_fit_signif$N_parasites))]


res_global_fit_signif$p_paco[which.min(c(res_global_fit_signif$N_hosts+res_global_fit_signif$N_parasites))]
res_global_fit_signif$p_parafit[which.min(c(res_global_fit_signif$N_hosts+res_global_fit_signif$N_parasites))]


network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";", header=TRUE)
tree_hosts <- read.tree(paste0("tree_hosts_", name, "_", seed, ".tre"))
tree_parasites <- read.tree(paste0("tree_parasites_", name, "_", seed, ".tre"))


network <- network[,which(colSums(network)>0)]
network <- network[which(rowSums(network)>0),]
tree_hosts <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[!tree_hosts$tip.label %in% rownames(network)])
tree_parasites <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[!tree_parasites$tip.label %in% colnames(network)])

assoc <- reshape2::melt(as.matrix(network))
assoc <- assoc[assoc$value>0,]

color_geo=c("#5dade2","#eb984e", "#58d68d")
names(color_geo) <- c("A", "B", "C")

pdf("../simulated_data/figure/cophylo_vicariance.pdf", width=4, height=4)
cophyloplot(tree_hosts, tree_parasites, assoc = assoc, use.edge.length = TRUE, show.tip.label = F, space=50, 
            length.line=0, lwd=1.25, col=color_geo[substring(assoc$Var1, 2, 2)])
dev.off()




# plot geography 

pdf("../simulated_data/figure/geo_values.pdf", width=3, height=4)
trait_info <- substring(tree_hosts$tip.label, 2, 2)  
names(trait_info) <-  tree_hosts$tip.label 
plot(tree_hosts, show.tip.label = FALSE, edge.width = 1.5)
tiplabels(pie=to.matrix(trait_info,sort(unique(trait_info))),piecol=c("#5dade2","#eb984e", "#58d68d"),cex=0.75)
add.scale.bar()

trait_info <- substring(tree_parasites$tip.label, 2, 2)  
names(trait_info) <-  tree_parasites$tip.label 
plot(tree_parasites, show.tip.label = FALSE, edge.width = 1.5)
tiplabels(pie=to.matrix(trait_info,sort(unique(trait_info))),piecol=c("#5dade2","#eb984e", "#58d68d"),cex=0.75)
add.scale.bar()
dev.off()





##### Step 3: Plot figure  - codiversification #####

rm(list=ls())

library(ape)
library(phytools)
library(mvMORPH)
library(paco)


setwd("/path/simulated_data_codiv/")


res_global_fit <- read.table(paste0("../global_fit/results_paco_parafit_codiv_all_nullmodel_1.csv"), sep=";", header=TRUE)


res_global_fit_signif <- res_global_fit[intersect(which(res_global_fit$p_parafit<0.05), which(res_global_fit$p_paco<0.05)),]

name=res_global_fit_signif$name[which.min(c(res_global_fit_signif$N_hosts+res_global_fit_signif$N_parasites))]
seed=res_global_fit_signif$seed[which.min(c(res_global_fit_signif$N_hosts+res_global_fit_signif$N_parasites))]

seed=res_global_fit_signif$seed[219]

pdf("figure/cophylo_traits.pdf", width=4, height=4)
for (seed in res_global_fit_signif$seed[res_global_fit_signif$name=="simul_codiv_1"]){
  
  network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";", header=TRUE)
  tree_hosts <- read.tree(paste0("tree_hosts_", name, "_", seed, ".tre"))
  tree_parasites <- read.tree(paste0("tree_parasites_", name, "_", seed, ".tre"))
  
  network <- network[,which(colSums(network)>0)]
  network <- network[which(rowSums(network)>0),]
  tree_hosts <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[!tree_hosts$tip.label %in% rownames(network)])
  tree_parasites <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[!tree_parasites$tip.label %in% colnames(network)])
  
  if ((Ntip(tree_hosts)+Ntip(tree_parasites))<40){
    
    assoc <- reshape2::melt(as.matrix(network))
    assoc <- assoc[assoc$value>0,]
    
    if (Ntip(tree_hosts)!=Ntip(tree_parasites)){
      print(seed)
      cophyloplot(tree_hosts, tree_parasites, assoc = assoc, use.edge.length = TRUE, show.tip.label = F, space=50, 
                  length.line=0, lwd=1.75, col="#626567")
    }}}
dev.off()


### plot trees

seed=231

network <- read.table(paste0("network_", name, "_", seed, ".csv"), sep=";", header=TRUE)
tree_hosts <- read.tree(paste0("tree_hosts_", name, "_", seed, ".tre"))
tree_parasites <- read.tree(paste0("tree_parasites_", name, "_", seed, ".tre"))

network <- network[,which(colSums(network)>0)]
network <- network[which(rowSums(network)>0),]
tree_hosts <- drop.tip(tree_hosts, tip=tree_hosts$tip.label[!tree_hosts$tip.label %in% rownames(network)])
tree_parasites <- drop.tip(tree_parasites, tip=tree_parasites$tip.label[!tree_parasites$tip.label %in% colnames(network)])

pdf("figure/trees.pdf", width=3, height=4)

plot(tree_hosts, show.tip.label = FALSE, edge.width = 1.5)
plot(tree_parasites, show.tip.label = FALSE, edge.width = 1.5)

tree_hosts_pagel <- tree_hosts

tree_hosts_pagel$edge.length[which(tree_hosts_pagel$edge[,2] %in% 1:Ntip(tree_hosts_pagel))] <- tree_hosts_pagel$edge.length[which(tree_hosts_pagel$edge[,2] %in% 1:Ntip(tree_hosts_pagel))] +10
plot(tree_hosts_pagel, show.tip.label = FALSE, edge.width = 1.5)

dev.off()












##### Step 4: Results eMPRess  #####


#### Results for trait matching  ####

rm(list=ls())

library(ggplot2)
library(ape)

setwd("/path/emPress/results/")
path_data="../../simulated_data/"
nsim=250

res_stats <- c()

for (costs in c("d1_t1_l1", "d4_t1_l1", "d2_t1_l2", "d4_t2_l1",  "d2_t3_l1")){
  list_pvalues <- c()
  print(costs)
  for (name in c("1", "2", "3", "4", "5", "6", "7", "8")){
    for (seed in 1:nsim){
      
      if (file.exists(paste0("res_",costs,"_simul_",name,"_",seed,".svg"))){
        res <- read.table(paste0("res_",costs,"_simul_",name,"_",seed,".svg"), comment.char = "", fill=TRUE, sep=";")
        
        res_pvalue <- res$V1[grep("p-value", res$V1)]
        res_pvalue <- gsub("    <!-- p-value = ", "", res_pvalue)
        res_pvalue <- as.numeric(gsub(" -->", "", res_pvalue))
        
        recon <- read.table(paste0("reconciliation_",costs,"_simul_",name,"_",seed,".csv"), sep=",")
        colnames(recon) <- c("parasite", "host", "event", "node_frequency", "event_frequency")
        
        if (!all(recon$node_frequency>=recon$event_frequency)){print("problem")} # makes sense
        
        
        # check the time consistency (check for back-in-time transfers)
        
        tree_hosts <- read.tree(paste0(path_data, "tree_hosts_only_int_simul_",name,"_",seed,".tre"))
        tree_parasites <- read.tree(paste0(path_data, "tree_parasites_only_int_simul_",name,"_",seed,".tre"))
        
        for (i in 1:nrow(recon)){
          if (length(grep("_",recon$parasite[i]))==1){
            recon$parasite[i] <- as.numeric(gsub("_p", "", recon$parasite[i]))+Ntip(tree_parasites)+1
          } else {
            recon$parasite[i] <- which(tree_parasites$tip.label==recon$parasite[i])
          }
          if (length(grep("_",recon$host[i]))==1){
            recon$host[i] <- as.numeric(gsub("_h", "", recon$host[i]))+Ntip(tree_hosts)+1
          } else {
            recon$host[i] <- which(tree_hosts$tip.label==recon$host[i])
          }
        }
        
        node_age_hosts <- node.depth.edgelength(tree_hosts)
        node_age_hosts <- max(node_age_hosts)-node_age_hosts
        names(node_age_hosts) <- 1:length(node_age_hosts)
        
        recon$age_min <- NA
        recon$age_max <- NA
        
        for (i in 1:nrow(recon)){
          
          if (recon$event[i]=="Contemporaneous"){
            recon$age_min[i] <- 0
            recon$age_max[i] <- 0
          }
          
          if (recon$event[i]=="Cospeciation"){
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            recon$age_max[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Transfer"){
            
            # the most ancient age of departure (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent age of arrival
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Loss"){
            
            # the most ancient age of loss (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent loss 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
          
          if (recon$event[i]=="Duplication"){
            
            # the most ancient age of duplication (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent duplication (ideally, it should be constrained with next nodes, but should not affect the measure of time consistency) 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
        }
        
        recon$diff_age <- recon$age_max-recon$age_min
        
        list_pvalues <- rbind(list_pvalues, c(name, seed, res_pvalue, table(recon$event)[c("Cospeciation", "Transfer", "Duplication", "Loss")], length(which(recon$diff_age<0)), length(which(recon$diff_age[recon$event=="Transfer"]<0)) ))
      }
    }
  }
  
  list_pvalues <- data.frame(list_pvalues)
  colnames(list_pvalues) <- c("name", "seed", "pvalue", "cospeciation", "transfer", "duplication", "loss", "time_consistency", "time_consistency_transfers")
  print(nrow(list_pvalues))
  list_pvalues$pvalue <- as.numeric(as.character(list_pvalues$pvalue))
  list_pvalues$cospeciation <- as.numeric(as.character(list_pvalues$cospeciation))
  list_pvalues$transfer <- as.numeric(as.character(list_pvalues$transfer))
  

  print(costs)
  print("permutations:")
  table <- table(list_pvalues$pvalue<0.05, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "permutations",table[nrow(table),]))
  
  
  print("time consistency:")
  list_pvalues$signif_time_consistency <- FALSE
  list_pvalues$signif_time_consistency[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency",table[nrow(table),]))
  
  print("number cospe/switches:")
  list_pvalues$signif_switches <- FALSE
  list_pvalues$signif_switches[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "cospe_switches",table[nrow(table),]))
  
  print("time consistency - number cospe/switches:")
  list_pvalues$signif_time_consistency_switches <- FALSE
  list_pvalues$signif_time_consistency_switches[intersect(intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0)), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency_switches",table[nrow(table),]))
  
  number_events_1 <- cbind(list_pvalues[,c(1,2,3,4)], rep("cospeciation", nrow(list_pvalues)))
  colnames(number_events_1) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_2 <- cbind(list_pvalues[,c(1,2,3,5)], rep("transfer", nrow(list_pvalues)))
  colnames(number_events_2) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_3 <- cbind(list_pvalues[,c(1,2,3,6)], rep("duplication", nrow(list_pvalues)))
  colnames(number_events_3) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_4 <- cbind(list_pvalues[,c(1,2,3,7)], rep("loss", nrow(list_pvalues)))
  colnames(number_events_4) <- c("name", "seed", "pvalue", "number_events", "event")
  
  number_events_all <- rbind(number_events_1,number_events_2, number_events_3, number_events_4)
  
  number_events_all <- data.frame(number_events_all)
  colnames(number_events_all) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_all$pvalue <- as.numeric(as.character(number_events_all$pvalue))
  number_events_all$number_events <- as.numeric(as.character(number_events_all$number_events))
  number_events_all$signif <- number_events_all$pvalue<0.05
  number_events_all$number_events[is.na(number_events_all$number_events)] <- 0
  number_events_all$name_seed <- paste0(number_events_all$name, "_", number_events_all$seed)
  
  number_events_all$event <- as.factor(number_events_all$event)
  number_events_all$event <- factor(number_events_all$event, levels = c("cospeciation", "transfer", "duplication", "loss"))
  
  pdf(paste0("../result_plots/emPress_number_events_", costs,".pdf"), width = 4.25, height = 2.5)
  
  print(ggplot(number_events_all[which(number_events_all$signif==TRUE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[which(number_events_all$signif==FALSE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d2b4de")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="1")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="2")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="3")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+scale_y_sqrt()+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="4")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  dev.off()
  
  write.table(list_pvalues, paste0("../result_plots/results_emPress_", costs,".csv"), sep=";", row.names=FALSE)
  
}

res_stats <- data.frame(res_stats)
colnames(res_stats) <- c("costs", "type", 1:8)
write.table(res_stats, paste0("../result_plots/results_stats_emPress_traits.csv"), sep=";", row.names=FALSE)






#### Results for trait matching - random bifurcations  ####


rm(list=ls())

library(ggplot2)
library(ape)

setwd("/path/emPress/results/")
path_data="../../simulated_data/"
nsim=250

res_stats <- c()

for (costs in c("d1_t1_l1", "d4_t1_l1", "d2_t1_l2", "d4_t2_l1",  "d2_t3_l1")){
  list_pvalues <- c()
  print(costs)
  for (name in c("1", "2", "3", "4", "5", "6", "7", "8")){
    for (seed in 1:nsim){
      
      if (file.exists(paste0("res_bifurcations_eMPRess_",costs,"_simul_",name,"_",seed,".svg"))){
        res <- read.table(paste0("res_bifurcations_eMPRess_",costs,"_simul_",name,"_",seed,".svg"), comment.char = "", fill=TRUE, sep=";")
        
        res_pvalue <- res$V1[grep("p-value", res$V1)]
        res_pvalue <- gsub("    <!-- p-value = ", "", res_pvalue)
        res_pvalue <- as.numeric(gsub(" -->", "", res_pvalue))
        
        recon <- read.table(paste0("reconciliation_bifurcations_eMPRess_",costs,"_simul_",name,"_",seed,".csv"), sep=",")
        colnames(recon) <- c("parasite", "host", "event", "node_frequency", "event_frequency")
        
        if (!all(recon$node_frequency>=recon$event_frequency)){print("problem")} # makes sense
        
        
        # check the time consistency (check for back-in-time transfers)
        
        tree_hosts <- read.tree(paste0(path_data, "tree_hosts_only_int_simul_",name,"_",seed,".tre"))
        tree_parasites <- read.tree(paste0(path_data, "tree_parasites_bifurcations_eMPRess_simul_",name,"_",seed,".tre"))
        
        for (i in 1:nrow(recon)){
          if (length(grep("_",recon$parasite[i]))==1){
            recon$parasite[i] <- as.numeric(gsub("_p", "", recon$parasite[i]))+Ntip(tree_parasites)+1
          } else {
            recon$parasite[i] <- which(tree_parasites$tip.label==recon$parasite[i])
          }
          if (length(grep("_",recon$host[i]))==1){
            recon$host[i] <- as.numeric(gsub("_h", "", recon$host[i]))+Ntip(tree_hosts)+1
          } else {
            recon$host[i] <- which(tree_hosts$tip.label==recon$host[i])
          }
        }
        
        node_age_hosts <- node.depth.edgelength(tree_hosts)
        node_age_hosts <- max(node_age_hosts)-node_age_hosts
        names(node_age_hosts) <- 1:length(node_age_hosts)
        
        recon$age_min <- NA
        recon$age_max <- NA
        
        for (i in 1:nrow(recon)){
          
          if (recon$event[i]=="Contemporaneous"){
            recon$age_min[i] <- 0
            recon$age_max[i] <- 0
          }
          
          if (recon$event[i]=="Cospeciation"){
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            recon$age_max[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Transfer"){
            
            # the most ancient age of departure (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent age of arrival
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Loss"){
            
            # the most ancient age of loss (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent loss 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
          
          if (recon$event[i]=="Duplication"){
            
            # the most ancient age of duplication (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent duplication (ideally, it should be constrained with next nodes, but should not affect the measure of time consistency) 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
        }
        
        recon$diff_age <- recon$age_max-recon$age_min
        
        list_pvalues <- rbind(list_pvalues, c(name, seed, res_pvalue, table(recon$event)[c("Cospeciation", "Transfer", "Duplication", "Loss")], length(which(recon$diff_age<0)), length(which(recon$diff_age[recon$event=="Transfer"]<0)) ))
      }
    }
  }
  
  list_pvalues <- data.frame(list_pvalues)
  colnames(list_pvalues) <- c("name", "seed", "pvalue", "cospeciation", "transfer", "duplication", "loss", "time_consistency", "time_consistency_transfers")
  print(nrow(list_pvalues))
  list_pvalues$pvalue <- as.numeric(as.character(list_pvalues$pvalue))
  list_pvalues$cospeciation <- as.numeric(as.character(list_pvalues$cospeciation))
  list_pvalues$transfer <- as.numeric(as.character(list_pvalues$transfer))
  
  print(costs)
  print("permutations:")
  table <- table(list_pvalues$pvalue<0.05, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "permutations",table[nrow(table),]))
  
  
  print("time consistency:")
  list_pvalues$signif_time_consistency <- FALSE
  list_pvalues$signif_time_consistency[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency",table[nrow(table),]))
  
  print("number cospe/switches:")
  list_pvalues$signif_switches <- FALSE
  list_pvalues$signif_switches[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "cospe_switches",table[nrow(table),]))
  
  print("time consistency - number cospe/switches:")
  list_pvalues$signif_time_consistency_switches <- FALSE
  list_pvalues$signif_time_consistency_switches[intersect(intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0)), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency_switches",table[nrow(table),]))
  
  number_events_1 <- cbind(list_pvalues[,c(1,2,3,4)], rep("cospeciation", nrow(list_pvalues)))
  colnames(number_events_1) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_2 <- cbind(list_pvalues[,c(1,2,3,5)], rep("transfer", nrow(list_pvalues)))
  colnames(number_events_2) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_3 <- cbind(list_pvalues[,c(1,2,3,6)], rep("duplication", nrow(list_pvalues)))
  colnames(number_events_3) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_4 <- cbind(list_pvalues[,c(1,2,3,7)], rep("loss", nrow(list_pvalues)))
  colnames(number_events_4) <- c("name", "seed", "pvalue", "number_events", "event")
  
  number_events_all <- rbind(number_events_1,number_events_2, number_events_3, number_events_4)
  
  number_events_all <- data.frame(number_events_all)
  colnames(number_events_all) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_all$pvalue <- as.numeric(as.character(number_events_all$pvalue))
  number_events_all$number_events <- as.numeric(as.character(number_events_all$number_events))
  number_events_all$signif <- number_events_all$pvalue<0.05
  number_events_all$number_events[is.na(number_events_all$number_events)] <- 0
  number_events_all$name_seed <- paste0(number_events_all$name, "_", number_events_all$seed)
  
  number_events_all$event <- as.factor(number_events_all$event)
  number_events_all$event <- factor(number_events_all$event, levels = c("cospeciation", "transfer", "duplication", "loss"))
  
  pdf(paste0("../result_plots/emPress_bifurcations_number_events_", costs,".pdf"), width = 4.25, height = 2.5)
  
  print(ggplot(number_events_all[which(number_events_all$signif==TRUE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[which(number_events_all$signif==FALSE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d2b4de")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="1")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="2")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="3")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+scale_y_sqrt()+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="4")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  dev.off()
  
  write.table(list_pvalues, paste0("../result_plots/results_bifurcations_eMPRess_", costs,".csv"), sep=";", row.names=FALSE)
  
  
}

res_stats <- data.frame(res_stats)
colnames(res_stats) <- c("costs", "type", 1:8)
write.table(res_stats, paste0("../result_plots/results_stats_emPress_trait_bifurcations.csv"), sep=";", row.names=FALSE)






#### Results for vicariance  ####

rm(list=ls())

library(ggplot2)
library(ape)

setwd("/path/emPress/results_geo/")
path_data="../../simulated_data_geo/"
nsim=250

res_stats <- c()

for (costs in c("d1_t1_l1", "d4_t1_l1", "d2_t1_l2", "d4_t2_l1",  "d2_t3_l1")){
  list_pvalues <- c()
  print(costs)
  for (name in c("geo_1", "geo_2", "geo_3", "geo_4", "geo_5", "geo_6", "geo_7", "geo_8")){
    for (seed in 1:nsim){
      
      if (file.exists(paste0("res_",costs,"_simul_",name,"_",seed,".svg"))){
        res <- read.table(paste0("res_",costs,"_simul_",name,"_",seed,".svg"), comment.char = "", fill=TRUE, sep=";")
        
        res_pvalue <- res$V1[grep("p-value", res$V1)]
        res_pvalue <- gsub("    <!-- p-value = ", "", res_pvalue)
        res_pvalue <- as.numeric(gsub(" -->", "", res_pvalue))
        
        recon <- read.table(paste0("reconciliation_",costs,"_simul_",name,"_",seed,".csv"), sep=",")
        colnames(recon) <- c("parasite", "host", "event", "node_frequency", "event_frequency")
        
        if (!all(recon$node_frequency>=recon$event_frequency)){print("problem")} # makes sense
        
        
        # check the time consistency (check for back-in-time transfers)
        
        tree_hosts <- read.tree(paste0(path_data, "tree_hosts_only_int_simul_",name,"_",seed,".tre"))
        tree_parasites <- read.tree(paste0(path_data, "tree_parasites_only_int_simul_",name,"_",seed,".tre"))
        
        for (i in 1:nrow(recon)){
          if (length(grep("_",recon$parasite[i]))==1){
            recon$parasite[i] <- as.numeric(gsub("_p", "", recon$parasite[i]))+Ntip(tree_parasites)+1
          } else {
            recon$parasite[i] <- which(tree_parasites$tip.label==recon$parasite[i])
          }
          if (length(grep("_",recon$host[i]))==1){
            recon$host[i] <- as.numeric(gsub("_h", "", recon$host[i]))+Ntip(tree_hosts)+1
          } else {
            recon$host[i] <- which(tree_hosts$tip.label==recon$host[i])
          }
        }
        
        node_age_hosts <- node.depth.edgelength(tree_hosts)
        node_age_hosts <- max(node_age_hosts)-node_age_hosts
        names(node_age_hosts) <- 1:length(node_age_hosts)
        
        recon$age_min <- NA
        recon$age_max <- NA
        
        for (i in 1:nrow(recon)){
          
          if (recon$event[i]=="Contemporaneous"){
            recon$age_min[i] <- 0
            recon$age_max[i] <- 0
          }
          
          if (recon$event[i]=="Cospeciation"){
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            recon$age_max[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Transfer"){
            
            # the most ancient age of departure (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent age of arrival
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Loss"){
            
            # the most ancient age of loss (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent loss 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
          
          if (recon$event[i]=="Duplication"){
            
            # the most ancient age of duplication (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent duplication (ideally, it should be constrained with next nodes, but should not affect the measure of time consistency) 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
        }
        
        recon$diff_age <- recon$age_max-recon$age_min
        
        list_pvalues <- rbind(list_pvalues, c(name, seed, res_pvalue, table(recon$event)[c("Cospeciation", "Transfer", "Duplication", "Loss")], length(which(recon$diff_age<0)), length(which(recon$diff_age[recon$event=="Transfer"]<0)) ))
      }
    }
  }
  
  list_pvalues <- data.frame(list_pvalues)
  colnames(list_pvalues) <- c("name", "seed", "pvalue", "cospeciation", "transfer", "duplication", "loss", "time_consistency", "time_consistency_transfers")
  print(nrow(list_pvalues))
  list_pvalues$pvalue <- as.numeric(as.character(list_pvalues$pvalue))
  list_pvalues$cospeciation <- as.numeric(as.character(list_pvalues$cospeciation))
  list_pvalues$transfer <- as.numeric(as.character(list_pvalues$transfer))
  
  print(costs)
  print("permutations:")
  table <- table(list_pvalues$pvalue<0.05, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "permutations",table[nrow(table),]))
  
  
  print("time consistency:")
  list_pvalues$signif_time_consistency <- FALSE
  list_pvalues$signif_time_consistency[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency",table[nrow(table),]))
  
  print("number cospe/switches:")
  list_pvalues$signif_switches <- FALSE
  list_pvalues$signif_switches[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "cospe_switches",table[nrow(table),]))
  
  print("time consistency - number cospe/switches:")
  list_pvalues$signif_time_consistency_switches <- FALSE
  list_pvalues$signif_time_consistency_switches[intersect(intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0)), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency_switches",table[nrow(table),]))
  
  number_events_1 <- cbind(list_pvalues[,c(1,2,3,4)], rep("cospeciation", nrow(list_pvalues)))
  colnames(number_events_1) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_2 <- cbind(list_pvalues[,c(1,2,3,5)], rep("transfer", nrow(list_pvalues)))
  colnames(number_events_2) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_3 <- cbind(list_pvalues[,c(1,2,3,6)], rep("duplication", nrow(list_pvalues)))
  colnames(number_events_3) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_4 <- cbind(list_pvalues[,c(1,2,3,7)], rep("loss", nrow(list_pvalues)))
  colnames(number_events_4) <- c("name", "seed", "pvalue", "number_events", "event")
  
  number_events_all <- rbind(number_events_1,number_events_2, number_events_3, number_events_4)
  
  number_events_all <- data.frame(number_events_all)
  colnames(number_events_all) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_all$pvalue <- as.numeric(as.character(number_events_all$pvalue))
  number_events_all$number_events <- as.numeric(as.character(number_events_all$number_events))
  number_events_all$signif <- number_events_all$pvalue<0.05
  number_events_all$number_events[is.na(number_events_all$number_events)] <- 0
  number_events_all$name_seed <- paste0(number_events_all$name, "_", number_events_all$seed)
  
  number_events_all$event <- as.factor(number_events_all$event)
  number_events_all$event <- factor(number_events_all$event, levels = c("cospeciation", "transfer", "duplication", "loss"))
  
  pdf(paste0("../result_plots/emPress_geo_number_events_", costs,".pdf"), width = 4.25, height = 2.5)
  
  print(ggplot(number_events_all[which(number_events_all$signif==TRUE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[which(number_events_all$signif==FALSE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d2b4de")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_1")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_2")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_3")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+scale_y_sqrt()+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_4")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  dev.off()
  
  write.table(list_pvalues, paste0("../result_plots/results_emPress_geo_", costs,".csv"), sep=";", row.names=FALSE)
  
}

res_stats <- data.frame(res_stats)
colnames(res_stats) <- c("costs", "type", 1:8)
write.table(res_stats, paste0("../result_plots/results_stats_emPress_vicariance.csv"), sep=";", row.names=FALSE)




#### Results for vicariance - random bifurcations  ####

rm(list=ls())

library(ggplot2)
library(ape)

setwd("/path/emPress/results_geo/")
path_data="../../simulated_data_geo/"
nsim=250

res_stats <- c()

for (costs in c("d1_t1_l1", "d4_t1_l1", "d2_t1_l2", "d4_t2_l1",  "d2_t3_l1")){
  list_pvalues <- c()
  print(costs)
  for (name in c("geo_1", "geo_2", "geo_3", "geo_4", "geo_5", "geo_6", "geo_7", "geo_8")){
    for (seed in 1:nsim){
      
      if (file.exists(paste0("res_bifurcations_eMPRess_",costs,"_simul_",name,"_",seed,".svg"))){
        res <- read.table(paste0("res_bifurcations_eMPRess_",costs,"_simul_",name,"_",seed,".svg"), comment.char = "", fill=TRUE, sep=";")
        
        res_pvalue <- res$V1[grep("p-value", res$V1)]
        res_pvalue <- gsub("    <!-- p-value = ", "", res_pvalue)
        res_pvalue <- as.numeric(gsub(" -->", "", res_pvalue))
        
        recon <- read.table(paste0("reconciliation_bifurcations_eMPRess_",costs,"_simul_",name,"_",seed,".csv"), sep=",")
        colnames(recon) <- c("parasite", "host", "event", "node_frequency", "event_frequency")
        
        if (!all(recon$node_frequency>=recon$event_frequency)){print("problem")} # makes sense
        
        
        # check the time consistency (check for back-in-time transfers)
        
        tree_hosts <- read.tree(paste0(path_data, "tree_hosts_only_int_simul_",name,"_",seed,".tre"))
        tree_parasites <- read.tree(paste0(path_data, "tree_parasites_bifurcations_eMPRess_simul_",name,"_",seed,".tre"))
        
        for (i in 1:nrow(recon)){
          if (length(grep("_",recon$parasite[i]))==1){
            recon$parasite[i] <- as.numeric(gsub("_p", "", recon$parasite[i]))+Ntip(tree_parasites)+1
          } else {
            recon$parasite[i] <- which(tree_parasites$tip.label==recon$parasite[i])
          }
          if (length(grep("_",recon$host[i]))==1){
            recon$host[i] <- as.numeric(gsub("_h", "", recon$host[i]))+Ntip(tree_hosts)+1
          } else {
            recon$host[i] <- which(tree_hosts$tip.label==recon$host[i])
          }
        }
        
        node_age_hosts <- node.depth.edgelength(tree_hosts)
        node_age_hosts <- max(node_age_hosts)-node_age_hosts
        names(node_age_hosts) <- 1:length(node_age_hosts)
        
        recon$age_min <- NA
        recon$age_max <- NA
        
        for (i in 1:nrow(recon)){
          
          if (recon$event[i]=="Contemporaneous"){
            recon$age_min[i] <- 0
            recon$age_max[i] <- 0
          }
          
          if (recon$event[i]=="Cospeciation"){
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            recon$age_max[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Transfer"){
            
            # the most ancient age of departure (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent age of arrival
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Loss"){
            
            # the most ancient age of loss (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent loss 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
          
          if (recon$event[i]=="Duplication"){
            
            # the most ancient age of duplication (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent duplication (ideally, it should be constrained with next nodes, but should not affect the measure of time consistency) 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
        }
        
        recon$diff_age <- recon$age_max-recon$age_min
        
        list_pvalues <- rbind(list_pvalues, c(name, seed, res_pvalue, table(recon$event)[c("Cospeciation", "Transfer", "Duplication", "Loss")], length(which(recon$diff_age<0)), length(which(recon$diff_age[recon$event=="Transfer"]<0)) ))
      }
    }
  }
  
  list_pvalues <- data.frame(list_pvalues)
  colnames(list_pvalues) <- c("name", "seed", "pvalue", "cospeciation", "transfer", "duplication", "loss", "time_consistency", "time_consistency_transfers")
  print(nrow(list_pvalues))
  list_pvalues$pvalue <- as.numeric(as.character(list_pvalues$pvalue))
  list_pvalues$cospeciation <- as.numeric(as.character(list_pvalues$cospeciation))
  list_pvalues$transfer <- as.numeric(as.character(list_pvalues$transfer))
  
  print(costs)
  print("permutations:")
  table <- table(list_pvalues$pvalue<0.05, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "permutations",table[nrow(table),]))
  
  print("time consistency:")
  list_pvalues$signif_time_consistency <- FALSE
  list_pvalues$signif_time_consistency[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency",table[nrow(table),]))
  
  print("number cospe/switches:")
  list_pvalues$signif_switches <- FALSE
  list_pvalues$signif_switches[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "cospe_switches",table[nrow(table),]))
  
  print("time consistency - number cospe/switches:")
  list_pvalues$signif_time_consistency_switches <- FALSE
  list_pvalues$signif_time_consistency_switches[intersect(intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0)), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency_switches",table[nrow(table),]))
  
  number_events_1 <- cbind(list_pvalues[,c(1,2,3,4)], rep("cospeciation", nrow(list_pvalues)))
  colnames(number_events_1) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_2 <- cbind(list_pvalues[,c(1,2,3,5)], rep("transfer", nrow(list_pvalues)))
  colnames(number_events_2) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_3 <- cbind(list_pvalues[,c(1,2,3,6)], rep("duplication", nrow(list_pvalues)))
  colnames(number_events_3) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_4 <- cbind(list_pvalues[,c(1,2,3,7)], rep("loss", nrow(list_pvalues)))
  colnames(number_events_4) <- c("name", "seed", "pvalue", "number_events", "event")
  
  number_events_all <- rbind(number_events_1,number_events_2, number_events_3, number_events_4)
  
  number_events_all <- data.frame(number_events_all)
  colnames(number_events_all) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_all$pvalue <- as.numeric(as.character(number_events_all$pvalue))
  number_events_all$number_events <- as.numeric(as.character(number_events_all$number_events))
  number_events_all$signif <- number_events_all$pvalue<0.05
  number_events_all$number_events[is.na(number_events_all$number_events)] <- 0
  number_events_all$name_seed <- paste0(number_events_all$name, "_", number_events_all$seed)
  
  number_events_all$event <- as.factor(number_events_all$event)
  number_events_all$event <- factor(number_events_all$event, levels = c("cospeciation", "transfer", "duplication", "loss"))
  
  pdf(paste0("../result_plots/emPress_bifurcations_geo_number_events_", costs,".pdf"), width = 4.25, height = 2.5)
  
  print(ggplot(number_events_all[which(number_events_all$signif==TRUE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[which(number_events_all$signif==FALSE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d2b4de")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_1")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_2")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_3")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+scale_y_sqrt()+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_4")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  dev.off()
  
  write.table(list_pvalues, paste0("../result_plots/results_bifurcations_emPress_geo_", costs,".csv"), sep=";", row.names=FALSE)
  
}

res_stats <- data.frame(res_stats)
colnames(res_stats) <- c("costs", "type", 1:8)
write.table(res_stats, paste0("../result_plots/results_stats_emPress_vicariance_bifurcation.csv"), sep=";", row.names=FALSE)




#### Results for codiv  ####

rm(list=ls())

library(ggplot2)
library(ape)

setwd("/path/emPress/results_codiv/")
path_data="../../simulated_data_codiv/"
nsim=250

res_stats <- c()

for (costs in c("d1_t1_l1", "d4_t1_l1", "d2_t1_l2", "d4_t2_l1",  "d2_t3_l1")){
  list_pvalues <- c()
  print(costs)
  for (name in c("codiv_1", "codiv_2", "codiv_3", "codiv_4")){
    for (seed in 1:nsim){
      
      if (file.exists(paste0("res_",costs,"_simul_",name,"_",seed,".svg"))){
        res <- read.table(paste0("res_",costs,"_simul_",name,"_",seed,".svg"), comment.char = "", fill=TRUE, sep=";")
        
        res_pvalue <- res$V1[grep("p-value", res$V1)]
        res_pvalue <- gsub("    <!-- p-value = ", "", res_pvalue)
        res_pvalue <- as.numeric(gsub(" -->", "", res_pvalue))
        
        recon <- read.table(paste0("reconciliation_",costs,"_simul_",name,"_",seed,".csv"), sep=",")
        colnames(recon) <- c("parasite", "host", "event", "node_frequency", "event_frequency")
        
        if (!all(recon$node_frequency>=recon$event_frequency)){print("problem")} # makes sense
        
        
        # check the time consistency (check for back-in-time transfers)
        
        tree_hosts <- read.tree(paste0(path_data, "tree_hosts_only_int_simul_",name,"_",seed,".tre"))
        tree_parasites <- read.tree(paste0(path_data, "tree_parasites_only_int_simul_",name,"_",seed,".tre"))
        
        for (i in 1:nrow(recon)){
          if (length(grep("_",recon$parasite[i]))==1){
            recon$parasite[i] <- as.numeric(gsub("_p", "", recon$parasite[i]))+Ntip(tree_parasites)+1
          } else {
            recon$parasite[i] <- which(tree_parasites$tip.label==recon$parasite[i])
          }
          if (length(grep("_",recon$host[i]))==1){
            recon$host[i] <- as.numeric(gsub("_h", "", recon$host[i]))+Ntip(tree_hosts)+1
          } else {
            recon$host[i] <- which(tree_hosts$tip.label==recon$host[i])
          }
        }
        
        node_age_hosts <- node.depth.edgelength(tree_hosts)
        node_age_hosts <- max(node_age_hosts)-node_age_hosts
        names(node_age_hosts) <- 1:length(node_age_hosts)
        
        recon$age_min <- NA
        recon$age_max <- NA
        
        for (i in 1:nrow(recon)){
          
          if (recon$event[i]=="Contemporaneous"){
            recon$age_min[i] <- 0
            recon$age_max[i] <- 0
          }
          
          if (recon$event[i]=="Cospeciation"){
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            recon$age_max[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Transfer"){
            
            # the most ancient age of departure (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent age of arrival
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
          }
          
          if (recon$event[i]=="Loss"){
            
            # the most ancient age of loss (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent loss 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
          
          if (recon$event[i]=="Duplication"){
            
            # the most ancient age of duplication (age of the previous node in the parasite tree)
            previous_node <- which(tree_parasites$edge[,2]==recon$parasite[i])
            if (length(previous_node)==1){
              recon$age_max[i] <- max(recon$age_max[which(recon$parasite==tree_parasites$edge[previous_node,1])])
            }else{
              if (length(which(tree_hosts$edge[,2]==recon$host[i]))==1){
                recon$age_max[i] <- node_age_hosts[tree_hosts$edge[which(tree_hosts$edge[,2]==recon$host[i]),1]]
              } else {
                recon$age_max[i] <- max(node_age_hosts)
              }
            }
            
            # the most recent duplication (ideally, it should be constrained with next nodes, but should not affect the measure of time consistency) 
            recon$age_min[i] <- node_age_hosts[recon$host[i]]
            
          }
        }
        
        recon$diff_age <- recon$age_max-recon$age_min
        
        list_pvalues <- rbind(list_pvalues, c(name, seed, res_pvalue, table(recon$event)[c("Cospeciation", "Transfer", "Duplication", "Loss")], length(which(recon$diff_age<0)), length(which(recon$diff_age[recon$event=="Transfer"]<0)) ))
      }
    }
  }
  
  list_pvalues <- data.frame(list_pvalues)
  colnames(list_pvalues) <- c("name", "seed", "pvalue", "cospeciation", "transfer", "duplication", "loss", "time_consistency", "time_consistency_transfers")
  print(nrow(list_pvalues))
  list_pvalues$pvalue <- as.numeric(as.character(list_pvalues$pvalue))
  list_pvalues$cospeciation <- as.numeric(as.character(list_pvalues$cospeciation))
  list_pvalues$transfer <- as.numeric(as.character(list_pvalues$transfer))
  list_pvalues$time_consistency_transfers <- as.numeric(as.character(list_pvalues$time_consistency_transfers))
  
  
  print(costs)
  print("permutations:")
  table <- table(list_pvalues$pvalue<0.05, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "permutations",table[nrow(table),]))
  
  
  print("time consistency:")
  list_pvalues$signif_time_consistency <- FALSE
  list_pvalues$signif_time_consistency[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency",table[nrow(table),]))
  
  print("number cospe/switches:")
  list_pvalues$signif_switches <- FALSE
  list_pvalues$signif_switches[intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "cospe_switches",table[nrow(table),]))
  
  print("time consistency - number cospe/switches:")
  list_pvalues$signif_time_consistency_switches <- FALSE
  list_pvalues$signif_time_consistency_switches[intersect(intersect(which(list_pvalues$pvalue<0.05), which(list_pvalues$time_consistency_transfers==0)), which(list_pvalues$cospeciation>list_pvalues$transfer))] <- TRUE
  table <- table(list_pvalues$signif_time_consistency_switches, list_pvalues$name)
  table <- print(t(round(t(table)/colSums(table),2)))
  res_stats <- rbind(res_stats, c(costs, "time_consistency_switches",table[nrow(table),]))
  
  number_events_1 <- cbind(list_pvalues[,c(1,2,3,4)], rep("cospeciation", nrow(list_pvalues)))
  colnames(number_events_1) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_2 <- cbind(list_pvalues[,c(1,2,3,5)], rep("transfer", nrow(list_pvalues)))
  colnames(number_events_2) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_3 <- cbind(list_pvalues[,c(1,2,3,6)], rep("duplication", nrow(list_pvalues)))
  colnames(number_events_3) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_4 <- cbind(list_pvalues[,c(1,2,3,7)], rep("loss", nrow(list_pvalues)))
  colnames(number_events_4) <- c("name", "seed", "pvalue", "number_events", "event")
  
  number_events_all <- rbind(number_events_1,number_events_2, number_events_3, number_events_4)
  
  number_events_all <- data.frame(number_events_all)
  colnames(number_events_all) <- c("name", "seed", "pvalue", "number_events", "event")
  number_events_all$pvalue <- as.numeric(as.character(number_events_all$pvalue))
  number_events_all$number_events <- as.numeric(as.character(number_events_all$number_events))
  number_events_all$signif <- number_events_all$pvalue<0.05
  number_events_all$number_events[is.na(number_events_all$number_events)] <- 0
  number_events_all$name_seed <- paste0(number_events_all$name, "_", number_events_all$seed)
  
  number_events_all$event <- as.factor(number_events_all$event)
  number_events_all$event <- factor(number_events_all$event, levels = c("cospeciation", "transfer", "duplication", "loss"))
  
  pdf(paste0("../result_plots/emPress_codiv_number_events_", costs,".pdf"), width = 4.25, height = 2.5)
  
  #scale_y_sqrt()+
  print(ggplot(number_events_all[which(number_events_all$signif==TRUE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[which(number_events_all$signif==FALSE),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d2b4de")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          #geom_boxplot(width=0.1,outlier.shape = NA)+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_1")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_2")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_3")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+scale_y_sqrt()+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  print(ggplot(number_events_all[intersect(which(number_events_all$signif==TRUE), which(number_events_all$name=="geo_4")),], aes(x=event, y=number_events,fill=event)) + xlab("")+ylab("Number of events")+
          geom_line(aes(group=name_seed),alpha = 0.5, color="#d7dbdd")+geom_hline(yintercept=0)+
          theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954","#2e86c1"))+theme(legend.position="none")+
          geom_boxplot(width=0.1,outlier.shape = NA)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
  
  dev.off()
  
  write.table(list_pvalues, paste0("../result_plots/results_emPress_codiv_", costs,".csv"), sep=";", row.names=FALSE)
  
}

res_stats <- data.frame(res_stats)
colnames(res_stats) <- c("costs", "type", 1:4)
write.table(res_stats, paste0("../result_plots/results_stats_emPress_codiv.csv"), sep=";", row.names=FALSE)




#####  Plot the ratio  of the number of events ######

rm(list=ls())

setwd("/path/emPress/result_plots/")

library(ggplot2)

type="all"
type_2 <- "trait matching"

res_all <- c()

for (type in c("", "geo_", "codiv_")){
    
    if (type=="") type_2 <- "trait matching"
    if (type=="geo_") type_2 <- "vicariance"
    if (type=="codiv_") type_2 <- "codiversification"
    
    
    print(type)
    print(type_2)
    
    res_empress <- read.table(paste0("results_emPress_",type,"d4_t1_l1.csv"), header=TRUE, sep=";")

    res_empress$type <- type_2
    res_all <- rbind(res_all, res_empress)
}

res_all$transfer[is.na(res_all$transfer)] <- 0

res_all$ratio <- res_all$transfer/res_all$cospeciation

res_all$signif <- res_all$pvalue<0.05

res_all$type_res_empress <- paste0(res_all$type,"_", res_all$signif)
res_all$type_res_empress <- as.factor(res_all$type_res_empress)
res_all$type_res_empress <- factor(res_all$type_res_empress, levels = levels(res_all$type_res_empress)[c(4,3,6,5,2,1)])

pdf("summary_Empress_stats.pdf", width=3, height=2.75)
ggplot(res_all, aes(y=ratio, x=type_res_empress, fill=type_res_empress))+ geom_hline(yintercept = 1)+geom_violin(trim=TRUE, draw_quantiles=0.5)+
  theme_bw()+xlab("")+ ylab("Ratio transfers vs. cospeciations")+scale_y_continuous(trans='sqrt', breaks=c(0, 0.25, 1, 2, 5, 10, 20, 30))+
  scale_fill_manual(values=c("#b03a2e","#1a5276","#b03a2e","#1a5276","#b03a2e","#1a5276"))+ theme(legend.position = "none")
dev.off()




#####  Plot the ratio  of the number of events - random bifurcations ######

rm(list=ls())

setwd("/path/emPress/result_plots/")

library(ggplot2)

type="all"
type_2 <- "trait matching"

res_all <- c()

for (type in c("", "geo_")){
  
  if (type=="") type_2 <- "trait matching"
  if (type=="geo_") type_2 <- "vicariance"
  
  print(type)
  print(type_2)
  
  res_empress <- read.table(paste0("results_bifurcations_emPress_",type,"d4_t1_l1.csv"), header=TRUE, sep=";")
  
  res_empress$type <- type_2
  res_all <- rbind(res_all, res_empress)
}

res_all$transfer[is.na(res_all$transfer)] <- 0

res_all$ratio <- res_all$transfer/res_all$cospeciation

res_all$signif <- res_all$pvalue<0.05

res_all$type_res_empress <- paste0(res_all$type,"_", res_all$signif)
res_all$type_res_empress <- as.factor(res_all$type_res_empress)
res_all$type_res_empress <- factor(res_all$type_res_empress, levels = levels(res_all$type_res_empress)[c(4,3,6,5,2,1)])

pdf("summary_Empress_stats_bifurcations.pdf", width=2.5, height=2.75)
ggplot(res_all, aes(y=ratio, x=type_res_empress, fill=type_res_empress))+ geom_hline(yintercept = 1)+geom_violin(trim=TRUE, draw_quantiles=0.5)+
  theme_bw()+xlab("")+ ylab("Ratio transfers vs. cospeciations")+scale_y_continuous(trans='sqrt', breaks=c(0, 0.25, 1, 2, 5, 10, 20, 30))+
  scale_fill_manual(values=c("#b03a2e","#1a5276","#b03a2e","#1a5276","#b03a2e","#1a5276"))+ theme(legend.position = "none")
dev.off()

