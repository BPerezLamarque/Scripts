#### Step 1: R script - Prepare the files for running empress on simulations ####

rm(list=ls())

library(ape)
library(phytools)
library(vegan)
library(mvMORPH)
library(reshape2)


name='simul_1'

for (name in c('simul_1', 'simul_2', 'simul_3', 'simul_4', 'simul_5')){
  
  
  source("../functions_biogeo.R")
  dir.create(file.path(getwd(), "data/random_biogeo_empress/"), showWarnings = FALSE)
  
  
  nb_rep <- 200 
  seed=1
  
  list_areas <- c("A","B","C","D","E")
  
  for (seed in 1:nb_rep){
    
    print(seed)
    
    tree_A= read.tree(paste0("data/tree_hosts_",name,"_",seed,".tre"))
    tree_B = read.tree(paste0("data/tree_parasites_",name,"_",seed,".tre"))
    
    # Hosts on columns and parasites on rows
    network = t(read.table(paste0("data/network_",name,"_",seed,".csv"),sep=";",header= TRUE))
    
    bgeo <- rep("A",ncol(network))
    for (biogeo in list_areas){
      bgeo[grep(biogeo,colnames(network))]=biogeo
    }
    table_biogeo = data.frame(species = colnames(network), biogeo = bgeo)
    
    # Prepare for empress (tree_B on rows first and tree_A on columns in second)
    list_links <- reshape2::melt(network)
    list_links <- list_links[list_links$value>0,]
    write.table(paste0(list_links$Var1,":", list_links$Var2), paste0("data/links_", name, "_", seed, ".txt"), col.names=F, row.names=F, quote=F)
    
    
    # Prepare the H0 biogeo for empress
    nperm=100
    for (perm in 1:nperm){
      set.seed(perm)
      rand_network <- network
      for (zone in unique(table_biogeo$biogeo)){
        w <- which(table_biogeo$biogeo==zone)
        s <- as.numeric(sample(as.character(w))) 
        rand_network[,w] <- as.matrix(rand_network[,s])
      }
      list_links <- reshape2::melt(rand_network)
      list_links <- list_links[list_links$value>0,]
      write.table(paste0(list_links$Var1,":", list_links$Var2), paste0("data/random_biogeo_empress/links_random_biogeo_", name, "_", seed, "_", perm,".txt"), col.names=F, row.names=F, quote=F)
    }
    
  }
}

#### Step 2: Bash script - Run eMPRess #####


#!/bin/bash

# Environment setup

module load python/3.7.10


simul="simul_1"


data_path="/PATH/"


d_rate=1
t_rate=1
l_rate=1


# H0 

for seed in $(shuf -i 1-200)
  do
    if test -f $data_path"/empress/res_d"$d_rate"_t"$t_rate"_l"$l_rate"_"$simul"_"$seed".svg"; then
      echo $seed
    else
      echo $seed
      python empress_cli.py p-value $data_path"/data/tree_hosts_"$simul"_"$seed".tre" $data_path"/data/tree_parasites_"$simul"_"$seed".tre" $data_path"/data/links_"$simul"_"$seed".txt" -d $d_rate -t $t_rate -l $l_rate --n 100 --outfile $data_path"/empress/res_d"$d_rate"_t"$t_rate"_l"$l_rate"_"$simul"_"$seed".svg"
    fi
    if test -f $data_path"/empress/reconciliation_d"$d_rate"_t"$t_rate"_l"$l_rate"_"$simul"_"$seed".csv"; then
      echo $seed
    else
      echo $seed
      python empress_cli.py reconcile $data_path"/data/tree_hosts_"$simul"_"$seed".tre" $data_path"/data/tree_parasites_"$simul"_"$seed".tre" $data_path"/data/links_"$simul"_"$seed".txt" -d $d_rate -t $t_rate -l $l_rate --csv $data_path"/empress/reconciliation_d"$d_rate"_t"$t_rate"_l"$l_rate"_"$simul"_"$seed".csv"
    fi
done



# H0 biogeo

for seed in $(shuf -i 1-200)
  do
    echo $seed
    for perm in {1..100}
    do
      if test -f $data_path"/empress/res_d"$d_rate"_t"$t_rate"_l"$l_rate"_"$simul"_"$seed"_biogeo_"$perm".svg"; then
      echo $seed
    else
      if test -f $data_path"/empress/cost_d"$d_rate"_t"$t_rate"_l"$l_rate"_"$simul"_"$seed"_biogeo_"$perm".csv"; then
        echo $seed
      else
        python empress_cli.py p-value $data_path"/data/tree_hosts_"$simul"_"$seed".tre" $data_path"/data/tree_parasites_"$simul"_"$seed".tre" $data_path"/data/random_biogeo_empress/links_random_biogeo_"$simul"_"$seed"_"$perm".txt" -d $d_rate -t $t_rate -l $l_rate --n 1 --outfile $data_path"/empress/res_d"$d_rate"_t"$t_rate"_l"$l_rate"_"$simul"_"$seed"_biogeo_"$perm".svg"
      fi
    fi
  done
done




#### Step 3: R script - Extract the MPR cost ####


rm(list=ls())

library(ape)
library(phytools)
library(vegan)
library(mvMORPH)
library(reshape2)
library(xml2)


extract_cost_svg <- function(svg_file){
  
  svg <- tryCatch(
    read_xml(svg_file),
    error = function(e) {
      message("Error reading ", svg_file, ": ", e$message)
      file.remove(svg_file)
      return(NULL)
    }
  )
  
  if (is.null(svg)){
    return(NA)
  }
  
  # 1. Extract the x-coordinate of the red vertical line
  paths <- xml_find_all(svg, "//*[local-name()='path']")
  
  red_paths <- paths[
    grepl("#ff0000",
          paste(xml_attr(paths, "style"), xml_attr(paths, "stroke")),
          ignore.case = TRUE)
  ]
  
  if (length(red_paths) == 0) {
    stop("No red path found in the SVG.")
  }
  
  # Extract x-coordinate from the first M command in the path
  red_d <- xml_attr(red_paths[[1]], "d")
  position <- as.numeric(
    sub("^\\s*M\\s*([-+0-9.eE]+).*", "\\1", red_d)
  )
  
  # 2. Extract x-axis tick values and their actual SVG x-coordinates
  ticks <- xml_find_all(
    svg,
    "//*[local-name()='g' and starts-with(@id, 'xtick_')]"
  )
  
  tick_values <- numeric()
  tick_positions <- numeric()
  
  for (tick in ticks) {
    
    # Tick label is stored as an XML comment (e.g. <!-- 157.5 -->)
    label_node <- xml_find_first(tick, ".//comment()")
    value <- as.numeric(trimws(xml_text(label_node)))
    
    # Actual tick position
    use_node <- xml_find_first(tick, ".//*[local-name()='use']")
    x <- as.numeric(xml_attr(use_node, "x"))
    
    if (!is.na(value) && !is.na(x)) {
      tick_values <- c(tick_values, value)
      tick_positions <- c(tick_positions, x)
    }
  }
  
  if (length(tick_values) < 2) {
    stop("Could not extract at least two x-axis ticks.")
  }
  
  # 3. Convert the red line position into the corresponding data value
  model <- lm(tick_values ~ tick_positions)
  
  total_cost_biogeo <- round(
    predict(model, newdata = data.frame(tick_positions = position))
  )
  return(total_cost_biogeo)
}


name='simul_1'


for (name in sample(c('simul_1', 'simul_2', 'simul_3', 'simul_4', 'simul_5'))){
  
  print(paste0("name: ", name))
  
  
  nb_rep <- 200 
  seed=1
  perm=1
  
  for (seed in sample(1:nb_rep)){
    
    print(paste0("seed: ", seed))
    
    # Cost for H0 
    
    svg_file <- paste0("empress/res_d1_t1_l1_", name, "_", seed,".svg")
    total_cost <- extract_cost_svg(svg_file)
    
    
    # Cost for H0 biogeo 
    for (perm in 1:100){
      if (!file.exists(paste0("empress/cost_d1_t1_l1_", name, "_", seed, "_biogeo_", perm,".csv"))){
        print(perm)
        if (file.exists(paste0("empress/res_d1_t1_l1_", name, "_", seed, "_biogeo_", perm,".svg"))){
          
          svg_file <- paste0("empress/res_d1_t1_l1_", name, "_", seed,"_biogeo_", perm, ".svg")
          
          total_cost_biogeo <- extract_cost_svg(svg_file)
          
          if (!is.na(total_cost_biogeo)){
            res <- c(name, seed, perm, total_cost, total_cost_biogeo)
            
            write.table(res, paste0("empress/cost_d1_t1_l1_", name, "_", seed, "_biogeo_", perm,".csv"), sep=";", row.names=FALSE, quote=FALSE, col.names=FALSE)
          }
          
        } else {
          print(paste0("Issue with not existing file: ", paste0("empress/res_d1_t1_l1_", name, "_", seed, "_biogeo_", perm,".svg")))
        }
      }
    }
  }
}








#### Step 4: R script - Make final file eMPRess ####


rm(list=ls())

library(ape)
library(phytools)
library(vegan)
library(mvMORPH)
library(reshape2)
library(xml2)



library(xml2)

extract_blue_distribution <- function(svg_file) {
  
  svg <- read_xml(svg_file)
  
  # Extract tick values and their actual SVG coordinates
  extract_ticks <- function(axis = "x") {
    
    tick_groups <- xml_find_all(
      svg,
      paste0(
        "//*[local-name()='g' and starts-with(@id, '",
        axis, "tick_')]"
      )
    )
    
    ticks <- lapply(tick_groups, function(tick) {
      
      # Tick label stored as an XML comment
      label_node <- xml_find_first(tick, ".//comment()")
      value <- as.numeric(trimws(xml_text(label_node)))
      
      # Actual tick coordinate
      use_node <- xml_find_first(
        tick, ".//*[local-name()='use']"
      )
      
      coord <- as.numeric(
        xml_attr(use_node, axis)
      )
      
      data.frame(value = value, coord = coord)
    })
    
    ticks <- do.call(rbind, ticks)
    ticks <- ticks[complete.cases(ticks), ]
    
    if (nrow(ticks) < 2) {
      stop("Could not extract at least two ", axis, "-axis ticks.")
    }
    
    ticks
  }
  
  x_ticks <- extract_ticks("x")
  y_ticks <- extract_ticks("y")
  
  # Linear relationships between SVG coordinates and axis values
  x_model <- lm(value ~ coord, data = x_ticks)
  y_model <- lm(value ~ coord, data = y_ticks)
  
  # Extract blue histogram bars
  paths <- xml_find_all(
    svg, "//*[local-name()='path']"
  )
  
  styles <- xml_attr(paths, "style")
  
  blue_paths <- paths[
    grepl("#0000ff", styles, ignore.case = TRUE)
  ]
  
  if (length(blue_paths) == 0) {
    stop("No blue histogram bars found.")
  }
  
  # Extract coordinates from each histogram rectangle
  bars <- lapply(blue_paths, function(path) {
    
    d <- xml_attr(path, "d")
    
    coords <- as.numeric(
      regmatches(
        d,
        gregexpr(
          "[-+]?(?:[0-9]*\\.)?[0-9]+(?:[eE][-+]?[0-9]+)?",
          d,
          perl = TRUE
        )
      )[[1]]
    )
    
    if (length(coords) < 8) {
      return(NULL)
    }
    
    # SVG path: M x0 y0 L x1 y0 L x1 ytop L x0 ytop
    x0 <- coords[1]
    baseline <- coords[2]
    x1 <- coords[3]
    y_top <- coords[6]
    
    # Convert SVG coordinates into data values
    x_center <- (x0 + x1) / 2
    
    x_value <- predict(
      x_model,
      newdata = data.frame(coord = x_center)
    )
    
    count <- predict(
      y_model,
      newdata = data.frame(coord = y_top)
    ) -
      predict(
        y_model,
        newdata = data.frame(coord = baseline)
      )
    
    data.frame(
      xmin = predict(x_model, newdata = data.frame(coord = x0)),
      xmax = predict(x_model, newdata = data.frame(coord = x1)),
      center = as.numeric(x_value),
      count = max(0, round(count))
    )
  })
  
  bars <- do.call(rbind, bars)
  
  # Remove empty classes if desired: here we retain them
  bars <- bars[order(bars$xmin), ]
  
  # Calculate weighted mean and sample standard deviation
  total_n <- sum(bars$count)
  
  if (total_n < 2) {
    stop("Not enough observations to calculate the standard deviation.")
  }
  
  mean_distribution <- weighted.mean(
    bars$center, bars$count
  )
  
  sd_distribution <- sqrt(
    sum(bars$count * (bars$center - mean_distribution)^2) /
      (total_n - 1)
  )
  
  list(
    distribution = bars,
    n = total_n,
    mean = mean_distribution,
    sd = sd_distribution
  )
}




name='simul_1'



for (name in c('simul_1', 'simul_2', 'simul_3', 'simul_4', 'simul_5')){
  
  res_all <- c()
  
  print(paste0("name: ", name))
  
  nb_rep <- 200
  seed=1
  perm=1
  
  for (seed in 1:nb_rep){
    
    print(paste0("seed: ", seed))
    
    file <- read.table(paste0("empress/res_d1_t1_l1_", name, "_", seed,".svg"), sep="\t", fill=TRUE)
    pvalue <- as.numeric(sapply(strsplit(split=" ", file$V1[grep("p-value", file$V1)]), "[[", 8))
    
    # get the full distribution 
    res <- extract_blue_distribution(paste0("empress/res_d1_t1_l1_", name, "_", seed,".svg"))
    mu_null <- res$mean
    sd_null <- res$sd
    
    
    # look at number of events 
    recon <- read.table(paste0("empress/reconciliation_d1_t1_l1_",name,"_",seed,".csv"), sep=",")
    colnames(recon) <- c("parasite", "host", "event", "node_frequency", "event_frequency")
    nb_events <- table(recon$event)[c("Cospeciation", "Transfer")]
    
    res <- c()
    
    # Cost for H0 biogeo 
    for (perm in 1:100){
      if (file.exists(paste0("empress/cost_d1_t1_l1_", name, "_", seed, "_biogeo_", perm,".csv"))){
        if (file.size(paste0("empress/cost_d1_t1_l1_", name, "_", seed, "_biogeo_", perm,".csv"))){
          res <- rbind(res, read.table(paste0("empress/cost_d1_t1_l1_", name, "_", seed, "_biogeo_", perm,".csv"))$V1)
        }
      }
    }
    
    res <- data.frame(res)
    colnames(res) <- c("name","seed","perm","cost","cost_random")
    res$cost <- as.numeric(res$cost)
    res$cost_random <- as.numeric(res$cost_random)
    cost <- res$cost[1]
    pvalue_biogeo <- length(which(res$cost_random<=res$cost))/nrow(res)
    
    mu_null_biogeo <- mean(res$cost_random, na.rm = TRUE)
    sd_null_biogeo <- sd(res$cost_random, na.rm = TRUE)
    
    res_all <- rbind(res_all, c(name, seed, cost, nb_events,
                                pvalue, pvalue_biogeo,
                                mu_null, sd_null,
                                mu_null_biogeo, sd_null_biogeo))
  }
  
  res_all <- data.frame(res_all)
  colnames(res_all) <- c("name", "seed", "cost", "number_cospeciations", "number_transfers",
                         "pvalue", "pvalue_biogeo", "mu_null", "sd_null", "mu_null_biogeo", "sd_null_biogeo")
  
  res_all$cost <- as.numeric(res_all$cost)
  res_all$number_cospeciations <- as.numeric(res_all$number_cospeciations)
  res_all$number_transfers <- as.numeric(res_all$number_transfers)
  res_all$pvalue <- as.numeric(res_all$pvalue)
  res_all$pvalue_biogeo <- as.numeric(res_all$pvalue_biogeo)
  res_all$mu_null <- as.numeric(res_all$mu_null)
  res_all$sd_null <- as.numeric(res_all$sd_null)
  res_all$mu_null_biogeo <- as.numeric(res_all$mu_null_biogeo)
  res_all$sd_null_biogeo <- as.numeric(res_all$sd_null_biogeo)
  
  
  write.table(res_all, paste0("results_",name,"_empress.csv"), sep=";", row.names=FALSE, quote=FALSE)
  
  
  # compute stats
  res_all$SES_null <- (res_all$cost - res_all$mu_null)/res_all$sd_null
  res_all$SES_null_biogeo <- (res_all$cost - res_all$mu_null_biogeo)/res_all$sd_null_biogeo
  
  # fraction of the original deviation explained by geographic structure.
  res_all$fraction_explained <- (res_all$mu_null_biogeo - res_all$mu_null) / (res_all$cost - res_all$mu_null)
  
  write.table(res_all, paste0("results_",name,"_empress.csv"), sep=";", row.names=FALSE, quote=FALSE)
  
}


