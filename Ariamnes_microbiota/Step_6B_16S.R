########################################################################################################
########################################################################################################
#### Step 6 B: Make the global OTU Table ####

rm(list=ls(all=TRUE)) 

############################################################################
####  Specify your path : #### 

param <- commandArgs(trailingOnly=TRUE)

path <- noquote(as.character(param[1]))
indexes <- as.numeric(param[2:length(param)])

print(indexes)

metadata <- read.table(paste(path,"metadata/metadata.csv",sep=''), header=T, sep=";")

setwd(paste(path,"OTU_table/",sep=''))


OTU_table <- read.table(paste(path,"OTU_table/OTU_list.txt",sep=''),sep="\t", header=F)
taxo <- OTU_table
colnames(OTU_table) <- c("OTU ")
n <- as.numeric(nrow(OTU_table))

for (k in indexes){
  print(k)
  table <- read.table(file=paste("OTU_S",k,"_16S_abundance.csv",sep=""), header=T, sep=";")
  m <- as.numeric(nrow(table))
  a <- noquote(as.character(metadata$Sample_ID[which(metadata$Name_microbiota==k)]))
  OTU_table[,a] <- 0.0
  OTU_table[,1] <- as.character(OTU_table[,1])
  table$OTU <- as.character(table$OTU)
  for(j in 1:m){ 
    i <- which(OTU_table$OTU==table$OTU[j])
    OTU_table[i,a] <- as.numeric(table[j,2])
    taxo[i,2] <- as.character(table[j,4])
  }
}

OTU_table_2 <- OTU_table

for(i in 1:ncol(OTU_table_2)){OTU_table_2[,i] <- noquote(as.character(paste(OTU_table[,i],"",sep=".0")))}

OTU_table_2[,1] <- OTU_table[,1]
colnames(OTU_table_2)[1]  <- "#OTU ID "
write.table(OTU_table_2, file=paste(path,"OTU_table/OTU_table_all.txt",sep=''),row.names=F, sep="\t")

OTU_table_3 <-data.frame(OTU_table_2, taxo[,c(2)])
colnames(OTU_table_3)[1]  <- "#OTU ID "
colnames(OTU_table_3)[ncol(OTU_table_3)]  <- "taxonomy"
write.table(OTU_table_3, file=paste(path,"OTU_table/OTU_table_taxonomy.txt",sep=''),row.names=F, sep="\t")


