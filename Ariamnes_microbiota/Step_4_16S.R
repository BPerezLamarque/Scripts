########################################################################################################
########################################################################################################
#### Step 4: Make the OTU Table for the negative controls : #### 

rm(list=ls(all=TRUE)) 

############################################################################
####  Load Data : #### 

param <- commandArgs(trailingOnly=TRUE)

path <- noquote(as.character(param[1]))
indexes <- as.numeric(param[2:length(param)])

metadata <- read.table(paste(path,"metadata/metadata.csv",sep=''),sep=";",header=T)

metadata$Number_sequences_16S <- NA
metadata$Number_OTU_16S <- NA

setwd(paste(path,"results_wNC/",sep=''))

for(k in indexes){  #negative control indexes
  print(k)
  table <- read.table(file=paste("OTU","_16s_taxonomy.txt",sep=""), header=F, sep="\t")
  count <- read.table(file=paste("OTU_S",k,"_16S_table.txt",sep=""), header=F, sep="\t")
  
  m <- as.numeric(ncol(count))
  
  metadata$Number_sequences_16S[which(metadata$Name_microbiota==k)] <- m-1  # Total number sequences
  print(m-1)
  n<- as.numeric(nrow(table)) #Number OTU in all the samples
  
  count[,c(2:m)] <- sapply(count[,c(2:m)],as.numeric) #convert columns content in numeric values (0 or 1)
  
  abund <- count[,1:2]
  w <- as.numeric(nrow(abund)) #Number OTU in the sample k
  
  metadata$Number_OTU_16S[which(metadata$Name_microbiota==k)] <- w-2
  print(w-2)
  abund[,2] <- apply(count[,2:m],1,sum)
  
  rm(count)
  
  for (j in 1:w){abund[j,3] <- as.character(table[which(as.character(table[,1])==as.character(abund[j,1])),2])}
  
  abund[,4] <- as.numeric(abund[,2])/m #relative abundance
  
  abund <- abund[order(abund[,2],decreasing=T),] #sorted by abundance
  row.names(abund)<-NULL
  
  table <-data.frame(abund[,c(1)], abund[,c(2,4)],abund[,c(3)])
  colnames(table) <- c("OTU","Absolute_Abundance","Relative_Abundance","Taxonomy")
  
  write.table(table, file=paste(path,"results_wNC/OTU_S",k,"_16S_abundance.csv",sep=""),row.names=F, sep=";")
  rm(table)
}

write.table(metadata, file=paste(path,"metadata/metadata_NC.csv",sep=''),row.names=F, sep=";")
write.table(metadata, file=paste(path,"metadata/metadata_NC_only.csv",sep=''),row.names=F, sep=";")
