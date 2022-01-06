# Scripts


## Pipelines for dealing with metabarcoding data

The folder  `Metabarcoding`  contains bash scripts to analyze metabarcoding data from Illumina paired-end sequencing. They can be used to cluster [OTUs at 97% (or any other thresholds)](https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/pipeline_OTU97.sh), [ASV](https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/pipeline_ASV.sh), or  [Swarm OTUs](https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/pipeline_OTU_Swarm.sh) ([Mahé et al., 2014](https://peerj.com/articles/593/); [Mahé et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4690345/pdf/peerj-03-1420.pdf)). They can be used for any marker gene (ITS, 16S, 18S...)

These scripts are derived from [Frédéric Mahé's metabarcoding pipeline](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline) and are mainly using [VSEARCH](https://github.com/torognes/vsearch). They are composed of 4 steps: Step 0 prepares a database for the taxonomic assignation, Step 1 merges the paired-end reads, Step 2 demultiplexes the sequences (i.e. assigns each read to its sample), and Step 3 performs the clustering. 

**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com




## Scripts for *Ariamnes* microbiota analyses 

The folder  `Ariamnes_microbiota`  contains the functions to reproduce the microbiota analyses from the manuscript: *A holobiont view of island biogeography: Unravelling patterns driving the nascent diversification of a Hawaiian spider and its microbial associates*, by Ellie E. Armstrong, Benoît Perez-Lamarque, Ke Bi, Cerise Chen, Leontine E. Becking, Jun Ying Lim, Tyler Linderoth, Henrik Krehenwinkel, Rosemary G. Gillespie, 2021, in Molecular Ecology, https://doi.org/10.1111/mec.16301.


**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com
