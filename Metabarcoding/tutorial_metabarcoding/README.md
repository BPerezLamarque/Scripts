
**Tutorial to analyze metabarcoding datasets using OTU clustering:**
====


This tutorial provides instructions on how to perform OTU clustering using VSEARCH to generate an OTU table from an amplicon metabarcoding dataset. It is part of the practical sessions of the Master course [Environmental genomics for microbial ecology](https://www.edu.bio.ens.psl.eu/spip.php?article271) (ENS, Université PSL).


The pipeline requires the installation of Python3, VSEARCH (v2), cutadapt, and FastQC. It also requires external scripts ([OTU_contingency_table.py](https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/OTU_contingency_table.py), [map2qiime.py](https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/map2qiime.py), and [make_stats.py](https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/make_stats.py)). The complete documentation of VSEARCH is available [here](https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch_manual.pdf). It is assumed that the dataset is **already demultiplexed** (see [here](https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/pipeline_OTU97.sh) for demultiplexing), with one R1 fastq file and R2 fastq file per sample. 


In this tutorial, we will analyze plant-associated root mycobiota from a disturbed community in La Réunion. The dataset includes root samples from both native and invasive plant species. Root-associated fungal communities have been characterized using metabarcoding of the 18S rRNA region (with the primer pair AMADf-AMDGr). The original and complete dataset is available in [Perez-Lamarque et al. (2022)](https://doi.org/10.1186/s40793-022-00434-0). 


**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com



# Contents:
**[Step 0: Prepare the database for taxonomic assignation](#step-0-prepare-the-database-for-taxonomic-assignation)**\
**[Step 1: Merge paired-ended reads and check the quality](#step-1-merge-paired-ended-reads-and-check-the-quality)**\
**[Step 2: OTU clustering at 97%](#step-2-otu-clustering-at-97)**

<br> 

## Set the working directory

All analyses will be performed in a given directory. Inside this directory, create a folder "raw_data/", which contains all the fastq files.  

Because this pipeline will generate large files, it is recommended to work in the /data/ folder, e.g. in /data/mag/bio21/login/EGME/praticals_mycorrhiza/
Warnings: This directory has to be empty before starting the analyses

```bash

path_analyses="WORKING_DIRECTORY" # path where all analyses will be done. 

cd $path_analyses

mkdir raw_data/

# copy all the fastq files into the folder raw data:
cp /XXXX/*fastq $path_analyses/raw_data/


```

<br> 

# Step 0: Prepare the database for taxonomic assignation


The goal of this step is to prepare the database that will be used for the taxonomic assignation of the OTUs. Here, the goal is the generate a database from the SILVA database,  which mostly contains sequences of the small subunit of ribosomal RNA (16S for prokaryotes and 18S for eukaryotes). We will focus here on eukaryotes using the two primers AMADf-AMDGr designed for root-associated fungi. 

After downloading the complete SILVA database, we use cutadapt to extract from the database the fungal reads that could have been amplified using the primer pair AMADf-AMDGr. Cutadapt gets the sequences matching both primers, trims forward & reverse primers to only keep the central amplicon region that could be present in our dataset, and store these reference amplicon sequences into a correct format. 


```bash

# Create a folder where the database will be stored 

mkdir database/
cd $path_analyses"/database/"


RELEASE=138.1 # version of SILVA to use 
URL="https://ftp.arb-silva.de/release_${RELEASE}/Exports"
INPUT="SILVA_${RELEASE}_SSURef_NR99_tax_silva.fasta.gz"

# Download the SILVA database
wget -c ${URL}/${INPUT}{,.md5} && md5 -r ${INPUT}.md5

# Define the name of the database (with the primer names)
OUTPUT="${INPUT/.fasta.gz/_AMADf_AMDGr.fasta}"
LOG="${INPUT/.fasta.gz/_AMADf_AMDGr.log}"


# Indicate the sequences of the primers used for metabarcoding
PRIMER_F="GGGAGGTAGTGACAATAAATAAC" # Forward primer (AMADf)
PRIMER_R="CCCAACTATCCCTATTAATCAT" #Reverse primer (AMDGr)

# Reverse complement the reverse primer
ANTI_PRIMER_R=RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )

# Paramerers of cutadapt
MIN_LENGTH=32 # minimum sequence length
MIN_F=$(( ${#PRIMER_F} / 2 )) # should match at least 50% of the forward primer
MIN_R=$(( ${#PRIMER_R} / 2 )) # should match at least 50% of the reverse primer
CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH}"

# Run cutadapt: 
zcat < "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
    ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
    ${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" | \
    sed '/^>/ s/;/|/g ; /^>/ s/ /_/g ; /^>/ s/_/ /1' > "${OUTPUT}"
    
```



Once the final database is ready, we can delete the original database:
```bash
rm ${INPUT}
```

<br> 


# Step 1: Merge paired-ended reads and check the quality


All analyses will now be performed in the folder "process".

```bash
mkdir $path_analyses/process/

cd $path_analyses/process/

# number of cores for the analyses
nb_cores=4

```


First, we will need to get the list of all the samples. This can be obtained using the following line: 

```bash

# List all samples ("fastq" format), remove the extension, and store the names

ls $path_analyses/raw_data/*.fastq | sed 's/\_1\.fastq$//' | sed 's/\_2\.fastq$//' | sed 's/\.fastq$//' | \
    sort -u > list_sample.txt


```

<br> 

## Step 1-A: Merge the paired-ended reads


We use a for loop to iterate over each line in list_sample.txt (i.e. one sample) and for each sample, we merge the R1 & R2 paired-ended reads:


```bash

mkdir $path_analyses/process/merged_reads/

for sample in $(cat list_sample.txt); do

    echo $sample
    
    INPUT_R1="../raw_data/"$sample"_1.fastq"
    INPUT_R2="../raw_data/$sample"_2.fastq"
    
    OUTPUT="merged_reads/MERGED_R1R2_"$sample".fastq"

    vsearch \
        --threads $nb_cores \
        --fastq_mergepairs ${INPUT_R1} \
        --reverse ${INPUT_R2} \
        --fastq_ascii 33 \
        --fastqout ${OUTPUT} \
        --quiet 2>> ${OUTPUT/.fastq/.log}
    
done

```

<br> 

## Step 1-B: Checking the quality with FastQC

Next, we use FastQC to check a for loop to iterate over each line in list_sample.txt (i.e. one sample), and for each sample, we merge the paired-ended reads:


```bash

mkdir $path_analyses/process/FastGC/

for sample in $(cat list_sample.txt); do

    echo $sample
    
    OUTPUT="merged_reads/MERGED_R1R2_"$sample".fastq"

    fastqc -t $nb_cores ${OUTPUT} -o FastQC/
    
done

```

Check the quality file generated for each sample and make sure that all the samples meet a certain quality. 

<br> 

## Step 1-C: Quality filtering 


In this step, we discard (i) the reads containing "N"s and (ii) the reads having >1 expected error (for a given read, the expected error is the sum of error probabilities for all its positions). 

The reads are then converted in fasta format and we use dereplication (i.e. identical reads are grouped together and we keep the information of their abundance in the header). During dereplication, the reads are renamed with a code name based on their sequence (sha1 encoding system); thus the same names are given to identical amplicons present in different samples.

```bash

for sample in $(cat list_sample.txt); do

    echo $sample
    
    INPUT="merged_reads/MERGED_R1R2_"$sample".fastq"
    OUTPUT="merged_reads/"$sample".fasta"
    
    vsearch \
    --fastq_filter $INPUT \
    --fastq_maxns 0 \
    --fastq_maxee 0.5 \
    --fastaout "temp_"$sample".fas"
    
    # Dereplicate at the sample level
    vsearch --quiet \
        --derep_fulllength "temp_"$sample".fas" \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel_sha1 \
        --output "${OUTPUT}"
        
    rm "temp_"$sample".fas" # remove the temporary file 
    
done

```

Note that the parameter "fastq_maxee" can be adjusted according to the average quality of the run (given by FastQC). 

 
Once all these steps have finished and the results are fine, we can remove the files that are no longer useful:

```bash

rm merged_reads/*fastq

```

<br> 

# Step 2: OTU clustering at 97%


All clustering analyses will also be performed in the folder "process". More details on this script are available [here](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline).

```bash
cd $path_analyses/process/

# number of cores for the analyses

nb_cores=4

# indicate the name of the database for taxonomic assignation (generated in Step 0):

database_taxonomy="SILVA_138.1_SSURef_NR99_tax_silva_AMADf_AMDGr.fasta"

```

<br> 

## Step 2-A: Dereplication of all the reads from all samples


```bash

cat $path_analyses/process/merged_reads/*.fas > reads_amplicon.fa

vsearch \
    --derep_fulllength reads_amplicon.fa \
    --sizein \
    --sizeout \
    --relabel_sha1 \
    --fasta_width 0 \
    --output reads_amplicon_derep.fa

rm reads_amplicon.fa

```

<br> 

## Step 2-B: Sort by size


We sort all the reads by size (from the most abundant to the least abundant). 

Here, we discard all the reads present in a single copy. Warning: We exclude singletons for increasing the speed of the clustering, but this step is not always recommended (use "minsize 1" instead")

```bash

vsearch -sortbysize reads_amplicon_derep.fa -output reads_amplicon_sorted.fa -minsize 2

```

<br> 

##  Step 2-C: OTU clustering at 97%

Now we perform the OTU clustering at 97%. For each OTU, we keep the centroid read as the representative sequence. 

```bash

vsearch -cluster_size  reads_amplicon_sorted.fa --id 0.97 --centroids reads_OTU97.fa --uc clusters_OTU97.uc --sizein --sizeout

```

Note, that 0.97 can be changed to 0.95, 0.99... to generate 99% or 95% OTUs instead. 

<br> 

The step of OTU clustering generates several outputs that will be necessary for building the OTU table (i.e. which OTUs are present in each sample). These files need proper formatting:

1) The fasta file "reads_OTU97.fa" contains the representative sequences of each OTU. It can be reformatted with one line per OTU sequence using the following command:

```bash

vsearch --fasta_width 0 --sortbysize reads_OTU97.fa --output reads_OTU97_final.fa

```

2) The step of OTU clustering also keeps track of the reads mapping to each OTU (file "clusters_OTU97.uc"). This file can be reformatted using the python script map2qiime.py:

```bash

python3 $path_scripts/map2qiime.py clusters_OTU97.uc > reads_mapped_OTU97.txt

```

3) From the fasta file "reads_OTU97_final.fa" we will also extract the abundance of each OTU (needed in the OTU table), with the following command:

```bash

python3 $path_scripts/make_stats.py reads_OTU97_final.fa > stats_file_OTU97.txt

```

<br> 


##  Step 2-D: Chimera checking

Next, we perform chimera filtering: we identify *de novo* the presence of chimeras among the OTUs and remove them. 

```bash

vsearch --uchime_denovo reads_OTU97_final.fa \
    --uchimeout reads_OTU97.uchime \
    --nonchimeras reads_OTU97_nonchimeras.fa

```

Non-chimeric OTUs are stored in the fasta file "reads_OTU97_nonchimeras.fa", while "reads_OTU97.uchime" keeps a summary of the chimera checking for each OTU. 

<br> 


##  Step 2-E: Taxonomic assignation 

We now assign a taxonomy to each OTU by blasting the OTUs against the SILVA database. 

```bash

vsearch --usearch_global reads_OTU97_nonchimeras.fa \
    --threads $nb_cores \
    --dbmask none \
    --qmask none \
    --rowlen 0 \
    --notrunclabels \
    --userfields query+id1+target \
    --maxaccepts 0 \
    --maxrejects 32 \
    --top_hits_only \
    --output_no_hits \
    --db "../database/"$database_taxonomy \
    --id 0.5 \
    --iddef 4 \
    --userout taxonomy_OTU97.txt

```

<br> 


##  Step 2-F:  Generate the OTU table

Finally, we can build the OTU table. 

```bash

STATS="stats_file_OTU97.txt"
OTUS="reads_mapped_OTU97.txt"
REPRESENTATIVES="reads_OTU97_final.fa"
UCHIME="reads_OTU97.uchime"
ASSIGNMENTS="taxonomy_OTU97.txt"
OTU_TABLE="OTU_table_OTU97.txt"


SCRIPT=$path_scripts"/OTU_contingency_table.py" 

python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${OTUS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    $path_analyses/process/merged_reads/*.fas > "${OTU_TABLE}"
    
```

<br> 

The generated OTU will likely be too large to be easily opened in R. We can already filter the OTU at this step to discard OTU, which won't be useful for statistical analyses. For instance, we discard the chimeric OTUs, the OTUs less than 200 bp, and the OTUs with an abundance lower than 10:

```bash

FILTERED="${OTU_TABLE/.txt/_filtered.txt}"
head -n 1 "${OTU_TABLE}" > "${FILTERED}"
cat "${OTU_TABLE}" | awk '$5 == "N" && $4 >= 200 && $2 >= 10' >> "${FILTERED}"

```

The final OTU table is then available in "OTU_table_OTU97_filtered.txt" and can be opened in R. 

<br> 


