#!/bin/bash

#####################################################################################################
######################  PERFORM OTU CLUSTERING FROM  METABARCODING DATASETS USING 97% OTUs ##########
#####################################################################################################


## You first need to install PYTHON3, VSEARCH, and CUTADAPT
## The python library numpy must also be installed (https://numpy.org/install/)
## see https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline for details


# You also need to download the python scripts OTU_contingency_table.py, map2qiime.py, and make_stats.py available in https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/

# Store them in a given folder
path_scripts="PATH_TOWARD_THE_SCRIPTS"



#####################################################################################################
######################  STEP0: PREPARE THE DATABASE FOR THE TAXONOMIC ASSIGNATION     ###############
#####################################################################################################



##### If you are working with the ITS2 region (e.g. with the primers ITS86F and ITS4)   ##############

# Download the UNITE database (https://unite.ut.ee)
# use 97% OTU of  Fungi + other Eukaryotes


path_analyses="YOUR_WORKING_DIRECTORY"
cd $path_analyses

mkdir database/
cd $path_analyses"/database/"


INPUT="UNITE_all_Euk_OTU97_1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.fasta" # UNITE database that you have previoulsy mannualt downloaded


# Define variables and output files for your primers
OUTPUT="UNITE_all_Euk_OTU97_ITS2_86F.fasta"
LOG="UNITE_all_Euk_OTU97_ITS2_86F.log"


# Precise the sequences of the primers you used for the metabarcoding
PRIMER_F="GTGAATCATCGAATCTTTGAA" # ITS86f
# Do not check for ITS4 because it is outside of most of the ITS sequences in the UNITE  database
MIN_LENGTH=200
MIN_F=$(( ${#PRIMER_F} / 2 ))

CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH}"

# Trim forward & reverse primers, format
cat "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
sed '/^>/ s/;/|/g ; /^>/ s/ /_/g ; /^>/ s/_/ /1' > "${OUTPUT}"




##### If you are working with the 16S or 18S region (e.g. with the primers 799F and 1193R of the 16S region)   ##############


path_analyses="YOUR_WORKING_DIRECTORY"
cd $path_analyses

mkdir database/
cd $path_analyses"/database/"


RELEASE=138.1
URL="https://ftp.arb-silva.de/release_${RELEASE}/Exports"
INPUT="SILVA_${RELEASE}_SSURef_NR99_tax_silva.fasta.gz"

# Download and check
wget -c ${URL}/${INPUT}{,.md5} && md5 -r ${INPUT}.md5

# Define variables and output files for your primers
OUTPUT="${INPUT/.fasta.gz/_799F_1193R.fasta}"
LOG="${INPUT/.fasta.gz/_799F_1193R.log}"


# Precise the sequences of the primers you used for the metabarcoding
PRIMER_F="AACMGGATTAGATACCCKG"
PRIMER_R="ACGTCATCCCCACCTTCC"
ANTI_PRIMER_R="GGAAGGTGGGGATGACGT"
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} / 2 ))
MIN_R=$(( ${#PRIMER_R} / 2 ))
CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH}"

# Trim forward & reverse primers, format
zcat < "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" | \
sed '/^>/ s/;/|/g ; /^>/ s/ /_/g ; /^>/ s/_/ /1' > "${OUTPUT}"

rm ${INPUT}




#####################################################################################################
######################     STEP1: MERGE R1 & R2 paired-ended sequences	       ######################
#####################################################################################################


# Step 1-A: Prepare the data


# Give the name of your marker gene
Marker="ITS2"
# or Marker="16S"


path_analyses="YOUR_WORKING_DIRECTORY"
cd $path_analyses

mkdir raw_data/
# Download the data and store them in the folder "raw_data"

mkdir process/
cd $path_analyses"/process/"

# Write the name of the paired-ended sequences
INPUT_R1="../raw_data/210920_SN234_A_L001_ATXJ-7_R1.fastq"  # bz2, gz and uncompressed fastq files are allowed
INPUT_R2="../raw_data/210920_SN234_A_L001_ATXJ-7_R2.fastq"


OUTPUT="Quality.encoding_"$Marker".log"

# Check the quality encoding (33 or 64?)
vsearch \
    --fastq_chars ${INPUT_R1} 2> ${OUTPUT/.fastq/.log}


# Step 1-B: Merge the paired-end reads

THREADS=2
ENCODING=33
OUTPUT="MERGED_R1R2_"$Marker".fastq"
OUTPUT_NOTMERGED="UNMERGED_R1R2_"$Marker


vsearch \
    --threads ${THREADS} \
    --fastq_mergepairs ${INPUT_R1} \
    --reverse ${INPUT_R2} \
    --fastq_ascii ${ENCODING} \
    --fastqout ${OUTPUT} \
    --fastq_allowmergestagger \
    --quiet 2>> ${OUTPUT/.fastq/.log}
    
    
# Step 1-C: Checking the quality with FASTQC

# FastQC is available on https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip)
mkdir FastQC/
/PATH_TO_FastQC/fastqc -t 4 ${OUTPUT} -o FastQC/




#####################################################################################################
######################  			STEP2: DEMULTIPLEX SEQUENCES	           ######################
#####################################################################################################


# Your 'mapping' file should contain at least 5 columns, the first one being the sample names (starting with '#'), and the following being, not necessary in this order but with the...
# ... exact names : "primerFw", "primerRev", "barcodeFw", "barcodeRev".
# see an example in https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/mapping_example_ITS.txt


# Step 2-A: Check the  mapping file


# Give the name of your marker gene
Marker="ITS2"
# or Marker="16S"

MAPPING="mapping_"$Marker".txt"


barcodeFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'barcodeFw' | cut -d: -f1) - 1 ))" # getting the index of column with forward barcode etc. "-1" is applied as array below counts from 0
primerFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'primerFw' | cut -d: -f1) -1 ))"
barcodeRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'barcodeRev' | cut -d: -f1) -1 ))"
primerRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'primerRev' | cut -d: -f1) -1 ))"


# Check that all the samples are present
while read -r line; do
    if [[ ! "$line" =~ ^#.* && ! "$line" =~ ^Sample.* ]]; then
        IFS=$'\t' read -r -a array <<< "$line"
        
        # Get sequences
        FwBarcode="${array[${barcodeFwColumnIdx}]}"
        FwPrimer="${array[${primerFwColumnIdx}]}"
        RevBarcode="${array[${barcodeRevColumnsIdx}]}"
        RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
        RevPrimer="${array[${primerRevColumnsIdx}]}"
        RevPrimerRC=$( echo "${RevPrimer}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
        SAMPLE_NAME="${array[0]}"
        
        # Some information
        echo "${SAMPLE_NAME} is being tested.."
        echo "Barcode Fw: ${FwBarcode}"
        echo "Primer Fw: ${FwPrimer}"
        echo "Primer Rev (RC): ${RevPrimerRC}"
        echo "Barcode Rev (RC): ${RevBarcodeRC}"
    fi
    
done < "${MAPPING}"


# Step 2-B: Perform the quality filtering


# Discard sequences containing Ns, remove if >2 error on average, convert to fasta
vsearch \
    --fastq_filter "MERGED_R1R2_"$Marker".fastq" \
    --fastq_maxns 0 \
    --fastq_maxee 2 \
    --fastaout "MERGED_R1R2_"$Marker".fasta"
    

 
# Step 2-C: Demultiplex the sequences
 
 
# Define binaries, temporary files and output files
INPUT="MERGED_R1R2_"$Marker".fasta"
INPUT_REVCOMP="${INPUT/.fasta/_RC.fasta}"


# Reverse complement fastq file
vsearch --quiet \
    --fastx_revcomp "${INPUT}" \
    --fastaout "${INPUT_REVCOMP}"


mkdir "Demultiplexed_data_"$Marker/

MIN_LENGTH=200


while read -r line; do
    if [[ ! "$line" =~ ^#.* && ! "$line" =~ ^Sample.* ]]; then
        IFS=$'\t' read -r -a array <<< "$line"

        # Get sequences
        FwBarcode="${array[${barcodeFwColumnIdx}]}"
        FwPrimer="${array[${primerFwColumnIdx}]}"
        RevBarcode="${array[${barcodeRevColumnsIdx}]}"
        RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
        RevPrimer="${array[${primerRevColumnsIdx}]}"
        RevPrimerRC=$( echo "${RevPrimer}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
        SAMPLE_NAME="${array[0]}"
        
        # Output file names
        LOG="Demultiplexed_data_"$Marker"/${SAMPLE_NAME}.log"
        FINAL_FASTA="Demultiplexed_data_"$Marker"/${SAMPLE_NAME}.fas"
        
        # Some information
        echo "${SAMPLE_NAME} is being processed.."
        echo "Barcode Fw: ${FwBarcode}"
        echo "Primer Fw: ${FwPrimer}"
        echo "Primer Rev (RC): ${RevPrimerRC}"
        echo "Barcode Rev (RC): ${RevBarcodeRC}"
        
        function trim_without_ambiguity {

            SEQTOT="${FwBarcode}${FwPrimer}${RevPrimerRC}${RevBarcodeRC}"
            MIN_MATCHED=${#SEQTOT}
            ERROR_RATE=0

            cat "${INPUT}" "${INPUT_REVCOMP}" | cutadapt -g "${FwBarcode}${FwPrimer}...${RevPrimerRC}${RevBarcodeRC}" --discard-untrimmed --minimum-length "${MIN_LENGTH}" -O ${MIN_MATCHED} -e "${ERROR_RATE}" - 2> "${LOG}" > "temp_"$Marker".fasta"
        }
        

        trim_without_ambiguity

        # Dereplicate at the study level
        vsearch --quiet \
            --derep_fulllength "temp_"$Marker".fasta" \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --relabel_sha1 \
            --output "${FINAL_FASTA}" 2>> "${LOG}"

    fi
done < "${MAPPING}"


# Note: the option sha1 (encoding system) is giving the same names to the identical amplicons across samples


# Remove the files that are no longer useful
rm -f "${INPUT}" "${INPUT_REVCOMP}" "temp_"$Marker".fastq" "temp_"$Marker".fasta"




###############################################################################################
###############################   STEP3:  CLUSTERING AT 97%      ##############################
###############################################################################################


## to see all the details: https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline

path_analyses="YOUR_WORKING_DIRECTORY"

cd $path_analyses"/process/"

Marker="ITS2"
database_taxonomy="UNITE_all_Euk_OTU97_ITS2_86F.fasta"

# or : Marker="16S
# database_taxonomy="SILVA_138.1_SSURef_NR99_tax_silva_799F_1193R.fasta"



# Step 3-A: Dereplicate of fasta file of every  individual sample  using VSEARCH (already done in Step 2)

# Important: for every sample, you must have a separate file containing all the metabarcoding reads from this sample that are **dereplicated with the option relabel_sha1** (as done in Step 2).



# Step 3-B: Dereplication of fasta file of all the reads using VSEARCH

PREFIX=""
cat $path_analyses/process/"Demultiplexed_data_"$Marker/"$PREFIX"*.fas > "reads_"$Marker".fa"


vsearch \
    --derep_fulllength "reads_"$Marker".fa" \
    --sizein \
    --sizeout \
    --relabel_sha1 \
    --fasta_width 0 \
    --output "reads_"$Marker"_derep.fa"


rm "reads_"$Marker".fa"


# Total number of reads in the dataset
grep ">" "reads_"$Marker"_derep.fa" | sed 's/>.*.;size=//g' | paste -sd+ - | bc



##  Step 3-C: Clustering at 97%

# Sort by size
vsearch -sortbysize "reads_"$Marker"_derep.fa" -output "reads_"$Marker"_sorted.fa" -minsize 1


## OTU clustering (at 97%)
vsearch -cluster_size  "reads_"$Marker"_sorted.fa" --id 0.97 --centroids "reads_"$Marker"_OTU97.fa" --uc "clusters_"$Marker"_OTU97.uc" --sizein --sizeout
# 0.97 can be changed to 0.95, 0.99...


# Make file that mapped reads to OTUs
python3 $path_scripts/map2qiime.py "clusters_"$Marker"_OTU97.uc" > "reads_"$Marker"_mapped_OTU97.txt" # python script present in https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/map2qiime.py


# One line per OTU sequence
vsearch --fasta_width 0 \
--sortbysize "reads_"$Marker"_OTU97.fa" \
--output "reads_"$Marker"_OTU97_final.fa"


# Make stats file
python3 $path_scripts/make_stats.py "reads_"$Marker"_OTU97_final.fa" > "stats_"$Marker"_OTU97.txt" # python script present in https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/make_stats.py



##  Step 3-D: Chimera checking

vsearch --uchime_denovo "reads_"$Marker"_OTU97_final.fa" --uchimeout "reads_"$Marker"_OTU97.uchime" --nonchimeras "reads_"$Marker"_OTU97_nonchimeras.fa"




##  Step 3-E: Assign taxonomy to OTUs using known reference sequences with blast


# Number of threads can be set to 1
vsearch --usearch_global "reads_"$Marker"_OTU97_nonchimeras.fa" \
--threads 2 \
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
--userout "taxonomy_"$Marker"_OTU97.txt"




##  Step 3-F:  Make the OTU table

STATS="stats_"$Marker"_OTU97.txt"
OTUS="reads_"$Marker"_mapped_OTU97.txt"
REPRESENTATIVES="reads_"$Marker"_OTU97_final.fa"
UCHIME="reads_"$Marker"_OTU97.uchime"
ASSIGNMENTS="taxonomy_"$Marker"_OTU97.txt"
OTU_TABLE="OTU_table_"$Marker"_OTU97.txt"


SCRIPT=$path_scripts"/OTU_contingency_table.py"  # this script is located https://github.com/BPerezLamarque/Scripts/blob/master/Metabarcoding/OTU_contingency_table.py

python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${OTUS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    $path_analyses/process/"Demultiplexed_data_"$Marker/"$PREFIX"*.fas > "${OTU_TABLE}"


# Filter per OTU size or spread, quality and chimeric status:
FILTERED="${OTU_TABLE/.txt/_filtered.txt}"
head -n 1 "${OTU_TABLE}" > "${FILTERED}"
cat "${OTU_TABLE}" | awk '$5 == "N" && $4 >= 200 && $2 >= 10 && $6 >= 1' >> "${FILTERED}"
# remove chimera, keep OTU of more than 200 bp, with an abundance of at least 10 reads, and spread of >=1


# Total number of reads in the OTU table
grep -v "OTU\tabundance" "${OTU_TABLE}" | cut -f2 | paste -sd+ - | bc

# Total number of reads in the filtered OTU table
grep -v "OTU\tabundance" "${FILTERED}" | cut -f2 | paste -sd+ - | bc






