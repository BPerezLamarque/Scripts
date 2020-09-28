#!/bin/sh
# NGS Analysis: Ariamnes microbiome 

# cd $path/script/
# chmod +x ./analyses_16s_100.sh
# ./analyses_16s_100.sh

######################################
echo "PATH : folder with 2 folders: raw_data/ with the forward and reverse fastq sequences, and database/ with the reference databases "

path=/Users/bperez/ownCloud/Recherche/Stages/Stage_M1/results_2019_97/
path_script=/Users/bperez/ownCloud/Recherche/Stages/Stage_M1/results_2019_97/script/

cd $path/

path_usearch=/Users/bperez/ownCloud/appli/uparsev11/usearch11.0.667_i86osx32
path_pipelines=/Users/bperez/ownCloud/appli/pipelines/scripts/




######################################
echo "Step 0: Prepare the metadata folder (/path/metadata/):"
# Make the metadata_populations.txt file with the columns "#SampleID	Island	Population"


######################################
echo "Step 0: Unzip your files:"

# Index of the samples, with and without negative controls

declare -a indexes=(41 43 46 50 54 56 76 110 116 68 73 86 105 106 112 114 51 55 84 99 117 120 45 62 66 71 77 78 88 98 100 102 104 107 118 48 49 53 58 64 67 72 81 109 113 119 42 44 69 79 87 97 108 111 115 47 52 59 60 61 65 70 74 80 82 83 85 101 57 63 103)
indexes_NC=`(seq 121 123)`
declare -a indexes_wNC=(41 43 46 50 54 56 76 110 116 68 73 86 105 106 112 114 51 55 84 99 117 120 45 62 66 71 77 78 88 98 100 102 104 107 118 48 49 53 58 64 67 72 81 109 113 119 42 44 69 79 87 97 108 111 115 47 52 59 60 61 65 70 74 80 82 83 85 101 57 63 103 121 122 123)



######################################
echo "Step 1: Apply the assembly, filtering, demultiplexing, dereplication and sorting for each sample:"

mkdir $path/results_wNC/
mkdir $path/OTU_table/
mkdir $path/diversity/

# same as 97%

######################################
echo "Step 2: Apply the assembly, filtering, demultiplexing, dereplication and sorting for all the samples together: taxonomic assignation"

echo "Dereplication"
source activate qiime1
$path_usearch -derep_fulllength $path/trait/reads_16s.fa -fastaout $path/trait/reads_16s_derep.fa -sizeout

echo "Abundance sort and discard singletons"
$path_usearch -sortbysize $path/trait/reads_16s_derep.fa -fastaout $path/trait/reads_16s_sorted.fa -minsize 2

echo "OTU Clustering de novo:"
source activate qiime1
$path_usearch -cluster_otus $path/trait/reads_16s_sorted.fa -otus $path/trait/reads_16s_OTU_clean.fa -relabel OTU -minsize 4

echo "Fasta Formatter:"
/Users/bperez/ownCloud/appli/fastx/fasta_formatter -i $path/trait/reads_16s_OTU_clean.fa -o $path/trait/reads_16s_OTU_format.fa

echo "Renamer:"
cd $path_pipelines/
perl bmp-otuName.pl -i $path/trait/reads_16s_OTU_format.fa -o $path/trait/16s_OTU.fa

cd $path


#echo "Taxonomy:"
source activate qiime1
assign_taxonomy.py -i $path/trait/16s_OTU.fa -o $path/trait/ -r /Users/bperez/ownCloud/Recherche/Stages/Stage_M1/results_2017_99/database/97_otus.fasta -t /Users/bperez/ownCloud/Recherche/Stages/Stage_M1/results_2017_99/database/97_otu_taxonomy.txt

cp $path/trait/16s_OTU.fa $path/results_wNC/OTU_16s.fa
cp $path/trait/16s_OTU_tax_assignments.txt $path/results_wNC/OTU_16s_taxonomy.txt


######################################
echo "Step 3: Prepare the OTU Table"

for k in "${indexes_wNC[@]}"
do
    echo "Map reads back to OTU database:"
    
    ####   97% threshold
    $path_usearch -usearch_global $path/results_wNC/reads_S"$k"_16s.fa -db $path/results_wNC/OTU_16s.fa -strand plus -id 0.97 -uc $path/results_wNC/S"$k"_16s_map.uc

    echo "Convert OTU Table:"
    echo "16S"
    python $path_pipelines/bmp-map2qiime.py $path/results_wNC/S"$k"_16s_map.uc > $path/results_wNC/OTU_S"$k"_16s_table_temp.txt
    rm $path/results_wNC/S"$k"_16s_map.uc
    touch $path/results_wNC/OTU_S"$k"_16s_table.biom
    rm $path/results_wNC/OTU_S"$k"_16s_table.biom
    make_otu_table.py -i $path/results_wNC/OTU_S"$k"_16s_table_temp.txt -t $path/trait/16s_OTU_tax_assignments.txt -o $path/results_wNC/OTU_S"$k"_16s_table.biom
    rm $path/results_wNC/OTU_S"$k"_16s_table_temp.txt
    touch $path/results_wNC/OTU_S"$k"_16s_table_summary.txt
    rm $path/results_wNC/OTU_S"$k"_16s_table_summary.txt
    biom summarize-table -i  $path/results_wNC/OTU_S"$k"_16s_table.biom -o $path/results_wNC/OTU_S"$k"_16s_table_summary.txt
    touch $path/results_wNC/OTU_S"$k"_16S_table.txt
    rm $path/results_wNC/OTU_S"$k"_16S_table.txt
    biom convert -i  $path/results_wNC/OTU_S"$k"_16s_table.biom  -o $path/results_wNC/OTU_S"$k"_16S_table.txt --to-tsv
done

# More than 90% matched  (up to 99%)  - some OTUS are lost in this step?


######################################
echo "Step 4: Make the OTU Table for the negative controls"
### Use the R script "Step_4_16S.R"

source deactivate
R CMD BATCH "--args $path $indexes_NC" $path_script/Step_4_16S.R


######################################
echo "Step 5: Remove the contaminants from the samples files"

cd $path/results_wNC/
# Select the main contaminant OTUs from your negative control
grep "Brachybacterium\|_acnes\|Veillonella\|Brevibacterium\|Corynebacteriaceae\|Burkholderia\|Unassigned|\Caldicellulosiruptor" $path/results_wNC/OTU_16s_taxonomy.txt | grep -o "OTU[0-9]*" | sort | uniq > contamination.txt

cd $path/results_wNC/
python3 $path_script/remove_contaminant.py


######################################
echo "Step 6: Make the OTU Table for the samples (one without taxonomy, the other with)"
### Use the “Step_6A_16S” R script

indexes_R=${indexes[@]}

R CMD BATCH "--args $path $indexes_R" $path_script/Step_6A_16S.R


echo "Make the list of OTU"
cd $path/OTU_table/

for i in ${indexes[@]}
do
    echo $i
	grep -o "OTU[0-9][0-9]*" OTU_S"$i"_16S_abundance.csv >> OTU2.txt
done
cat OTU2.txt | sort |uniq > OTU_list.txt
rm OTU2.txt

# 525 OTUs (non contaminated)



echo "Make the OTU Table in txt file"
### Use the Step 6 B) on R script "Creation_OTU_table.R"

cd $path/OTU_table/

indexes_R=${indexes[@]}

R CMD BATCH "--args $path $indexes_R" $path_script/Step_6B_16S.R


sed 's/"//g' OTU_table_all.txt > OTU_table_all_2.txt
cat OTU_table_all_2.txt > OTU_table_all.txt
rm OTU_table_all_2.txt

sed 's/"//g' OTU_table_taxonomy.txt > OTU_table_taxonomy_2.txt
cat OTU_table_taxonomy_2.txt > OTU_table_taxonomy.txt
rm OTU_table_taxonomy_2.txt


echo "Convert the OTU Table in biom file"
cd $path/OTU_table/
source activate qiime1
touch OTU_table_all.biom
rm OTU_table_all.biom
biom convert -i OTU_table_all.txt -o OTU_table_all.biom --table-type="OTU table" --to-hdf5
touch OTU_table_taxonomy.biom
rm OTU_table_taxonomy.biom
biom convert -i OTU_table_taxonomy.txt -o OTU_table_taxonomy.biom --table-type="OTU table" --to-hdf5  --process-obs-metadata taxonomy



echo "Class endosymbionts vs gut content" 
cd $path/OTU_table/

sed '1 s/_/-/g' OTU_table_all.txt > OTU_table_all_order.txt
sed '1 s/_/-/g' OTU_table_taxonomy.txt | sed '1 s/taxo...c.2../taxonomy/g' > OTU_table_taxonomy_order.txt

grep "Rickettsia\|Rickettsiella\|Wolbachia\|Canditatus\|Cardinium" OTU_table_taxonomy.txt | grep -o "OTU[0-9]*" > endosymbiontes.txt

## 9  OTUs correspond to endosymbionts

cp OTU_table_all.txt test.txt
while read p; do
    grep -v $p test.txt > test2.txt
    rm test.txt
    cat test2.txt > test.txt
    rm test2.txt
done <endosymbiontes.txt
cp test.txt OTU_table_gut.txt
rm test.txt

cp OTU_table_taxonomy.txt test.txt
while read p; do
    grep -v $p test.txt > test2.txt
    rm test.txt
    cat test2.txt > test.txt
    rm test2.txt
done <endosymbiontes.txt
cp test.txt OTU_table_gut_taxonomy.txt
rm test.txt

cp OTU_table_all.txt test.txt
while read p; do
	touch test2.txt
    grep $p test.txt >> test2.txt
done <endosymbiontes.txt
cp test2.txt OTU_table_endo.txt
rm test.txt
rm test2.txt

cp OTU_table_taxonomy.txt test.txt
while read p; do
	touch test2.txt
    grep $p test.txt >> test2.txt
done <endosymbiontes.txt
cp test2.txt OTU_table_endo_taxonomy.txt
rm test.txt
rm test2.txt
# Remove mitochondria # 2 OTUs
grep "mito" OTU_table_taxonomy.txt | grep -o "OTU[0-9]*" > mitochondria.txt
cp OTU_table_endo_taxonomy.txt test.txt
while read p; do
    grep -v $p test.txt > test2.txt
    rm test.txt
    cat test2.txt > test.txt
    rm test2.txt
done <mitochondria.txt
cp test.txt OTU_table_endo_taxonomy.txt
rm test.txt

cp OTU_table_endo.txt test.txt
while read p; do
    grep -v $p test.txt > test2.txt
    rm test.txt
    cat test2.txt > test.txt
    rm test2.txt
done <mitochondria.txt
cp test.txt OTU_table_endo.txt
rm test.txt

#### Add header (automatic)
head -n 1 OTU_table_gut_taxonomy.txt > OTU_table_endo_taxonomy_2.txt
cat OTU_table_endo_taxonomy.txt >> OTU_table_endo_taxonomy_2.txt

cat OTU_table_endo_taxonomy_2.txt > OTU_table_endo_taxonomy.txt
rm OTU_table_endo_taxonomy_2.txt

head -n 1 OTU_table_gut.txt > OTU_table_endo_2.txt
cat OTU_table_endo.txt >> OTU_table_endo_2.txt

cat OTU_table_endo_2.txt > OTU_table_endo.txt
rm OTU_table_endo_2.txt

echo "Convert the OTU Table in biom file"
cd $path/OTU_table/
source activate qiime1
touch OTU_table_gut.biom
rm OTU_table_gut.biom
biom convert -i OTU_table_gut.txt -o OTU_table_gut.biom --table-type="OTU table" --to-hdf5

touch OTU_table_gut_taxonomy.biom
rm OTU_table_gut_taxonomy.biom
biom convert -i OTU_table_gut_taxonomy.txt -o OTU_table_gut_taxonomy.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy

touch OTU_table_endo.biom
rm OTU_table_endo.biom
biom convert -i OTU_table_endo.txt -o OTU_table_endo.biom --table-type="OTU table" --to-hdf5

touch OTU_table_endo_taxonomy.biom
rm OTU_table_endo_taxonomy.biom
biom convert -i OTU_table_endo_taxonomy.txt -o OTU_table_endo_taxonomy.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy




######################################
echo "Step 7: Make a tree with all the OTU sequences"
cd $path/OTU_table/
cp $path/results_wNC/OTU_16s.fa $path/OTU_table/
source activate qiime1
#align_seqs.py -i OTU_16s.fa -m muscle -p 0.1 -o rep_set_align
align_seqs.py -i OTU_16s.fa -o rep_set_align
filter_alignment.py -i rep_set_align/OTU_16s_aligned.fasta -s -g 0.9 -o filtered_alignment
make_phylogeny.py -i filtered_alignment/OTU_16s_aligned_pfiltered.fasta -o OTU_16S.tre



######################################
echo "Step 8: Rarefaction"
source activate qiime1
cd $path/OTU_table/

# from m to x sequences by sample (by s), n repetitions (calculate mean)
multiple_rarefactions.py -i OTU_table_all.biom -m 400 -x 8000 -s 400 -n 20 -o ../diversity/rarefaction/
multiple_rarefactions.py -i OTU_table_endo.biom -m 10 -x 1000 -s 100 -n 20 -o ../diversity/rarefaction_endo/
multiple_rarefactions.py -i OTU_table_gut.biom -m 1000 -x 4000 -s 500 -n 20 -o ../diversity/rarefaction_gut/




######################################
echo "Step 9: Alpha-diversity"

echo "Use rarefied data with a tree"
source activate qiime1
cd $path/diversity/
alpha_diversity.py -i rarefaction/ -o alpha_div/ -m observed_species,shannon,chao1,PD_whole_tree -t $path/OTU_table/OTU_16S.tre
alpha_diversity.py -i rarefaction_endo/ -o alpha_div_endo/ -m observed_species,shannon,chao1,PD_whole_tree -t $path/OTU_table/OTU_16S.tre
alpha_diversity.py -i rarefaction_gut/ -o alpha_div_gut/ -m observed_species,shannon,chao1,PD_whole_tree -t $path/OTU_table/OTU_16S.tre



echo "Visualisation Alpha Diversity"
cd $path/diversity/
collate_alpha.py -i alpha_div/ -o alpha_diversity/
collate_alpha.py -i alpha_div_endo/ -o alpha_diversity_endo/
collate_alpha.py -i alpha_div_gut/ -o alpha_diversity_gut/
make_rarefaction_plots.py -i alpha_diversity/ -m $path/metadata/metadata_populations.txt -o alpha_diversity/plot
make_rarefaction_plots.py -i alpha_diversity_endo/ -m $path/metadata/metadata_populations.txt -o alpha_diversity_endo/plot
make_rarefaction_plots.py -i alpha_diversity_gut/ -m $path/metadata/metadata_populations.txt -o alpha_diversity_gut/plot



######################################
echo "Step 10: Beta-diversity"

echo "Use rarefied data with a tree"
source activate qiime1
cd $path/diversity/
beta_diversity.py -i rarefaction/ -o beta_diversity/ -m bray_curtis,unweighted_unifrac,unweighted_unifrac_full_tree,weighted_normalized_unifrac,weighted_unifrac -t $path/OTU_table/OTU_16S.tre

echo "Principal coordinate analysis"
cd $path/diversity/beta_diversity/

mkdir weighted_unifrac_rarefaction_3200/
mkdir unweighted_unifrac_rarefaction_3200/
mkdir weighted_normalized_unifrac_rarefaction_3200/
mkdir bray_curtis_rarefaction_3200
mkdir unweighted_unifrac_full_tree_rarefaction_3200

mv weighted_unifrac_rarefaction_3200* weighted_unifrac_rarefaction_3200/
mv unweighted_unifrac_rarefaction_3200* unweighted_unifrac_rarefaction_3200/
mv weighted_normalized_unifrac_rarefaction_3200* weighted_normalized_unifrac_rarefaction_3200/
mv bray_curtis_rarefaction_3200* bray_curtis_rarefaction_3200/
mv unweighted_unifrac_full_tree_rarefaction_3200* unweighted_unifrac_full_tree_rarefaction_3200/

principal_coordinates.py -i weighted_unifrac_rarefaction_3200 -o weighted_unifrac_3200_pcoa
principal_coordinates.py -i unweighted_unifrac_rarefaction_3200 -o unweighted_unifrac_3200_pcoa
principal_coordinates.py -i weighted_normalized_unifrac_rarefaction_3200 -o weighted_normalized_unifrac_3200_pcoa
principal_coordinates.py -i bray_curtis_rarefaction_3200 -o bray_curtis_3200_pcoa
principal_coordinates.py -i unweighted_unifrac_full_tree_rarefaction_3200 -o unweighted_unifrac_full_tree_3200_pcoa

make_2d_plots.py -i $path/diversity/beta_diversity/weighted_unifrac_3200_pcoa -m $path/metadata/metadata_populations.txt
make_2d_plots.py -i $path/diversity/beta_diversity/unweighted_unifrac_3200_pcoa -m $path/metadata/metadata_populations.txt
make_2d_plots.py -i $path/diversity/beta_diversity/weighted_normalized_unifrac_3200_pcoa -m $path/metadata/metadata_populations.txt
make_2d_plots.py -i $path/diversity/beta_diversity/bray_curtis_3200_pcoa -m $path/metadata/metadata_populations.txt
make_2d_plots.py -i $path/diversity/beta_diversity/unweighted_unifrac_full_tree_3200_pcoa -m $path/metadata/metadata_populations.txt


echo "Dissimilarity matrix"
source activate qiime1
cd $path/diversity/beta_diversity/
dissimilarity_mtx_stats.py -i weighted_unifrac_rarefaction_3200/ -o stats_weighted_unifrac_3200
dissimilarity_mtx_stats.py -i bray_curtis_rarefaction_3200/ -o bray_curtis_unifrac_3200


make_distance_boxplots.py -d bray_curtis_rarefaction_3200/bray_curtis_rarefaction_3200_12.txt -f Island -o boxplot_bray_curtis_island_3200 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d bray_curtis_rarefaction_3200/bray_curtis_rarefaction_3200_12.txt -f Population -o boxplot_bray_curtis_pop_weighted_unifrac_3200 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d unweighted_unifrac_rarefaction_3200/unweighted_unifrac_rarefaction_3200_12.txt -f Island -o boxplot_unweighted_unifrac_3200 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d weighted_normalized_unifrac_rarefaction_3200/weighted_normalized_unifrac_rarefaction_3200_12.txt -f Island -o boxplot_weighted_normalized_unifrac_3200 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d unweighted_unifrac_full_tree_rarefaction_3200/unweighted_unifrac_full_tree_rarefaction_3200_12.txt -f Island -o boxplot_unweighted_unifrac_full_tree_3200 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d weighted_unifrac_rarefaction_3200/weighted_unifrac_rarefaction_3200_12.txt -f Island -o boxplot_weighted_unifrac_3200 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d weighted_unifrac_rarefaction_3200/weighted_unifrac_rarefaction_3200_12.txt -f Population -o boxplot_weighted_pop_unifrac_3200 -m $path/metadata/metadata_populations.txt



echo "Beta-diversity  endo "

cd $path/diversity/
beta_diversity.py -i rarefaction_endo/ -o beta_diversity_endo/ -m bray_curtis,weighted_unifrac -t $path/OTU_table/OTU_16S.tre
cd $path/diversity/beta_diversity_endo/
mkdir weighted_unifrac_rarefaction_110/
mkdir bray_curtis_rarefaction_110
mv weighted_unifrac_rarefaction_110* weighted_unifrac_rarefaction_110/
mv bray_curtis_rarefaction_110* bray_curtis_rarefaction_110/
principal_coordinates.py -i weighted_unifrac_rarefaction_110 -o weighted_unifrac_110_pcoa
principal_coordinates.py -i bray_curtis_rarefaction_110 -o bray_curtis_110_pcoa
make_2d_plots.py -i $path/diversity/beta_diversity_endo/weighted_unifrac_110_pcoa -m $path/metadata/metadata_populations.txt
make_2d_plots.py -i $path/diversity/beta_diversity_endo/bray_curtis_110_pcoa -m $path/metadata/metadata_populations.txt
dissimilarity_mtx_stats.py -i weighted_unifrac_rarefaction_110/ -o stats_weighted_unifrac_110
make_distance_boxplots.py -d bray_curtis_rarefaction_110/bray_curtis_rarefaction_110_12.txt -f Island -o boxplot_bray_curtis_island_110 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d bray_curtis_rarefaction_110/bray_curtis_rarefaction_110_12.txt -f Population -o boxplot_bray_curtis_pop_weighted_unifrac_110 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d weighted_unifrac_rarefaction_110/weighted_unifrac_rarefaction_110_12.txt -f Island -o boxplot_weighted_unifrac_110 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d weighted_unifrac_rarefaction_110/weighted_unifrac_rarefaction_110_12.txt -f Population -o boxplot_weighted_pop_unifrac_110 -m $path/metadata/metadata_populations.txt



echo "Beta-diversity  gut "

cd $path/diversity/
beta_diversity.py -i rarefaction_gut/ -o beta_diversity_gut/ -m bray_curtis,weighted_unifrac -t $path/OTU_table/OTU_16S.tre
cd $path/diversity/beta_diversity_gut/
mkdir weighted_unifrac_rarefaction_3000/
mkdir bray_curtis_rarefaction_3000
mv weighted_unifrac_rarefaction_3000* weighted_unifrac_rarefaction_3000/
mv bray_curtis_rarefaction_3000* bray_curtis_rarefaction_3000/
principal_coordinates.py -i weighted_unifrac_rarefaction_3000 -o weighted_unifrac_3000_pcoa
principal_coordinates.py -i bray_curtis_rarefaction_3000 -o bray_curtis_3000_pcoa
make_2d_plots.py -i $path/diversity/beta_diversity_gut/weighted_unifrac_3000_pcoa -m $path/metadata/metadata_populations.txt
make_2d_plots.py -i $path/diversity/beta_diversity_gut/bray_curtis_3000_pcoa -m $path/metadata/metadata_populations.txt
dissimilarity_mtx_stats.py -i weighted_unifrac_rarefaction_3000/ -o stats_weighted_unifrac_3000
make_distance_boxplots.py -d bray_curtis_rarefaction_3000/bray_curtis_rarefaction_3000_12.txt -f Island -o boxplot_bray_curtis_island_3000 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d bray_curtis_rarefaction_3000/bray_curtis_rarefaction_3000_12.txt -f Population -o boxplot_bray_curtis_pop_weighted_unifrac_3000 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d weighted_unifrac_rarefaction_3000/weighted_unifrac_rarefaction_3000_12.txt -f Island -o boxplot_weighted_unifrac_3000 -m $path/metadata/metadata_populations.txt
make_distance_boxplots.py -d weighted_unifrac_rarefaction_3000/weighted_unifrac_rarefaction_3000_12.txt -f Population -o boxplot_weighted_pop_unifrac_3000 -m $path/metadata/metadata_populations.txt



echo "Plot the abundance chart"
source activate qiime1
cd $path/OTU_table/
summarize_taxa_through_plots.py -f -o taxa_all_taxonomy -i OTU_table_taxonomy.biom -m $path/metadata/metadata_populations.txt -c Population
summarize_taxa_through_plots.py -f -o taxa_endo_taxonomy -i OTU_table_endo_taxonomy.biom -m $path/metadata/metadata_populations.txt -c Population
summarize_taxa_through_plots.py -f -o taxa_gut_taxonomy -i OTU_table_gut_taxonomy.biom -m $path/metadata/metadata_populations.txt -c Population 

summarize_taxa_through_plots.py -f -o taxa_all_taxonomy_order_sort -i OTU_table_taxonomy.biom -m $path/metadata/metadata_populations.txt -s
summarize_taxa_through_plots.py -f -o taxa_endo_taxonomy_order_sort -i OTU_table_endo_taxonomy.biom -m $path/metadata/metadata_populations.txt -s
summarize_taxa_through_plots.py -f -o taxa_gut_taxonomy_order_sort -i OTU_table_gut_taxonomy.biom -m $path/metadata/metadata_populations.txt -s


