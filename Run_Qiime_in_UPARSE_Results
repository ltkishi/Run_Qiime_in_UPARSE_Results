#!/bin/bash

## This is a bash script to process Qiime in Uparse results from amplicon sequences sequenced with MiSeq or Hiseq. 

### 

clear

echo "-------------------------------------------------------------------"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo " "
echo " "
echo " adapted workflow by Ramiro Logares https://doi.org/10.5281/zenodo.259579 and DOI:10.1111/1462-2920.12250"
echo " adapted by Luciano Takeshi Kishi  September  2016                                  "
echo " Works with Usearch v8x | 32bits [Another version for vsearch v1.5x] "
echo " Qiime version 1.9.1+dfsg-1biolinux4"
echo " "
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "--------------------------------------------------------------------"

##### Steps #####
#
# Run UPARSE_workflow_97clust_v1.5_usearch_no_singletons_pear_BayesHammer.sh of Ramiro Logares
# Make mapping.txt file in the results directory 
# Run Run_Qiime_in_UPARSE_Results.sh in zsh qiime
# Open the results in core_diversity folder
#
################# 

##### Databases #####
rdp_fas="/data/db/mothur/rdp/trainset14_032015.rdp.fasta"
rdp_tax="/data/db/mothur/rdp/trainset14_032015.rdp.tax"

silva="/data/db/mothur/silva/silva.nr_v119.align"

#####  #####

## Checking the bmp-map2qiime.py in the $PATH (www.brmicrobiome.org) 
bmp_map=$(which bmp-map2qiime.py)

## Edit line 91 for the correct parameter, for example: group.




date


echo "Classifying sequences with mothur, 80% confidence"

assign_taxonomy.py -i otus97_repset_clean.fa -m mothur -c 0.8 -o output_taxonomy -r ${rdp_fas} -t ${rdp_tax}


echo "Aligning sequences with Muscle using silva.nr_v119.align"

align_seqs.py -m muscle -i otus97_repset_clean.fa -o rep_set_align -t ${silva} 


echo "Filtering alignment ... "

filter_alignment.py -i rep_set_align/otus97_repset_clean_aligned.fasta --remove_outliers -o filtered_alignment


echo "Building the Tree ... "

make_phylogeny.py -i filtered_alignment/otus97_repset_clean_aligned_pfiltered.fasta -o rep_set.tre


echo "Building otu_table.txt ..."

cp map97.uc map.uc

perl -i -pe 's/.seqnum/_/g' map.uc 

python ${bmp_map} > otu_table.txt

echo "Building OTU table ... "

make_otu_table.py -i otu_table.txt -t output_taxonomy/otus97_repset_clean_tax_assignments.txt -o otu_table.biom


echo "Building table indexes ... "

alpha_diversity.py -i otu_table.biom -m chao1,ace,observed_otus,equitability,simpson,shannon,observed_species -o adiv_chao1_pd.txt -t rep_set.tre


echo "Buildind OTU table groups ... "

biom summarize-table -i otu_table.biom -o results_biom_table.txt


echo "Building Core Diversity Analysis ... "

cat results_biom_table.txt | perl -e 'while(<>){if(/Min: (\d+)/){system(`core_diversity_analyses.py -i otu_table.biom -o core_diversity -m mapping.txt -t rep_set.tre -e $1`)}}'


echo "Analyzing significant groups ... "

group_significance.py -i otu_table.biom -m mapping.txt -c group -o group_significance.txt -s ANOVA


date

echo "Finished ... "
