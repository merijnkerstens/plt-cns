#!/bin/bash

#MEME RUNS DICOTS
#Author: Yvet Boele
#Goal: Script searches for conserved elements across promoters of all PLT clades of dicots species using the tool MEME
#Output: per PLT clade a list of conserved motifs 
#Script should be run in command line.

directory="/lustre/BIF/nobackup/boele028/SA_2023/promoters_PLT/111124_20kb-cds_plt-chopped/prepping_promoters_dicots"
# List of input FASTA files

# Loop through each input file to run a zoops MEME
for fasta_file in ${directory}/*_truncated_promoters_minlen50.fasta; do
    # # Extract file name without extension
    file_name=$(basename -- "${fasta_file}")
    file_name_no_ext="${file_name%.*}"
    # Run MEME with the specified parameters
    meme ${fasta_file} -dna -oc "./zoops_dicot_${file_name_no_ext}" -mod zoops -nmotifs 15 -minw 15 -maxw 100 -objfun classic -minsites 15 -markov_order 0 

    # Print a separator for better readability
    echo "----------------------------------------"
done

#meme ${directory}/PLT3-7_truncated_promoters_minlen50.fasta -dna -oc "./zoops_dicot_30motifs_PLT3-7_truncated_promoters_minlen50" -mod zoops -nmotifs 30 -minw 15 -maxw 100 -objfun classic -minsites 15 -markov_order 0
