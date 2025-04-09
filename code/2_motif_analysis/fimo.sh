#!/bin/bash

#FIMO RUNS OF MEME MOTIFS (***meme.txt) AGAINST FULL 20KB PROMOTERS
#Author: Yvet Boele
#Goal: Script searches for previously detected motifs across promoters of all PLT clades using the tool FIMO
#Output: a folder with all FIMO output per motif x promoter type combination 
#Script should be run in command line

promoter_directory=/lustre/BIF/nobackup/boele028/SA_2023/promoters_PLT/111124_20kb-cds_plt-chopped/fimo/promoters
for motif in ./motifs/*/*/*meme.txt; do
    for promoter in $promoter_directory/*.fasta; do 
        motif_name=$(basename ${motif%_meme.txt})
        promoter_name=$(basename $promoter | sed 's/\_full_promoters_minlen50.fasta//')
	output_directory=${motif_name}_against_${promoter_name}
	echo $motif_name
	echo $promoter_name
	echo $output_directory
        mkdir ./fimo_output/$output_directory
    	fimo --oc ./fimo_output/$output_directory $motif $promoter
    done
done

#for fi in ./fimo_output/*/fimo.tsv; 
#do filename=$(echo $fi | sed 's/\.\/fimo_output\///' | sed 's/\/fimo//'); 
#mv $fi ./ALL_FIMO_TSVs/$filename; 
#done
