"""
TITLE
Authors: Yvet Boele

does this and this

Output:

Script should be run as <python>
requires python3.9
"""

import os
import subprocess
import csv
import pandas as pd
import argparse


def plt_type(filename):
    # Extract species and PLT from the filename
    plt_type = filename.split('_')[6].split('/')[1]
    return plt_type

def parse_and_grep_all_homologs(filename, gffs):
    with open(filename, 'r') as file:
        elements = [line.strip() for line in file]
    # Prepare a list to store the results
    results = []
    gff_directory_pattern = gffs 
    #elements = ['Oeu052052.1', 'MBI01_g39230_MAGBIO', 'AT1G01030']
    #elements = ['Nitab4.5_0005002g0010']
    # Grep each element in the specified directory and capture the output
    for element in elements:
        # Execute grep command
        command = f"grep '{element}' {gff_directory_pattern}"

        grep_result = subprocess.run(
            command,
            text=True,
            capture_output=True,
            shell=True
        )
        
        # Check if there is any output
        if grep_result.stdout:
            # get species
            # Split the output into lines and each line into tab-separated entries
            lines = grep_result.stdout.strip().split('\n')
            for line in lines:
                grep_file, grep_result = line.split('.gff3:')
                species = get_species_from_grep(grep_file)
                entries = grep_result.split('\t')
                entries[8] = entries[8].split(';')
                entries.append(species)
                results.append(entries)
    return results

def get_species_from_grep(grep_filename):
    short_species = grep_filename.split('.')[-1]
    return short_species
    
def gene_id(parsed_grep_output):
    all_genes = []
    for element in parsed_grep_output:
        if element[2] == 'gene':  # Check if the type is 'gene'
            # Split the 9th column by semicolons to get individual attributes
            attributes = element[8]
            
            # Iterate through the attributes to find the one starting with 'gene_id='
            gene_id = None
            for attr in attributes:
                if attr.startswith('gene_id='):
                    gene_id = attr.split('=')[1]  # Get the part after 'gene_id='
                    break  # Stop once we find the gene_id
            
            # Only proceed if a gene_id was found
            if gene_id:
                species = element[-1]  # Extract species from the last column
                gene_tuple = (species, gene_id)  # Create a tuple (species, gene)
                if gene_tuple not in all_genes:  # Ensure uniqueness
                    all_genes.append(gene_tuple)
            else:
                print('problem!')
    return all_genes


def get_CDS_per_gene(parsed_grep_output, gene, species):
    #filtered_list = [entry for entry in parsed_grep_output if entry[2] == 'CDS' and entry[8][0].startswith(gene)]
    #filtered_list = [entry for entry in parsed_grep_output if entry[2] == 'CDS' and gene in entry[8]]
    filtered_list = [
    entry for entry in parsed_grep_output
    if entry[2] == 'CDS' and any(gene in attr for attr in entry[8])
    ]
    df_fl = pd.DataFrame(filtered_list)
    df_fl.iloc[:, 3] = pd.to_numeric(df_fl.iloc[:, 3])
    df_fl.iloc[:, 4] = pd.to_numeric(df_fl.iloc[:, 4])
    chrom = df_fl.iloc[0,0]
    strand = df_fl.iloc[0,6]
    cds_1 = df_fl.iloc[:,3].min()
    cds_2 = df_fl.iloc[:,4].max()
    bedfile_entry = (chrom, cds_1, cds_2, gene.removesuffix("gene_id="), '.', strand)
    output_file = f"{species}.{gene.removeprefix('gene_id=')}.bed" 
    with open(output_file, 'w',) as bedfile:
        bedfile.write('\t'.join(map(str, bedfile_entry))+ '\n')
    return output_file

def get_promoter_bed_per_gene(cds_bed):
    promoter_bed = f"{cds_bed.removesuffix('.bed')}.promoter.bed"
    species = cds_bed.split('.')[0]
    cmd = f"bedtools flank -i {cds_bed} -g {base}PLAZA5.0-files-dicots/genomes/{species}.bed -l 20000 -r 0 -s > {promoter_bed}"
    result = subprocess.run(cmd, shell=True)
    return promoter_bed

def update_promoter_bed_upstream_gene(cds_promoter):
    species = cds_promoter.split('.')[0]
    cmd = f"bedtools intersect -a {base}PLAZA5.0-files-dicots/gffs/{species}/annotation.selected_transcript.all_features.{species}.gff3.tmp -b {cds_promoter}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    results = result.stdout.split('\n')
    real_results = []
    for elem in results:
        if elem:
            real_result = elem.split('\t')
            #don't match to self
            if not cds_promoter.removesuffix('.promoter.bed').removeprefix(species).removeprefix('.') in real_result[8]:
                real_results.append(real_result)
    outfile = f"{cds_promoter.removesuffix('.promoter.bed')}.truncated_promoter.bed"
    with open(cds_promoter, 'r') as f:
        line = f.readline().strip().split('\t')
        with open(outfile, 'w') as outputfile:
            if real_results:
                updated_coordinates = update_coordinates(line, real_results)
                tsv_output = csv.writer(outputfile, delimiter = '\t')
                tsv_output.writerow(updated_coordinates)
            else:
                tsv_output = csv.writer(outputfile, delimiter = '\t')
                tsv_output.writerow(line)

def update_promoter_bed_upstream_PLT(cds_promoter):
    species = cds_promoter.split('.')[0]
    all_plts = get_all_plts()
    cmd = f"bedtools intersect -a {base}PLAZA5.0-files-dicots/gffs/{species}/annotation.selected_transcript.all_features.{species}.gff3.tmp -b {cds_promoter}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    results = result.stdout.split('\n')
    upstream_plt = []
    outfile = f"{cds_promoter.removesuffix('.promoter.bed')}.full_except_PLT_promoter.bed"
    for elem in results:
        if elem:
            real_result = elem.split('\t')
            if not cds_promoter.removesuffix('.promoter.bed').removeprefix(species).removeprefix('.') in real_result[8]:
                for plt in all_plts:
                    if plt in real_result[8]:
                        upstream_plt.append(real_result)
    with open(cds_promoter, 'r') as f:
        line = f.readline().strip().split('\t')
        with open(outfile, 'w') as outputfile:
            if upstream_plt:
                updated_coordinates = update_coordinates(line, upstream_plt)
                tsv_output = csv.writer(outputfile, delimiter = '\t')
                tsv_output.writerow(updated_coordinates)
            else:
                tsv_output = csv.writer(outputfile, delimiter = '\t')
                tsv_output.writerow(line)

def update_coordinates(coordinates, overlapping_genes):
    left_coordinates = []
    right_coordinates = []
    bedfile_strand = coordinates[5]
    for line in overlapping_genes:
        gene_start = line[3]
        gene_end = line[4]
        if bedfile_strand == '+':
            right_coordinates.append(int(gene_end))
        if bedfile_strand == '-':
            left_coordinates.append(int(gene_start))
    if bedfile_strand == '+':
        closest_upstream = max(right_coordinates)
        coordinates[1] = closest_upstream 
    if bedfile_strand == '-':
        closest_upstream = min(left_coordinates)
        coordinates[2] = closest_upstream
    return coordinates

def get_all_plts():
    list_files = [f for f in os.listdir("/lustre/BIF/nobackup/boele028/SA_2023/promoters_PLT/111124_20kb-cds_plt-chopped/prepping_promoters_dicots/") if f.endswith('.list')]
    all_plts = []
    for file in list_files:
        with open(f"{base}prepping_promoters_dicots/{file}", 'r') as f:
            line = f.readlines()
            for elem in line:
                all_plts.append(elem.strip())
    return all_plts

def main(inputfile, gffs):
    plt_clade = plt_type(inputfile)
    hl = parse_and_grep_all_homologs(inputfile, gffs)
    all_genes = gene_id(hl)
    total_genes = len(all_genes)  
    
    for i, (species, gene) in enumerate(all_genes, start=1):
        cds = get_CDS_per_gene(hl, gene, species)
        promoter = get_promoter_bed_per_gene(cds)
        update_promoter_bed_upstream_gene(promoter)
        update_promoter_bed_upstream_PLT(promoter)
        
        # Print progress message
        print(f"[{i}/{total_genes}] {gene} is done")
    
if __name__ == "__main__":
    #PLT filename is an argument when running the script (has to be the full length path)
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help="Path to the input file containing gene data (e.g., homologs list).")
    args = parser.parse_args()
    
    gffs = '/lustre/BIF/nobackup/boele028/SA_2023/promoters_PLT/111124_20kb-cds_plt-chopped/PLAZA5.0-files-dicots/gffs/*/annotation.selected_transcript.all_features.*.gff3'
    base = '/lustre/BIF/nobackup/boele028/SA_2023/promoters_PLT/111124_20kb-cds_plt-chopped/'
    main(args.filename, gffs)
