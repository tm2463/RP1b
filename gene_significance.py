"""
This program takes the export SNP output from Mauve and returns the top 10 genes sorted by SNP density
File handling is the same as in Mauve.py
The GenBank files for which the alignment was created are needed for the GenBank parser
"""

import numpy as np
import pandas as pd
from collections import Counter
from Bio import SeqIO


# file handling; remove redundant columns from snp file; store SNPs and their locations
with open('mauve_snp', 'r') as f:
    file = np.loadtxt(f, delimiter='\t', dtype=str)
    file = np.delete(file, 0, 0)
    columns = file.shape[1]
    rm_list = [n for n in range(columns) if n % 3 != 0]
    data = pd.DataFrame(file).drop(rm_list, axis=1)
    snps = np.array(data[0], dtype=str)
    data = pd.DataFrame(data.drop(0, axis=1))

# assign SNP location to correct genome
snp_filter = np.zeros(shape=data.shape, dtype=int)
for snp, x in zip(snps, range(data.shape[0])):
    count = Counter(snp).most_common(1)[0][0]
    for base, y in zip(snp, range(data.shape[1])):
        if base == count:
            snp_filter[x, y] += 0
        else:
            snp_filter[x, y] += 1
data = data.mul(snp_filter)
data.replace('', 0, inplace=True)
data = np.array(data, dtype=int)
data_swap = np.swapaxes(data, 1, 0)


# parse genbank files and extract relevant features
def genbank_parser():
    file_order = ['SPARK_516', 'SPARK_1954', 'SPARK_1059', 'SPARK_2540', 'SPARK_1507', 'SPARK_1852']
    locus = []
    end_pos = []
    full_gene_name = []
    translation = []
    start_pos = []
    for name in file_order:
        with open('genbank_files/' + name + '.gbk', 'r') as gbk:
            genbank_recs = SeqIO.parse(gbk, "genbank")
            genbank = next(genbank_recs, None)
            locus_seq = []
            end_seq = []
            full_gene_name_seq = []
            translation_seq = []
            start_seq = []
            for item in genbank.features:
                if item.type == "CDS":
                    locus_seq.append(item.qualifiers.get("locus_tag", ["Unknown"])[0])
                    end_seq.append(int(item.location.end))
                    full_gene_name_seq.append(item.qualifiers.get("product", ["Unknown"])[0])
                    translation_seq.append(item.qualifiers.get("translation", ["Unknown"])[0])
                    start_seq.append(int(item.location.start))
        locus.append(locus_seq)
        end_pos.append(end_seq)
        full_gene_name.append(full_gene_name_seq)
        translation.append(translation_seq)
        start_pos.append(start_seq)
    return locus, end_pos, full_gene_name, translation, start_pos


# create dictionary for gene names and snp density
def gene_snp_density():
    gene_info = genbank_parser()
    snp_array = data_swap
    start_pos = gene_info[4]
    end_pos = gene_info[1]
    locus = gene_info[0]
    density_dict = {}
    for n in range(len(start_pos)):
        genome_snp_array = snp_array[n]
        genome_start_pos = start_pos[n]
        genome_end_pos = end_pos[n]
        genome_locus = locus[n]
        genome_snp_array.sort()
        for x in range(len(genome_start_pos)):
            gene_start = genome_start_pos[x]
            gene_end = genome_end_pos[x]
            gene_locus = genome_locus[x]
            gene_length = gene_end - gene_start
            start = np.searchsorted(genome_snp_array, gene_start, side="left")
            end = np.searchsorted(genome_snp_array, gene_end, side="right")
            snp_count = end - start
            gene_snp_density = snp_count / gene_length
            density_dict[(n, str(gene_locus))] = float(gene_snp_density)
    return density_dict


# dictionaries to translate gene locus to gene name or protein sequence
def translate_dicts():
    data = genbank_parser()
    locus = data[0]
    name = data[2]
    translation = data[3]
    my_dict = {}
    protein_dict = {}
    for x, y in zip(locus, name):
        for i, j in zip(x, y):
            my_dict[str(i)] = j
    for a, b in zip(locus, translation):
        for c, d in zip(a, b):
            protein_dict[str(c)] = d
    return my_dict, protein_dict


# write file with top genes
def sort_dict(n):
    gene_snp_density_data = gene_snp_density()
    density_dict = gene_snp_density_data
    genome_dict = {}
    for (genome, locus), density in density_dict.items():
        if genome not in genome_dict:
            genome_dict[genome] = []
        genome_dict[genome].append((locus, density))
    top_genes = {}
    for genome, values in genome_dict.items():
        sorted_values = sorted(values, key=lambda x: x[1], reverse=True)[:n]
        top_genes[genome] = sorted_values
    translate = translate_dicts()
    with open('top_10_genes_2.txt', 'w') as file:
        for genome, values in top_genes.items():
            for gene_locus, density in values:
                file.write(f"{genome}\t{translate[0][gene_locus]}\t{density:.6f}\t{translate[1][gene_locus]}")
                file.write('\n')


sort_dict(10)
