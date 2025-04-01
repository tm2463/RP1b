"""
The program extracts sequences from fasta files based on user defined start and end genes
The purpose of this program is to aid exploratory analysis but extracting sections identified as containing recombination
The output is a multifasta file that can be aligned with mauve
SNPs can then be exported and processed through Mauve.py
"""

from Bio import SeqIO


def genbank_parser():
    file_order = ['SPARK_516', 'SPARK_1954', 'SPARK_1059', 'SPARK_2540', 'SPARK_1507', 'SPARK_1852']
    gene_name = []
    start_pos = []
    end_pos = []
    for name in file_order:
        with open('genbank_files/' + name + '.gbk', 'r') as gbk:
            genbank_recs = SeqIO.parse(gbk, "genbank")
            genbank = next(genbank_recs, None)
            gene_name_seq = []
            start_seq = []
            end_seq = []
            for item in genbank.features:
                if item.type == "CDS":
                    if item.qualifiers.get("gene", ["Unknown"])[0] != 'Unknown':
                        gene_name_seq.append(item.qualifiers.get("gene", ["Unknown"])[0])
                    else:
                        gene_name_seq.append(item.qualifiers.get("locus_tag", ["Unknown"])[0])
                    start_seq.append(int(item.location.start))
                    end_seq.append(int(item.location.end))
        gene_name.append(gene_name_seq)
        start_pos.append(start_seq)
        end_pos.append(end_seq)
    return gene_name, start_pos, end_pos


def fasta_parser():
    file_order = ['SPARK_516_C1_v2', 'SPARK_1954_C1_V2', 'SPARK_1059_C1', 'SPARK_2540_C2_V2', 'SPARK_1507_C1', 'SPARK_1852_C1']
    sequences = []
    for name in file_order:
        with open('fasta_files/' + name + '.fasta', 'r') as f:
            file = f.read().split('\n')[1:]
            file = ''.join([x for x in file])
            sequences.append(file)
    return sequences


def get_seqs(file_name, first, last):
    file_dict = {0: 'SPARK_516_C1_V2', 1: 'SPARK_1954_C1_V2', 2: 'SPARK_1059_C1', 3: 'SPARK_2540_C2_V2', 4: 'SPARK_1507_C1', 5: 'SPARK_1852_C1'}
    data = genbank_parser()
    for genome in range(6):
        sequences = fasta_parser()[genome]
        gene_name_list = data[0][genome]
        start_list = data[1][genome]
        end_list = data[2][genome]
        start_index = gene_name_list.index(first)
        end_index = gene_name_list.index(last)
        start_coord = start_list[start_index]
        end_coord = end_list[end_index]
        with open('sequence_sections/' + file_name + '.txt', 'a') as txt:
            txt.write('>' + file_dict[genome] + '\n')
            txt.write(sequences[start_coord:end_coord] + '\n')


# 3 inputs for the function: Outfile file name, the first gene name and the last gene name
# in order to work the first and last genes must be common among all genomes
get_seqs('recomb_check_band_1', 'barA', 'dld')
