''' Format GENCODE gff3 file
Description: https://medium.com/intothegenomics/annotate-genes-and-genomic-coordinates-using-python-9259efa6ffc2

### Tested with python 3.6
Usage:
    python gencode_format.py <gencode annotated gff3 file>
'''

import pandas as pd
import numpy as np
import pysam
import sys
import os


## Read the gff3 file into a pandas dataframe
gencode_gff3 = sys.argv[1]
print (f"Reading the GENCODE input file: {gencode_gff3}")
gencode = pd.read_table(gencode_gff3, comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])

## Extract Genes in the gff3 file “feature = gene”
print ("Extract Genes...")
gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)

## Extract gene_name, gene_type, gene_status, level of each gene
print ("Extract gene_name, gene_type, gene_status, level of each gene...")
def gene_info(x):
    g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
    g_type = list(filter(lambda x: 'gene_type' in x,  x.split(";")))[0].split("=")[1]
    g_status = list(filter(lambda x: 'gene_status' in x,  x.split(";")))[0].split("=")[1]
    g_leve = int(list(filter(lambda x: 'level' in x,  x.split(";")))[0].split("=")[1])
    return (g_name, g_type, g_status, g_leve)

gencode_genes["gene_name"], gencode_genes["gene_type"], gencode_genes["gene_status"], gencode_genes["gene_level"] = zip(*gencode_genes.attribute.apply(lambda x: gene_info(x)))

## Extract all known protein_coding genes
print ("Extract all known protein_coding genes...")
gencode_genes = gencode_genes[gencode_genes['gene_status'] == 'KNOWN'].reset_index().drop('index', axis=1)
gencode_genes = gencode_genes[gencode_genes['gene_type'] == 'protein_coding'].reset_index().drop('index', axis=1)

## Remove duplicates — Prioritize verified and manually annotated loci over automatically annotated loci
print ("Remove duplicates...")
gencode_genes = gencode_genes.sort_values(['gene_level', 'seqname'], ascending=True).drop_duplicates('gene_name', keep='first').reset_index().drop('index', axis=1)

## Generating the indexed GENCODE gene source dataset
print ("Generating the TABIX indexed file...")
gencode_genes.to_csv("gencode.v19.annotation.gff3_all_known_genes.txt", header=False, sep="\t", index=False)
os.system("cut -f 1,2,3,5 gencode.v19.annotation.gff3_all_known_genes.txt | sortBed -i > gencode.v19.annotation.gff3_all_known_genes.txt.sorted.formatted.bed")
os.system("bgzip gencode.v19.annotation.gff3_all_known_genes.txt.sorted.formatted.bed")
os.system("tabix -p bed gencode.v19.annotation.gff3_all_known_genes.txt.sorted.formatted.bed.gz")
