''' Find genes in a set of genomic coordinates
Description: https://medium.com/intothegenomics/annotate-genes-and-genomic-coordinates-using-python-9259efa6ffc2

### Tested with python 3.6
Usage:
    python <source file from gencode_format.py script> <bed file>
    bed file: Accepts 6 column bed file, if you have a bed file in a different format, please change the 46th line accordingly.
'''

import pandas as pd
import numpy as np
import pysam
import sys
import os

## Calculate the overlap
def overlap(q_st, q_end, res_st, res_end):
    o  = min(q_end, res_end)-max(q_st, res_st)
    return o

## Find genes
def gencode_all_known_genes(a, tb):
    genes = []

    try:
        for region in tb.fetch(a['chr'], int(a['start']), int(a['end'])):
            if region:
                r = region.split('\t')
                overlap_len = overlap(int(a['start']), int(a['end']), int(r[1]), int(r[2]))
                ret_val = '{}({})'.format(r[3], np.round(overlap_len/float(int(a['end'])-int(a['start']))*100, 2)) ### Percentage of the input interval that overlap with the gene
                genes.append(ret_val)

        if len(genes)>0:
            return ";".join(genes)
        else:
            return "NA(0)"
    except ValueError:
        return "NA(0)"

## Read TABIX annotated file GENCODE gene file...
print ("Read TABIX annotated input file...")
gencode_v19 = pysam.TabixFile(sys.argv[1])

## Read input bed file
print ("Read input bed file...")
df = pd.read_table(sys.argv[2], names=['chr', 'start', 'end', 'name', 'score', 'strand'])

print ("Looking for genes...")
df['genes'] = df.apply(lambda x: gencode_all_known_genes(x[['chr', 'start', 'end']], gencode_v19), axis=1)

## Remove all the intervals that do not overlap with genes
df = df[df['genes'] != "NA(0)"].reset_index(drop=True)

print ("Transforming the dataset so that each row will contain a single gene...")
new_rows = []
for i,r in df.iterrows():
    g_list = r['genes'].split(";")
    for g in g_list:
        g = g.replace(" ","")
        new_rows.append(np.append(r[['chr', 'start', 'end', 'name', 'score', 'strand', 'genes']].values, g))

df_perGene = pd.DataFrame()
df_perGene = df_perGene.append(pd.DataFrame(new_rows, columns=['chr', 'start', 'end', 'name', 'score', 'strand', 'genes', 'gene_ID'])).reset_index().drop('index', axis=1)

df_perGene['gene_name'] = df_perGene['gene_ID'].apply(lambda x: x.split("(")[0])
df_perGene['gene_coverage'] = df_perGene['gene_ID'].apply(lambda x: x.split("(")[1].replace(")", ""))
df_perGene = df_perGene.drop(["genes", "gene_ID"], axis=1)

df_perGene.to_csv("genes.txt", header=True, index=False, sep="\t")
