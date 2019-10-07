# 10/6/19
# Dave Granas
# Ribomethseq end-counts 

import pysam
import collections
import pandas as pd
import os
import numpy as np
import sys

def get_rrna_abbreviations():
    # TODO
    pass

def get_end_count(bamfile):
    '''
    Given: bamfile for a particular sample
    Gets 5' end of R1 (upstream) reads and 3' of R2 (downstream) reads 
    '''
    pass
    # to do

if len(sys.argv) != 3:
    raise SystemError('usage: python ribomethseq.py [bamfile] [output folder]')

bamfile, output_folder = sys.argv[1:]

# if folder exists, make sure it is empty
if os.path.exists(output_folder):
    if os.listdir(output_folder):
        raise SystemError(f'{output_folder} is not empty folder')
else: # make folder if doesn't exist already
    os.mkdir(output_folder)

g5p, g3p = get_end_count(bamfile)

for gene, val5, val3 in zip(g5p.keys(), g5p.values(), g3p.values()):
    print(gene)
    df = pd.DataFrame()
    df['5-prime'] = pd.Series(g5p[gene])
    df['3-prime'] = pd.Series(g3p[gene])
    df = df.fillna(0)
    df = df.iloc[1:]
    max_index = max(df.index)
    new_index = pd.Index(np.arange(1,max_index+1))
    df = df.reindex(new_index)

    df.to_csv(os.path.join(output_folder, f'{gene}.csv'))