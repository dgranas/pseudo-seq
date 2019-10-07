# 10/6/19
# Dave Granas
# Ribomethseq end-counts 

import pysam
import collections
import pandas as pd
import os
import numpy as np
import sys

def get_rrna_abbreviations(bamfile):
    '''
    make a dict to convert long reference names from bam file
    into shorter abbreviations for the corresponding rRNA genes
    '''
    # get the reference names from a bam file
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    ref_names = samfile.references

    # make a dict to convert the long name to a short abbreviation (if aligned to rRNA)
    rrnas = ('5-8S','18S','28S')
    ref_name_to_abbrev = {}
    for ref_name in ref_names:
        found_rrnas = [rrna for rrna in rrnas if rrna in ref_name]
        if len(found_rrnas) > 1:
            raise SystemError(f'Found multiple rRNA names in {ref_name}')
        if not found_rrnas:
            print(f'Warning: did not find expected rRNA name in {ref_name}')
            ref_name_to_abbrev[ref_name] = ref_name
        else:
            ref_name_to_abbrev[ref_name] = found_rrnas[0]

    return ref_name_to_abbrev

def get_end_count(bamfile):
    '''
    Given: bamfile for a particular sample
    Gets 5' end of R1 (upstream) reads and 3' of R2 (downstream) reads 
    '''
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    reads = samfile.fetch()
    
    print(bamfile)
    genes = set()
    gene_to_5p = {}
    gene_to_3p = {}
    
    for read in reads:
        gene = read.reference_name
        if gene not in genes:
            genes.add(gene)
            gene_to_5p[gene] = collections.defaultdict(int)
            gene_to_3p[gene] = collections.defaultdict(int)
            
        start_pos = read.reference_start + 1 # add 1 since pysam is 0-based
        '''
        from Birkedal et al:
        The nucleotides at the 3′ read‐ends directly depend on their 2′‐OH function, 
        whereas the nucleotides at the 5′ read‐ends result from the 2′‐OH function of 
        their 5′ neighbors. Thus, the 5′ read‐ends are all shifted one position upstream 
        in the data treatment such that reads corresponding to a methylated nucleotide 
        will align in the two datasets.
        '''
        start_pos -= 1 # shift start 1bp upstream due to fragmentation chemistry mentioned above
        stop_pos = read.reference_end
        
        if read.is_read1 and not read.is_reverse:
            gene_to_5p[gene][start_pos] += 1
        elif read.is_read2 and read.is_reverse:
            gene_to_3p[gene][stop_pos] += 1
            
    print(genes)
    return gene_to_5p, gene_to_3p

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
