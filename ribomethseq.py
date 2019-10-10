# 10/6/19
# Dave Granas
# Ribomethseq end-counts 

import pysam
import collections
import pandas as pd
import os
import numpy as np
import sys

def a_score(s):
    '''
    Given a series of numbers containing left-flank - n_i - right-flank
    Calculate the A score from Birkedal et al
    Supplementary page 6
    '''
    s.reset_index(drop=True, inplace=True)

    i = len(s)//2
    n = s[i]

    l = s[:i]
    r = s[i+1:]
    num = 2*n+1
    den = abs(np.mean(l) - np.std(l,ddof=1))/2 + n + abs(np.mean(r) - np.std(r,ddof=1))/2 + 1
    a = max(1-num/den, 0)
    return a

def b_score(s):
    '''
    Given a series of numbers containing left-flank - n_i - right-flank
    Calculate B score from Birkedal et al
    Supplementary page 6
    '''
    # TO DO
    pass
    
def c_score(s):
    '''
    Given a series of numbers containing left-flank - n_i - right-flank
    calculate the C score from Birkedal et al
    Supplementary page 7
    '''
    # TO DO
    pass

def sliding_window(s, flank_size, score_function):
    '''
    Given a series of total read counts for each position in gene
    Step through 1 bp at a time
    For each window call a score function (for scores A, B, or C)
    Note that we add zeroes before and after to allow edge calculations
    Returns a list of scores
    '''
    # TO DO
    pass

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

# get a dict to convert long reference names to abbreviations
ref_name_to_abbrev = get_rrna_abbreviations(bamfile)

g5p, g3p = get_end_count(bamfile)

# store end counts in pandas dataframe, save as csv files for each gene
for long_gene, val5, val3 in zip(g5p.keys(), g5p.values(), g3p.values()):
    gene = ref_name_to_abbrev[long_gene]
    print(gene)
    df1 = pd.DataFrame.from_dict(val5,orient='index',columns=['5-prime'])
    df2 = pd.DataFrame.from_dict(val3,orient='index',columns=['3-prime'])
    max_index = max(max(df1.index), max(df2.index))
    new_index = pd.Index(np.arange(1,max_index+1))
    df1 = df1.reindex(new_index)
    df2 = df2.reindex(new_index)
    df1 = df1.fillna(0)
    df2 = df2.fillna(0)
    df = pd.concat([df1, df2],axis=1)
    df = df[['5-prime','3-prime']]
    df['total'] = df['5-prime'] + df['3-prime']
    for score_name, score_func in zip(('a_score','b_score','c_score'), (a_score, b_score, c_score)):
        df[score_name] = sliding_window(df['total'], flank_size, score_func)
    print(len(df))

    df.to_csv(os.path.join(output_folder, f'{gene}.csv'))
