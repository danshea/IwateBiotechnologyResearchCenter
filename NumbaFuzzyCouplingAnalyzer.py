#!/usr/bin/env python

################################################################################
#
# Name: CouplingAnalyzer.py
# Date: 2019-10-04
# Author: Dan Shea
# Description: Check all SNPs for coupling of genotype in all founders
#
################################################################################

import pandas as pd
import numpy as np
import os.path
import os
import re
from collections import OrderedDict
import gzip
from numba import jit
import time
import sys

def main():
    samples = ['N01','N03','N04','N05','N06','N07','N08','N09','N10','N11','N12','N13','N14','N16','N17','N18','N19','N20','N21','N22',]
    filenames = ['N01_imputed_SNP_filtered_skip10000.vcf.gz','N03_imputed_SNP_filtered_skip10000.vcf.gz','N04_imputed_SNP_filtered_skip10000.vcf.gz',
                 'N05_imputed_SNP_filtered_skip10000.vcf.gz','N06_imputed_SNP_filtered_skip10000.vcf.gz','N07_imputed_SNP_filtered_skip10000.vcf.gz',
                 'N08_imputed_SNP_filtered_skip10000.vcf.gz','N09_imputed_SNP_filtered_skip10000.vcf.gz','N10_imputed_SNP_filtered_skip10000.vcf.gz',
                 'N11_imputed_SNP_filtered_skip10000.vcf.gz','N12_imputed_SNP_filtered_skip10000.vcf.gz','N13_imputed_SNP_filtered_skip10000.vcf.gz',
                 'N14_imputed_SNP_filtered_skip10000.vcf.gz','N16_imputed_SNP_filtered_skip10000.vcf.gz','N17_imputed_SNP_filtered_skip10000.vcf.gz',
                 'N18_imputed_SNP_filtered_skip10000.vcf.gz','N19_imputed_SNP_filtered_skip10000.vcf.gz','N20_imputed_SNP_filtered_skip10000.vcf.gz',
                 'N21_imputed_SNP_filtered_skip10000.vcf.gz','N22_imputed_SNP_filtered_skip10000.vcf.gz',]
    
    dfs = OrderedDict()
    for key, f in zip(samples, filenames):
        with gzip.open(f, 'rt') as fh:
            dfs[key] = pd.read_csv(fh, sep='\t', header=None, comment='#')
    
    for key in samples:
        dfs[key].drop(columns=[2,5,6,7,8], inplace=True)
    
    for key in samples:
        dfs[key].rename(columns={0: 'CHROM', 1: 'POS', 3:'REF', 4:'ALT', 9:'HITOMEBORE', 10:'FOUNDER'}, inplace=True)
    
    for key in samples:
        mapping = {k: 'RIL_{}'.format(v) for k, v in zip(dfs[key].columns[6:], range(0,len(dfs[key].columns[6:])))}
        dfs[key].rename(columns=mapping, inplace=True)
    
    def parse_gt(x):
        match = re.match('([01]/[01])', x)
        if match is not None:
            return match.group()
        else:
            return np.NaN
    
    for key in samples:
        tmp = dfs[key].iloc[:, 4:]
        for c in tmp.columns:
            tmp[c] = tmp[c].apply(parse_gt)
        dfs[key] = pd.concat([dfs[key].iloc[:, 0:4], tmp], axis=1)
    
    recoded = OrderedDict()
    for key in samples:
        tmp_df = list()
        for nt in dfs[key].itertuples(index=False, name=None):
            tmp_row = list()
            if nt[4] == '0/0' and nt[5] == '1/1':
                for ril in nt[4:]:
                    if ril == '0/0':
                        tmp_row.append(0)
                    elif ril == '1/1':
                        tmp_row.append(1)
                    elif (ril == '0/1') or (ril == '1/0'):
                        tmp_row.append(0.5)
                    else:
                        tmp_row.append(0.5)
            else:
                for ril in nt[4:]:
                    if ril == '0/0':
                        tmp_row.append(1)
                    elif ril == '1/1':
                        tmp_row.append(0)
                    elif (ril == '0/1') or (ril == '1/0'):
                        tmp_row.append(0.5)
                    else:
                        tmp_row.append(0.5)
            tmp_df.append(tmp_row)
        recoded[key] = pd.DataFrame(tmp_df, columns=dfs[key].columns[4:])
    
    output_dfs = OrderedDict()
    for key in samples:
        output_dfs[key] = pd.merge(dfs[key].iloc[:,0:4], recoded[key], how='inner', left_index=True, right_index=True)
    
    @jit(nopython=True)
    def comparator(a, b):
        tmp_sum = np.array([0.0], dtype=np.float64)
        tmp_sum = a + (-1 * b)
        tmp_min = np.min(tmp_sum)
        tmp_max = np.max(tmp_sum)
        if (tmp_min == -1.0) or (tmp_max == 1.0):
            return False
        else:
            return True

    results_df = OrderedDict()
    for key in samples:
        results_df[key] = list()
        now = time.localtime()
        print('{}.{}.{} {}:{}:{} Started processing {}'.format(now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec, key), flush=True)
        for i in range(0, output_dfs[key].shape[0]-1):
            for j in range(i+1, output_dfs[key].shape[0]):
                if output_dfs[key].iloc[i, 0] != output_dfs[key].iloc[j, 0]:
                    tmp_a = np.array(output_dfs[key].iloc[i, 6:], dtype=np.float64)
                    tmp_b = np.array(output_dfs[key].iloc[j, 6:], dtype=np.float64)
                    flag = comparator(tmp_a, tmp_b)
                    if flag:
                        tmp_row = list(output_dfs[key].iloc[i, :])
                        tmp_row.extend(list(output_dfs[key].iloc[j, :]))
                        results_df[key].append(tmp_row)
        results_df[key] = pd.DataFrame(results_df[key])
        now = time.localtime()
        print('{}.{}.{} {}:{}:{} Finished processing {}'.format(now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec, key), flush=True)
        results_df[key].to_csv('{}_NumbaFuzzyCouplingAnalyzer.tsv.gz'.format(key), sep='\t', header=False, index=False)
    
    sys.exit(0)

if __name__ == '__main__':
    main()

