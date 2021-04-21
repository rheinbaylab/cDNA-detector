from collections import Counter
import pysam
import pandas as pd
import numpy as np
import re
import os
import sys
import collections
import scipy
from scipy import stats
import statsmodels


def f_get_consensus_sequence(list_seq):
    ## get all possible consensus seq in each position
    consensus_seq = str();
    list_seq.sort(key = len)
    if len(list_seq)==0:
        consensus_seq = '';
    for i in range(0,len(list_seq[-1])):
        tmp = [x[0] for x in list_seq if x]
        tmp_count = [ tmp.count(n) for n in list('ACGT') ]
        tmp_sum = sum(tmp_count)+ 0.01
        tmp_freq = [float(x)/tmp_sum for x in tmp_count]
        list_seq = [x[1:] for x in list_seq if x]
        if max(tmp_freq)>=0.6:
            consensus_seq = consensus_seq +  ['ACGT'][0][tmp_freq.index(max(tmp_freq))]
        else:
            break
    return(consensus_seq)


def f_consensus_seq_count_exact(consense_seq, list_detect_seq):
    ## count consensus exact match sequences.
    seqlist_count = collections.Counter()
    tmp_cons_seqlist = [consense_seq for detect_seq in list_detect_seq if consense_seq.startswith(detect_seq) and detect_seq!=""]
    seqlist_count = collections.Counter(tmp_cons_seqlist)
    if tmp_cons_seqlist == []:
        seqlist_count = collections.Counter({consense_seq:0})
    return(seqlist_count)




def f_sequence_pair_compare_n(g_str,c_str,type):
    i=0
    if type == "start":
        g_str = g_str[::-1]
        c_str = c_str[::-1]
    for i in range(min(len(g_str),len(c_str))):
        if g_str[i] == c_str[i]:
            i = i+1
        else:
            break
    return(i)
