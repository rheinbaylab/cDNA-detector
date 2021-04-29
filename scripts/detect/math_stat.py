from scipy.special import betaln, gammaln, logsumexp
from scipy.stats import binom
import math
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
np.seterr(divide = 'ignore') 
# np.seterr(divide = 'warn') 



def beta_binomial_significance(nc,nt,Nc,Nt):
    # nc = 129
    # nt = 191
    # Nc = 46979
    # Nt = 2926481
    if nt == 0:
            return np.nan
    if nc == 0:
            return 1
    else:
        k = nc - 1
        n = nt
        a = Nc
        b = Nt - Nc + 1
        pvalue = 1 - stats.betabinom.cdf(k, n, a, b, loc=0)
        return pvalue


def f_combin_p(list_p):
    # get combined pvalue for dependent situations
    g_p = np.array(list_p)
    g_p = np.nan_to_num(g_p, nan = 1)
    pv_exon = stats.hmean(g_p)
    return(pv_exon)