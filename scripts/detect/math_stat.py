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



def beta_binomial_significance(test_alt_count, test_total, bg_alt, bg_total):
        if test_total ==0:
                return np.nan
        if test_alt_count ==0:
                return 1
        else:
                return  max(1-math.exp(calculate_beta_binomial_logcdf(test_alt_count-1, test_total, bg_alt+1, bg_total-bg_alt+1)), 0)

#### THE FOLLOWING TWO FUNCTIONS WERE WRITTEN BY SOMEONE ELSE FOR RECPSEG
def calculate_log_n_choose_k(n, k):
        """
        This function computes the logarithm of the binomial coefficient.
        :param n: non-negative integer
        :param k: non-negative integer
        :return: log(n choose k)
        """
        return gammaln(n+1) - gammaln(n-k+1) - gammaln(k+1)

def calculate_beta_binomial_logpdf(k, n, alpha, beta):
        """
        This function computes the logarithm of the likelihood by the beta-binomial distribution. In the beta-binomial
        distribution, the likelihood of success is unknown and it follows a beta distribution.
        :param k: non-negative integer used in the binomial
        :param n: non-negative integer used in the binomial
        :param alpha: positive real used to model prior
        :param beta: positive real used to model prior
        :return:
        """
        return betaln(k+alpha, n-k+beta) + calculate_log_n_choose_k(n, k) - betaln(alpha, beta)

###### END BORROWED FROM RECAPSEG

def calculate_beta_binomial_logcdf(n1, N1, alpha, beta):
        ## Note that python intervals are left-closed and right-open, thus n1+1
        return logsumexp([(calculate_beta_binomial_logpdf(k, N1, alpha, beta)) for k in range(0,n1+1)])




# def f_combin_p(list_p):
#     # get combine pvalue
#     g_p = np.array(list_p)
#     pv_exon = stats.combine_pvalues(np.nan_to_num(g_p, nan = 1))[1]
#     return(pv_exon)
    
# def f_combin_p(list_p):
#     # get combine pvalue
#     g_p = np.array(list_p)
#     pv_exon = np.nanmin(np.nan_to_num(g_p, nan = 1))
#     return(pv_exon)


def f_combin_p(list_p):
    # get combined pvalue for dependent situations
    g_p = np.array(list_p)
    g_p = np.nan_to_num(g_p, nan = 1)
    pv_exon = stats.hmean(g_p)
    return(pv_exon)