#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'
from collections import namedtuple
import numpy as np
from scipy import stats
from scipy.stats import find_repeats,distributions,ttest_1samp
from statsmodels.stats import multitest


WilcoxonResult = namedtuple('WilcoxonResult', ('statistic', 'pvalue'))
def wilcoxon_greater(x, y, zero_method ="wilcox", correction = False):
	"""
	data if x is larger than y, single-sided.
	"""

	if np.allclose(x,y,equal_nan=True):
		return WilcoxonResult(np.nan, np.nan)
	"""
	shamelessly stolen from scipy
	"""
	if len(x) < 10 and not (np.allclose(x,x[0]) and np.allclose(y,y[0])):
		#sample size too small, using the ttest
		t_statistic,t_pvalue = ttest_1samp(x-y,popmean=0)
		if np.mean(x-y) >0:
			t_pvalue /= 2.0
		else:
			t_pvalue = 1 - t_pvalue / 2.0
		return WilcoxonResult(t_statistic,t_pvalue)

	if zero_method not in ["wilcox", "pratt", "zsplit"]:
		raise ValueError("Zero method should be either 'wilcox' "
						 "or 'pratt' or 'zsplit'")
	if y is None:
		d = np.asarray(x)
	else:
		x, y = map(np.asarray, (x, y))
		if len(x) != len(y):
			raise ValueError('Unequal N in wilcoxon.  Aborting.')
		d = x - y
		d[(d==0) & (x+y!=0)] = -1 #penalty for equal value


	if zero_method == "wilcox":
		# Keep all non-zero differences
		d = np.compress(np.not_equal(d, 0), d, axis=-1)

	count = len(d)
	# if count < 10:
	# 	warnings.warn("Warning: sample size too small for normal approximation.")

	r = stats.rankdata(abs(d))
	r_plus = np.sum((d > 0) * r, axis=0)
	r_minus = np.sum((d < 0) * r, axis=0)

	if zero_method == "zsplit":
		r_zero = np.sum((d == 0) * r, axis=0)
		r_plus += r_zero / 2.
		r_minus += r_zero / 2.

	T = min(r_plus, r_minus)
	mn = count * (count + 1.) * 0.25
	se = count * (count + 1.) * (2. * count + 1.)

	if zero_method == "pratt":
		r = r[d != 0]

	replist, repnum = find_repeats(r)
	if repnum.size != 0:
		# Correction for repeated elements.
		se -= 0.5 * (repnum * (repnum * repnum - 1)).sum()

	se = np.sqrt(se / 24)
	correction = 0.5 * int(bool(correction)) * np.sign(T - mn)
	z = (T - mn - correction) / se
	if r_plus > r_minus:
		prob = distributions.norm.sf(abs(z))
	else:
		prob = 1-distributions.norm.sf(abs(z))

	return WilcoxonResult(T, prob)

# def combine_pvals(pvalues, method="stouffer"):
# 	"""
# 	:param pvs
# 	:return: combined pvalue
# 	"""

# 	pvs = pvalues[~np.isnan(pvalues)]
# 	if pvs.size != 2:
# 		comb_pv = np.nan
# 	else:
# 		comb_pv = stats.combine_pvalues(pvalues,method=method)[1]

# 	return comb_pv


def cor_mat(x,y):
	"""
	create a k x k symmetric matrix that reflects the correlation structure among the tests.
	"""
	if x.sum() == 0 or y.sum() == 0:
		R = np.diag([1,1])
	else:
		R = np.corrcoef(x,y)
	return R

def meff(R, method ="galwey"):

	if method not in ["nyholt","liji","gao","galwey"]:
		raise ValueError("method must be nyholt, liji, gao, galwey")
	
	#get eigenvalues of 'R' matrix
	evs = np.linalg.eigvalsh(R)

	if method == "nyholt":
		k = evs.size
		m = 1 + (k-1) * (1 - np.var(evs,ddof=1) / k)
	elif method == "liji":
		evs_abs = np.abs(evs)
		tmp = evs_abs >= 1
		m = np.sum(tmp.astype(int) + (evs_abs - np.floor(evs_abs)))

	elif method == "gao":
		C = 0.995
		evs_rev_sort = -np.sort(-evs)
		m = np.where(np.cumsum(evs_rev_sort) / np.sum(evs) > C)[0][0]
	elif method == "galwey":
		if np.any(evs<0):
			evs[evs <0 ] = 0
		
		m = np.power(np.sum(np.sqrt(evs)), 2) / np.sum(evs)
		m = np.floor(m)

	return m

def combine_pvals(pvalues, method="stouffer", adjust="none", R="none"):
	"""
	see more details at: https://search.r-project.org/CRAN/refmans/poolr/html/stouffer.html
	"""

	if adjust not in ["nyholt","liji","gao","galwey","none"]:
		raise ValueError("method must be nyholt, liji, gao, galwey, or none")

	if adjust != "none" and isinstance(R,str) and R == "none":
		raise ValueError("R must be specified if adjust is set to something other than 'none'")

	
	pvs = pvalues[~np.isnan(pvalues)]
	if pvs.size != 2:
		return np.nan
	
	k = pvs.size
	
	if method == "stouffer":
		statistic = np.sum(stats.norm.isf(pvs)) / np.sqrt(pvs.size)			
		if adjust in ["nyholt", "liji", "gao", "galwey"]:
			m = meff(R,method=adjust)
			statistic = statistic * np.sqrt( m / k)
		elif adjust == "generalized":
			statistic = statistic * np.sqrt(k) / np.sqrt(np.sum(R))
		comb_pv = stats.norm.sf(statistic)
	else:
		raise ValueError("Currently, only stouffer is supported!")
	return comb_pv

def pvals_adjust(comb_pvs,method="bonferroni"):
	return multitest.multipletests(comb_pvs,method=method)[1]

	
