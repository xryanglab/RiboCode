#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

from Bio.Seq import translate
from sys import stderr
def translation(seq,table=1,cds=True):
	"""
	translate the DNA to protein sequence using the translation table
	table = 1, is the standard table, ref: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
	"""
	if len(seq) % 3 != 0:
		stderr.write("Warning: sequence is not divisible by 3")
		seq = seq[:-(len(seq) % 3)]
	return translate(seq,table,cds=cds)
