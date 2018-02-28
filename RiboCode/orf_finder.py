#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

import re

def searchCodon(seq,START_CODON,ALTERNATIVE_START_CODON_LIST,STOP_CODON_LIST):
	"""
	search all start codon and stop codons for each frame
	:param seq:
	:return: start_codon_dict, stop_codon_dict, the key is frame.
	"""

	seq = seq.upper()

	start_regx = re.compile('|'.join(START_CODON))
	stop_regx = re.compile('|'.join(STOP_CODON_LIST))

	start_idx = [m.start(0) for m in re.finditer(start_regx,seq)]
	stop_idx = [m.start(0) for m in re.finditer(stop_regx,seq)]
	if ALTERNATIVE_START_CODON_LIST:
		alt_start_regx = re.compile('|'.join(ALTERNATIVE_START_CODON_LIST))
		alt_start_idx = [m.start(0) for m in re.finditer(alt_start_regx,seq)]
	else:
		alt_start_idx = False
	return start_idx,stop_idx,alt_start_idx

def orf_find(start_idx,stop_idx,alt_start_idx,MIN_AA_LENGTH):
	"""
	extract ORFs
	"""
	commonstop_dict = {}
	for f in (0,1,2):
		inframe_stops = [x for x in stop_idx if x%3==f]
		inframe_starts = [x for x in start_idx if x%3==f]
		if alt_start_idx:
			inframe_alt_starts = [x for x in alt_start_idx if x%3==f]
		else:
			inframe_alt_starts = None

		last_stop = -1
		for j in inframe_stops:
			tmp_starts = []
			alt_flag = 1
			if inframe_starts:
				for i in inframe_starts:
					if i < j and i > last_stop and (j-i) >= MIN_AA_LENGTH*3:
						alt_flag = 0
						tmp_starts.append(i)
			if alt_flag and inframe_alt_starts:
				for i in inframe_alt_starts:
					if i < j and i > last_stop and (j-i) >= MIN_AA_LENGTH*3:
						tmp_starts.append(i)
			last_stop = j
			if tmp_starts:
				commonstop_dict[j] = tmp_starts
	return commonstop_dict

def orf_finder(seq,START_CODON,ALTERNATIVE_START_CODON_LIST,STOP_CODON_LIST,MIN_AA_LENGTH):
	start_idx,stop_idx,alt_start_idx = searchCodon(seq,START_CODON,ALTERNATIVE_START_CODON_LIST,STOP_CODON_LIST)
	orf_dict = orf_find(start_idx,stop_idx,alt_start_idx,MIN_AA_LENGTH)
	return orf_dict
