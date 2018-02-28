#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from builtins import map, zip, range,object
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'
"""
loading and preparing the data
"""
import sys
class LoadConfig(object):
	"""
	Read the ribosome profiling alignment information from config file,
	See the sample file in the data folder.
	"""

	def __init__(self,configFile):
		self.filename = configFile
		self.configList = []
		self._parsing()

	def _parsing(self):
		i = 0
		with open(self.filename) as fin:
			samplenames = []
			for line in fin:
				if line.startswith("#"):
					continue
				if not line.strip():
					continue
				samplename, bamfile, stranded, plen, psite = line.strip().split()
				if "-" in plen:
					p_s,p_e = plen.split("-")
					plen = list(range(int(p_s),int(p_e) + 1))
				else:
					plen = plen.strip().split(",")
				psite = psite.strip().split(",")
				plen = list(map(int,plen))
				psite = list(map(int,psite))
				if len(plen) != len(psite):
					sys.stderr.write("Error, pleas check you config file\n, \
					                  the number of read lengths and P-site offsets is different! \
	  				                  Use commas to separate read lengths and P-site offsets,\
					                  e.g. 28,29 11,12\n")
					sys.exit()
				if stranded not in ["yes", "reverse"]:
					sys.stderr.write("Error, pleas check your config file\n, \
					                  the stranded should be yes or reverse.\
					                  reverse means reversed strand interpretation\n")
					sys.exit()
				else:
					stranded = True if stranded == "yes" else False
				if samplename in samplenames:
					sys.stderr.wirte("Error, pls check you config file\n, bam file name is duplicated: %s.\n" % samplename)
					sys.exit()
				samplenames.append(samplename)
				self.configList.append(dict(zip(["samplename","filepath","stranded","psites_dict"],
				                                 [samplename,bamfile,stranded,dict(zip(plen,psite))])) )
				i += 1
		if i == 0:
			sys.stderr.write("Error, can not determine the P-site locations." +
			                 "Please check your input parameters and try again.")
			sys.exit()
