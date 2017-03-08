#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

import pysam
import numpy as np
import h5py
from itertools import izip
from collections import defaultdict
from prepare_transcripts import *

def write_psites(tpsites,psites_number,filename):
	with h5py.File(filename,"w") as fout:
		ds = h5py.special_dtype(vlen=str)
		dt = h5py.special_dtype(vlen=np.dtype("int32"))
		fout.create_dataset("transcript_ids",data=tpsites.keys(),dtype=ds)
		fout.create_dataset("p_sites",data=tpsites.values(),dtype=dt, compression="gzip")
		fout.create_dataset("psites_number",data=psites_number,dtype="int32")
	return None

def load_psites(filename):
	with h5py.File(filename,"r") as fin:
		k = fin["transcript_ids"][:]
		v = fin["p_sites"][:]
		tpsites = dict(izip(k,v))
		psites_number = fin["psites_number"].value
	return tpsites,psites_number

def read_bam(configData):
	"""
	read RPF bam file, calculate psites reads and reads number
	"""
	name = configData["samplename"]
	bamFile = configData["filepath"]
	stranded = configData["stranded"]
	psites_dict = configData["psites_dict"]
	transcript_dict = configData["transcript_dict"]
	if os.path.exists(name + "_psites.hd5"):
		tpsites,psites_number = load_psites(name + "_psites.hd5" )
	else:
		if not os.path.exists(bamFile):
			raise IOError("bam file does not exist: %s" % bamFile)
		sys.stdout.write("\tReading bam file: %s......\n" % bamFile)

		tpsites = defaultdict()
		total_psites = set()
		# init
		for tid in transcript_dict.iterkeys():
			tpsites[tid] = np.zeros(transcript_dict[tid].length,dtype="int32")

		tracks = pysam.AlignmentFile(bamFile)
		if tid not in tracks.references:
			sys.stderr.write("Error, the references in bam is different from transcriptome annotation, \n" +
			                 "you should input the transcriptome BAM file.")
			sys.exit()

		for r in tracks.fetch(until_eof=True):
			if r.is_unmapped:
				continue
			if r.is_reverse == stranded:
				continue
			tid = r.reference_name

			if r.query_length in psites_dict:
				t_psite = r.reference_start + psites_dict[r.query_length]
				if t_psite > transcript_dict[tid].length:
					continue
				else:
					total_psites.add(r.query_name)
					tpsites[tid][t_psite] += 1
		psites_number = len(total_psites)
		del total_psites
		write_psites(tpsites,psites_number, name + "_psites.hd5")
		sys.stdout.write("Finished reading bam file!\n")

	return tpsites,psites_number

def parallel_read(configList,thread_num):
	from multiprocessing import Pool
	mypool = Pool(thread_num)
	tpsites_list = mypool.map(read_bam,configList)
	mypool.close()
	mypool.join()
	transcript_dict = configList[0]["transcript_dict"]
	tpsites_sum = defaultdict()
	total_psites_number = 0
	# init
	for tid in transcript_dict.iterkeys():
		tpsites_sum[tid] = np.zeros(transcript_dict[tid].length,dtype="int32")
	for tpsites,psites_number in tpsites_list:
		total_psites_number += psites_number
		for tid in transcript_dict.iterkeys():
			tpsites_sum[tid] += tpsites[tid]

	return tpsites_sum,total_psites_number

def psites_count(configList,transcript_dict,thread_num=1):
	if len(configList) == 1:
		configList[0]["transcript_dict"] = transcript_dict
		tpsites_sum,total_psites_number = read_bam(configList[0])
	else:
		if thread_num == 1:
			tpsites_sum = defaultdict()
			total_psites_number = 0
			# init
			for tid in transcript_dict.iterkeys():
				tpsites_sum[tid] = np.zeros(transcript_dict[tid].length,dtype="int32")

			for configData in configList:
				configData["transcript_dict"] = transcript_dict
				tpsites,psites_number = read_bam(configData)
				total_psites_number += psites_number
				for tid in transcript_dict.iterkeys():
					tpsites_sum[tid] += tpsites[tid]
		else:
			for configData in configList:
				configData["transcript_dict"] = transcript_dict
			tpsites_sum,total_psites_number = parallel_read(configList,thread_num)

	return tpsites_sum,total_psites_number

if __name__ == "__main__":
	from loadconfig import LoadConfig
	config = LoadConfig("config.txt")
	psites_count(config.configList,transcript_dict,1)
