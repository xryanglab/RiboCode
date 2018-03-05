#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from builtins import map,range
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

from collections import Counter,defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.backends.backend_pdf import PdfPages
from .prepare_transcripts import *
from .detectORF import extract_frame,test_frame
from .parsing_opts import parsing_metaplots
import os
import sys

def filter_transcript(gene_dict,transcript_dict):
	"""
	For each coding gene, only one principal transcript is selected for meta-analysis
	"""
	filter_tids = []
	for gobj in gene_dict.values():
		if gobj.chrom in ("chrM","chrMT","M","MT"):
			continue
		select_tid = None
		if gobj.principal_transcripts:
			if len(gobj.principal_transcripts) == 1:
				tid = gobj.principal_transcripts[0]
				if transcript_dict[tid].startcodon and transcript_dict[tid].stopcodon:
					select_tid = tid
			else:
				for tid in gobj.principal_transcripts:
					if transcript_dict[tid].startcodon and transcript_dict[tid].stopcodon:
						if select_tid:
							if transcript_dict[tid].length > transcript_dict[select_tid].length:
								select_tid = tid
						else:
							select_tid = tid
		elif gobj.transcripts:
			# if only transcript level1
			levels = set([transcript_dict[tid].attr.get("level",None) for tid in gobj.transcripts])
			if len(levels) > 1:
				level = list(sorted(levels))[0]
			else:
				level = list(levels)[0]
			for tid in gobj.transcripts:
				if transcript_dict[tid].attr.get("level",None):
					if transcript_dict[tid].attr["level"] == level:
						if transcript_dict[tid].startcodon and transcript_dict[tid].stopcodon:
							if select_tid:
								if transcript_dict[tid].length > transcript_dict[select_tid].length:
									select_tid = tid
							else:
								select_tid = tid
				else:
					if transcript_dict[tid].startcodon and transcript_dict[tid].stopcodon:
						if select_tid:
							if transcript_dict[tid].length > transcript_dict[select_tid].length:
								select_tid = tid
						else:
							select_tid = tid
		if select_tid:
			filter_tids.append(select_tid)

	return set(filter_tids)

def readTranscriptBam(bamFile,filter_tids,transcript_dict,stranded,minLength,maxLength):

	tracks = pysam.AlignmentFile(bamFile)
	if tracks.references[0] not in transcript_dict:
		sys.stderr.write("Error, the references in bam are different from transcriptome annotation, \n" +
		                 "you should input the transcriptome BAM file.")
		sys.exit()

	distance_to_start_count = defaultdict(lambda: np.zeros(101,dtype=int))
	distance_to_stop_count = defaultdict(lambda: np.zeros(101,dtype=int))
	length_counter = Counter()
	for r in tracks:
		if  r.is_unmapped:
			continue
		if r.is_reverse == stranded:
			continue

		tid = r.reference_name
		if tid not in filter_tids:
			continue

		length_counter[r.query_length] += 1
		if minLength <= r.query_length <= maxLength:
			distance_to_start = r.reference_start - transcript_dict[tid].startcodon.start
			distance_to_stop = r.reference_start - transcript_dict[tid].stopcodon.end
			if abs(distance_to_start) <=50:
				distance_to_start_count[r.query_length][50 + distance_to_start] += 1
			if abs(distance_to_stop) <=50:
				distance_to_stop_count[r.query_length][50 + distance_to_stop] += 1

	return (distance_to_start_count,distance_to_stop_count,length_counter)


def distancePlot(distance_to_start_count,distance_to_stop_count,pre_psite_dict,length_counter,outname):
	length_set = set(list(distance_to_start_count.keys()) + list(distance_to_stop_count.keys()))
	total_reads = sum(length_counter.values())
	with PdfPages(outname + ".pdf") as pdf:
		x = np.arange(-50,51,dtype=int)
		colors = np.tile(["b","g","r"], 34)
		for l in sorted(length_set):
			#plt.figure(figsize=(5,3))
			if l not in pre_psite_dict:
				xticks = [-40,-20,0,20,40]
			else:
				xticks = sorted([-40,-20,0,20,40] + [pre_psite_dict[l] -50])
			perct = '{:.2%}'.format(length_counter[l] / total_reads)
			fig,(ax1,ax2) = plt.subplots(nrows=2,ncols=1)
			y1 = distance_to_start_count[l]
			y2 = distance_to_stop_count[l]
			ax1.vlines(x,ymin=np.zeros(101),ymax=y1,colors=colors[:-1])
			ax1.tick_params(axis='x',which="both",top="off",direction='out')
			ax1.set_xticks(xticks)
			ax1.set_xlim((-50,50))
			ax1.set_xlabel("Distance (nt)")
			ax1.set_ylabel("Alignments")

			ax1.set_title("({} nt reads,proportion:{})".format(l,perct) + "\n Distance 5'- start codons")

			ax2.vlines(x,ymin=np.zeros(101),ymax=y2,colors=colors[:-1])
			ax2.tick_params(axis='x',which="both",top="off",direction='out')
			ax2.set_xticks(xticks)
			ax2.set_xlim((-50,50))
			ax2.set_xlabel("Distance (nt)")
			ax2.set_ylabel("Alignments")
			ax2.set_title("Distance 5'- stop codons")

			fig.tight_layout()
			pdf.savefig(fig)
			plt.close()

	return None

def lengthDistribution(length_counter,outname):
	w,h = plt.figaspect(0.4)
	plt.figure(figsize=(w,h))
	x = sorted(length_counter.keys())
	y = [length_counter[i] for i in x]
	plt.bar(x,y,width=0.95,edgecolor="white",align="center",color="#297FFF")
	plt.savefig(outname + "_readlength_distribution.pdf")
	plt.close()

def _write_to_file(distance_dict,filename):
	with open(filename,"w") as fout:
		for k,v in distance_dict.items():
			fout.write("%i:\t%s\n" % (k,"\t".join(map(str,v))))

def _predict_psite(rlength,rdensity,frame0_percent_cutoff,pv1_cutoff,pv2_cutoff):
	predict_psite = None
	others = [] #f0sum,f1sum,f2sum,f0_percent,pv1,pv2,pv
	tmp_f0sum = 0
	for i in range(4,min(rlength,20)):
		idx = 50 - i
		f0,f1,f2 = extract_frame(rdensity[idx:idx+51])
		if f0.sum() < 10:
			continue
		f0_percent = f0.sum() / (f0.sum() + f1.sum() + f2.sum())
		pv1,pv2,pv = test_frame(f0, f1, f2)
		if f0_percent < frame0_percent_cutoff:
			continue
		if (pv1.pvalue < pv1_cutoff) and (pv2.pvalue < pv2_cutoff):
			if f0.sum() > tmp_f0sum:
				predict_psite = idx
				tmp_f0sum = f0.sum()
				others=[f0.sum(),f1.sum(),f2.sum(),f0_percent,pv1.pvalue,pv2.pvalue,pv]
	return predict_psite,others


def meta_analysis(gene_dict,transcript_dict,args):
	# filter transcripts
	filter_tids = filter_transcript(gene_dict,transcript_dict)
	if len(filter_tids) == 0:
		sys.stderr.write("Oops, no start codons and stop codons are found in annotation files.\n" +
		                 "If no any start codons and stop codons are annotated in GTF file \n" +
		                 ",please skip this step and create a config file to specify the P-site based on other evidence." )
		sys.exit()
	# read bam file
	distance_to_start_count,distance_to_stop_count,length_counter = readTranscriptBam(
		args.rpf_mapping_file,filter_tids,transcript_dict,args.stranded,args.minLength,args.maxLength)
	# predefine the psite
	pre_psite_dict = {}
	total_reads = sum(length_counter.values())
	fout = open(args.outname + "_pre_config.txt", "w")
	fout.write("#read_length\tproportion(per total mapped reads)\tpredicted_psite\tf0_sum\tf1_sum\tf2_sum\tf0_percent\tpvalue1\tpvalue2\tpvalue_combined\n")
	for l,d in distance_to_start_count.items():
		if d.sum() < 10:
			continue
		mask_max_psite = d[:50].argmax()
		predict_psite,others = _predict_psite(l,d,args.frame0_percent,args.pvalue1_cutoff,args.pvalue2_cutoff)
		if predict_psite:
			if mask_max_psite != predict_psite:
				sys.stderr.write("Warning:The predicted P-site location(%i) for length %i is not the highest peak(%i),\
				                 please confirm according metagene plots.\n" % (50-predict_psite,l,50-mask_max_psite))
			read_percent = '{:.2%}'.format(length_counter[l] / total_reads)
			pre_psite_dict[l] = predict_psite
			f0sum,f1sum,f2sum,f0_percent,pv1,pv2,pv = others
			fout.write("# " + "\t".join(map(str,[l,read_percent,50-predict_psite,f0sum,f1sum,f2sum,
			                                     '{:.2%}'.format(f0_percent),pv1,pv2,pv])) + "\n")

	#print the psite lines
	fout.write("\n")
	if pre_psite_dict:
		fout.write("# " + "\t".join(["SampleName","AlignmentFile","Stranded(yes/reverse)","P-siteReadLength","P-siteLocations"]) + "\n")
		stranded = "yes" if args.stranded is True else "reverse"
		pre_psite_len = list(map(str,sorted(pre_psite_dict.keys())))
		pre_psite_loc = list(map(str,[-pre_psite_dict[i]+50 for i in sorted(pre_psite_dict.keys())]))
		sampleName = os.path.splitext(os.path.basename(args.rpf_mapping_file))[0]
		fout.write("\t".join(map(str,[sampleName,args.rpf_mapping_file,stranded,",".join(pre_psite_len),",".join(pre_psite_loc)])) + "\n")
		fout.close()
		distancePlot(distance_to_start_count,distance_to_stop_count,pre_psite_dict,length_counter,args.outname)
		#lengthDistribution(length_counter,args.outname)
	else:
		distancePlot(distance_to_start_count,distance_to_stop_count,pre_psite_dict,length_counter,args.outname)
		sys.stderr.write("No obviously periodicity were detected from alignment reads in annotated start codons,\n" +
		                 "it could be due to poor quality sequencing.\n" +
		                 "Please check the metagene plots and try again by lowering the value of frame0_percent")


def main():
	verboseprint("Create metaplot file and predict the P-site locations ...")
	# load gene dict and transcript dict
	args = parsing_metaplots()
	gene_dict,transcript_dict = load_transcripts_pickle(os.path.join(args.annot_dir,"transcripts.pickle"))
	meta_analysis(gene_dict,transcript_dict,args)
	verboseprint("Complete prediction of the P-site locations")


if __name__ == "__main__":
	main()
