#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

from collections import Counter,defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.backends.backend_pdf import PdfPages
from prepare_transcripts import *
import os
import sys

def filter_transcript(gene_dict,transcript_dict):
	"""
	For each coding gene, only one principal transcript is selected for meta-analysis
	"""
	filter_tids = []
	for gobj in gene_dict.values():
		if gobj.gene_type != "protein_coding":
			continue
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
		else:
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
	if next(iter(filter_tids)) not in tracks.references:
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


def distancePlot(distance_to_start_count,distance_to_stop_count,length_counter,outname):
	length_set = set(distance_to_start_count.keys() + distance_to_stop_count.keys())
	total_reads = sum(length_counter.values())
	with PdfPages(outname + ".pdf") as pdf:
		x = np.arange(-50,51,dtype=int)
		colors = np.tile(["b","g","r"], 34)
		for l in sorted(length_set):
			#plt.figure(figsize=(5,3))
			perct = '{:.2%}'.format(float(length_counter[l]) / total_reads)
			fig,(ax1,ax2) = plt.subplots(nrows=2,ncols=1)
			y1 = distance_to_start_count[l]
			y2 = distance_to_stop_count[l]
			ax1.vlines(x,ymin=np.zeros(101),ymax=y1,colors=colors[:-1])
			ax1.tick_params(axis='x',which="both",top="off",direction='out')
			ax1.set_xticks([-40,-20,-12,0,20,40])
			ax1.set_xlim((-50,50))
			ax1.set_xlabel("Distance (nt)")
			ax1.set_ylabel("Alignments")

			ax1.set_title("({} nt reads,proportion:{})".format(l,perct) + "\n Distance 5'- start codons")

			ax2.vlines(x,ymin=np.zeros(101),ymax=y2,colors=colors[:-1])
			ax2.tick_params(axis='x',which="both",top="off",direction='out')
			ax2.set_xticks([-40,-20,-12,0,20,40])
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
		for k,v in distance_dict.iteritems():
			fout.write("%i:\t%s\n" % (k,"\t".join(map(str,v))))



def main():
	import parsing_opts
	args = parsing_opts.parsing_metaplots()
	# load gene dict and transcript dict
	gene_dict,transcript_dict = load_transcripts_pickle(os.path.join(args.annot_dir,"transcripts.pickle"))
	# filter transcripts
	filter_tids = filter_transcript(gene_dict,transcript_dict)
	# read bam file
	distance_to_start_count,distance_to_stop_count,length_counter = readTranscriptBam(
		args.ribo_bam,filter_tids,transcript_dict,args.stranded,args.minLength,args.maxLength)
	distancePlot(distance_to_start_count,distance_to_stop_count,length_counter,args.outname)
	#lengthDistribution(length_counter,args.outname)

if __name__ == "__main__":
	main()
