#!/usr/bin/env python
from __future__ import division
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

"""
One-step script for running RiboCode program.
Usage:
	RiboCode_onestep.py -g <gtf.txt> -f <genome.fa> -r <ribo_toTranscript.bam> -o <result.txt>
Details:
	Type "RiboCode_onestep.py -h" to get the help information.
"""

from prepare_transcripts import *
def main():
	"""
	Master function to call different functionalities of RiboCode
	"""
	import sys
	import os

	from parsing_opts import parsing_ribo_onestep
	args = parsing_ribo_onestep()
	#1. prepare transcript annotation

	verboseprint("Preparing annotation files ...")
	if not os.path.exists("annot"):
		try:
			os.mkdir("annot")
		except OSError as e:
			raise e
	gene_dict, transcript_dict = processTranscripts(args.genomeFasta,args.gtfFile,"annot")
	verboseprint("The step of preparing transcript annotation has been completed!")

	#2. determine p-sites locations
	import metaplots
	from detectORF import extract_frame,test_frame
	min_readLength = 24
	max_readLength = 36
	frame0_percent = 0.65
	pvalue1_cutoff = pvalue2_cutoff = 0.001
	# filter transcripts
	filter_tids = metaplots.filter_transcript(gene_dict,transcript_dict)
	# read bam file
	distance_to_start_count,distance_to_stop_count,length_counter = metaplots.readTranscriptBam(
		args.rpf_mapping_file,filter_tids,transcript_dict,args.stranded,min_readLength,max_readLength)
	# predefine the psite
	pre_psite_dict = {}
	total_reads = sum(length_counter.values())
	if os.path.exists("config.txt"):
		sys.stderr.write("Warning, the config.txt exists, it will be overwritten !\n")
	fout = open("config.txt", "w")
	fout.write("#read_length\tproportion\tpredicted_psite\tnumber_of_codons_chosen\tf0_sum\tf1_sum\t_f2_sum\tf0_percent\tpvalue1\tpvalue2\tpvalue_combined\n")
	for l,d in distance_to_start_count.iteritems():
		if d.sum() < 10:
			continue
		pre_psite = d[:50].argmax()
		f0,f1,f2 = extract_frame(d[pre_psite:])
		if f0.sum() < 10:
			continue
		pv1,pv2,pv = test_frame(f0, f1, f2)
		f0_percent = f0.sum() / (f0.sum() + f1.sum() + f2.sum())
		if f0_percent < frame0_percent:
			continue
		if (pv1.pvalue < pvalue1_cutoff) and (pv2.pvalue < pvalue2_cutoff):
			read_percent = '{:.2%}'.format(length_counter[l] / total_reads)
			if length_counter[l] / total_reads >= 0.05:
				pre_psite_dict[l] = pre_psite
			num_of_codons = len(f0)
			fout.write("# " + "\t".join(map(str,[l,read_percent,-pre_psite+50,num_of_codons,f0.sum(),f1.sum(),f2.sum(),
			                                     '{:.2%}'.format(f0_percent),pv1.pvalue,pv2.pvalue,pv])) + "\n")
	#print the psite lines
	stranded = "yes" if args.stranded is True else "reverse"
	fout.write("\n")
	fout.write("# " + "\t".join(["SampleName","AlignmentFile","Stranded(yes/reverse)","P-siteReadLength","P-siteLocations"]) + "\n")
	pre_psite_len = map(str,sorted(pre_psite_dict.keys()))
	pre_psite_loc = map(str,[-pre_psite_dict[i]+50 for i in sorted(pre_psite_dict.keys())])
	sampleName = os.path.splitext(os.path.basename(args.rpf_mapping_file))[0]
	fout.write("\t".join(map(str,[sampleName,args.rpf_mapping_file,stranded,",".join(pre_psite_len),",".join(pre_psite_loc)])) + "\n")
	fout.close()
	metaplots.distancePlot(distance_to_start_count,distance_to_stop_count,pre_psite_dict,length_counter,"config.txt")

	#3 detectORF
	#  reading the bam file
	import process_bam
	from loadconfig import LoadConfig
	# read the config file
	configIn = LoadConfig("config.txt")

	tpsites_sum, total_psites_number = process_bam.psites_count(configIn.configList,transcript_dict,thread_num=1)
	import detectORF
	if args.longest_orf == "yes":
		longest_orf = True
	else:
		longest_orf = False

	if args.start_codon:
		START_CODON = args.start_codon.strip().split(",")
	if args.alternative_start_codons:
		ALTERNATIVE_START_CODON_LIST = args.alternative_start_codons.strip().split(",")
	else:
		ALTERNATIVE_START_CODON_LIST = None
	if args.stop_codon:
		STOP_CODON_LIST = args.stop_codon.strip().split(",")

	output_gtf = args.output_gtf
	output_bed = args.output_bed
	detectORF.main(gene_dict=gene_dict, transcript_dict=transcript_dict, annot_dir = "annot",
	               tpsites_sum=tpsites_sum, total_psites_number=total_psites_number,
	               pval_cutoff = args.pval_cutoff, only_longest_orf=longest_orf, START_CODON=START_CODON,
	               ALTERNATIVE_START_CODON_LIST=ALTERNATIVE_START_CODON_LIST, STOP_CODON_LIST=STOP_CODON_LIST,
	               MIN_AA_LENGTH=args.min_AA_length, outname=args.output_name,
	               output_gtf=output_gtf, output_bed=output_bed)

if __name__ == "__main__":
	verboseprint("Detecting ORFs ...")
	main()
	verboseprint("Finished !")
