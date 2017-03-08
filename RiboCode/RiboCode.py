#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

"""
Master script for running RiboCode program.
Usage:
	RiboCode.py -e <config_table.txt> -o <result.txt>
Details:
	Type "python RiboCode.py -h" to get the help information.
"""
import os
def main():
	"""
	Master function to call different functionalities of RiboCode
	"""
	import sys

	from parsing_opts import parsing_ribo
	args = parsing_ribo()
	from loadconfig import LoadConfig
	# read the config file
	configIn = LoadConfig(args.config_file)

	# check if annotation directory is exists.
	if not os.path.exists(args.annot_dir):
		sys.stderr.write("Error, the annotation directory not exists, pls run prepare_transcript.py first!\n")
		sys.exit()
	else:
		from prepare_transcripts import load_transcripts_pickle
		gene_dict, transcript_dict = load_transcripts_pickle(os.path.join(args.annot_dir,"transcripts.pickle"))

	#  reading the bam file
	import process_bam
	tpsites_sum, total_psites_number = process_bam.psites_count(configIn.configList,transcript_dict)
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
	detectORF.orf_detect(gene_dict,transcript_dict,args.annot_dir,tpsites_sum,total_psites_number,args.pval_cutoff,
	                     longest_orf,START_CODON,ALTERNATIVE_START_CODON_LIST,STOP_CODON_LIST,args.min_AA_length,
	                     args.output_file)
if __name__ == "__main__":
	main()
