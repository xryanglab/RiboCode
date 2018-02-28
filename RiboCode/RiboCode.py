#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

"""
Master script for running RiboCode program.
Usage:
	RiboCode.py -e <config_table.txt> -o <result.txt>
Details:
	Type "RiboCode.py -h" to get the help information.
"""
import os
from .prepare_transcripts import *
def main():
	"""
	Master function to call different functionalities of RiboCode
	"""
	import sys

	from .parsing_opts import parsing_ribo
	args = parsing_ribo()
	from .loadconfig import LoadConfig
	# read the config file
	configIn = LoadConfig(args.config_file)

	# check if annotation directory is exists.
	if not os.path.exists(args.annot_dir):
		sys.stderr.write("Error, the annotation directory not exists, pls run prepare_transcript.py first!\n")
		sys.exit()
	else:
		from .prepare_transcripts import load_transcripts_pickle
		gene_dict, transcript_dict = load_transcripts_pickle(os.path.join(args.annot_dir,"transcripts.pickle"))

	#  reading the bam file
	from . import process_bam
	tpsites_sum, total_psites_number = process_bam.psites_count(configIn.configList,transcript_dict,thread_num=1)
	from . import detectORF
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
	detectORF.main(gene_dict=gene_dict, transcript_dict=transcript_dict, annot_dir = args.annot_dir,
	               tpsites_sum=tpsites_sum, total_psites_number=total_psites_number,
	               pval_cutoff = args.pval_cutoff, only_longest_orf=longest_orf, START_CODON=START_CODON,
	               ALTERNATIVE_START_CODON_LIST=ALTERNATIVE_START_CODON_LIST, STOP_CODON_LIST=STOP_CODON_LIST,
	               MIN_AA_LENGTH=args.min_AA_length, outname=args.output_name,
	               output_gtf=output_gtf, output_bed=output_bed)

if __name__ == "__main__":
	verboseprint("Detecting ORFs ...")
	main()
	verboseprint("Finished !")
