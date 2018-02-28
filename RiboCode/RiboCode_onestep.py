#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

"""
One-step script for running RiboCode program.
Usage:
	RiboCode_onestep.py -g <gtf.txt> -f <genome.fa> -r <ribo_toTranscript.bam> -o <result.txt>
Details:
	Type "RiboCode_onestep.py -h" to get the help information.
"""

from .prepare_transcripts import *
import argparse


def main():
	"""
	Master function to call different functionalities of RiboCode
	"""
	import sys
	import os

	from .parsing_opts import parsing_ribo_onestep
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
	from . import metaplots
	args_metaplot = argparse.ArgumentParser(description="Storing the arguments for metaplots analysis")
	args_metaplot.rpf_mapping_file = args.rpf_mapping_file
	args_metaplot.stranded = args.stranded
	args_metaplot.minLength = args.minLength
	args_metaplot.maxLength = args.maxLength
	args_metaplot.frame0_percent = args.frame0_percent
	args_metaplot.pvalue1_cutoff = args_metaplot.pvalue2_cutoff = 0.001
	args_metaplot.outname = "metaplots"
	metaplots.meta_analysis(gene_dict,transcript_dict,args_metaplot)

	#3 detectORF
	#  reading the bam file
	from . import process_bam
	from .loadconfig import LoadConfig
	# read the config file
	configIn = LoadConfig(args_metaplot.outname + "_pre_config.txt")
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
