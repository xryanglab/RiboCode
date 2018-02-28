#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from builtins import str, zip, map, range
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

from .test_func import wilcoxon_greater,combine_pvals
from .prepare_transcripts import *
from .translate_dna import translation
from .orf_finder import orf_finder
from collections import OrderedDict,defaultdict
from .pie_plot import pie
import numpy as np
import os
import sys

def extract_frame(orf_psite_array):
	"""
	extract frame0 , frame1, frame2 vector
	"""
	if orf_psite_array.size % 3 != 0:
		shiftn = orf_psite_array.size % 3
		orf_psite_array2 = orf_psite_array[:-shiftn]
		f0 = orf_psite_array2[0:orf_psite_array2.size:3]
		f1 = orf_psite_array2[1:orf_psite_array2.size:3]
		f2 = orf_psite_array2[2:orf_psite_array2.size:3]
	else:
		f0 = orf_psite_array[0:orf_psite_array.size:3]
		f1 = orf_psite_array[1:orf_psite_array.size:3]
		f2 = orf_psite_array[2:orf_psite_array.size:3]
	return f0,f1,f2

def test_frame(f0, f1, f2):
	"""
	data if f0>f1, f0>f2
	"""

	pv1 = wilcoxon_greater(f0, f1)
	pv2 = wilcoxon_greater(f0, f2)
	pv = combine_pvals(np.array([pv1.pvalue, pv2.pvalue]))
	return pv1,pv2,pv

def percentage_format(value,decimal=2):
	"""
	Format number as percentage
	"""
	formatStr = "{:." + str(decimal) + "%}"
	return formatStr.format(value)
	#return "{:.2%}".format(value)

def cal_coverage(inArray):
	nozeros = np.flatnonzero(inArray).size
	return nozeros / inArray.size

def cal_RPKM(input_counts,input_length,total_counts):
	"""
	calculate the RPKM value
	"""
	return ((10**9) * input_counts) / (total_counts * input_length)

def classfy_orf(orfiv, cdsiv):
	"""
	define the orf type according the orf position related to cds region
	"""
	if not cdsiv:
		orftype = "novel"
	else:
		orftype = "annotated"
		if orfiv.end == cdsiv.end:
			orftype = "annotated"
		elif orfiv.end < cdsiv.start:
			orftype = "uORF"
		elif orfiv.start > cdsiv.end:
			orftype = "dORF"
		elif orfiv.start < cdsiv.start and orfiv.end < cdsiv.end and orfiv.end > cdsiv.start:
			orftype = "Overlap_uORF"
		elif cdsiv.start < orfiv.start < cdsiv.end and orfiv.end > cdsiv.end:
			orftype = "Overlap_dORF"
		elif cdsiv.start < orfiv.start < orfiv.end < cdsiv.end and (cdsiv.end-orfiv.end) % 3 !=0:
			orftype = "internal"
	return orftype

def start_check(commonstop_dict, only_longest_orf, tpsites, pval_cutoff):
	"""
	automatically determine the start codon
	"""
	candicate_orf_list = []
	for stop,start_list in commonstop_dict.items():
		# if np.sum(tpsites[start_list[0]:stop:3]) < 10:
		# 	continue
		start_num = len(start_list)
		start = None
		if only_longest_orf or start_num == 1:
			start = start_list[0]
			candicate_orf_list.append(Interval(start,stop,"+"))
			continue
		for start_idx in range(start_num):
			cur_start = start_list[start_idx]
			if start_idx == start_num -1:
				next_start = stop
			else:
				next_start = start_list[start_idx + 1]

			downstream_region_psites = tpsites[cur_start:next_start]
			f0_downstream,f1_downstream,f2_downstream = extract_frame(downstream_region_psites)
			if f0_downstream.max() == 0:
				continue
			if np.flatnonzero(f0_downstream).size >= 10:
				pv1,pv2,pv = test_frame(f0_downstream,f1_downstream,f2_downstream)
				if pv < pval_cutoff:
					start = cur_start
					break
				else:
					continue
			else:
				d1_pos = d1_neg = 0
				d2_pos = d2_neg = 0
				for i in range(f0_downstream.size):
					if f0_downstream[i] == f1_downstream[i] == 0:
						pass
					elif f0_downstream[i] > f1_downstream[i]:
						d1_pos += 1
					else:
						d1_neg += 1

					if f0_downstream[i] == f2_downstream[i] == 0:
						pass
					elif f0_downstream[i] > f2_downstream[i]:
						d2_pos += 1
					else:
						d2_neg += 1
				if d1_pos + d2_pos >  d1_neg + d2_neg:
					start = cur_start
					break
				else:
					continue
		if start is not None:
			candicate_orf_list.append(Interval(start,stop,"+"))
	return candicate_orf_list

def ORF_record(only_ATG=True):
	keys = ["ORF_ID","ORF_type","transcript_id","transcript_type","gene_id","gene_name","gene_type","chrom","strand",
	        "ORF_length","ORF_tstart","ORF_tstop","ORF_gstart","ORF_gstop","annotated_tstart","annotated_tstop",
	        "annotated_gstart","annotated_gstop","Psites_sum_frame0","Psites_sum_frame1","Psites_sum_frame2",
	        "Psites_coverage_frame0","Psites_coverage_frame1","Psites_coverage_frame2","Psites_frame0_RPKM",
	        "pval_frame0_vs_frame1","pval_frame0_vs_frame2","pval_combined","AAseq","orf_iv"]
	if only_ATG is False:
		keys.insert(9,"start_codon")
	r = OrderedDict(zip(keys,[None]*len(keys)))
	return r

def write_result(orf_results,outname):

	header = "\t".join(list(orf_results[0].keys())[:-1])
	with open(outname + '.txt',"w") as fout:
		fout.write(header + "\n")
		for v in orf_results:
			fout.write("\t".join(map(str,list(v.values())[:-1])) + "\n")

def write_to_gtf(gene_dict, transcript_dict, orf_results, collapsed_orf_idx, outname):
	"""
	Here is a brief description of the GFF fields:
	1. seqname - The name of the sequence. Must be a chromosome or scaffold.
	2. source - The program that generated this feature.
	3. feature - The name of this type of feature. Some examples of standard feature types are "CDS" "start_codon" "stop_codon" and "exon"li>
	4. start - The starting position of the feature in the sequence. The first base is numbered 1.
	5. end - The ending position of the feature (inclusive).
	6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ":.":.
	7. strand - Valid entries include "+", "-", or "." (for don't know/don't care).
	8. frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be ".".
	9. group - All lines with the same group are linked together into a single item.
	"""

	with open(outname + ".gtf","w") as fout1, open(outname + "_collapsed.gtf","w") as fout2:
		for idx,orf in enumerate(orf_results):
			orf_id = orf["ORF_ID"]
			orf_type = orf["ORF_type"]
			source = "RiboCode"
			strand = orf["strand"]
			orf_gstart = orf["ORF_gstart"] - 1
			orf_gstop = orf["ORF_gstop"]
			if strand == "-":
				orf_gstart, orf_gstop = orf_gstop, orf_gstart
				orf_gstart -= 1
				orf_gstop += 1

			gobj = gene_dict[orf["gene_id"]]
			tobj = transcript_dict[orf["transcript_id"]]
			orf_iv = orf["orf_iv"]
			exon_ivs = transcript_iv_transform(tobj, orf_iv)
			start_codon_ivs = transcript_iv_transform(tobj,Interval(orf_iv.start,orf_iv.start+3,"+"))
			stop_codon_ivs = transcript_iv_transform(tobj,Interval(orf_iv.end,orf_iv.end + 3,"+"))

			flag = 0 # for collapsed file
			if idx in collapsed_orf_idx:
				flag = 1

			#ORF
			tmp = [gobj.chrom,source,"ORF",orf_gstart + 1,orf_gstop,".",strand,"."]
			fields = 'orf_id "%s"; orf_type "%s"; gene_id "%s"; gene_name "%s"; transcript_id "%s"' % \
			         (orf_id, orf_type, gobj.gene_id, gobj.gene_name, tobj.transcript_id)
			fout1.write("\t".join(map(str,tmp)) + "\t" + fields + "\n")
			if flag:
				fout2.write("\t".join(map(str,tmp)) + "\t" + fields + "\n")
			#start_codon
			for i in start_codon_ivs:
				tmp = [gobj.chrom,source,"start_codon",i.start+1,i.end,".",strand,"."]
				fields = 'orf_id "%s"; orf_type "%s"; gene_id "%s"; gene_name "%s"; transcript_id "%s"' % \
				         (orf_id, orf_type, gobj.gene_id, gobj.gene_name, tobj.transcript_id)
				fout1.write("\t".join(map(str,tmp)) + "\t" + fields + "\n")
				if flag:
					fout2.write("\t".join(map(str,tmp)) + "\t" + fields + "\n")
			#stop_codon
			for i in stop_codon_ivs:
				tmp = [gobj.chrom,source,"stop_codon",i.start+1,i.end,".",strand,"."]
				fields = 'orf_id "%s"; orf_type "%s"; gene_id "%s"; gene_name "%s"; transcript_id "%s"' % \
				         (orf_id, orf_type, gobj.gene_id, gobj.gene_name, tobj.transcript_id)
				fout1.write("\t".join(map(str,tmp)) + "\t" + fields + "\n")
				if flag:
					fout2.write("\t".join(map(str,tmp)) + "\t" + fields + "\n")
			#exons
			for j,i in enumerate(exon_ivs):
				tmp = [gobj.chrom,source,"exon",i.start + 1,i.end,".",strand,"."]
				fields = 'orf_id "%s"; orf_type "%s"; exon_number %i; gene_id "%s"; gene_name "%s"; transcript_id "%s"' % \
				         (orf_id, orf_type, j, gobj.gene_id, gobj.gene_name, tobj.transcript_id)
				fout1.write("\t".join(map(str,tmp)) + "\t" + fields + "\n")
				if flag:
					fout2.write("\t".join(map(str,tmp)) + "\t" + fields + "\n")

def write_to_bed(gene_dict, transcript_dict, orf_results, collapsed_orf_idx, outname):
	# The first three fields in each feature line are required:
	#
	# 1. chrom - name of the chromosome or scaffold. Any valid seq_region_name can be used, and chromosome names can be
	#    given with or without the 'chr' prefix.
	# 2. chromStart - Start position of the feature in standard chromosomal coordinates (i.e. first base is 0).
	# 3. chromEnd - End position of the feature in standard chromosomal coordinates

	# Nine additional fields are optional. Note that columns cannot be empty - lower-numbered fields must always be
	# populated if higher-numbered ones are used.
	#
	# 4. name - Label to be displayed under the feature, if turned on in "Configure this page".
	# 5. score - A score between 0 and 1000. See track lines, below, for ways to configure the display style of scored data.
	# 6. strand - defined as + (forward) or - (reverse).
	# 7. thickStart - coordinate at which to start drawing the feature as a solid rectangle
	# 8. thickEnd - coordinate at which to stop drawing the feature as a solid rectangle
	# 9. itemRgb - an RGB colour value (e.g. 0,0,255). Only used if there is a track line with the value of itemRgb set
	#    to "on" (case-insensitive).
	# 10. blockCount - the number of sub-elements (e.g. exons) within the feature
	# 11. blockSizes - the size of these sub-elements
	# 12. blockStarts - the start coordinate of each sub-element
	with open(outname + ".bed","w") as fout1, open(outname + "_collapsed.bed","w") as fout2:
		for idx,orf in enumerate(orf_results):
			orf_id = orf["ORF_ID"]
			orf_type = orf["ORF_type"]
			strand = orf["strand"]
			orf_gstart = orf["ORF_gstart"] - 1
			orf_gstop = orf["ORF_gstop"]
			if strand == "-":
				orf_gstart, orf_gstop = orf_gstop, orf_gstart
				orf_gstart -= 1
				orf_gstop += 1

			gobj = gene_dict[orf["gene_id"]]
			tobj = transcript_dict[orf["transcript_id"]]
			orf_iv = orf["orf_iv"]
			exon_ivs = transcript_iv_transform(tobj, orf_iv)

			flag = 0 # for collapsed file
			if idx in collapsed_orf_idx:
				flag = 1

			name = ";".join([orf_id,gobj.gene_name,tobj.transcript_id,orf_type])
			for i in exon_ivs:
				tmp = [gobj.chrom,i.start,i.end,name,".",strand]
				fout1.write("\t".join(map(str,tmp)) + "\n")
				if flag:
					fout2.write("\t".join(map(str,tmp)) + "\n")

def main(gene_dict, transcript_dict, annot_dir, tpsites_sum, total_psites_number, pval_cutoff, only_longest_orf,
         START_CODON, ALTERNATIVE_START_CODON_LIST, STOP_CODON_LIST, MIN_AA_LENGTH, outname, output_gtf, output_bed):

	PSITE_SUM_CUTOFF = F0_NONZEROS = 5
	transcript_seq = GenomeSeq(os.path.join(annot_dir,"transcripts_sequence.fa"))
	only_ccds = True
	if len(START_CODON) == 1 and (not ALTERNATIVE_START_CODON_LIST):
		only_ATG = True
	else:
		only_ATG = False

	orf_results = []
	transcript_orf_hash = defaultdict(list)

	tid_num = len(transcript_dict)
	counts = 0
	transcript_reads_count = {}
	for i,tid in enumerate(transcript_dict.keys()):
		if i % 1000 == 0:
			sys.stderr.write(percentage_format(i / tid_num, decimal=0) + " has finished! \r")

		tobj = transcript_dict[tid]
		tpsites = tpsites_sum[tid]
		gobj = gene_dict[tobj.gene_id]
		transcript_reads_count[tid] = tpsites.sum()
		if tpsites.sum() <PSITE_SUM_CUTOFF:
			continue

		tseq = transcript_seq.get_seq(tobj.transcript_id)
		commonstop_dict = orf_finder(tseq,START_CODON,ALTERNATIVE_START_CODON_LIST,STOP_CODON_LIST,MIN_AA_LENGTH)
		candicate_orfivs = start_check(commonstop_dict, only_longest_orf, tpsites, pval_cutoff)
		for orf_iv in candicate_orfivs:
			if orf_iv.length < MIN_AA_LENGTH * 3:
				continue
			orf_psites_array = tpsites[orf_iv.start:orf_iv.end]
			f0,f1,f2 = extract_frame(orf_psites_array)
			if f0.sum() < PSITE_SUM_CUTOFF:
				continue
			if np.flatnonzero(f0).size < F0_NONZEROS:
				continue
			pv1,pv2,pv = test_frame(f0,f1,f2)
			if pv <= pval_cutoff:
				orf_dict = ORF_record(only_ATG)
				orf_dict["orf_iv"] = orf_iv
				orf_dict["transcript_id"] = tobj.transcript_id
				orf_dict["transcript_type"] = tobj.transcript_type
				orf_dict["gene_id"] = tobj.gene_id
				orf_dict["gene_name"] = tobj.gene_name
				orf_dict["gene_type"] = gobj.gene_type
				orf_dict["chrom"] = tobj.chrom
				orf_dict["strand"] = tobj.genomic_iv.strand
				orf_dict["ORF_tstart"] = orf_iv.start + 1 # convert to 1-base
				orf_dict["ORF_tstop"] = orf_iv.end + 3    # include stop-codon
				orf_gstart = transcript_pos_transform(tobj, orf_iv.start) + 1 # convert to 1-base
				orf_gstop = transcript_pos_transform(tobj, orf_iv.end + 2) + 1 # include stop-codon
				orf_dict["ORF_gstart"] = orf_gstart
				orf_dict["ORF_gstop"] = orf_gstop
				if tobj.startcodon is not None:
					annotated_tstart = tobj.startcodon.start + 1
					annotated_gstart = transcript_pos_transform(tobj, tobj.startcodon.start) + 1
				else:
					annotated_tstart = None
					annotated_gstart = None
				if tobj.stopcodon is not None:
					annotated_tstop = tobj.stopcodon.end
					annotated_gstop = transcript_pos_transform(tobj, tobj.stopcodon.end - 1) + 1
				else:
					annotated_tstop = None
					annotated_gstop = None
				orf_dict["annotated_tstart"] = annotated_tstart
				orf_dict["annotated_tstop"] = annotated_tstop
				orf_dict["annotated_gstart"] = annotated_gstart
				orf_dict["annotated_gstop"] = annotated_gstop
				orf_dict["ORF_ID"] = "%s_%s_%s_%i" % (tobj.gene_id,orf_gstart,orf_gstop,orf_iv.length/3)
				orf_dict["ORF_length"] = orf_iv.length
				orf_dict["Psites_sum_frame0"] = f0.sum()
				orf_dict["Psites_sum_frame1"] = f1.sum()
				orf_dict["Psites_sum_frame2"] = f2.sum()
				orf_dict["Psites_coverage_frame0"] = percentage_format(cal_coverage(f0))
				orf_dict["Psites_coverage_frame1"] = percentage_format(cal_coverage(f1))
				orf_dict["Psites_coverage_frame2"] = percentage_format(cal_coverage(f2))
				orf_dict["pval_frame0_vs_frame1"] = pv1.pvalue
				orf_dict["pval_frame0_vs_frame2"] = pv2.pvalue
				orf_dict["pval_combined"] = pv
				orf_dict["Psites_frame0_RPKM"] = cal_RPKM(orf_dict["Psites_sum_frame0"],orf_iv.length,total_psites_number)
				orf_seq = transcript_seq.get_seq(tobj.transcript_id,orf_iv.start,orf_iv.end,"+")
				orf_aaseq = translation(orf_seq,cds=False)
				orf_dict["AAseq"] = orf_aaseq
				orf_dict["ORF_type"] = classfy_orf(orf_iv,tobj.cds)
				if not only_ATG:
					orf_dict["start_codon"] = orf_seq[:3]
				orf_results.append(orf_dict)
				transcript_orf_hash[tid].append(counts)
				counts += 1

	#combine ORFs from transcripts
	#For each stop codon of a gene, the most upstream start is left
	sys.stdout.write("\n\tWriting the results to file .....\n")
	write_result(orf_results,outname)
	fout = open(outname + "_collapsed.txt","w")
	header = "\t".join(list(orf_results[0].keys())[:-1])
	fout.write(header + "\n")
	ORFs_category_dict = OrderedDict()
	for k in ["annotated","uORF","dORF","Overlapped","novel_PCGs","novel_NonPCGs"]:
		ORFs_category_dict[k] = set()

	collapsed_orf_idx = set()
	for gobj in gene_dict.values():
		orf_transcript_dict = {} #key is orf_gstop, value is tid,orf_gstart,orf_f0_sum,orf index
		ccds_tids = []

		if only_ccds:
			# only CCDS transcript or all transcript.
			for tid in gobj.transcripts:
				if "CCDS" in transcript_dict[tid].attr.get("tag",[]):
					ccds_tids.append(tid)
			if ccds_tids:
				tids = ccds_tids
			else:
				tids = gobj.transcripts
		else:
			tids = gobj.transcripts

		for tid in tids:
			for i in transcript_orf_hash[tid]:
				orf = orf_results[i]
				orf_gstop = orf["ORF_gstop"]
				orf_gstart = orf["ORF_gstart"]

				if orf_gstop in orf_transcript_dict:
					if gobj.genomic_iv.strand == "+" and orf_gstart < orf_transcript_dict[orf_gstop][1]:
						orf_transcript_dict[orf_gstop] = (tid,orf_gstart,orf["Psites_sum_frame0"],i)
					elif gobj.genomic_iv.strand == "-" and orf_gstart > orf_transcript_dict[orf_gstop][1]:
						orf_transcript_dict[orf_gstop]= (tid,orf_gstart,orf["Psites_sum_frame0"],i)
					elif (orf_gstart == orf_transcript_dict[orf_gstop][1]) and \
							(orf["Psites_sum_frame0"] > orf_transcript_dict[orf_gstop][2]):
						orf_transcript_dict[orf_gstop]= (tid,orf_gstart,orf["Psites_sum_frame0"],i)
					else:
						continue
				else:
					orf_transcript_dict[orf_gstop] = (tid,orf_gstart,orf["Psites_sum_frame0"],i)

		for v in orf_transcript_dict.values():
			_,_,_,orf_idx = v
			orf = orf_results[orf_idx]
			collapsed_orf_idx.add(orf_idx)
			if orf["ORF_type"] == "annotated":
				ORFs_category_dict["annotated"].add(gobj.gene_id)
			elif orf["ORF_type"] == "uORF":
				ORFs_category_dict["uORF"].add(gobj.gene_id)
			elif orf["ORF_type"] == "dORF":
				ORFs_category_dict["dORF"].add(gobj.gene_id)
			elif orf["ORF_type"] in ["Overlap_uORF","Overlap_dORF","internal"]:
				ORFs_category_dict["Overlapped"].add(gobj.gene_id)
			elif orf["ORF_type"] == "novel" and gobj.gene_type == "protein_coding":
				ORFs_category_dict["novel_PCGs"].add(gobj.gene_id)
			else:
				ORFs_category_dict["novel_NonPCGs"].add(gobj.gene_id)

			fout.write("\t".join(map(str,list(orf.values())[:-1])) + "\n")

	fout.close()
	ORFs_category_dict2 = OrderedDict()
	for k in ORFs_category_dict.keys():
		num = len(ORFs_category_dict[k])
		if num >0: ORFs_category_dict2[k] = num
	try:
		pie(ORFs_category_dict2,outname)
	except:
		pass
	if output_gtf:
		write_to_gtf(gene_dict, transcript_dict, orf_results, collapsed_orf_idx, outname)
	if output_bed:
		write_to_bed(gene_dict, transcript_dict, orf_results, collapsed_orf_idx, outname)
	verboseprint("Finished!")
	return None
