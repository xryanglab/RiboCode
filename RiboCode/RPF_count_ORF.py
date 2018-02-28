#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

import HTSeq
import pysam
import sys

def readGTF(gtfFile):
	gtf = HTSeq.GFF_Reader(gtfFile)
	start_codon_sites = {}
	stop_codon_sites = {}
	counts = {}

	ORF_features = HTSeq.GenomicArrayOfSets("auto",stranded="yes")
	# for i,f in enumerate(gtf):
	for f in gtf:
		# i += 1
		# if i % 10000 == 0:
		# 	sys.stderr.write("%d GFF lines processed.\r" % i)
		orf_id = f.attr["orf_id"]
		if f.type == "exon":
			ORF_features[f.iv] += orf_id
			counts[orf_id] = 0
		elif f.type == "start_codon":
			start_codon_sites[orf_id] = f.iv.start_d
		elif f.type == "stop_codon":
			stop_codon_sites[orf_id] = f.iv.end_d
	return start_codon_sites,stop_codon_sites,ORF_features,counts


def invert_strand(iv):
	iv2 = iv.copy()
	if iv2.strand == "+":
		iv2.strand = "-"
	elif iv2.strand == "-":
		iv2.strand = "+"
	else:
		raise ValueError("Illegal strand")
	return iv2


def count_reads(start_codon_sites,stop_codon_sites,ORF_features,counts,map_file,stranded,min_quality,count_mode,
                first_exclude_codons,last_exclude_codons,min_read,max_read,exclude_min_ORF):

	lowqual = 0
	notaligned = 0
	nonunique = 0
	too_short = 0
	too_long = 0
	min_read_string = "__too_short(<%i)" % min_read
	max_read_string = "__too_long(<%i)" % max_read
	first_exclude_nt = first_exclude_codons * 3
	last_exclude_nt = last_exclude_codons * 3

	pysam_fh = pysam.AlignmentFile(map_file)
	is_bam = pysam_fh.is_bam
	pysam_fh.close()
	if is_bam:
		tracks = HTSeq.BAM_Reader(map_file)
	else:
		tracks = HTSeq.SAM_Reader(map_file)
	# for i,r in enumerate(tracks):
	for r in tracks:
		# if i % 100000 == 0:
		# 	sys.stderr.write("%d alignment record processed.\r" % i)
		if not r.aligned:
			notaligned += 1
			continue
		try:
			if r.optional_field("NH") >1:
				nonunique += 1
				continue
		except KeyError:
			pass
		if r.aQual < min_quality:
			lowqual += 1
			continue
		read_len = len(r.read.seq)
		if read_len < min_read:
			too_short += 1
			continue
		if read_len > max_read:
			too_long += 1
			continue
		if stranded == "yes":
			iv_seq = (co.ref_iv for co in r.cigar if co.type =="M" and co.size >0)
		else:
			iv_seq = (invert_strand(co.ref_iv) for co in r.cigar if co.type=="M" and co.size>0)

		try:
			if count_mode == "intersection-strict":
				fs = None
				for iv in iv_seq:
					for iv2,fs2 in ORF_features[iv].steps():
						if fs is None:
							fs = fs2.copy()
						else:
							fs = fs.intersection(fs2)
			elif count_mode == "union":
				fs = set()
				for iv in iv_seq:
					for iv2, fs2 in ORF_features[iv].steps():
						fs = fs.union(fs2)
			if fs is None or len(fs) == 0:
				continue
			elif len(fs) > 1:
				continue
			else:
				orf_id = list(fs)[0]
				if read_len < exclude_min_ORF:
					counts[orf_id] +=1
					continue
				try:
					if abs(start_codon_sites[orf_id] - r.iv.start_d) < first_exclude_nt:
						continue
					elif abs(r.iv.end_d - stop_codon_sites[orf_id]) < last_exclude_nt:
						continue
					else:
						counts[orf_id] += 1
				except:
					counts[orf_id] += 1
		except:
			sys.stderr.write("Error occurred when processing mapping file in line:%s\n" % r.get_sam_line())
	counts["__too_low_quality"] += lowqual
	counts["__not_aligned"] += notaligned
	counts[min_read_string] += too_short
	counts[max_read_string] += too_long
	counts["__alignment_not_unique"] += nonunique

	return counts

def main():
	from .parsing_opts import parsing_ORF_count
	from .prepare_transcripts import verboseprint
	args = parsing_ORF_count()
	verboseprint("Reading the GTF file ...")

	start_codon_sites,stop_codon_sites,ORF_features,counts = readGTF(args.gtf_file)
	counts["__too_low_quality"] = 0
	counts["__not_aligned"] = 0
	counts["__too_short(<%i)" % args.min_read] = 0
	counts["__too_long(<%i)" % args.max_read] = 0
	counts["__alignment_not_unique"] = 0
	verboseprint("Reading the mapping file ...")

	for f in args.rpf_mapping_file:
		counts = count_reads(start_codon_sites,stop_codon_sites,ORF_features,counts,f.strip(),args.stranded,
		                     args.min_quality,args.count_mode,args.first_exclude_codons,args.last_exclude_codons,
		                     args.min_read,args.max_read,args.exclude_min_ORF)

	if args.output_file == "-":
		fout = sys.stdout
	else:
		fout = open(args.output_file,"w")

	for i in sorted(counts.keys()):
		fout.write("%s\t%d\n" % (i,counts[i]))

	fout.close()
	verboseprint("Finished!")
if __name__ == "__main__":
	main()
