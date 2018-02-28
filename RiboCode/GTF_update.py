#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from builtins import str, map, range

# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

from .prepare_transcripts import *
from collections import defaultdict
from sys import stderr
# perform various operations to generate the standard GTF file

def GroupGeneSubsets(gtf_file):
	"""
	Group gtf lines into subset that are the same gene
	Group is based on the ID attribute:

	:param gtf_data:
	:return:
	"""
	subsets = defaultdict(list)
	sourted_gset_keys = []
	sourted_gset_keys_set = set()
	with open(gtf_file) as fin:

		for line in fin:
			if line[0] =="#" or (not line.strip()):
				continue

			field_dict = parsing_line(line)
			field_dict["line"] = line
			gid=field_dict["attr"]["gene_id"]
			if "gene_name" not in field_dict["attr"]:
				field_dict["attr"]["gene_name"] = gid
			chrom = field_dict["chrom"]
			strand = field_dict["iv"].strand
			gid0 = chrom + ":"+gid+":"+strand
			subsets[gid0].append(field_dict)
			if gid0 not in sourted_gset_keys_set:
				sourted_gset_keys_set.add(gid0)
				sourted_gset_keys.append(gid0)
	#select duplicated gene id
	gene_id_list = [i.split(":")[1] for i in subsets.keys()]
	gene_id_uniq = set(gene_id_list)
	duplicated_gene_id = [i for i in gene_id_uniq if gene_id_list.count(i) >1]
	duplicated_gene_id = set(duplicated_gene_id)
	for i in subsets.keys():
		chrom,gid,strand = i.split(":")
		if gid in duplicated_gene_id:
			stderr.write("Warnning, gene_id is duplicated: %s , will be renamed as: %s\n" % (gid, i))
			for j in range(len(subsets[i])):
				subsets[i][j]["attr"]["gene_id"] = i
	return subsets,sourted_gset_keys

def TranscriptFeature(gene_set,gattr):
	tid_set = [i["attr"]["transcript_id"] for i in gene_set]
	tid_set=set(tid_set)
	chrom = gene_set[0]["chrom"]
	for tid in tid_set:
		t_subsets = [i for i in gene_set if i["attr"]["transcript_id"] == tid]
		ivs_set = [i["iv"] for i in t_subsets]

		strand = ivs_set[0].strand
		if strand == "+":
			ivs_set.sort(key=lambda x: x.start, reverse=False)
			tstart = ivs_set[0].start+1
			tend = ivs_set[-1].end
		else:
			ivs_set.sort(key=lambda x: x.start, reverse=True)
			tstart = ivs_set[-1].start + 1
			tend = ivs_set[0].end
		source = t_subsets[0]["source"]
		ttype = t_subsets[0]["attr"].get("transcript_type",None) or t_subsets[0]["attr"].get("transcript_biotype",None)

		tattr = gattr + ' transcript_id "%s";' % tid
		if ttype:
			tattr += ' transcript_type "%s";' % ttype
		out_list = [chrom,source,"transcript",tstart,tend,".",strand,".",tattr]
		print("\t".join(map(str,out_list)))

		for data in t_subsets:
			print(data["line"],end='')

def AddGeneFeature(subsets,sourted_gset_keys,add_transcript=True):
	for g in sourted_gset_keys:
		chrom,gid,strand =  g.split(":")
		ivs_set = [i["iv"] for i in subsets[g]]
		if strand == "+":
			ivs_set.sort(key=lambda x: x.start, reverse=False)
			gstart = ivs_set[0].start+1
			gend = ivs_set[-1].end
		else:
			ivs_set.sort(key=lambda x: x.start, reverse=True)
			gstart = ivs_set[-1].start + 1
			gend = ivs_set[0].end

		source = subsets[g][0]["source"]
		gname = subsets[g][0]["attr"].get("gene_name",subsets[g][0]["attr"]["gene_id"])
		gtype = subsets[g][0]["attr"].get("gene_type",None) or subsets[g][0]["attr"].get("gene_biotype",None)

		gattr = 'gene_id "%s"; gene_name "%s";' % (gid,gname)
		if gtype:
			gattr += ' gene_type "%s";' % gtype
		out_list = [chrom,source,"gene",gstart,gend,".",strand,".",gattr]
		print("\t".join(map(str,out_list)))
		if add_transcript:
			TranscriptFeature(subsets[g],gattr)
		else:
			for data in subsets[g]:
				print(data["line"],end='')

def main():
	from .parsing_opts import parsing_gtf_update
	args = parsing_gtf_update()
	gset,sourted_gset_keys = GroupGeneSubsets(args.gtfFile)
	#if the transcript feature exists.
	add_transcript = True
	for i in gset[sourted_gset_keys[1]]:
		if i["feature"] == "transcript":
			add_transcript = False
			break

	AddGeneFeature(gset,sourted_gset_keys,add_transcript)

if __name__ == "__main__":
	main()
