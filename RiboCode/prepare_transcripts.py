#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from builtins import range,object
# -*- coding:UTF-8 -*-
"""
preparing the transcripts annotation files.
"""
import sys
if sys.version_info.major ==2:
	import cPickle as pickle
else:
	import pickle
	from sys import intern
import os
from pyfasta import Fasta
from time import strftime

class ParsingError(Exception):
	pass

class Interval(object):
	""" basic interval """
	def __init__(self,start,end,strand):
		if end < start:
			raise ValueError('Error: end smaller than start !')

		self.start = start
		self.end = end
		self.strand = strand
		if strand == "+":
			self.start_d = start
			self.end_d = end
		elif strand == "-":
			self.start_d = end - 1
			self.end_d = start - 1
		else:
			raise ValueError('Error: strand is neither "+" nor "-" !')
		self.length = end - start

	def contain(self,start,end):
		if end < start:
			raise ValueError('Error: end smaller than start !')
		if self.start <= start and self.end >= end:
			return True
		else:
			return False

	def contain_iv(self,iv):
		return self.contain(iv.start,iv.end)

	def is_overlapped_with(self, iv):
		if self.start <= iv.start <= self.end or self.start <= iv.end <= self.end:
			return True
		else:
			return False

	def __str__(self):
		return "iv:%i-%i" % (self.start,self.end)

	__repr__ = __str__

def Interval_from_directional(start_d,end_d,strand="-"):
	if strand == "-":
		return Interval(end_d+1,start_d+1,"-")
	else:
		return Interval(start_d,end_d,"+")

def parsing_attr(attr_string):
	attr_dict = {}
	for i in attr_string.split(";"):
		if i:
			k,v = i.strip().split(" ",1)
			k = intern(k.strip())
			v = intern(v.strip('"'))
			if k == "tag":
				attr_dict.setdefault(k,[]).append(v)
			else:
				attr_dict[k] = v
	return attr_dict

def parsing_line(line):
	"""
        # GTF columns:
        # 1) chrom: str ("1", "X", "chrX", etc...)
        # 2) source : str (not used)
        # 3) feature : str ("gene", "transcript", &c)
        # 4) start : int
        # 5) end : int
        # 6) score : float or "." (not used)
        # 7) strand : "+", "-", or "."
        # 8) frame : 0, 1, 2 or "." (not used)
        # 9) attribute : key-value pairs separated by semicolons
	"""
	fields = line.strip().split("\t")
	chrom,source,feature,start,end,score,strand,frame,attr = fields
	if len(fields) != 9:
		raise ParsingError("Wrong number of fields %d (expected 9)" % len(fields))
	#convert to 0-based.
	iv = Interval(int(start) - 1, int(end), strand)
	field_dict = {"chrom": intern(chrom),"source":source,"feature": intern(feature),"iv":iv,"attr":parsing_attr(attr)}
	return field_dict

class Gene(object):
	""" Gene object """

	def __init__(self,field_dict):
		self.chrom = field_dict["chrom"]
		self.genomic_iv = field_dict["iv"]
		self.attr = field_dict["attr"]
		self.gene_id = self.attr["gene_id"]
		self.gene_name = self.attr.get("gene_name",self.attr.get("gene_id",None))
		self.gene_type = self.attr.get("gene_type",self.attr.get("gene_biotype",None))
		self.transcripts = []
		self.principal_transcripts = []

	def __str__(self):
		return "GeneID:%s;GeneName:%s;GeneType:%s" % (self.gene_id,self.gene_name,self.gene_type)

	__repr__ = __str__

class Transcript(object):
	""" Transcript object """

	def __init__(self,field_dict):
		self.chrom = field_dict["chrom"]
		self.genomic_iv = field_dict["iv"]
		self.attr = field_dict["attr"]
		self.gene_id = self.attr["gene_id"]
		self.gene_name = self.attr.get("gene_name",self.attr.get("gene_id",None))
		self.transcript_type = self.attr.get("transcript_type",self.attr.get("transcript_biotype",None))
		self.transcript_id = self.attr["transcript_id"]
		self.genomic_exons = []
		self.genomic_cds = []
		self.genomic_startcodon = []
		self.genomic_stopcodon = []
		self.length = 0

		self.cds = None
		self.startcodon = None
		self.stopcodon = None

	def add_feature(self,field_dict):
		feature = field_dict["feature"]
		iv = field_dict["iv"]
		if feature == "exon":
			self.genomic_exons.append(iv)
			self.length += iv.length
		elif feature == "CDS":
			self.genomic_cds.append(iv)
		elif feature == "start_codon":
			self.genomic_startcodon.append(iv)
		elif feature == "stop_codon":
			self.genomic_stopcodon.append(iv)
		else:
			raise ValueError("Error: the feature is not recognized: %s" % feature)

	def __str__(self):
		return "TranscriptID:%s;TranscriptType:%s" % (self.transcript_id,self.transcript_type)

	__repr__ = __str__

def readGTF(filename):
	if not os.path.exists(filename):
		raise IOError("\tGTF file does not exist: %s" % filename)
	sys.stdout.write("\tReading the GTF file: %s .......\n" % filename)
	gene_dict = {}
	transcript_dict = {}
	with open(filename) as fin:
		for i,line in enumerate(fin):
			if line[0] == "#" or (not line.strip()):
				continue
			field_dict = parsing_line(line)
			if field_dict["feature"] == "gene":
				gobj = Gene(field_dict)
				gene_dict[gobj.gene_id] = gobj
			elif field_dict["feature"] == "transcript":
				tobj = Transcript(field_dict)
				transcript_dict[tobj.transcript_id] = tobj
			elif field_dict["feature"] in ["exon","CDS","start_codon","stop_codon"]:
				tid = field_dict["attr"]["transcript_id"]
				try:
					transcript_dict[tid].add_feature(field_dict)
				except KeyError:
					raise ParsingError("Error in line %i. The annotation in GTF file should be three-level hierarchy of \
					                    gene => transcript => exon (or CDS)" % i)
			else:
				pass
	return gene_dict,transcript_dict

def get_chrom(name):
	if " " in name:
		return name.split()[0]
	elif "|" in name:
		return name.split("|")
	else:
		return name

class GenomeSeq(object):
	""" genomic sequence"""

	def __init__(self,filename):
		self.filename = filename
		self.fh = Fasta(filename, key_fn = get_chrom)
	def get_seq(self,chrom,start=0,end=False,strand="+"):
		if end is False:
			end = len(self.fh[chrom])
		return self.fh.sequence({"chr":chrom, "start":start, "stop": end, "strand":strand}, one_based=False)


def genomic_iv_transform(genomic_exons_sorted, genomic_ivs_sorted):
	"""
	transform the genomic interval (of cds, start and stop codons) to inner coordinate in transcript,
	transcript_exons must be sorted
	return 0-based iv
	"""

	length = 0
	innerIV = []
	for i in genomic_exons_sorted:
		for j in genomic_ivs_sorted:
			if i.contain_iv(j):
				start = length + abs(j.start_d - i.start_d)
				end = start + (j.end - j.start)
				if innerIV and start == innerIV[-1].end:
					innerIV[-1] = Interval(innerIV[-1].start, end, "+")
				else:
					innerIV.append(Interval(start,end,"+"))
		length += i.end - i.start
	if len(innerIV) == 0:
		raise ValueError("\tCan't transform the genomic interval, please check!")
	return innerIV

def transcript_pos_transform(tobj, transcript_pos):
	"""
	transform the transcript position to genomic position
	0-based
	if strand is "-", return position of strand.
	"""
	length = 0
	strand = tobj.genomic_iv.strand
	genomic_pos = -1
	for i in tobj.genomic_exons:
		if transcript_pos < length + i.length:
			if strand == "+":
				genomic_pos = transcript_pos - length + i.start
			else:
				genomic_pos = i.start_d - (transcript_pos - length)
			break
		length += i.length
	return genomic_pos

def transcript_iv_transform(tobj, transcript_iv):
	"""
	transform the transcript interval to genomic interval, zero-based.
	"""
	strand = tobj.genomic_iv.strand
	genomic_start = transcript_pos_transform(tobj,transcript_iv.start)
	genomic_end = transcript_pos_transform(tobj,transcript_iv.end-1)
	if strand == "+":
		genomic_end += 1
	else:
		genomic_end -= 1
	exons_bound = [genomic_start]


	for i in tobj.genomic_exons:
		if strand == "+":
			if genomic_start < i.start < genomic_end:
				exons_bound.append(i.start)
			if genomic_start < i.end < genomic_end:
				exons_bound.append(i.end)
		else:
			if genomic_end < i.start_d < genomic_start:
				exons_bound.append(i.start_d)
			if genomic_end < i.end_d < genomic_start:
				exons_bound.append(i.end_d)

	exons_bound.append(genomic_end)
	if len(exons_bound) % 2 != 0:
		print("error when transform the transcript interval to genomic!")
	exons_ivs = []
	for i in range(0,len(exons_bound),2):
		exons_ivs.append(Interval_from_directional(exons_bound[i],exons_bound[i+1],strand))

	return exons_ivs

def load_transcripts_pickle(pickle_file):
	if os.path.exists(pickle_file):
		sys.stdout.write("\tLoading transcripts.pickle ...\n")
		with open(pickle_file,"rb") as fin:
			gene_dict, transcript_dict = pickle.load(fin)
	else:
		raise IOError("\tError, %s file not found\n" % pickle_file)

	return gene_dict,transcript_dict

def processTranscripts(genomeFasta,gtfFile,out_dir):
	"""
	transform the genomic position into the transcript position
	"""
	pickle_file = os.path.join(out_dir,"transcripts.pickle")
	if os.path.exists(pickle_file):
		gene_dict, transcript_dict = load_transcripts_pickle(pickle_file)
		return gene_dict, transcript_dict
	else:
		gene_dict,transcript_dict = readGTF(gtfFile)
		if not os.path.exists(genomeFasta):
			raise IOError("\tError, the genomic fasta file not found: %s" % genomeFasta)
		genomic_seq = GenomeSeq(genomeFasta)
		transcript_seq_file = open(os.path.join(out_dir,"transcripts_sequence.fa"),"w") #title: transcriptid length
		transcript_cds_file = open(os.path.join(out_dir,"transcripts_cds.txt"),"w") # tid, cds_start, cds_end

		sys.stderr.write("\tProcess the transcripts ....\n")
		for tobj in transcript_dict.values():
			# store the transcript id in gene object
			if tobj.transcript_id not in gene_dict[tobj.gene_id].transcripts:
				gene_dict[tobj.gene_id].transcripts.append(tobj.transcript_id)
				if "appris_principal_1" in tobj.attr.get("tag",[]) or "appris_principal" in tobj.attr.get("tag", []):
					gene_dict[tobj.gene_id].principal_transcripts.append(tobj.transcript_id)
			# sort the exons,cds,startcodon,stopcodon according their position on genomic
			if tobj.genomic_iv.strand == "+":
				sorted_reverse = False
			else:
				sorted_reverse = True

			tobj.genomic_exons.sort(key=lambda x: x.start, reverse=sorted_reverse)
			# convert the CDS genomic coordination into transcript coordination
			if tobj.genomic_cds:
				tobj.genomic_cds.sort(key=lambda x: x.start, reverse=sorted_reverse)
				tobj.cds = genomic_iv_transform(tobj.genomic_exons, tobj.genomic_cds)
				if len(tobj.cds) > 1:
					sys.stderr.write("Warning: the CDS is discontinuous," +
					                  "only first region is used, %s\n" % tobj.transcript_id)
				tobj.cds = tobj.cds[0]
				transcript_cds_file.write("%s\t%i\t%i\n" % (tobj.transcript_id,tobj.cds.start +1,tobj.cds.end + 3))

			# start and stop codon
			if tobj.genomic_startcodon:
				tobj.genomic_startcodon.sort(key=lambda x: x.start, reverse=sorted_reverse)
				tobj.startcodon = genomic_iv_transform(tobj.genomic_exons,tobj.genomic_startcodon)
				if len(tobj.startcodon) >1:
					sys.stderr.write("Warning: the start codon is discontinuous," +
					                "only first region is used, %s\n" % tobj.transcript_id)
				tobj.startcodon = tobj.startcodon[0]
			elif tobj.cds:
				tobj.startcodon = Interval(tobj.cds.start, tobj.cds.start + 3, "+")
			else:
				pass

			if tobj.genomic_stopcodon:
				tobj.genomic_stopcodon.sort(key=lambda x: x.start, reverse=sorted_reverse)
				tobj.stopcodon = genomic_iv_transform(tobj.genomic_exons,tobj.genomic_stopcodon)
				if len(tobj.stopcodon) >1:
					sys.stderr.write("Warning: the stop codon is discontinuous," +
					                "only first region is used, %s\n" % tobj.transcript_id)
				tobj.stopcodon = tobj.stopcodon[0]
			elif tobj.cds:
				tobj.stopcodon = Interval(tobj.cds.end, tobj.cds.end + 3, "+")
			else:
				pass

			transcript_seq = [genomic_seq.get_seq(tobj.chrom,exon_iv.start,exon_iv.end,exon_iv.strand) for exon_iv in tobj.genomic_exons]
			transcript_seq = "".join(transcript_seq)
			# write the transcript sequence to file
			transcript_seq_file.write(">%s %i\n%s\n" % (tobj.transcript_id,tobj.length,transcript_seq))

		transcript_seq_file.close()
		transcript_cds_file.close()

		## dump the pickle
		sys.stderr.write("\tSaving the transcripts.pickle\n")
		with open(pickle_file,"wb") as fout:
			pickle.dump([gene_dict,transcript_dict],fout,protocol=pickle.HIGHEST_PROTOCOL)

		return gene_dict,transcript_dict

def verboseprint(printstring):
	sys.stderr.write('[%s] %s\n' % (strftime("%Y-%m-%d %H:%M:%S"), printstring))
	sys.stderr.flush()


def main():
	from .parsing_opts import parsing_transcript
	args = parsing_transcript()
	verboseprint("Preparing annotation files ...")
	processTranscripts(args.genomeFasta,args.gtfFile,args.out_dir)
	verboseprint("The step of preparing transcript annotation has been completed.")
if __name__ == "__main__":
	main()
