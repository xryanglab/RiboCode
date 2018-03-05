#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import h5py
import os

def generate_colors(length,shiftn=0):
	"""
	generate r,b,g colors
	"""
	colors0 = ["red","blue","green"]
	colors1 = ["green","red","blue"]
	colors2 = ["blue","green","red"]
	rep_times = int(np.ceil(length/3))
	if shiftn == 0:
		colors = np.tile(colors0,rep_times)
	elif shiftn==1:
		colors = np.tile(colors1,rep_times)
	else:
		colors = np.tile(colors2,rep_times)
	return colors[:length]

def plot_ORF(ax,tpsites,orf_start):
	colors = generate_colors(tpsites.size,orf_start%3)
	ax.vlines(np.arange(tpsites.size),ymin=0,ymax=tpsites,colors=colors)
	ax.tick_params(axis='x',which="both",top="off",direction='out',labelsize=15)
	ax.tick_params(axis='y',which="both",labelsize=15)
	ax.set_ylabel("P-site reads density",fontsize=18)
	ax.set_xlim((0,tpsites.size))

def plot_annotation(ax,tlength,start,stop,label,color):
	width = 0.15
	ax.set_xlim((0,tlength))
	ax.fill((start-1,stop,stop,start-1),
	        (1+width/2,1+width/2,1-width/2,1-width/2),
			color=color,lw=0.5,zorder=20)
	ax.axhline(1,color="gray",lw=0.5)
	ax.set_frame_on(False)
	ax.xaxis.set_ticks_position("none")
	ax.yaxis.set_ticks_position("none")
	ax.set_xticks([])
	ax.set_yticks([1])
	ax.set_yticklabels([label],fontsize=18)
	ax.set_ylim((1-width/2,1+width/2))

def read_psites_array(filename,transcript_id):
	with h5py.File(filename,"r") as fin:
		k = fin["transcript_ids"][:]
		idx = np.where(k == transcript_id)[0][0]
		v = fin["p_sites"][idx]
	return v

def plot_main(cds_start,cds_end,psites_array,orf_tstart,orf_tstop,outname):
	"""
	the main plot function
	"""
	plt.figure(figsize=(8,4))
	if cds_start is not None:
		gs = gridspec.GridSpec(3,1,height_ratios=[10,1,1],hspace=0.6,left=0.2,right=0.95)
	else:
		gs = gridspec.GridSpec(2,1,height_ratios=[11,1],hspace=0.6,left=0.2,right=0.95)

	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])
	plot_ORF(ax1,psites_array,orf_tstart)
	plot_annotation(ax2,psites_array.size,orf_tstart,orf_tstop,"Predicted","#3994FF")

	if cds_start is not None:
		ax3 = plt.subplot(gs[2])
		plot_annotation(ax3,psites_array.size,cds_start,cds_end,"Annotated","#006DD5")
	# plt.tight_layout()
	plt.savefig(outname + ".pdf")

def main():
	from .parsing_opts import parsing_plot_orf_density
	from .loadconfig import LoadConfig
	args = parsing_plot_orf_density()
	transcripts_cds_file = os.path.join(args.annot_dir,"transcripts_cds.txt")
	cds_start = None
	cds_end = None
	with open(transcripts_cds_file) as fin:
		for line in fin:
			if line.startswith(args.transcript_id):
				_,s,e = line.strip().split("\t")
				cds_start = int(s)
				cds_end = int(e)
				break
	configIn = LoadConfig(args.config_file)
	samples = [i["samplename"] for i in configIn.configList]
	for i,s in enumerate(samples):
		if i == 0:
			psites_array = read_psites_array(s + "_psites.hd5",args.transcript_id)
		else:
			psites_array += read_psites_array(s + "_psites.hd5",args.transcript_id)
	if not args.outname:
		args.outname = "%s_%i_%i" % (args.transcript_id,args.orf_tstart,args.orf_tstop)
	plot_main(cds_start,cds_end,psites_array,args.orf_tstart,args.orf_tstop,args.outname)
