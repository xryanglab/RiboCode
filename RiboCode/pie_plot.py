#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from builtins import str, range
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'
"""
===============
Basic pie chart
===============
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def percentage_format(value,decimal=2):
	"""
	Format number as percentage
	"""
	formatStr = "{:." + str(decimal) + "%}"
	return formatStr.format(value)
	#return "{:.2%}".format(value)

def get_colors(NUM_COLORS):
	cm = plt.get_cmap('gist_rainbow')
	colors = []
	for i in range(NUM_COLORS):
		colors.append(cm(1.*i/NUM_COLORS))  # color will now be an RGBA tuple
	return colors

def pie(orderDict, outname):
	sizes = list(orderDict.values())
	total = sum(sizes)
	labels = []
	for i in orderDict.keys():
		perct = orderDict[i] / float(total)
		labels.append("%s: %i (%s)" % (i,int(orderDict[i]),percentage_format(perct,1)))

	colors = get_colors(len(labels))
	plt.figure(figsize=(8,6))
	patches, text = plt.pie(sizes, startangle=90,colors=colors)
	plt.legend(patches,labels,loc=4,bbox_to_anchor=(1, 0),fontsize=15)
	plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
	plt.tight_layout()
	plt.savefig(outname + "_ORFs_category.pdf")
	plt.close()
