#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'
from RiboCode import RiboCode

print "RiboCode",dir(RiboCode)

def test_prepare_transcripts():
	RiboCode.prepare_transcripts.main()

def test_metaplots():
	RiboCode.metaplots.main()

def test_plot_orf_density():
	RiboCode.plot_orf_density.main()

def test_ORFcount():
	RiboCode.RPF_count_ORF.main()
