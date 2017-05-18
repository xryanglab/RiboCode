#!/usr/bin/env python
# -*- coding:UTF-8 -*-
__author__ = 'Zhengtao Xiao'
from RiboCode import prepare_transcripts,metaplots,plot_orf_density

def test_prepare_transcripts():
	prepare_transcripts.main()

def test_metaplots():
	metaplots.main()

def test_plot_orf_density():
	plot_orf_density.main()
