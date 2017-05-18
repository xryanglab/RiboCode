#!/usr/bin/env bash

pushd .

WORKDIR=`pwd`

#1
echo
echo "test prepare_transcripts"
echo

prepare_transcripts -h

if [ $? != 0 ]; then
	exit 1
fi

#2
echo
echo "test metaplots"
echo

metaplots -h

if [ $? != 0 ]; then
	exit 1
fi

#3
echo
echo "test plot_orf_density"
echo

plot_orf_density -h

if [ $? != 0 ]; then
	exit 1
fi
