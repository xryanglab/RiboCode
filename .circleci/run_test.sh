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

#4
echo
echo "test ORFcount"
echo

ORFcount -h

if [ $? != 0 ]; then
	exit 1
fi

#5
echo
echo "test RiboCode"
echo

RiboCode -h

if [ $? != 0 ]; then
	exit 1
fi

#6
echo
echo "test RiboCode_onestep"
echo

RiboCode_onestep -h

if [ $? != 0 ]; then
	exit 1
fi

#7
echo
echo "test GTFupdate"
echo

GTFupdate -h

if [ $? != 0 ]; then
	exit 1
fi
