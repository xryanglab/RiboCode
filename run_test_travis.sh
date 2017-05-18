#!/usr/bin/env bash

pushd .

WORKDIR=`pwd`

echo
echo "test prepare_transcripts"
echo

prepare_transcripts -h

if [ $? != 0 ]; then
	exit 1
fi
