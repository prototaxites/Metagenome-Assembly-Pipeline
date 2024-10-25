#!/usr/bin/env bash

bindir=$1
outfile=$2
extension=$3

rm -rf $outfile
for file in $bindir/*.$extension
do
	bin=`basename $file .$extension`
	contigs=`grep $'^>' $file | cut -d'>' -f2`
	for contig in $contigs
	do
		echo -e "$contig\t$bin" >> $outfile
	done
done