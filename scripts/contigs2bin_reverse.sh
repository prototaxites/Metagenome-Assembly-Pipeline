#!/usr/bin/env bash

contig2bin_file=$1
assembly=$2
outdir=$3

mkdir -p $outdir

bins=`cut -f2 $contig2bin_file | sort -u`
for bin in $bins
do
	grep $'\t'$bin$'$' $contig2bin_file \
		| cut -f1 \
		> $outdir/tmp.list
	select_fasta_by_list.pl \
		-i $assembly \
		-l $outdir/tmp.list \
		-o $outdir/$bin.fa
done
rm -rf $outdir/tmp.list