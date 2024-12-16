#!/bin/bash

set -e

module load seqkit/2.8.2--h9ee0642_0

contig2bin=$1
contigs=$2
outdir=$3

awk '{print $2}' ${contig2bin} | sort -u | while read bin
do
	seqkit grep -f <(grep -w ${bin} ${contig2bin} | awk '{ print $1 }') ${contigs} > $3/${bin}.fa
done
