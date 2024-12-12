#!/usr/bin/env bash

bindir=$1
outfile=$2
extension=$3

awk -v ext=$extension \
        'BEGIN { OFS = "\t" }
        BEGINFILE {
                cmd=sprintf("basename %s .%s", FILENAME, ext)
                cmd | getline bin
        }
        /^>/ {
                sub(/>/, "", $1)
                print $1,bin
        }' $bindir/*.$extension > $outfile
