#!/usr/bin/env bash

# bin_scrubbing.sh

bindir=$1
extension=$2
outdir=$3
contig_info=$4

mkdir -p $outdir
cd $outdir
if ! test -d $outdir/checkm.raw || ! test -f $outdir/checkm.raw/results.tsv || [ `cat $outdir/checkm.raw/results.tsv | wc -l` -eq 0 ]
then
	rm -rf $outdir/checkm.raw
	bsub -n16 -R"span[hosts=1]" \
		-o $outdir/checkm_raw.o \
		-e $outdir/checkm_raw.e \
		-M50000 \
		-R 'select[mem>50000] rusage[mem=50000]' \
			"singularity exec -B /lustre,/nfs \
				$LOCAL_IMAGES/checkm.sif \
				checkm lineage_wf \
					-t 16 -x $extension \
					--tab_table \
					-f $outdir/checkm.raw/results.tsv \
					$bindir \
					$outdir/checkm.raw"
fi
# Filter contigs in bins by size and circularity
## Keep all contigs >= 1MB
## Keep all circular contigs
## If no circular contigs >= 100KB...
##		- Keep contigs >= 100kb
rm -rf $outdir/output_bins
mkdir -p $outdir/output_bins
for bin_file in $bindir/*.$extension
do
	bin=`basename $bin_file .$extension`
	if ! test -f $outdir/output_bins/$bin.$extension || [ `head $outdir/output_bins/$bin.$extension  | wc -l` -eq 0 ]
	then
		rm -rf $outdir/tmp.list
		touch $outdir/tmp.list
		contigs=`grep $'^>' $bin_file | cut -d'>' -f2`
		large_circ="FALSE"
		for contig in $contigs
		do
			if [ `grep $'^'$contig$'\t' $contig_info | cut -f4` == "Y" ]
			then
				echo $contig >> $outdir/tmp.list
				if [ `grep $'^'$contig$'\t' $contig_info | cut -f2` -ge 100000 ]
				then
					large_circ="TRUE"
				fi
			elif [ `grep $'^'$contig$'\t' $contig_info | cut -f2` -ge 1000000 ]
			then
				echo $contig >> $outdir/tmp.list
			fi
		done
		if [ "$large_circ" == "FALSE" ]
		then
			fastalength $bin_file | awk '{if($1 >= 100000) {print $2}}' >> $outdir/tmp.list
		fi
		if [ `cat $outdir/tmp.list | wc -l` -gt 0 ] && [ `cat $outdir/tmp.list | wc -l` -lt `echo $contigs | wc -w` ]
		then
			bin=`basename $bin_file .$extension`
			echo $bin
			select_fasta_by_list.pl \
				-l $outdir/tmp.list \
				-i $bin_file \
				-o $outdir/output_bins/$bin.$extension
		fi
	fi
done

if ! test -d $outdir/checkm.scrubbed || ! test -f $outdir/checkm.scrubbed/results.tsv || [ `cat $outdir/checkm.scrubbed/results.tsv | wc -l` -eq 0 ]
then
	rm -rf $outdir/checkm.scrubbed
	bsub -n16 -R"span[hosts=1]" \
		-o $outdir/checkm_scrubbed.o \
		-e $outdir/checkm_scrubbed.e \
		-M50000 \
		-R 'select[mem>50000] rusage[mem=50000]' \
			"singularity exec -B /lustre,/nfs \
				$LOCAL_IMAGES/checkm.sif \
					checkm lineage_wf \
						-t 16 \
						-x $extension \
						--tab_table \
						-f $outdir/checkm.scrubbed/results.tsv \
						$outdir/output_bins \
						$outdir/checkm.scrubbed"
fi

while ! test -d $outdir/checkm.raw || ! test -f $outdir/checkm.raw/results.tsv || ! test -d $outdir/checkm.scrubbed || ! test -f $outdir/checkm.scrubbed/results.tsv
do
	sleep 5m
done

# Compare our cleaned bins with the original bins
## If cleaning the bin changed completness score, use old bin
for bin_file in $bindir/*.$extension
do
	bin=`basename $bin_file .$extension`
	if ! test -f $outdir/output_bins/$bin.$extension
	then
		cp $bin_file $outdir/output_bins/
	else
		new=`grep $'^'$bin$'\t' $outdir/checkm.scrubbed/results.tsv | awk '{print $13}'`
		old=`grep $'^'$bin$'\t' $outdir/checkm.raw/results.tsv | awk '{print $13}'`
		if [ "$old" != "$new" ]
		then
			cp $bin_file $outdir/output_bins/
		elif [ `grep $'^>' $outdir/output_bins/$bin.$extension | wc -l` -lt `grep $'^>' $bindir/$bin.$extension | wc -l` ]
		then
			mv $outdir/output_bins/$bin.$extension $outdir/output_bins/${bin}_scrubbed.$extension
		fi
	fi
done
contigs2bin.sh \
	$outdir/output_bins \
	$outdir/contigs2bin.tsv \
	$extension


