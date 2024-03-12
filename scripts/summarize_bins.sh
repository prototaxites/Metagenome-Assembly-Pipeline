#!/usr/bin/env bash

asm_fofn=$1
outfasta_dir=$2
outfile=$3

mkdir -p $outfasta_dir
rm -rf $outfile
while read line
do
	program=`echo $line | cut -d' ' -f1`
	asm_dir=`echo $line | cut -d' ' -f2`
	tolid=`basename $asm_dir | cut -d'.' -f1`
	date=`basename $asm_dir | cut -d'.' -f3`
	if `echo $program | grep -q "hifiasm"`
	then
		if test -f $asm_dir/$tolid.p_ctg.fa.gz
		then
			cat $asm_dir/$tolid.p_ctg.fa.gz \
				> $outfasta_dir/$program.$date.fa.gz
		elif test -f $asm_dir/$tolid.p_ctg.fa
		then
			gzip -c $asm_dir/$tolid.p_ctg.fa \
				> $outfasta_dir/$program.$date.fa.gz
		elif test -f $asm_dir/$tolid.p_ctg.gfa
		then
			gfa2fasta.sh $asm_dir/$tolid.p_ctg.gfa \
				| gzip -c - \
				> $outfasta_dir/$program.$date.fa.gz
		else
			echo "Could not find primary contig file in $asm_dir."
		fi
		if test -f $asm_dir/$tolid.a_ctg.fa.gz
		then
			cat $asm_dir/$tolid.a_ctg.fa.gz \
				>> $outfasta_dir/$program.$date.fa.gz
		elif test -f $asm_dir/$tolid.a_ctg.fa
		then
			gzip -c $asm_dir/$tolid.a_ctg.fa \
				>> $outfasta_dir/$program.$date.fa.gz
		elif test -f $asm_dir/$tolid.a_ctg.gfa
		then
			gfa2fasta.sh $asm_dir/$tolid.a_ctg.gfa \
				| gzip -c - \
				>> $outfasta_dir/$program.$date.fa.gz
		else
			echo "Could not find alternate contig file in $asm_dir."
		fi
	elif [ "$program" == "meta-mdbg" ]
	then
		if test -f $asm_dir/contigs.fasta.gz
		then
			cp $asm_dir/contigs.fasta.gz \
				$outfasta_dir/$program.$date.fa.gz
		elif test -f $asm_dir/contigs.fasta
		then
			gzip -c $asm_dir/contigs.fasta \
				> $outfasta_dir/$program.$date.fa.gz
		else
			echo "Cound not find contig file in $asm_dir"
		fi
	elif `echo $program | grep -q "flye"`
	then
		if test -f $asm_dir/assembly.fasta.gz
		then
			cp $asm_dir/assembly.fa.gz \
				$outfasta_dir/$program.$date.fa.gz
		elif test -f $asm_dir/assembly.fasta
		then
			gzip -c $asm_dir/assembly.fasta \
				> $outfasta_dir/$program.$date.fa.gz
		else
			echo "Cound not find contig file in $asm_dir"
		fi
	fi
	if test -d $asm_dir/mags
	then
		for group_dir in $asm_dir/mags/*
		do
			group=`basename $group_dir`
			for binner_dir in $group_dir/*
			do
				if test -d $binner_dir/output_bins
				then
					binner=`basename $binner_dir`
					refiner="NA"
					for file in $binner_dir/output_bins/*
					do
						if test -f $file && `grep -q $'^>' $file`
						then
							extension=`echo $file | awk -F',' '{print $NF}'`
							bin=`basename $file .$extension`
							contigs=`grep $'^>' $file | cut -d'>' -f2`
							for contig in $contigs
							do
								echo -e "$program.$date,$group,$binner,$refiner,$contig,$bin" >> $outfile
							done
						fi
					done
				fi
				if test -d $binner_dir/dastool && test -d $binner_dir/dastool/output_bins
				then
					binner=`basename $binner_dir`
					refiner="dastool"
					for file in $binner_dir/output_bins/*
					do
						if test -f $file && `grep -q $'^>' $file`
						then
							extension=`echo $file | awk -F'.' '{print $NF}'`
							bin=`basename $file .$extension`
							contigs=`grep $'^>' $file | cut -d'>' -f2`
							for contig in $contigs
							do
								echo -e "$program.$date,$group,$binner,$refiner,$contig,$bin" >> $outfile
							done
						fi
					done
				fi
			done
		done
	fi
done < $asm_fofn
				
	
	