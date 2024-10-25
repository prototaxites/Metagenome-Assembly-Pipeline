#!/usr/bin/env bash

CHECKM_DATA_PATH=
GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release207_v2

usage="SYNOPSIS: map_runner.sh [OPTIONS]\n\
OPTIONS:\n\
    -D | --assembly_directories Tab-delimited file containing assembly directoy info\n\
                                  (Assemblers = [hifiasm, hifiasm-meta, hifiasm-hic, flye]) (REQUIRED).\n\
                                    Format --> 'assembler    assembly_directory'\n\
    -p | --prefix               Assembly name prefix\n\
    -A | --assembly_haps        Comma separated list of hap assemblies to bin [primary,alternate_primary].\n\
                                    Options:\n\
                                        alternate\n\
                                        primary\n\
                                        alternate_primary (bin alternate and primary assemblies\n\
                                                               separately and combine with refinement program(s))\n\
                                        all (bin combined alternate and primary assemblies)\n\
    -r | --reads                Long reads (REQUIRED)\n\	
    -H | --hic_directory        Directory containing HiC cram files (if running MetaTOR or bin3C).\n\
    -c | --hic_reads            Comma separated paths to HiC reads fastqs (can be gzipped). Overrides option (-H).\n\
    -B | --binning_programs     Comma separated list of binning programs to run. [metabat2,maxbin2]\n\
                                    Options:\n\
                                        metabat2\n\
                                        maxbin2\n\
                                        metator (requires HiC)\n\
                                        bin3c (requires HiC)\n\
    -R | --refining_programs.   Comma separated list of bin refinement programs [dastool,magscot].\n\
                                    Options:\n\
                                        none\n\
                                        dastool\n\
                                        magscot\n\
    -O | --outdir               Output directory [./map_outputs]\n\
    -C | --checkm_db            CheckM database (see https://github.com/Ecogenomics/CheckM/wiki/Installation)\n\
	                                [$CHECKM_DATA_PATH]\n\
    -G | --gtdbtk_db            GTDB-TK database (see https://ecogenomics.github.io/GTDBTk/installing/index.html)\n\
                                     [$GTDBTK_DB]\n\
    -h | --help                 Print help message\n\n\
AUTHOR:\n\
Noah Gettle 2023"

outdir="./map_outputs"
assembly_haps="primary alternate_primary"
binning_programs=""
refining_programs="dastool,magscot"
assembly_directories=""
reads=""
hic_directory=""
hic_reads=""
threads=16
while [ "$1" != "" ]
do
	case $1 in
		-D | --assembly_directories)
			if ! `beginswith "-" $2`
			then
				shift
				assembly_directories=$1
				if ! test -f $assembly_directories
				then
					echo -e "ERROR: Could not find assembly directory list $assembly_directories.\n\n$usage"
					exit 1
				elif [ `head $assembly_directories | wc -l` -eq 0 ]
				then
					echo -e "ERROR: $assembly_directories file appears to be empty.\n\n$usage"
					exit 1
				fi
				for assembler in `cut -f1 $assembly_directories`
				do
					if [ "$assembler" != "hifiasm" ] && [ "$assembler" != "hifiasm-meta" ] && [ "$assembler" != "hifiasm-hic" ] && [ "$assembler" != "flye" ] && [ "$assembler" != "metaflye" ] && [ "$assembler" != "meta-mdbg" ]
					then
						echo -e "ERROR: Assembler $assembler in $assembly_directories is not a valid value.\n\n$usage"
						exit 1
					fi
				done
				for dir in `cut -f2 $assembly_directories`
				do
					if ! test -d $dir
					then
						echo -e "ERROR: $dir in $assembly_directories does not exist.\n\n$usage"
						exit 1
					fi
				done
			fi;;
		-p | --prefix)
			if ! `beginswith "-" $2`
			then
				shift
				prefix=$1
			fi;;
		-A | --assembly_haps)
			if ! `beginswith "-" $2`
			then
				shift
				assembly_haps=`echo $1 | sed "s|,| |g"`
				for hap in $assembly_haps
				do
					if [ "$hap" != "primary" ] && [ "$hap" != "alternate" ] && [ "$hap" != "alternate_primary" ] && [ "$hap" != "all" ]
					then
						echo -e "ERROR: $hap is not a valid value for -A.\n\n$usage"
						exit 1
					fi
				done
			fi;;
		-r | --reads)
			if ! `beginswith "-" $2`
			then
				shift
				reads=$1
				if ! test -f $reads
				then
					echo -e "ERROR: Could not find reads file $reads.\n\n$usage"
					exit 1
				fi
			fi;;
		-H | --hic_directory)
			if ! `beginswith "-" $2`
			then
				shift
				hic_directory=$1
				if ! test -d $hic_directory
				then
					echo -e "ERROR: HiC directory $hic_directory does not exist.\n\n$usage"
					exit 1
				fi
			fi;;
		-c | --hic_reads)
			if ! `beginswith "-" $2`
			then
				shift
				hic_reads=$1
				for file in `echo $hic_reads | sed "s|,| |g"`
				do
					if ! test -f $file
					then
						echo -e "ERROR: Could not find HiC reads $file.\n\n$usage"
						exit 1
					fi
				done
			fi;;
		-B | --binning_programs)
			if ! `beginswith "-" $2`
			then
				shift
				binning_programs=$1
				for binner in `echo $binning_programs | sed "s|,| |g"`
				do
					if [ "$binner" != "bin3c" ] && [ "$binner" != "metator" ] && [ "$binner" != "metabat2" ] && [ "$binner" != "maxbin2" ]
					then
						echo -e "ERROR: $binner is not a valid value for -B.\n\n$usage"
						exit 1
					fi
				done
			fi;;
		-R | --refining_programs)
			if ! `beginswith "-" $2`
			then
				shift
				refining_programs=$1
				for refiner in `echo $refining_programs | sed "s|,| |g"`
				do
					if [ "$refiner" != "dastool" ] && [ "$refiner" != "magscot" ]
					then
						echo -e "ERROR: $refiner is not a valid value for -R.\n\n$usage"
						exit 1
					fi
				done
			fi;;
		-O | --outdir)
			if ! `beginswith "-" $2`
			then
				shift
				outdir=$1
				mkdir -p $outdir
			fi;;
		-C | --checkm_db)
			if ! `beginswith "-" $2`
			then
				shift
				checkm_db=$1
				if ! test -d $checkm_db
				then
					echo -e "ERROR: Could not find CheckM database directory $checkm_db.\n\n$usage"
					exit 1
				else
					export CHECKM_DATA_PATH=$checkm_db
				fi
			fi;;
		-G | --gtdbtk_db)
			if ! `beginswith "-" $2`
			then
				shift
				gtdbtk_db=$1
				if ! test -d $gtdbtk_db
				then
					echo -e "ERROR: Could not find GTDB-Tk database directory $cgtdbtk_db.\n\n$usage"
					exit 1
				else
					export GTDBTK_DB=$gtdbtk_db
				fi
			fi;;
		-h | --help)
			echo -e $usage
			exit 0;;
		*)
			echo -e "Option $1 not valid\n\n$usage"
			exit 1
	esac
	shift
done

###################################################################
##### Sanity Check ######
#########################

# Check for missing required values
if [ "$assembly_directories" == "" ]
then
	echo -e "ERROR: No assembly directory list provided.\n\n$usage"
	exit 1
fi
if [ "$reads" == "" ]
then
	echo -e "ERROR: No long read file provided.\n\n$usage"
	exit 1
fi
if [ "$binning_programs" == "" ]
then
	binning_programs="metabat2,maxbin2"
	if [ "$hic_directory" != "" ] || [ "$hic_reads" != "" ]
	then
		binning_programs="$binning_programs,bin3c,metator"
	fi
fi
if [ "$hic_directory" == "" ] && [ "$hic_reads" == "" ]
then
	if `echo $binning_programs | grep -q "bin3c"` || `echo $binning_programs | grep -q "metator"`
	then
		echo "ERROR: Some of the binning tools requested require HiC data.\n\n$usage"
		exit 1
	fi
fi

alternate_primary=FALSE
if `echo $assembly_haps | grep -q "alternate_primary"`
then
	alternate_primary=TRUE
	assembly_haps=`echo $assembly_haps | sed "s|alternate_primary||g"`
	if ! `echo $assembly_haps | grep -q "alternate"`
	then
		assembly_haps="$assembly_haps alternate"
	fi
fi
###################################################################
# Convert HiC reads from cram to fastq if necessary
mkdir -p $outdir

if [ "$hic_directory" != "" ] && [ "$hic_reads" == "" ]
then
	if ! test -f $outdir/hic.1.fastq && ! test -f $outdir/hic.1.fastq.gz
	then
		for cram in $hic_directory/*.cram
		do
			samtools fastq -n -f0x40 -F0xB00 $cram >> $outdir/hic.1.fastq
			samtools fastq -n -f0x80 -F0xB00 $cram >> $outdir/hic.2.fastq
		done
	fi
	if ! test -f $outdir/hic.1.fastq.gz
	then
		gzip $outdir/hic.1.fastq
		gzip $outdir/hic.2.fastq
	fi
	hic_reads="$outdir/hic.1.fastq.gz,$outdir/hic.2.fastq.gz"
fi


# Check and prep inputs
while read line
do
	assembler=`echo $line | cut -d' ' -f1`
	main_dir=`echo $line | cut -d' ' -f2`
	dir_name=`basename $main_dir`
	assembler_parts=$assembly_haps
	# Record the locations of unzipped assemblies for each part
	rm -rf $outdir/$dir_name.assemblies.txt
	contig_info_arg=""
	rm -rf $outdir/c1.txt
	rm -rf $outdir/c2.txt
	echo -e "PREPARING INPUTS FOR RUN: $assembler"
	if [ "$assembler" == "metaflye" ] || [ "$assembler" == "flye" ]
	then
		assembler_parts="all"
		if ! test -f $main_dir/assembly.fasta
		then
			if test -f $main_dir/assembly.fasta.gz
			then
				gunzip -c $main_dir/assembly.fasta.gz \
					> $outdir/$dir_name.tmp.fa
				echo -e "all\t$outdir/$dir_name.tmp.fa" >> $outdir/$dir_name.assemblies.txt
			else
				echo -e "ERROR: Could not find expected assembly file 'assembly.fasta[.gz]' in $main_dir.\n\n"
				exit 1
			fi
		else
			echo -e "all\t$main_dir/assembly.fasta" >> $outdir/$dir_name.assemblies.txt
		fi
		if test -f $main_dir/contig_info.tsv
		then
			contig_info_input="-c $main_dir/contig_info.tsv"
		fi
	elif [ "$assembler" == "meta-mdbg" ]
	then
		assembler_parts="all"
		if ! test -f $main_dir/contigs.fasta
		then
			if test -f $main_dir/contigs.fasta.gz
			then
				gunzip -c $main_dir/contigs.fasta.gz \
					> $outdir/$dir_name.tmp.fa
				asm=$outdir/$dir_name.tmp.fa
				echo -e "all\t$outdir/$dir_name.tmp.fa" >> $outdir/$dir_name.assemblies.txt
			else
				echo -e "ERROR: Could not find expected assembly file 'contigs.fasta[.gz]' in $main_dir.\n\n"
				exit 1
			fi
		else
			asm=$main_dir/contigs.fasta
			echo -e "all\t$main_dir/contigs.fasta" >> $outdir/$dir_name.assemblies.txt
		fi
		if ! test -f $main_dir/contig_info.tsv
		then
			fastalength $asm \
				| awk '{OFS="\t"; print $2"\t"$1"\tNA\t"$2}' \
				| sed "s|\tctg.*c$|\tY|g" | sed "s|\tctg.*l$|\tN|g" \
				>  $main_dir/contig_info.tsv
		fi
		contig_info_input="-c $main_dir/contig_info.tsv"
	elif [ "$assembler" == "hifiasm-hic" ]
	then
		rm -rf $outdir/$dir_name.contig_info.tsv
		if `echo $assembler_parts | grep -q 'primary'` || `echo $assembler_parts | grep -q 'all'`
		then
			if ! test -f $main_dir/$prefix.hap1.p_ctg.fa
			then
				if ! test -f $main_dir/$prefix.hap1.p_ctg.fa.gz
				then
					if ! test -f $main_dir/$prefix.hap1.p_ctg.gfa
					then
						echo -e "ERROR: Could not find expected assembly file '$prefix.hap1.p_ctg[.fa][.gz][.gfa]' in $main_dir"
						exit 1
					else
						gfa2fasta.sh $main_dir/$prefix.hap1.p_ctg.gfa \
							> $outdir/$dir_name.primary.tmp.fa
						echo -e "primary\t$outdir/$dir_name.primary.tmp.fa" >> $outdir/$dir_name.assemblies.txt
					fi
				else
					gunzip -c $main_dir/$prefix.hap1.p_ctg.fa.gz \
						>  $outdir/$dir_name.primary.tmp.fa
					echo -e "primary\t$outdir/$dir_name.primary.tmp.fa" >> $outdir/$dir_name.assemblies.txt
				fi
			else
				echo -e "primary\t$main_dir/$prefix.hap1.p_ctg.fa" >> $outdir/$dir_name.assemblies.txt
			fi
			if test -f $main_dir/$prefix.hap1.p_ctg.gfa
			then
				grep 'rd:i' $main_dir/$prefix.hap1.p_ctg.gfa | cut -f2,4,5 \
					| sed "s|LN:i:||g" | sed "s|rd:i:||g" \
					>> $outdir/c1.txt
			fi
		fi
		if `echo $assembler_parts | grep -q 'alternate'` || `echo $assembler_parts | grep -q 'all'`
		then
			if ! test -f $main_dir/$prefix.hap2.p_ctg.fa
			then
				if ! test -f $main_dir/$prefix.hap2.p_ctg.fa.gz
				then
					if ! test -f $main_dir/$prefix.hap2.p_ctg.gfa
					then
						echo -e "ERROR: Could not find expected assembly file '$prefix.hap2.p_ctg[.fa][.gz][.gfa]' in $main_dir"
						exit 1
					else
						gfa2fasta.sh $main_dir/$prefix.hap2.p_ctg.gfa \
							> $outdir/$dir_name.alternate.tmp.fa
						echo -e "alternate\t$outdir/$dir_name.alternate.tmp.fa" >> $outdir/$dir_name.assemblies.txt
					fi
				else
					gunzip -c $main_dir/$prefix.hap1.p_ctg.fa.gz \
						>  $outdir/$dir_name.alternate.tmp.fa
					echo -e "alternate\t$outdir/$dir_name.alternate.tmp.fa" >> $outdir/$dir_name.assemblies.txt
				fi
			else
				echo -e "alternate\t$main_dir/$prefix.hap2.p_ctg.fa" >> $outdir/$dir_name.assemblies.txt
			fi
			if test -f $main_dir/$prefix.hap2.p_ctg.gfa
			then
				grep 'rd:i' $main_dir/$prefix.hap2.p_ctg.gfa | cut -f2,4,5 \
					| sed "s|LN:i:||g" | sed "s|rd:i:||g" \
					>> $outdir/c1.txt
			fi
		fi
		if test -f  $outdir/c1.txt
		then
			cut -f1  $outdir/c1.txt \
				| awk '{print substr($0,length,1)}' \
				| sed "s|l|N|g" | sed "s|c|Y|g" \
				> $outdir/c2.txt
			paste $outdir/c1.txt $outdir/c2.txt > $outdir/$dir_name.contig_info.tsv
			rm -rf $outdir/c1.txt $outdir/c2.txt
			contig_info_input="-c $outdir/$dir_name.contig_info.tsv"
		fi
	else
		if `echo $assembler_parts | grep -q 'primary'` || `echo $assembler_parts | grep -q 'all'`
		then
			if ! test -f $main_dir/$prefix.p_ctg.fa
			then
				if ! test -f $main_dir/$prefix.p_ctg.fa.gz
				then
					if ! test -f $main_dir/$prefix.p_ctg.gfa
					then
						echo -e "ERROR: Could not find expected assembly file '$prefix.p_ctg[.fa][.gz][.gfa]' in $main_dir"
						exit 1
					else
						gfa2fasta.sh $main_dir/$prefix.p_ctg.gfa | sed "s|ctg|ptg|g" \
							> $outdir/$dir_name.primary.tmp.fa
						echo -e "primary\t$outdir/$dir_name.primary.tmp.fa" >> $outdir/$dir_name.assemblies.txt
					fi
				else
					gunzip -c $main_dir/$prefix.p_ctg.fa.gz | sed "s|ctg|ptg|g" \
						>  $outdir/$dir_name.primary.tmp.fa
					echo -e "primary\t$outdir/$dir_name.primary.tmp.fa" >> $outdir/$dir_name.assemblies.txt
				fi
			else
				if [ "$assembler" == "hifiasm-meta" ]
				then
					sed "s|ctg|ptg|g" $main_dir/$prefix.p_ctg.fa \
						> $outdir/$dir_name.primary.tmp.fa
					echo -e "primary\t$outdir/$dir_name.primary.tmp.fa" >> $outdir/$dir_name.assemblies.txt
				else
					echo -e "primary\t$main_dir/$prefix.p_ctg.fa" >> $outdir/$dir_name.assemblies.txt
				fi
			fi
			if test -f $main_dir/$prefix.p_ctg.gfa
			then
				grep 'rd:i' $main_dir/$prefix.p_ctg.gfa | cut -f2,4,5 \
					| sed "s|LN:i:||g" | sed "s|rd:i:||g" | sed "s|ctg|ptg|g" \
					> $outdir/c1.txt
				grep 'dp:f' $main_dir/$prefix.p_ctg.gfa | cut -f2,4,5 \
					| sed "s|LN:i:||g" | sed "s|dp:f:||g" | sed "s|ctg|ptg|g" \
					>> $outdir/c1.txt
			fi
		fi
		if `echo $assembler_parts | grep -q 'alternate'` || `echo $assembler_parts | grep -q 'all'`
		then
			if ! test -f $main_dir/$prefix.a_ctg.fa
			then
				if ! test -f $main_dir/$prefix.a_ctg.fa.gz
				then
					if ! test -f $main_dir/$prefix.a_ctg.gfa
					then
						echo -e "ERROR: Could not find expected assembly file '$prefix.a_ctg[.fa][.gz][.gfa]' in $main_dir"
						exit 1
					else
						gfa2fasta.sh $main_dir/$prefix.a_ctg.gfa | sed "s|ctg|atg|g" \
							> $outdir/$dir_name.alternate.tmp.fa
						echo -e "alternate\t$outdir/$dir_name.alternate.tmp.fa" >> $outdir/$dir_name.assemblies.txt
					fi
				else
					gunzip -c $main_dir/$prefix.a_ctg.fa.gz | sed "s|ctg|atg|g" \
						>  $outdir/$dir_name.alternate.tmp.fa
					echo -e "alternate\t$outdir/$dir_name.alternate.tmp.fa" >> $outdir/$dir_name.assemblies.txt
				fi
			else
				if [ "$assembler" == "hifiasm-meta" ]
				then
					sed "s|ctg|atg|g" $main_dir/$prefix.a_ctg.fa \
						> $outdir/$dir_name.alternate.tmp.fa
					echo -e "alternate\t$outdir/$dir_name.alternate.tmp.fa" >> $outdir/$dir_name.assemblies.txt
				else
					echo -e "alternate\t$main_dir/$prefix.a_ctg.fa" >> $outdir/$dir_name.assemblies.txt
				fi
			fi
			if test -f $main_dir/$prefix.a_ctg.gfa
			then
				grep 'rd:i' $main_dir/$prefix.a_ctg.gfa | cut -f2,4,5 \
					| sed "s|LN:i:||g" | sed "s|rd:i:||g" | sed "s|ctg|atg|g" \
					>> $outdir/c1.txt
				grep 'dp:f' $main_dir/$prefix.a_ctg.gfa | cut -f2,4,5 \
					| sed "s|LN:i:||g" | sed "s|dp:f:||g" | sed "s|ctg|atg|g" \
					>> $outdir/c1.txt
			fi
		fi
		if test -f $outdir/c1.txt
		then
			cut -f1  $outdir/c1.txt \
				| awk '{print substr($0,length,1)}' \
				| sed "s|l|N|g" | sed "s|c|Y|g" \
				> $outdir/c2.txt
			paste $outdir/c1.txt $outdir/c2.txt > $outdir/$dir_name.contig_info.tsv
			rm -rf $outdir/c1.txt $outdir/c2.txt
			contig_info_input="-c $outdir/$dir_name.contig_info.tsv"
		fi
	fi
	assemblies=`cut -f2 $outdir/$dir_name.assemblies.txt`
	rm -rf $outdir/$dir_name.all.tmp.fa
	for assembly in $assemblies
	do
		cat $assembly >> $outdir/$dir_name.all.tmp.fa
	done
	if `echo $binning_programs | grep -q "metabat2"` || `echo $binning_programs | grep -q "maxbin2"`
	then
		echo "line 511"
		if ! test -d $main_dir/mags/all || ! test -f $main_dir/mags/all/depth.tsv
		then
			echo "line 514"
			mkdir -p $main_dir/mags/all
			assembly=$outdir/$dir_name.all.tmp.fa
			echo "#### Mapping reads to assembly ####"
			if ! test -f $assembly.fai
			then
				samtools faidx $assembly
			fi
			bsub -n16 -q long -R"span[hosts=1]" \
				-o $outdir/minimap.$dir_name.%J.o \
				-e $outdir/minimap.$dir_name.%J.e \
				-M100000 \
				-R "select[mem>100000] rusage[mem=100000]" \
					"singularity exec -B /lustre,/nfs \
						$LOCAL_IMAGES/minimap2.sif \
							minimap2 \
								-a -t $threads \
								-x map-hifi \
								$assembly $reads \
								| samtools view -T $assembly -b - \
								| samtools sort -T $assembly -@ $threads - \
								> $main_dir/mags/all/$prefix.$assembler.all.bam \
					&& singularity exec -B /lustre,/nfs \
						$LOCAL_IMAGES/metabat.sif \
							jgi_summarize_bam_contig_depths \
								--outputDepth $main_dir/mags/all/depth.tsv \
								$main_dir/mags/all/$prefix.$assembler.all.bam"
		fi
	fi
	if `echo $refining_programs | grep -q 'magscot'`
	then
		hmm_dir=$main_dir/mags/all/hmm
		mkdir -p $hmm_dir
		if ! test -f $hmm_dir/input.hmm
		then
			assembly=$outdir/$dir_name.all.tmp.fa
			bsub -n8 -q long -R"span[hosts=1]" \
				-o $outdir/hmm.$dir_name.%J.o \
				-e $outdir/hmm.$dir_name.%J.e \
				-M20000 \
				-R "select[mem>20000] rusage[mem=20000]" \
					"singularity exec -B /lustre,/nfs \
						$LOCAL_IMAGES/prodigal.sif \
							prodigal \
								-i $assembly \
								-p meta \
								-a $hmm_dir/prodigal.faa \
								-d $hmm_dir/prodigal.ffn \
								-o $hmm_dir/tmpfile \
					&& singularity exec -B /lustre,/nfs \
						$LOCAL_IMAGES/hmmer.sif \
							hmmsearch \
								-o $hmm_dir/hmm.tigr.out \
								--tblout $hmm_dir/hmm.tigr.hit.out \
								--noali --notextw --cut_nc --cpu $threads \
								$GTDBTK_DB/hmm/gtdbtk_rel207_tigrfam.hmm \
								$hmm_dir/prodigal.faa \
					&& singularity exec -B /lustre,/nfs \
						$LOCAL_IMAGES/hmmer.sif \
							hmmsearch \
								-o $hmm_dir/hmm.pfam.out \
								--tblout $hmm_dir/hmm.pfam.hit.out \
								--noali --notextw --cut_nc --cpu $threads \
								$GTDBTK_DB/hmm/gtdbtk_rel207_Pfam-A.hmm \
								$hmm_dir/prodigal.faa"
		fi
	fi
done < $assembly_directories


if [ "$hic_reads" != "" ]
then
	input_hic="-H $hic_reads"
fi
while read line
do
	assembler=`echo $line | cut -d' ' -f1`
	main_dir=`echo $line | cut -d' ' -f2`
	dir_name=`basename $main_dir`
	counter=0
	while ! test -f $main_dir/mags/all/depth.tsv && [ $counter -le 1440 ]
	do
		if [ $counter -eq 0 ]
		then
			date=`date`
			echo "$date : Waiting for completion of read mapping for $dir_name"
		fi
		sleep 1m
		counter=`expr $count + 1`
	done
	if `echo $refining_programs | grep -q 'magscot'`
	then
		counter=0
		hmm_dir=$main_dir/mags/all/hmm
		while (! test -f $hmm_dir/hmm.tigr.hit.out || ! test -f $hmm_dir/hmm.pfam.hit.out) && [ $counter -le 1440 ]
		do
			if [ $counter -eq 0 ]
			then
				date=`date`
				echo "$date : Waiting for completion of hmm creation for $dir_name"
			fi
			sleep 1m
			counter=`expr $counter + 1`
		done
		if ! test -f $hmm_dir/input.hmm || [ `head $hmm_dir/input.hmm | wc -l` -eq 0 ]
		then
			cat $hmm_dir/hmm.tigr.hit.out \
				| grep -v "^#" \
				| awk '{print $1"\t"$3"\t"$5}' \
				> $hmm_dir/tigr
			cat $hmm_dir/hmm.pfam.hit.out \
				| grep -v "^#" \
				| awk '{print $1"\t"$4"\t"$5}' \
				> $hmm_dir/pfam
			cat $hmm_dir/pfam $hmm_dir/tigr > $hmm_dir/input.hmm
		fi
	fi
	assembler_parts=$assembly_haps
	contig_info_arg=""
	echo -e "SETTING UP RUN: $assembler"
	if [ "$assembler" == "metaflye" ] || [ "$assembler" == "flye" ] || [ "$assembler" == "meta-mdbg" ]
	then
		assembler_parts="all"
		if test -f $main_dir/contig_info.tsv
		then
			contig_info_input="-c $main_dir/contig_info.tsv"
		fi
	elif `echo $assembler | grep -q 'hifiasm'`
	then
		contig_info_input="-c $outdir/$dir_name.contig_info.tsv"
	fi
	assembler_parts=`echo $assembler_parts | sed "s|alternate_primary||g"`
	# Start the binning! #
	for part in $assembler_parts
	do
		echo -e "\t\t$part"
		mkdir -p $main_dir/mags/$part
		# rm -rf $main_dir/mags/$part/complete
		if ! test -f $main_dir/mags/$part/complete
		then
			if ! test -f $main_dir/mags/$part/depth.tsv
			then
				cp $main_dir/mags/all/depth.tsv $main_dir/mags/$part/depth.tsv
			fi
			if `echo $refining_programs | grep -q 'magscot'` && ! test -d $main_dir/mags/$part/hmm
			then
				ln -fs $main_dir/mags/all/hmm $main_dir/mags/$part/hmm
			fi
			echo "Running"
			rm -rf $main_dir/mags/$part/logs/*
			if [ "$part" == "all" ]
			then
				assembly=$outdir/$dir_name.all.tmp.fa
			else
				assembly=`grep $'^'$part$'\t' $outdir/$dir_name.assemblies.txt | cut -f2`
			fi
			mkdir -p  $main_dir/mags/$part/logs
			rm -rf  $main_dir/mags/$part/logs/*
			# prefix_input=${dir_name}_${part}
			prefix_input=$prefix.$assembler.$part
			mem=1000
			# rm -rf $main_dir/mags/$part/complete
			bsub -n16 -q long -R"span[hosts=1]" \
				-o $main_dir/mags/$part/logs/metagenome_assembly_pipeline.%J.o \
				-e $main_dir/mags/$part/logs/metagenome_assembly_pipeline.%J.e \
				-M$mem \
				-R "select[mem>$mem] rusage[mem=$mem]" \
					"metagenome_assembly_pipeline.sh -s \
						-P $binning_programs \
						-R $refining_programs \
						-a $assembly \
						$contig_info_input \
						-p $prefix_input \
						-r $reads \
						$input_hic \
						-o $main_dir/mags/$part \
						-t $threads \
						-E $main_dir/mags/$part/errors.txt \
						-D $main_dir/mags/$part/bin_stats.csv"
					# metagenome_assembly_pipeline.sh -s \
					# 	-P $binning_programs \
					# 	-R $refining_programs \
					# 	-a $assembly \
					# 	$contig_info_input \
					# 	-p $prefix_input \
					# 	-r $reads \
					# 	$input_hic  \
					# 	-o $main_dir/mags/$part \
					# 	-t $threads \
					# 	-E $main_dir/mags/$part/errors.txt \
					# 	-D $main_dir/mags/$part/bin_stats.csv
		fi
	done
done < $assembly_directories

# while read line
# do
# 	assembler=`echo $line | cut -d' ' -f1`
# 	main_dir=`echo $line | cut -d' ' -f2`
# 	dir_name=`basename $main_dir`
# 	assembler_parts=`echo $assembly_haps | sed "s|alternate_primary||g"`
# 	if [ "$assembler" == "flye" ] || [ "$assembler" == "metaflye" ] || [ "$assembler" == "meta-mdbg" ]
# 	then
# 		assembler_parts="all"
# 	fi
# 	for part in $assembler_parts
# 	do
# 		counter=0
# 		while ! test -f $main_dir/mags/$part/complete && [ $counter -le 1440 ]
# 		do
# 			if [ $counter -eq 0 ]
# 			then
# 				date=`date`
# 				echo "$date : Waiting for completion of $dir_name - $part"
# 			fi
# 			sleep 1m
# 			counter=`expr $counter + 1`
# 		done
# 		if test -f $main_dir/mags/$part/bin_stats.csv
# 		then
# 			if ! test -f $outdir/bin_stats.csv
# 			then
# 				cp $main_dir/mags/$part/bin_stats.csv $outdir/bin_stats.csv
# 			else
# 				tail -n +2 $main_dir/mags/$part/bin_stats.csv >> $outdir/bin_stats.csv
# 			fi
# 		fi
# 	done
# done < $assembly_directories
#
# if [ "$alternate_primary" == "TRUE" ]
# then
# 	while read line
# 	do
# 		assembler=`echo $line | cut -d' ' -f1`
# 		main_parts="alternate primary"
# 		main_dir=`echo $line | cut -d' ' -f2`
# 		mkdir -p $main_dir/mags/alternate_primary
# 		touch $main_dir/mags/alternate_primary/error.txt
# 		sed -i "/magscot/d" $main_dir/mags/alternate_primary/error.txt
# 		dir_name=`basename $main_dir`
# 		contig_info_input="-c $outdir/$dir_name.contig_info.tsv"
# 		prefix_input=$prefix.$assembler.alternate_primary
# 		maps=""
# 		rm -rf $main_dir/mags/alternate_primary/bin_stats.csv
# 		if [ "$assembler" != "metaflye" ] && [ "$assembler" != "flye" ]
# 		then
# 			echo -e "RUNNING ALTERNATE-PRIMARY BINNING FOR $assembler"
# 			for part in $main_parts
# 			do
# 				for program in `echo $binning_programs | sed "s|,| |g"`
# 				do
# 					if test -d $main_dir/mags/$part/$program
# 					then
# 						if test -f $main_dir/mags/$part/$program/contigs2bin.tsv && [ `head $main_dir/mags/$part/$program/contigs2bin.tsv | wc -l` -gt 0 ]
# 						then
# 							if [ "$maps" == "" ]
# 							then
# 								maps=$main_dir/mags/$part/$program/contigs2bin.tsv
# 							else
# 								maps="$maps,$main_dir/mags/$part/$program/contigs2bin.tsv"
# 							fi
# 						fi
# 					fi
# 				done
# 			done
# 			if ! test -f $outdir/$dir_name.all.tmp.fa
# 			then
# 				assemblies=`cut -f2 $outdir/$dir_name.assemblies.txt`
# 				for assembly in $assemblies
# 				do
# 					cat $assembly >> $outdir/$dir_name.all.tmp.fa
# 				done
# 			fi
# 			assembly=$outdir/$dir_name.all.tmp.fa
# 			# if `echo $refining_programs | grep -q 'magscot'` && ! test -d $main_dir/mags/alternate_primary/hmm
# 			# then
# 			# 	mkdir -p $main_dir/mags/alternate_primary/hmm
# 			# 	ln -fs $main_dir/mags/all/hmm $main_dir/mags/alternate_primary
# 			# fi
# 			for refiner in `echo $refining_programs | sed "s|,| |g"`
# 			do
# 				mkdir -p $main_dir/mags/alternate_primary/${refiner}_raw
# 				rm -rf $main_dir/mags/alternate_primary/${refiner}_raw/binstats_complete
# 				rm -rf $main_dir/mags/alternate_primary/${refiner}_raw/gtdb/output/complete
# 				if ! test -f $main_dir/mags/alternate_primary/${refiner}_raw/binstats_complete
# 				then
# 					mkdir -p $main_dir/mags/alternate_primary/logs
# 					if test -d $main_dir/mags/alternate_primary/${refiner}_raw/gtdb/genomes
# 					then
# 						cp -R $main_dir/mags/alternate_primary/${refiner}_raw/gtdb/genomes $main_dir/mags/alternate_primary/${refiner}_raw/output_bins
# 					fi
# 					# rm -rf $main_dir/mags/alternate_primary/${refiner}_raw/binstats_complete
# 					echo "RUNNING $refiner for alternate_primary : $prefix_input"
# 					rm -rf $main_dir/mags/alternate_primary/${refiner}_raw/prokka
# 					bsub -n$threads -q long -R"span[hosts=1]" \
# 						-o $main_dir/mags/alternate_primary/logs/${refiner}_raw.%J.o \
# 						-e $main_dir/mags/alternate_primary/logs/${refiner}_raw.%J.e \
# 						-M50000 \
# 						-R 'select[mem>50000] rusage[mem=50000]' \
# 							"${refiner}_wrapper.sh -e \
# 								-H $main_dir/mags/all/hmm \
# 								-a $assembly \
# 								-i $maps \
# 								$contig_info_input \
# 								-d $main_dir/mags/all/depth.tsv \
# 								-b ${refiner}.raw \
# 								-p $prefix_input \
# 								-o $main_dir/mags/alternate_primary/${refiner}_raw \
# 								-t $threads \
# 								-E $main_dir/mags/alternate_primary/error.txt \
# 								-D $main_dir/mags/alternate_primary/bin_stats.csv \
# 								-l $main_dir/mags/alternate_primary/logs"
# 					# ${refiner}_wrapper.sh -e \
# 					# 	-H $main_dir/mags/all/hmm \
# 					# 	-a $assembly \
# 					# 	-i $maps \
# 					# 	$contig_info_input \
# 					# 	-d $main_dir/mags/all/depth.tsv \
# 					# 	-b ${refiner}.raw \
# 					# 	-p $prefix.$assembler.alternate_primary \
# 					# 	-o $main_dir/mags/alternate_primary/${refiner}_raw \
# 					# 	-t $threads \
# 					# 	-E $main_dir/mags/alternate_primary/error.txt \
# 					# 	-D $main_dir/mags/alternate_primary/bin_stats.csv \
# 					# 	-l $main_dir/mags/alternate_primary/logs
# 				fi
# 			done
# 		fi
# 	done < $assembly_directories
# 	while read line
# 	do
# 		assembler=`echo $line | cut -d' ' -f1`
# 		main_dir=`echo $line | cut -d' ' -f2`
# 		dir_name=`basename $main_dir`
# 		assembler_parts=`echo $assembly_haps | sed "s|alternate_primary||g"`
# 		if [ "$assembler" != "metaflye" ] && [ "$assembler" != "flye" ] && [ "$assembler" != "meta-mdbg" ]
# 		then
# 			for refiner in `echo $refining_programs | sed "s|,| |g"`
# 			do
# 				counter=0
# 				while ! test -f $main_dir/mags/alternate_primary/${refiner}_raw/binstats_complete && [ $counter -le 1440 ] && ! test -f $main_dir/mags/alternate_primary/error.txt && ! `grep -q "$prefix.$assembler.alternate_primary"$'\t'$refiner".raw"$'\t'$refiner $main_dir/mags/alternate_primary/error.txt`
# 				do
# 					if [ $counter -eq 0 ]
# 					then
# 						date=`date`
# 						echo "$date : Waiting for completion of $dir_name - alternate_primary - $refiner"
# 					fi
# 					sleep 1m
# 					counter=`expr $counter + 1`
# 				done
# 			done
# 			if test -f $main_dir/mags/alternate_primary/bin_stats.csv
# 			then
# 				if ! test -f $outdir/bin_stats.csv
# 				then
# 					cp $main_dir/mags/alternate_primary/bin_stats.csv $outdir/bin_stats.csv
# 				else
# 					tail -n +2 $main_dir/mags/alternate_primary/bin_stats.csv >> $outdir/bin_stats.csv
# 				fi
# 			fi
# 		fi
# 	done < $assembly_directories
# fi
#
# sed -i "s|archaea,bacteria|archaea;bacteria|g" $outdir/bin_stats.csv
# sed -i "s|archaea,euk|archaea;euk|g" $outdir/bin_stats.csv
# sed -i "s|bacteria,archaea|bacteria;archaea|g" $outdir/bin_stats.csv
# sed -i "s|bacteria,euk|bacteria;euk|g" $outdir/bin_stats.csv
# sed -i "s|euk,archaea|euk;archaea|g" $outdir/bin_stats.csv
# sed -i "s|euk,bacteria|euk;bacteria|g" $outdir/bin_stats.csv
#
#
# score_mag_outputs.py -p $outdir/final $outdir/bin_stats.csv
#
# rm -rf $outdir/*tmp*
# touch $outdir/complete
#