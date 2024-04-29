#!/usr/bin/env bash

# BINSTATS #
usage="SYNOPSIS: binstats.sh [OPTIONS]\n\
OPTIONS:\n\
	-a | --assembly\t\tAssembly (can be gzipped) [REQUIRED]\n\
	-i | --input_bins\t\tDirectory containing bins in fasta format [\$outdir/gtdb/genomes]\n\
	-c | --contig_info\t\tContig information in format:\n\
						\t\t\"CONTIG_NAME\tLENGTH\tCOVERAGE\tIS_CIRCULAR[YN]\"\n\
	-d | --depth\t\tTab-delimited contig depth file in format...\n\
			   \t\t\t\"contigName contigLen totalAvgDepth {PREFIX}.bam {PREFIX}.bam-var\"\n\
	-e | --extension\t\tBin fasta file extension [fa]\n\
	-b | --binning_program\t\tBinning program [NA]\n\
	-R | --refinement_tool\t\tTool used to refine bins [NA]\n\
	-p | --prefix\t\t\tFile prefix [NA]\n\
	-G | --gtdbtk\t\t\tRun GTDB-TK [FALSE]\n\
	-o | --outdir\t\t\tOutput directory [./binstats]\n\
	-t | --threads\t\t\tNumber of threads per job [1]\n\
	-E | --error_out\t\tFile to output error information [\$outdir/error.tsv]\n\
	-D | --data_out\t\tFile to output bin statistics [\$outdir/bin_info.tsv]\n\
	-C | --checkm_db\t\tCheckM database (see https://github.com/Ecogenomics/CheckM/wiki/Installation) [\$CHECKM_DATA_PATH]\n\
	-G | --gtdbtk_db\t\tgGTDB-TK database (see https://ecogenomics.github.io/GTDBTk/installing/index.html) [\$GTDBTK_DB]\n\
	-h | --help\t\t\tPrint help message\n\n\
AUTHOR:\n\
Noah Gettle 2023"

threads=1
outdir="./binstats"
extension=fa
refinement_tool=NA
binning_program=NA
prefix=NA
gtdbtk=FALSE
while [ "$1" != "" ]
do
	case $1 in
		-i | --input_bins)
			if ! `beginswith "-" $2`
			then
				shift
				input_bins=$1
				if ! test -d $input_bins
				then
					echo -e "Could not find bin directory $input_bins.\n"
					exit 1
				fi
			fi;;
		-a | --assembly)
			if ! `beginswith "-" $2`
			then
				shift
				assembly=$1
				if ! test -f $assembly
				then
					echo -e "Could not find assembly file $assembly.\n"
					exit 1
				fi
			fi;;
		-c | --contig_info)
			if ! `beginswith "-" $2`
			then
				shift
				contig_info=$1
				if ! test -f $contig_info
				then
					echo -e "Could not find contig information file $contig_info.\n"
					exit 1
				fi
			fi;;
		-d | --depth)
			if ! `beginswith "-" $2`
			then
				shift
				depth=$1
				if ! test -f $depth
				then
					echo -e "Could not find contig depth file $depth.\n"
					exit 1
				fi
			fi;;
		-G | --gtdbtk)
			gtdbtk=TRUE;;
		-e | --extension)
			if ! `beginswith "-" $2`
			then
				shift
				extension=$1
			fi;;
		-b | --binning_program)
			if ! `beginswith "-" $2`
			then
				shift
				binning_program=$1
			fi;;
		-R | --refinement_tool)
			if ! `beginswith "-" $2`
			then
				shift
				refinement_tool=$1
			fi;;
		-p | --prefix)
			if ! `beginswith "-" $2`
			then
				shift
				prefix=$1
			fi;;
		-o | --outdir)
			if ! `beginswith "-" $2`
			then
				shift
				outdir=$1
			fi;;
		-t | --threads)
			if ! `beginswith "-" $2`
			then
				shift
				threads=$1
			fi;;
		-E | --error_out)
			if ! `beginswith "-" $2`
			then
				shift
				error_out=$1
				export 	ERROR_OUT=$error_out
			fi;;	
		-D | --data_out)
			if ! `beginswith "-" $2`
			then
				shift
				data_out=$1
				export 	DATA_OUT=$data_out
			fi;;
		-C | --checkm_db)
			if ! `beginswith "-" $2`
			then
				shift
				checkm_db=$1
				if ! test -d $checkm_db
				then
					echo -e "Could not find CheckM database directory $checkm_db.\n"
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
					echo -e "Could not find GTDB-Tk database directory $cgtdbtk_db.\n"
					exit 1
				else
					export GTDBTK_DB=$gtdbtk_db
				fi
			fi;;
		-l | --logout)
			if ! `beginswith "-" $2`
			then
				shift
				logout=$1
			fi;;
		-h | --help)
			echo -e $usage
			exit 0;;
		*)
			echo -e "Option $1 not valid\m$usage"
			exit 1
	esac
	shift
done

###################################################################
##### Sanity Check ######
#########################

if [ "$logout" == "" ]
then
	logout=$outdir/log
fi

if [ "$DATA_OUT" == "" ]
then
	export DATA_OUT=$outdir/bin_info.tsv
fi

if [ "$ERROR_OUT" == "" ]
then
	export ERROR_OUT=$outdir/error.tsv
fi

if [ "$CHECKM_DATA_PATH" == "" ]
then
	export CHECKM_DATA_PATH=/lustre/scratch124/tol/projects/darwin/users/ng13/checkm_db
fi

if [ "$assembly" == "" ]
then
	echo -e "#### ERROR --- BINSTATS WRAPPER: Assembly file required. ####"
	if test -d $logout && test -f $logout/mag_pipe_progress.log
	then
		echo -e "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): Assembly file not provided. ####" >> $logout/mag_pipe_progress.log
		echo -e "$prefix\t$binning_program\t$refinement_tool\tbinstats" >> $ERROR_OUT
	fi
	exit 1
fi
if [ "$input_bins" == "" ]
then
	input_bins=$outdir/output_bins
fi

if [ `ls $input_bins/*.$extension | wc -l` -eq 0 ]
then
	echo "#### ERROR --- BINSTATS WRAPPER: No bins with extension $extension found in $input_bins. ####"
	if test -d $logout && test -f $logout/mag_pipe_progress.log
	then
		echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): No bins with extension $extension found in $input_bins. ####" >>  $logout/mag_pipe_progress.log
		echo -e "$prefix\t$binning_program\t$refinement_tool\tbinstats" >> $ERROR_OUT
	fi
	exit 1
fi


rm -rf $outdir/binstats_complete

###################################################################
##### GTDB-TK ######
####################
if [ "$gtdbtk" == "TRUE" ]
then
	mkdir -p $outdir/gtdb/output
	# rm -rf $outdir/gtdb/output/complete
	if ! test -f $outdir/gtdb/output/complete
	then
		rm -rf $logout/gtdb-$binning_program-$refinement_tool.*
		if ! `ls $outdir/gtdb/output | grep -q $'gtdbtk.bac.*.summary.tsv'` && ! `ls $outdir/gtdb/output | grep -q $'gtdbtk.ar.*.summary.tsv'`
		then
			bsub -n$threads -q long -R"span[hosts=1]" \
				-o $logout/gtdb-$binning_program-$refinement_tool.o \
				-e $logout/gtdb-$binning_program-$refinement_tool.e \
				-M100000 \
				-R 'select[mem>100000] rusage[mem=100000]' \
					"gtdbtk_wrapper.sh \
						-d $input_bins \
						-x $extension \
						-o $outdir/gtdb/output"
		else
			gtdbtk_wrapper.sh \
				-d $input_bins \
				-x $extension \
				-o $outdir/gtdb/output
		fi
	fi
fi
###################################################################
##### CheckM ######
###################
if ! test -f $outdir/checkm/results.tsv || [ `cat $outdir/checkm/results.tsv | wc -l` -le 1 ]
then
	echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Generating bin quality metrics with CheckM ####" >> $logout/mag_pipe_progress.log
	echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Generating bin quality metrics with CheckM ####"
	checkm lineage_wf \
		-t $threads -x $extension \
		--tab_table \
		-f $outdir/checkm/results.tsv \
		$input_bins \
		$outdir/checkm
fi
if ! test -f $outdir/checkm/results.tsv || [ `cat $outdir/checkm/results.tsv | wc -l` -eq 1 ]
then
	echo -e "$prefix\t$binning_program\t$refinement_tool\tcheckm" >> $ERROR_OUT
	echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): CheckM produced no output. ####" >> $logout/mag_pipe_progress.log
	echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): CheckM produced no output. ####"
	exit 1
fi

###################################################################
##### CheckM SSU ######
#######################

if ! test -f $outdir/checkm_ssu/ssu_summary.tsv || [ `cat $outdir/checkm_ssu/ssu_summary.tsv | wc -l` -le 1 ]
then
	echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Identifying SSUs in binned assemblies with CheckM ####" >> $logout/mag_pipe_progress.log
	echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Identifying SSUs in binned assemblies with CheckM ####"
	checkm ssu_finder \
		-t $threads -x $extension \
		$assembly \
		$input_bins \
		$outdir/checkm_ssu
fi
if ! test -f $outdir/checkm_ssu/ssu_summary.tsv || [ `cat $outdir/checkm_ssu/ssu_summary.tsv | wc -l` -le 1 ]
then
	echo -e "$prefix\t$binning_program\t$refinement_tool\tcheckm_ssu" >> $ERROR_OUT
	echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): CheckM SSU Finder produced no outputs. ####" >> $logout/mag_pipe_progress.log
	echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): CheckM SSU Finder produced no outputs. ####"
	exit 1
fi

if ! test -f $outdir/checkm_ssu/ssu_bin_summary.tsv || [ `tail -n +2 $outdir/checkm_ssu/ssu_bin_summary.tsv | wc -l` -ne `ls $input_bins/*.$extension | wc -l` ]
then
	echo -e "bin_id\tssu_count\thmms\tssu_coords" > $outdir/checkm_ssu/ssu_bin_summary.tsv
	for file in $input_bins/*.$extension
	do
		bin=`basename $file .$extension`
		count=0
		seqs="NA"
		hmms="NA"
		grep $'^'$bin$'\t' $outdir/checkm_ssu/ssu_summary.tsv \
			> $outdir/checkm_ssu/tmp.tsv
		while read line
		do
			count=`expr $count + 1`
			tig=`echo $line | cut -d' ' -f2`
			hmm=`echo $line | cut -d' ' -f3`
			start=`echo $line | cut -d' ' -f5`
			end=`echo $line | cut -d' ' -f6`
			comp=`echo $line | cut -d' ' -f7`
			strand="+"
			if [ "$comp" == "False" ]
			then
				strand="-"
			fi
			coords="$tig:$start-$end:$strand"
			if [ "$seqs" == "NA" ]
			then
				seqs="$coords"
				hmms=$hmm
			else
				seqs="$seqs;$coords"
				if ! `echo $hmms | grep -q $hmm`
				then
					hmms="$hmms;$hmm"
				fi
			fi 
		done < $outdir/checkm_ssu/tmp.tsv
		echo -e "$bin\t$count\t$hmms\t$seqs" >> $outdir/checkm_ssu/ssu_bin_summary.tsv
	done
fi

###################################################################
##### PROKKA ######
###################
# rm -rf $outdir/prokka
mkdir -p $outdir/prokka
if ! test -f $outdir/prokka/annot_info.tsv || [ `tail -n +2 $outdir/prokka/annot_info.tsv | wc -l` -ne `ls $input_bins/*.$extension | wc -l` ]
then
	echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): PROKKA bin annotation ####" >> $logout/mag_pipe_progress.log
	echo "#### BINSTATS WRAPPER ($binning_program): PROKKA bin annotation ####"
	rm -rf $outdir/prokka/annot_info.tsv
	echo -e "bin\ttotal_trnas\tunique_trnas\trrna_23s\trrna_16s\trrna_5s" > $outdir/prokka/annot_info.tsv
	for file in $input_bins/*.$extension
	do
		bin=`basename $file .$extension`
		contamination=`grep $bin$'\t' $outdir/checkm/results.tsv \
			| cut -f13`
		# domain=`grep $bin$'\t' $outdir/gtdb/output/gtdbtk.summary.tsv \
		# 	| cut -f2 \
		# 	| cut -d';' -f1 \
		# 	| sed "s|d__||g" | \
		# 	sed "s|Unclassified_||g"`
		clade=`grep $bin$'\t' $outdir/checkm/results.tsv | cut -f2`
		if [ "$clade" != 'Root' ]
		then
			if `echo $clade | grep -q 'Archaea'`
			then
				domain="Archaea"
			else
				domain="Bacteria"
			fi
		elif `grep -q $'^'$bin$'\t' $outdir/checkm_ssu/ssu_summary.tsv`
		then
			domain=`grep $'^'$bin$'\t' $outdir/checkm_ssu/ssu_summary.tsv | head -1 | cut -f3`
		else
			domain="Root"
		fi
		if ([ "$domain" == "Archaea" ] || [ "$domain" == "Bacteria" ]) && [ `echo "$contamination <= 10.0" | bc -l` -eq 1 ]
		then
			rm -rf $outdir/prokka/temp
			mkdir -p $outdir/prokka/temp
			cd $outdir/prokka/temp
			prokka \
				--rfam \
				--cpus 0 \
				--noanno \
				--rawproduct \
				--kingdom $domain \
				$file
			total_trna=`grep "tRNA" $outdir/prokka/temp/PROKKA_*/PROKKA_*.txt | cut -d' ' -f2`
			if [ "$total_trna" == "" ]
			then
				total_trna=0
			fi
			unique_trna=`grep "tRNA" $outdir/prokka/temp/PROKKA_*/PROKKA_*.tsv \
				| awk '{print $NF}' | awk -F'[-(]' '{print $2}' | sort -u | wc -l`
			rrna_23s="N"
			rrna_16s="N"
			rrna_5s="N"
			if `grep -q "23S ribosomal RNA" $outdir/prokka/temp/PROKKA_*/PROKKA_*.tsv`
			then
				rrna_23s="Y"
			fi
			if `grep -q "16S ribosomal RNA" $outdir/prokka/temp/PROKKA_*/PROKKA_*.tsv`
			then
				rrna_16s="Y"
			fi
			if `grep -q "5S ribosomal RNA" $outdir/prokka/temp/PROKKA_*/PROKKA_*.tsv`
			then
				rrna_5s="Y"
			fi
		else
			total_trna="NA"
			unique_trna="NA"
			rrna_23s="NA"
			rrna_16s="NA"
			rrna_5s="NA"
		fi
		echo -e "$bin\t$total_trna\t$unique_trna\t$rrna_23s\t$rrna_16s\t$rrna_5s" \
			>> $outdir/prokka/annot_info.tsv
	done
fi
if ! test -f $outdir/prokka/annot_info.tsv || [ `tail -n +2 $outdir/prokka/annot_info.tsv | wc -l` -ne `ls $input_bins/*.$extension | wc -l` ]
then
	echo -e "$prefix\t$binning_program\t$refinement_tool\tprokka" >> $ERROR_OUT
	echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): Prokka produced an error or no output ####" >> $logout/mag_pipe_progress.log
	echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): Prokka produced an error or no output ####"
	exit 1
fi



###################################################################
##### dRep ######
#################
mkdir -p $outdir/drep/prefilter
cd $outdir
head -1 $outdir/checkm/results.tsv \
	| sed "s|Bin Id|genome|g" | sed "s|C|c|g" | sed "s|\t|,|g" \
	| cut -d',' -f1,12,13 >  $outdir/drep/checkm.csv
awk -F'\t' -v ext=$extension '{if($13 <= 10){OFS=","; $1 = $1"."ext; print $1,$12,$13}}' $outdir/checkm/results.tsv \
	>> $outdir/drep/checkm.csv

filtered_bins=`tail -n +2 $outdir/drep/checkm.csv \
	| cut -d',' -f1`
rm -rf $outdir/drep/prefilter/*
for bin in $filtered_bins
do
	# if [ `grep -v $'^>' $outdir/gtdb/genomes/$bin | wc -c` -ge 50000 ]
	# then
	cp $input_bins/$bin \
		$outdir/drep/prefilter/
	# fi
done


if [ `ls $outdir/drep/prefilter | wc -l` -gt 1 ] && [ `ls $outdir/drep/prefilter/ | grep $extension$'$' | wc -l` -gt 1 ]
then
	if (! test -d $outdir/drep/dereplicated_genomes || [ `ls $outdir/drep/prefilter/ | grep $extension$'$' | wc -l` -eq 0 ]) && (! test -f $ERROR_OUT || ! `grep -q "$prefix\t$binning_program\t$refinement_tool\tdrep" $ERROR_OUT`)
	then
		echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Dereplicating bins with DREP ####" >> $logout/mag_pipe_progress.log
		echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Dereplicating bins with DREP ####"
		singularity exec -B /lustre,/nfs \
			$LOCAL_IMAGES/drep.sif \
				dRep dereplicate drep \
					-p 30  \
					-g $outdir/drep/prefilter/*.$extension \
					--genomeInfo $outdir/drep/checkm.csv
		if ! test -d $outdir/drep/dereplicated_genomes || [ `ls $outdir/drep/dereplicated_genomes | wc -l` -eq 0 ] || [ `ls $outdir/drep/prefilter/ | grep $extension$'$' | wc -l` -eq 0 ]
		then
			echo -e "$prefix\t$binning_program\t$refinement_tool\tdrep" >> $ERROR_OUT
			echo "#### WARNING --- BINSTATS WRAPPER ($binning_program-$refinement_tool): dRep produced no output ####" >> $logout/mag_pipe_progress.log
			echo "#### WARNING --- BINSTATS WRAPPER ($binning_program-$refinement_tool): dRep produced no output ####"
			mv $outdir/drep/prefilter/*.$extension $outdir/drep/dereplicated_genomes
			# exit 1
		fi
	fi
fi

###################################################################
##### Combine Data Output ######
################################
if [ "$gtdbtk" == "TRUE" ]
then
	# Wait for GTDB-TK to complete
	n=0
	while ! test -f $outdir/gtdb/output/complete && (! test -f $logout/gtdb-$binning_program-$refinement_tool.o || ! `grep -q 'ERROR' $logout/gtdb-$binning_program-$refinement_tool.o`) && [ $n -le 288 ]
	do
		if [ $n -eq 0 ]
		then
			date=`date`
			echo "$date : Waiting for completion of GTDB - $prefix-$binning_program-$refinement_tool"
		fi
		sleep 1m
		n=`expr $n + 1`
	done
	if ! test -f $outdir/gtdb/output/complete && ! test -f $outdir/gtdb/output/complete && ! `grep -q 'ERROR' $logout/gtdb-$binning_program-$refinement_tool`
	then
		echo -e "$prefix\t$binning_program\t$refinement_tool\tgtdbtk" >> $ERROR_OUT
		echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): GTDB-TK did not finish within 24 hrs. Exiting ####" >> $logout/mag_pipe_progress.log
		echo "#### ERROR --- BINSTATS WRAPPER ($binning_program-$refinement_tool): GTDB-TK did not finish within 24 hrs. Exiting ####"
		exit 1
	fi
	gtdbtk_file=$outdir/gtdb/output/gtdbtk.summary.csv
fi
# Collection
checkm_file=$outdir/checkm/results.tsv
ssu_file=$outdir/checkm_ssu/ssu_bin_summary.tsv
prokka_file=$outdir/prokka/annot_info.tsv
checkm_file_header=`head -1 $checkm_file | cut -f2- | sed "s|\t|,|g"`
ssu_file_header=`head -1 $ssu_file | cut -f2- | sed "s|\t|,|g"`
prokka_file_header=`head -1 $prokka_file | cut -f2- | sed "s|\t|,|g"`
gtdbtk_file_header="classification,fastani_reference,fastani_reference_radius,fastani_taxonomy,fastani_ani,fastani_af,closest_placement_reference,closest_placement_radius,closest_placement_taxonomy,closest_placement_ani,closest_placement_af,pplacer_taxonomy,classification_method,note,other_related_references(genome_id;species_name;radius;ANI;AF),msa_percent,translation_table,red_value,warnings,ncbi_classification,ncbi_taxon,taxon_id"

if ! test -f $DATA_OUT
then
	echo -e "name,binning_program,refining_program,bin_id,quality,size,contigs,circular,mean_coverage,$prokka_file_header,$checkm_file_header,drep,$ssu_file_header,$gtdbtk_file_header" \
		> $DATA_OUT
fi

echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Outputting bin statistics ####" >> $logout/mag_pipe_progress.log
echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Outputting bin statistics ####"
gtdbtk_info=`python -c 'print("NA," * 22)' | sed "s|,$||g"`
for bin in $input_bins/*.$extension
do
	bin_id=`basename $bin .$extension`
	# echo $bin
	drep="FAILED"
	if test -d $outdir/drep/dereplicated_genomes && `ls $outdir/drep/dereplicated_genomes | grep -q $bin_id.$extension`
	then
		drep="PASSED"
	fi
	checkm_info=`grep $bin_id$'\t' $checkm_file | cut -f2- | sed "s|\t|,|g"`
	ssu_info=`grep $bin_id$'\t' $ssu_file | cut -f2- | sed "s|\t|,|g"`
	prokka_info=`grep $bin_id$'\t' $prokka_file | cut -f2- | sed "s|\t|,|g"`
	if [ "$gtdbtk" == "TRUE" ] && test -f $gtdbtk_file && [ `head $gtdbtk_file | wc -l` -gt 1 ]
	then
		gtdbtk_info=`grep $bin_id',' $gtdbtk_file | cut -d',' -f2-`
	fi
	completeness=`echo $checkm_info | cut -d',' -f11`
	contamination=`echo $checkm_info | cut -d',' -f12`
	size=`grep -v $'^>' $bin | wc -c`
	contigs=`grep $'^>' $bin | cut -d'>' -f2`
	rm -rf $outdir/contig_info.tmp
	if [ "$depth" != "" ]
	then
		rm -rf $outdir/depth.tmp
		for contig in $contigs
		do
			grep $contig$'\t' $depth >> $outdir/depth.tmp
		done
		mean_coverage=`awk '{total += ($2 * $3); bps += $2} END{print total/bps}' $outdir/depth.tmp`
	else
		mean_coverage=NA
	fi
	for contig in $contigs
	do
		grep $contig$'\t' $contig_info >>  $outdir/contig_info.tmp
	done
	circular=`grep $'\tY$' $outdir/contig_info.tmp | wc -l`
	contig_count=`cat $outdir/contig_info.tmp | wc -l`
	quality=LOW
	if [ `echo "$contamination <= 10.0" | bc -l` -eq 1 ] && [ `echo "$completeness >= 50.0" | bc -l` -eq 1 ]
	then
		quality=MEDIUM
		if [ `echo "$contamination <= 5.0" | bc -l` -eq 1 ] && [ `echo $prokka_info | cut -d',' -f3` == "Y" ] && [ `echo $prokka_info | cut -d',' -f4` == "Y" ] && [ `echo $prokka_info | cut -d',' -f5` == "Y" ] && [ `echo $prokka_info | cut -d',' -f2` != "NA" ] && [ `echo $prokka_info | cut -d',' -f2` -ge 18 ] 
		then
			if ([ `echo "$completeness >= 50.0" | bc -l` -eq 1 ] && [ $contig_count == $circular ]) || [ `echo "$completeness >= 90.0" | bc -l` -eq 1 ]
			then
				quality=HIGH
			fi
		fi
	fi
	echo -e "$prefix,$binning_program,$refinement_tool,$bin_id,$quality,$size,$contig_count,$circular,$mean_coverage,$prokka_info,$checkm_info,$drep,$ssu_info,$gtdbtk_info" \
		>> $DATA_OUT
done
echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Complete ####" >> $logout/mag_pipe_progress.log
echo "#### BINSTATS WRAPPER ($binning_program-$refinement_tool): Complete ####"
touch $outdir/binstats_complete