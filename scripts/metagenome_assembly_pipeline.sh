#!/usr/bin/env bash


usage="SYNOPSIS: metagenome_assembly_pipeline.sh [OPTIONS]\n\
OPTIONS:\n\
	-a | --assembly\t\tAssembly (can be gzipped) [REQUIRED]\n\
	-c | --contig_info\t\tContig information in format \"CONTIG_NAME\tLENGTH\tCOVERAGE\tIS_CIRCULAR[YN]\"\n\
	-p | --prefix\t\tFile prefix [assembly basename]\n\
	-P | --programs\t\tBinning programs [metabat2,metator,bin3c,maxbin2]\n\
		\t\tNote: If more than one program listed outputs will be combined using DASTool\n\
		\t\t\tOptions...\n\
		\t\t\t\tmetabat2\n\
		\t\t\t\tmetator\n\
		\t\t\t\tbin3c\n\
		\t\t\t\tmaxbin2\n\
	-R | --refining_programs\t\tComma separated list of bin refinement programs [dastool,magscot].\n\
		\t\t\tOptions:\n\
		\t\t\t\tnone\n\
		\t\t\t\tdastool\n\
		\t\t\t\tmagscot\n\
	-s | --combine_refined\t\tCombine and refine refined output from each binning program in additon to raw output\n\
	-r | --reads\t\tLong reads [Required]\n\
	-H | --hic\t\tHiC reads \n\
	-o | --outdir\t\tOutput directory [./metagenome_pipeline]\n\
	-t | --threadst\tNumber of threads per job [1]\n\
	-E | --error_out\t\tFile to output error information [\$outdir/error.tsv]\n\
	-D | --data_out\t\tFile to output bin statistics [\$outdir/bin_info.tsv]\n\
	-C | --checkm_db\tCheckM database (see https://github.com/Ecogenomics/CheckM/wiki/Installation) [$CHECKM_DATA_PATH]\n\
	-G | --gtdbtk_db\tgGTDB-TK database (see https://ecogenomics.github.io/GTDBTk/installing/index.html) [$GTDBTK_DB]\n\
	-f | --force\t\tOverwrite previous data\n\
	-h | --help\t\tPrint help message\n\n\
AUTHOR:\n\
Noah Gettle 2023"

outdir="./metagenome_pipeline"
combine_refined="FALSE"
threads=1
while [ "$1" != "" ]
do
	case $1 in
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
		-p | --prefix)
			if ! `beginswith "-" $2`
			then
				shift
				prefix=$1
			fi;;
		-P | --programs)
			if ! `beginswith "-" $2`
			then
				shift
				programs=$1
			fi;;
		-r | --reads)
			if ! `beginswith "-" $2`
			then
				shift
				reads=$1
				if ! test -f $reads
				then
					echo -e "Could not find long read file $reads.\n"
					exit 1
				fi
			fi;;
		-H | --hic)
			if ! `beginswith "-" $2`
			then
				shift
				hic=$1
				for h in `echo $hic | sed "s|,| |g"`
				do
					if ! test -f $h
					then
						echo -e "Could not find hic read file $h.\n"
						exit 1
					fi
				done
			fi;;
		-o | --outdir)
			if ! `beginswith "-" $2`
			then
				shift
				outdir=$1
			fi;;
		-R | --refining_programs)
			if ! `beginswith "-" $2`
			then
				shift
				refining_programs=`echo $1 | sed "s|,| |g"`
				for refiner in $refining_programs
				do
					if [ "$refiner" != "dastool" ] && [ "$refiner" != "magscot" ]
					then
						echo -e "ERROR: $refiner is not a valid value for -R.\n\n$usage"
						exit 1
					fi
				done
			fi;;
		-s | --combine_refined)
			combine_refined="TRUE";;
		-t | --threads)
			if ! `beginswith "-" $2`
			then
				shift
				threads=$1
			fi;;
		-D | --data_out)
			if ! `beginswith "-" $2`
			then
				shift
				data_out=$1
				export 	DATA_OUT=$data_out
			fi;;
		-E | --error_out)
			if ! `beginswith "-" $2`
			then
				shift
				error_out=$1
				export 	ERROR_OUT=$error_out
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
		-h | --help)
			echo -e $usage
			exit 0;;
		*)
			echo -e "Option $1 not valid\n$usage"
			exit 1
	esac
	shift
done

###################################################################
##### Sanity Check ######
#########################

if [ "$assembly" == "" ]
then
	echo -e "Assembly file required.\n$usage"
	exit 1
fi

if [ "$prefix" == "" ]
then
	prefix=`basename $assembly .gz`
	prefix=`basename $prefix .fa`
	prefix=`basename $prefix .fasta`
fi

if [ "$reads" == "" ]
then
	echo -e "Long read file required.\n$usage"
	exit 1
fi

if [ "$programs" == "" ]
then
	if [ "$hic" != "" ]
	then
		programs="metabat2,metator,bin3c,maxbin2"
	else
		programs="metabat2,maxbin2"
	fi
fi
echo "Running with binning programs $programs"

if (`echo $programs | grep -q "bin3c"` || `echo $programs | grep -q "metator"`) && [ "$hic" == "" ]
then
	echo -e "HiC data is required to run bin3C and Metator but was not provided"
	exit 1
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

if [ "$GTDBTK_DB" == "" ]
then
	export GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release207_v2
fi

###################################################################
#####  ######
#########################
logout=$outdir/logs
mkdir -p $logout
rm -rf $logout/*
rm -rf $logout/mag_pipe.log
rm -rf $logout/mag_pipe_progress.log

# Preprocessing
depth=""
if `echo $programs | grep -q "metabat2"` || `echo $programs | grep -q "maxbin2"`
then
	if ! test -f $outdir/depth.tsv && [ "$depth" == "" ]
	then
		if ! test -f $outdir/$prefix.bam
		then
			echo "#### MAIN: Mapping reads to assembly ####" >> $logout/mag_pipe_progress.log
			echo "#### MAIN: Mapping reads to assembly ####"
			if ! test -f $assembly.fai
			then
				samtools faidx $assembly
			fi
			minimap2 \
				-a -t $threads \
				-x map-hifi \
				$assembly $reads \
				| samtools view -T $assembly -b - \
				| samtools sort -T $assembly -@ $threads - \
				> $outdir/$prefix.bam 
		fi
		echo "#### MAIN: Calculating contig read depth ####" >> $logout/mag_pipe_progress.log
		echo "#### MAIN: Calculating contig read depth ####"
		singularity exec -B /lustre,/nfs \
			$LOCAL_IMAGES/metabat.sif \
				jgi_summarize_bam_contig_depths \
					--outputDepth $outdir/depth.tsv \
					$outdir/$prefix.bam
	fi
fi

if ! test -f $outdir/depth.tsv && [ "$depth" == "" ]
then
	echo "#### WARNING -- MAIN: Creation of depth file failed and none was provided ####" >> $logout/mag_pipe_progress.log
elif test -f $outdir/depth.tsv && [ "$depth" == "" ]
then
	depth=$outdir/depth.tsv
fi

rm -rf $data_out
# rm -rf $error_out
touch $error_out
sed -i "/magscot/d" $error_out

rm -rf $outdir/bin_stats.csv
programs=`echo $programs | sed "s|,| |g"`
mkdir -p $outdir/logs
for program in $programs
do
	if test -d $outdir/$program/raw && ! test -d $outdir/$program/output_bins && test -d $outdir/$program/raw/gtdb/genomes
	then
		cp -R $outdir/$program/raw/gtdb/genomes $outdir/$program/output_bins
	fi
	if test -d $outdir/$program/dastool && ! test -d $outdir/$program/dastool/output_bins && test -d $outdir/$program/dastool/gtdb/genomes
	then
		cp -R $outdir/$program/dastool/gtdb/genomes $outdir/$program/dastool/output_bins
	fi
	rm -rf $outdir/$program/raw/binstats_complete
	# rm -rf $outdir/$program/dastool/binstats_complete
	# rm -rf $outdir/$program/dastool/gtdb/output/complete
	# rm -rf $outdir/$program/raw/prokka
	# rm -rf $outdir/$program/dastool/prokka
	# rm -rf $outdir/$program/raw/drep
	# rm -rf $outdir/$program/dastool/drep
	if (! test -f $outdir/$program/${program}_complete || \
		! test -f $outdir/$program/dastool/dastool_complete || \
		! test -f $outdir/$program/dastool/binstats_complete) && \
		! `grep -q $'^'$prefix$'\t'$program$'$' $error_out` && \
		! `grep -q $'^'$prefix$'\t'$program$'\tdastool$' $error_out`
	then
		echo "#### RUNNING WRAPPER FOR $program ####" >> $logout/mag_pipe_progress.log
		echo "#### RUNNING WRAPPER FOR $program ####"
		bsub -n$threads -q long -R"span[hosts=1]" \
			-o $outdir/logs/${program}.%J.o \
			-e $outdir/logs/${program}.%J.e \
			-M300000 \
			-R 'select[mem>300000] rusage[mem=300000]' \
				"${program}_wrapper.sh -R \
					-a $assembly \
					-c $contig_info \
					-p $prefix \
					-r $reads \
					-d $depth \
					-H $hic \
					-o $outdir/$program \
					-t $threads \
					-E $error_out \
					-D $data_out \
					-l $logout"
		# ${program}_wrapper.sh -R \
		# 					-a $assembly \
		# 					-c $contig_info \
		# 					-p $prefix \
		# 					-r $reads \
		# 					-d $depth \
		# 					-H $hic \
		# 					-o $outdir/$program \
		# 					-t $threads \
		# 					-E $error_out \
		# 					-D $data_out \
		# 					-l $logout
	fi
done

binning_program=`echo $programs | sed "s| |_|g"`
if ! `echo $refining_programs | grep -q 'none'`
then
	for refiner in $refining_programs
	do
		echo -e "#### RUNNING $refiner WRAPPER for combined raw outputs ####"
		if test -d $outdir/${refiner}_raw/gtdb/genomes
		then
			cp -R $outdir/${refiner}_raw/gtdb/genomes $outdir/${refiner}_raw/output_bins
		fi
		maps=""
		# rm -rf $outdir/${refiner}_raw/prokka
		# rm -rf $outdir/${refiner}_raw/drep
		rm -rf $outdir/${refiner}_raw/binstats_complete
		# rm -rf $outdir/${refiner}_raw/gtdb/output/complete
		# rm -rf $outdir/${refiner}_raw/${refiner}_complete
		if (! test -d $outdir/${refiner}_raw || ! test -f $outdir/${refiner}_raw/binstats_complete) &&  ! `grep -q $'^'$prefix$'\t.*.raw\t'$refiner $ERROR_OUT`
		then
			maps=""
			if ! test -f $outdir/${refiner}_raw/${refiner}_complete || ! test -f $outdir/${refiner}_raw/binstats_compete
			then
				for program in $programs
				do
					counter=0
					while ! test -f $outdir/$program/${program}_complete && ! `grep -q $'^'$prefix$'\t'$program$'$' $ERROR_OUT` 
					do
						date=`date`
						if [ $counter -eq 0 ]
						then
							echo "$date : MAIN: Waiting for completion of $program ####" >> $logout/mag_pipe_progress.log
							echo "$date : MAIN: Waiting for completion of $program ####"
						fi
						sleep 1m
						counter=`expr $counter + 1`
					done
					if ! test -f $outdir/$program/contigs2bin.tsv || [ `head $outdir/$program/contigs2bin.tsv | wc -l` -eq 0 ]
					then
						echo -e "#### WARNING -- MAIN: Binning with $program either errored out or failed to produce the expected outputs. ####\
							\n####\tProceeding with DASTOOL refinement without these results. ####" >> $logout/mag_pipe_progress.log
						programs=`echo $programs | sed "s|$program||g"`
						binning_program=`echo $binning_program \
							| sed "s|${program}_||g" | sed "s|_$program||g"`
					else
						if [ "$maps" == "" ]
						then
							maps=$outdir/$program/contigs2bin.tsv
						else
							maps="$maps,$outdir/$program/contigs2bin.tsv"
						fi
					fi
				done
			fi
			if [ "$maps" != "" ]
			then
				echo "#### MAIN: Running $refiner refinement on raw binned results ####" >> $logout/mag_pipe_progress.log
				echo "#### MAIN: Running $refiner refinement on raw binned results ####"
				mkdir -p $outdir/${refiner}_raw
				bsub -n$threads -q long -R"span[hosts=1]" \
					-o $outdir/logs/${refiner}_raw.%J.o \
					-e $outdir/logs/${refiner}_raw.%J.e \
					-M100000 \
					-R 'select[mem>100000] rusage[mem=100000]' \
						"${refiner}_wrapper.sh -e \
							-H $outdir/hmm \
							-a $assembly \
							-i $maps \
							-c $contig_info \
							-d $depth \
							-b $binning_program.raw \
							-p $prefix \
				 			-o $outdir/${refiner}_raw \
							-t $threads \
							-E $error_out \
							-D $data_out \
							-C $checkm_db \
							-G $gtdbtk_db \
							-l $logout"
				# ${refiner}_wrapper.sh -e \
				# 							-H $outdir/hmm \
				# 							-a $assembly \
				# 							-i $maps \
				# 							-c $contig_info \
				# 							-d $depth \
				# 							-b $binning_program.raw \
				# 							-p $prefix \
				# 				 			-o $outdir/${refiner}_raw \
				# 							-t $threads \
				# 							-E $error_out \
				# 							-D $data_out \
				# 							-C $checkm_db \
				# 							-G $gtdbtk_db \
				# 							-l $logout
			else
				echo -e "$prefix\t$refiner.raw" >> $ERROR_OUT
			fi
		fi
	done
fi

binning_program=`echo $programs | sed "s| |_|g"`
if [ "$combine_refined" == "TRUE" ]
then
	for refiner in $refining_programs
	do
		if test -d $outdir/${refiner}_filtered/gtdb/genomes
		then
			cp -R $outdir/${refiner}_filtered/gtdb/genomes $outdir/${refiner}_filtered/output_bins
		fi
		maps=""
		# rm -rf $outdir/${refiner}_filtered/prokka
		# rm -rf $outdir/${refiner}_filtered/drep
		rm -rf $outdir/${refiner}_filtered/binstats_complete
		# rm -rf $outdir/${refiner}_filtered/gtdb/output/complete
		# rm -rf $outdir/${refiner}_filtered/${refiner}_complete
		if (! test -d $outdir/${refiner}_filtered || ! test -f $outdir/${refiner}_filtered/binstats_complete) &&  ! `grep -q $'^'$prefix$'\t.*.filtered\t'$refiner $ERROR_OUT`
		then
			maps=""
			if ! test -f $outdir/${refiner}_filtered/${refiner}_complete || ! test -f $outdir/${refiner}_filtered/binstats_compete
			then
				for program in $programs
				do
					counter=0
					while ! test -f $outdir/$program/dastool/dastool_complete && ! `grep -q $'^'$prefix$'\t'$program$'$' $ERROR_OUT` && ! `grep -q $'^'$prefix$'\t'$program$'\tdastool$' $ERROR_OUT`
					do
						date=`date`
						if [ $counter -eq 0 ]
						then
							echo "$date : MAIN: Waiting for completion of refinement of $program outputs  ####" >> $logout/mag_pipe_progress.log
							echo "$date : MAIN: Waiting for completion of refinement of $program outputs  ####"
						fi
						sleep 1m
						counter=`expr $counter + 1`
					done
					if ! test -f $outdir/$program/dastool/contigs2bin.tsv || [ `head $outdir/$program/dastool/contigs2bin.tsv | wc -l` -eq 0 ]
					then
						echo -e "#### WARNING: DASTOOL refinement of binned results from $program either errored out or failed to produce the expected outputs. ####\
							\n####\tProceeding with DASTOOL refinement without these results. ####" >> $logout/mag_pipe_progress.log
						programs=`echo $programs | sed "s|$program||g"`
						binning_program=`echo $binning_program \
							| sed "s|${program}_||g" | sed "s|_$program||g"`
					else
						if [ "$maps" == "" ]
						then
							maps=$outdir/$program/dastool/contigs2bin.tsv
						else
							maps="$maps,$outdir/$program/dastool/contigs2bin.tsv"
						fi
					fi
				done
			fi
			if [ "$maps" != "" ]
			then
				# if [ "$refiner" == "magscot" ]
				# then
				# 	while ! test -f $outdir/${refiner}_raw/${refiner}_complete
				# 	do
				# 		date=`date`
				# 		echo "$date : MAIN: Waiting for ${refiner} refinement of raw bins to complete  ####" >> $logout/mag_pipe_progress.log
				# 		echo "$date : MAIN: Waiting for ${refiner} refinement of raw bins to complete  ####"
				# 		sleep 5m
				# 	done
				# fi
				echo "#### MAIN: Running ${refiner} refinement on refined binned results ####" >> $logout/mag_pipe_progress.log
				echo "#### MAIN: Running ${refiner} refinement on refined binned results ####"
				mkdir -p $outdir/${refiner}_filtered
				bsub -n$threads -q long -R"span[hosts=1]" \
					-o $outdir/logs/${refiner}_filtered.%J.o \
					-e $outdir/logs/${refiner}_filtered.%J.e \
					-M100000 \
					-R 'select[mem>100000] rusage[mem=100000]' \
						"${refiner}_wrapper.sh -e \
							-H $outdir/hmm \
							-a $assembly \
							-i $maps \
							-c $contig_info \
							-d $depth \
							-b $binning_program.filtered \
							-p $prefix \
				 			-o $outdir/${refiner}_filtered \
							-t $threads \
							-E $error_out \
							-D $data_out \
							-C $checkm_db \
							-G $gtdbtk_db \
							-l $logout"
				# ${refiner}_wrapper.sh -e \
				# 							-H $outdir/hmm \
				# 							-a $assembly \
				# 							-i $maps \
				# 							-c $contig_info \
				# 							-d $depth \
				# 							-b $binning_program.filtered \
				# 							-p $prefix \
				# 				 			-o $outdir/${refiner}_filtered \
				# 							-t $threads \
				# 							-E $error_out \
				# 							-D $data_out \
				# 							-C $checkm_db \
				# 							-G $gtdbtk_db \
				# 							-l $logout
			else
				echo -e "$prefix\t$refiner.filtered" >> $ERROR_OUT
			fi
		fi
	done
fi

# Wait for all programs to complete #
for program in $programs
do
	count=0
	while ! test -f $outdir/$program/raw/binstats_complete && ! test -f $outdir/$program/dastool/binstats_complete && ! `grep -q $'^'$prefix$'\t'$program$'$' $ERROR_OUT` && ! `grep -q $'^'prefix\t$program\traw\t*$'$' $ERROR_OUT` 
	do
		if [ $count -eq 0 ]
		then
			date=`date`
			echo -e "$date : MAIN: Waiting on $program. ####"
		fi
		sleep 1m
		count=`expr $count + 1`
	done
done

if ! `echo $refining_programs | grep -q 'none'`
then
	for refiner in $refining_programs
	do
		count=0
		while ! test -f $outdir/${refiner}_raw/binstats_complete && ! `grep -q $'^'$prefix$'\t.*.raw\t'$refiner $ERROR_OUT` && ! `grep -q $'^'$prefix$'\t'$refiner$'.raw$' $ERROR_OUT` 
		do
			if [ $count -eq 0 ]
			then
				date=`date`
				echo -e "$date : MAIN: Waiting on ${refiner}_raw. ####"
			fi
			sleep 1m
			count=`expr $count + 1`
		done
	done
fi

if [ "$combine_refined" == "TRUE" ]
then
	for refiner in $refining_programs
	do
		count=0
		while ! test -f $outdir/${refiner}_filtered/binstats_complete && ! `grep -q $'^'$prefix$'\t.*.filtered\t'$refiner $ERROR_OUT` && ! `grep -q $'^'$prefix$'\t'$refiner$'.filtered$' $ERROR_OUT`
		do
			if [ $count -eq 0 ]
			then
				date=`date`
				echo -e "$date : MAIN: Waiting on ${refiner}_filtered. ####"
			fi
			sleep 1m
			count=`expr $count + 1`
		done
	done
fi

touch $outdir/complete
