#!/usr/bin/env bash

# BIN3C WRAPPER #
usage="SYNOPSIS: dastool_wrapper.sh [OPTIONS]\n\
OPTIONS:\n\
	-a | --assembly\t\tAssembly (can be gzipped) [REQUIRED]\n\
	-i | --input_maps\t\tFiles mapping contigs to bins [REQUIRED]\n\
	-c | --contig_info\t\tContig information in format \"CONTIG_NAME\tLENGTH\tCOVERAGE\tIS_CIRCULAR[YN]\"\n\
	-d | --depth\t\tTab-delimited contig depth file in format...\n\
			   \t\t\t\"contigName contigLen totalAvgDepth {PREFIX}.bam {PREFIX}.bam-var\"\n\
	-b | --binning_program\t\tBinning program [NA]\n\
	-p | --prefix\t\tFile prefix [NA]\n\
	-o | --outdir\t\tOutput directory [./bin3c]\n\
	-e | --eval\t\tEvaluate output bins\n\
	-t | --threadst\tNumber of threads per job [1]\n\
	-E | --error_out\t\tFile to output error information [\$outdir/error.tsv]\n\
	-D | --data_out\t\tFile to output bin statistics [\$outdir/bin_info.tsv]\n\
	-C | --checkm_db\tCheckM database (see https://github.com/Ecogenomics/CheckM/wiki/Installation) [\$CHECKM_DATA_PATH]\n\
	-G | --gtdbtk_db\tgGTDB-TK database (see https://ecogenomics.github.io/GTDBTk/installing/index.html) [\$GTDBTK_DB]\n\
	-h | --help\t\tPrint help message\n\n\
AUTHOR:\n\
Noah Gettle 2023"

outdir="./dastool"
threads=1
binning_program=NA
prefix=NA
eval=FALSE
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
		-i | --input_maps)
			if ! `beginswith "-" $2`
			then
				shift
				input_maps=$1
				for map in `echo $input_maps | sed "s|,| |g"`
				do
					if ! test -f $map
					then
						echo -e "Could not find contig-bin file $map.\n"
						exit 1
					fi
				done
			fi;;
		-b | --binning_program)
			if ! `beginswith "-" $2`
			then
				shift
				binning_program=$1
			fi;;
		-p | --prefix)
			if ! `beginswith "-" $2`
			then
				shift
				prefix=$1
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
		-o | --outdir)
			if ! `beginswith "-" $2`
			then
				shift
				outdir=$1
			fi;;
		-e | --eval)
			eval=TRUE;;
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
		-H | --hmm_dir)
			if  ! `beginswith "-" $2`
			then
				shift
				hmm_dir=$1
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
if [ "$GTDBTK_DB" == "" ]
then
	export GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release207_v2
fi


if [ "$assembly" == "" ]
then
	echo -e "#### ERROR --- DASTOOL WRAPPER: Assembly file required. ####"
	if test -d $logout && test -f $logout/mag_pipe_progress.log
	then
		echo -e "#### ERROR --- DASTOOL WRAPPER ($binning_program): Assembly file not provided. ####" >> $logout/mag_pipe_progress.log
		echo -e "$prefix\t$binning_program\tdastool" >> $ERROR_OUT
	fi
	exit 1
fi

if [ "$input_maps" == "" ]
then
	echo -e "#### ERROR --- DASTOOL WRAPPER: No input mapping files provided. ####"
	if test -d $logout && test -f $logout/mag_pipe_progress.log
	then
		echo -e "#### ERROR --- DASTOOL WRAPPER ($binning_program): No input mapping files provided. ####" >> $logout/mag_pipe_progress.log
		echo -e "$prefix\t$binning_program\tdastool" >> $ERROR_OUT
	fi
	exit 1
fi

mkdir -p $logout

#################################################################################
#### DAS_tools ####
###################
mkdir -p $outdir
rm -rf $outdir/dastool_complete
if ! test -d $outdir/output_bins || [ `ls $outdir/output_bins/*.fa | wc -l` -eq 0 ]
then
	echo "####  DASTOOL WRAPPER ($binning_program): Running DASTOOL ####" >> $logout/mag_pipe_progress.log
	echo "####  DASTOOL WRAPPER ($binning_program): Running DASTOOL ####"
	singularity exec --bind /lustre:/lustre \
		$LOCAL_IMAGES/das_tool.sif \
			DAS_Tool \
				-t $threads \
				-i $input_maps \
				-c $assembly \
				-o $outdir/dastool \
				--write_bins
	if ! test -d $outdir/dastool_DASTool_bins || [ `ls $outdir/dastool_DASTool_bins/*.fa | wc -l` -eq 0 ]
	then
		echo -e "$prefix\t$binning_program\tdastool" >> $ERROR_OUT
		echo "#### ERROR --- DASTOOL WRAPPER ($binning_program): No output produced by dastool" >> $logout/mag_pipe_progress.log
		echo "#### ERROR --- DASTOOL WRAPPER ($binning_program): No output produced by dastool"
		exit 1
	fi
	mkdir -p $outdir/output_bins
	mv $outdir/dastool_DASTool_bins/*.fa \
		$outdir/output_bins/
fi
contigs2bin.sh \
	$outdir/output_bins \
	$outdir/contigs2bin.tsv \
	fa

if [ "$eval" == "TRUE" ]
then
	if ! test -f $outdir/binstats_complete
	then
		echo "#### DASTOOL WRAPPER ($binning_program): Running Binstat Wrappers for DASTOOL refined bins ####" >> $logout/mag_pipe_progress.log
		echo "#### DASTOOL WRAPPER ($binning_program): Running Binstat Wrappers for DASTOOL refined bins ####"
		# bsub -n$threads -q long -R"span[hosts=1]" \
		# 	-o $logout/$binning_program.dastool.binstats.%J.o \
		# 	-e $logout/$binning_program.dastool.binstats.%J.e \
		# 	-M100000 \
		# 	-R 'select[mem>100000] rusage[mem=100000]' \
		# 		"binstats.sh -G \
		# 			-a $assembly \
		# 			-o $outdir \
		# 			-i $outdir/output_bins \
		# 			-c $contig_info \
		# 			-d $depth \
		# 			-e fa \
		# 			-b $binning_program \
		# 			-R dastool \
		# 			-p $prefix \
		# 			-t $threads \
		# 			-E $error_out \
		# 			-D $data_out \
		# 			-C $checkm_db \
		# 			-G $gtdbtk_db \
		# 			-l $logout"
		binstats.sh -G \
							-a $assembly \
							-o $outdir \
							-i $outdir/output_bins \
							-c $contig_info \
							-d $depth \
							-e fa \
							-b $binning_program \
							-R dastool \
							-p $prefix \
							-t $threads \
							-E $error_out \
							-D $data_out \
							-C $checkm_db \
							-G $gtdbtk_db \
							-l $logout
	fi
fi
echo "#### DASTOOL WRAPPER ($binning_program): Complete ####" >> $logout/mag_pipe_progress.log
echo "#### DASTOOL WRAPPER ($binning_program): Complete ####"
touch $outdir/dastool_complete