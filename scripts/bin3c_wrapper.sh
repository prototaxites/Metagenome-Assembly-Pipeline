#!/usr/bin/env bash

# BIN3C WRAPPER #
usage="SYNOPSIS: bin3c_wrapper.sh [OPTIONS]\n\
OPTIONS:\n\
	-a | --assembly\tAssembly (can be gzipped) [REQUIRED]\n\
	-c | --contig_info\tTab-delimited contig information in format...\n\
			   \t\t\t\"CONTIG_NAME LENGTH COVERAGE IS_CIRCULAR[YN]\"\n\
	-p | --prefix\t\tFile prefix [assembly basename]\n\
	-r | --reads\t\tLong reads \n\
	-H | --hic\t\tHiC reads [REQUIRED]\n\
	-d | --depth\t\tTab-delimited contig depth file in format...\n\
			   \t\t\t\"contigName contigLen totalAvgDepth {PREFIX}.bam {PREFIX}.bam-var\"\n\
	-o | --outdir\t\tOutput directory [./bin3c]\n\
	-e | --eval\t\tEvaluate output bins\n\
	-R | --refine\t\tRefine bins using DASTOOL\n\
	-t | --threads\t\tNumber of threads per job [1]\n\
	-E | --error_out\tFile to output error information [\$outdir/error.tsv]\n\
	-D | --data_out\tFile to output bin statistics [\$outdir/bin_info.tsv]\n\
	-C | --checkm_db\tCheckM database [\$CHECKM_DATA_PATH]\n\
				\t\t\t(see https://github.com/Ecogenomics/CheckM/wiki/Installation)\n\
	-G | --gtdbtk_db\tgGTDB-TK database [\$GTDBTK_DB]\n\
				\t\t\t(see https://ecogenomics.github.io/GTDBTk/installing/index.html)\n\
	-h | --help\t\tPrint help message\n\n\
AUTHOR:\n\
Noah Gettle 2023"

outdir="./bin3c"
threads=1
eval=FALSE
refine=FALSE
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
		-H | --hic)
			if ! `beginswith "-" $2`
			then
				shift
				hic=$1
			fi;;
		-o | --outdir)
			if ! `beginswith "-" $2`
			then
				shift
				outdir=$1
			fi;;
		-e | --eval)
			eval=TRUE;;
		-R | --refine)
			refine=TRUE;;
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


if [ "$prefix" == "" ]
then
	prefix=`basename $assembly .gz`
	prefix=`basename $prefix .fa`
	prefix=`basename $prefix .fasta`
fi

if [ "$DATA_OUT" == "" ]
then
	export DATA_OUT=$outdir/bin_info.tsv
fi

if [ "$ERROR_OUT" == "" ]
then
	export ERROR_OUT=$outdir/error.tsv
fi

if [ "$assembly" == "" ]
then
	echo -e "#### ERROR --- BIN3C WRAPPER: Assembly file required. ####\n$usage"
	if test -d $logout && test -f $logout/mag_pipe_progress.log
	then
		echo -e "#### ERROR --- BIN3C WRAPPER: Assembly file not provided. ####" >> $logout/mag_pipe_progress.log
		echo -e "$prefix\tbin3c" >> $ERROR_OUT
	fi
	exit 1
fi


if [ "$CHECKM_DATA_PATH" == "" ]
then
	export CHECKM_DATA_PATH=/lustre/scratch124/tol/projects/darwin/users/ng13/checkm_db
fi
if [ "$GTDBTK_DB" == "" ]
then
	export GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release207_v2
fi

if [ "$hic" == "" ]
then
	echo -e "#### ERROR --- BIN3C_WRAPPER: HiC reads required. ####\n$usage"
	if test -d $logout && test -f $logout/mag_pipe_progress.log
	then
		echo -e "#### ERROR --- BIN3C WRAPPER: HiC reads not provided. ####" >> $logout/mag_pipe_progress.log
		echo -e "$prefix\tbin3c" >> $ERROR_OUT
	fi
	exit 1
fi
for h in `echo $hic | sed "s|,| |g"`
do
	if ! test -f $h
	then
		echo -e "#### ERROR --- BIN3C_WRAPPER: Could not find hic read file $h. ####\n$usage"
		if test -d $logout && test -f $logout/mag_pipe_progress.log
		then
			echo -e "#### ERROR --- BIN3C WRAPPER: Could not find hic read file $h. ####" >> $logout/mag_pipe_progress.log
			echo -e "$prefix\tbin3c" >> $ERROR_OUT
		fi
		exit 1
	fi
done
hic1=`echo $hic | cut -d',' -f1`
hic2=`echo $hic | cut -d',' -f2`

###################################################################
mkdir -p $logout

# Bin3C
mkdir -p $outdir
rm -rf $outdir/bin3c_complete
cd $outdir

if ! test -d $outdir/output_bins || [ `ls $outdir/output_bins/*.fna | wc -l` == 0 ]
then
	echo "#### BIN3C WRAPPER: Indexing assembly ####" >> $logout/mag_pipe_progress.log
	echo "#### BIN3C WRAPPER: Indexing assembly ####"
	singularity exec -B /lustre,/nfs \
		/software/tola/images/bwa-mem2-2.2.1.sif \
			bwa-mem2 index $assembly
	if ! test -f $outdir/hic2ctg.bam
	then
		echo "#### BIN3C WRAPPER: Mapping HiC reads to assembly ####" >> $logout/mag_pipe_progress.log
		echo "#### BIN3C WRAPPER: Mapping HiC reads to assembly ####"
		singularity exec -B /lustre,/nfs \
			/software/tola/images/bwa-mem2-2.2.1.sif \
				bwa-mem2 mem \
					-t 16 -5SP \
					$assembly \
					$hic1 \
					$hic2 \
					| samtools view -bS - \
					| samtools sort -@ $threads -n - \
					> $outdir/hic2ctg.bam
	fi
	if ! test -f $outdir/bin3c_out/contact_map.p.gz
	then
		echo "#### BIN3C WRAPPER: Creating contact map ####" >> $logout/mag_pipe_progress.log
		echo "#### BIN3C WRAPPER: Creating contact map ####"
		rm -rf $outdir/bin3c_out
		singularity run --bind /lustre/:/lustre/ \
			$LOCAL_IMAGES/bin3c.sif \
				bin3C mkmap -v \
					-e DpnII -e HinfI -e DdeI -e MseI \
					$assembly \
					$outdir/hic2ctg.bam \
					$outdir/bin3c_out
	fi

	if ! test -d $outdir/bin3c_clust/fasta || [ `ls $outdir/bin3c_clust/fasta/*.fna | wc -l` == 0 ]
	then
		echo "#### BIN3C WRAPPER: Binning assembly ####" >> $logout/mag_pipe_progress.log
		echo "#### BIN3C WRAPPER: Binning assembly ####"
		rm -rf $outdir/bin3c_clust
		singularity run --bind /lustre/:/lustre/ \
			$LOCAL_IMAGES/bin3c.sif \
				bin3C cluster -v \
					--fasta $assembly \
					$outdir/bin3c_out/contact_map.p.gz \
					$outdir/bin3c_clust
		if [ `ls $outdir/bin3c_clust/fasta/*.fna | wc -l` == 0 ]
		then
			echo -e "$prefix\tbin3c" >> $ERROR_OUT
			echo -e "#### ERROR --- BIN3C WRAPPER: No output produced by bin3c ####" >>  $logout/mag_pipe_progress.log
			echo -e "#### ERROR --- BIN3C WRAPPER: No output produced by bin3c ####"
			exit 1
		fi
	fi
	mkdir -p $outdir/output_bins
	for fasta in $outdir/bin3c_clust/fasta/*.fna
	do
		bin=`basename $fasta .fna`
		rm -rf $outdir/output_bins/$bin.fna
		awk -F'[ :]' '{if($1 ~ /^>/){print ">"$3} else {print}}' $fasta \
			>> $outdir/output_bins/$bin.fna
	done
else
	echo "#### BIN3C WRAPPER: Previously binned fastas found in $outdir/raw/gtdb/genomes ####" >> $logout/mag_pipe_progress.log
	echo "#### BIN3C WRAPPER: Previously binned fastas found in $outdir/raw/gtdb/genomes ####"
fi
contigs2bin.sh \
	$outdir/output_bins \
	$outdir/contigs2bin.tsv \
	fna
touch $outdir/bin3c_complete

###################################################################

extension=fna
binning_program=bin3c
refinement_tool=raw

if [ "$eval" == "TRUE" ]
then
	if ! test -d $outdir/raw || ! test -f $outdir/raw/binstats_complete
	then
		echo "#### BIN3C WRAPPER: Running Binstat Wrapper for binned assemblies ####" >> $logout/mag_pipe_progress.log
		echo "#### BIN3C WRAPPER: Running Binstat Wrapper for binned assemblies ####"
		bsub -n$threads -q long -R"span[hosts=1]" \
			-o $logout/$binning_program.$refinement_tool.binstats.%J.o \
			-e $logout/$binning_program.$refinement_tool.binstats.%J.e \
			-M100000 \
			-R 'select[mem>100000] rusage[mem=100000]' \
				"binstats.sh -G \
					-a $assembly \
					-o $outdir \
					-i $outdir/output_bins \
					-c $contig_info \
					-d $depth \
					-e $extension \
					-b $binning_program \
					-R $refinement_tool \
					-p $prefix \
					-t $threads \
					-E $error_out \
					-D $data_out \
					-C $checkm_db \
					-G $gtdbtk_db \
					-l $logout"
	fi
fi
#################################################################################
#### DAS_tools ####
###################
if [ "$refine" == "TRUE" ]
then
	mkdir -p $outdir/dastool
	if ! test -f $outdir/dastool/dastool_complete || ! test -f $outdir/dastool/binstats_complete
	then
		echo "#### BIN3C WRAPPER: Running DASTOOL wrapper on binned assemblies ####" >> $logout/mag_pipe_progress.log
		echo "#### BIN3C WRAPPER: Running DASTOOL wrapper on binned assemblies ####"
		# bsub -n$threads -q long -R"span[hosts=1]" \
		# 	-o $logout/$binning_program.dastool.%J.o \
		# 	-e $logout/$binning_program.dastool.%J.e \
		# 	-M50000 \
		# 	-R 'select[mem>50000] rusage[mem=50000]' \
		# 		"dastool_wrapper.sh -e \
		# 			-a $assembly \
		# 			-i $outdir/contigs2bin.tsv \
		# 			-c $contig_info \
		# 			-d $depth \
		# 			-b $binning_program \
		# 			-p $prefix \
		#  			-o $outdir/dastool \
		# 			-t $threads \
		# 			-E $error_out \
		# 			-D $data_out \
		# 			-C $checkm_db \
		# 			-G $gtdbtk_db \
		# 			-l $logout"
		dastool_wrapper.sh -e \
							-a $assembly \
							-i $outdir/contigs2bin.tsv \
							-c $contig_info \
							-d $depth \
							-b $binning_program \
							-p $prefix \
				 			-o $outdir/dastool \
							-t $threads \
							-E $error_out \
							-D $data_out \
							-C $checkm_db \
							-G $gtdbtk_db \
							-l $logout
	fi
fi
echo "#### BIN3C WRAPPER: Complete ####" >> $logout/mag_pipe_progress.log
echo "#### BIN3C WRAPPER: Complete ####"
touch $outdir/bin3c_complete