#!/usr/bin/env bash

# MAGSCOT WRAPPER #
usage="SYNOPSIS: magscot_wrapper.sh [OPTIONS]\n\
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
	-H | --hmm_dir\tDirectory containing prodigal predictions (if doesn't exist it will be created)[\$outdir/hmm]\n\
	-E | --error_out\t\tFile to output error information [\$outdir/error.tsv]\n\
	-D | --data_out\t\tFile to output bin statistics [\$outdir/bin_info.tsv]\n\
	-C | --checkm_db\tCheckM database (see https://github.com/Ecogenomics/CheckM/wiki/Installation) [\$CHECKM_DATA_PATH]\n\
	-G | --gtdbtk_db\tgGTDB-TK database (see https://ecogenomics.github.io/GTDBTk/installing/index.html) [\$GTDBTK_DB]\n\
	-h | --help\t\tPrint help message\n\n\
AUTHOR:\n\
Noah Gettle 2023"

outdir="./magscot"
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
				input_maps=`echo $1 | sed "s|,| |g"`
				for map in $input_maps
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
		-H | --hmm_dir)
			if  ! `beginswith "-" $2`
			then
				shift
				hmm_dir=$1
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
mkdir -p $outdir
if [ "$logout" == "" ]
then
	logout=$outdir/log
fi

if [ "$hmm_dir" == "" ]
then
	hmm_dir=$outdir/hmm
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
	echo -e "#### ERROR --- MAGSCOT WRAPPER: Assembly file required. ####"
	if test -d $logout && test -f $logout/mag_pipe_progress.log
	then
		echo -e "#### ERROR --- MAGSCOT WRAPPER ($binning_program): Assembly file not provided. ####" >> $logout/mag_pipe_progress.log
		echo -e "$prefix\t$binning_program\tmagscot" >> $ERROR_OUT
	fi
	exit 1
fi

if [ "$input_maps" == "" ]
then
	echo -e "#### ERROR --- MAGSCOT WRAPPER: No input mapping files provided. ####"
	if test -d $logout && test -f $logout/mag_pipe_progress.log
	then
		echo -e "#### ERROR --- MAGSCOT WRAPPER ($binning_program): No input mapping files provided. ####" >> $logout/mag_pipe_progress.log
		echo -e "$prefix\t$binning_program\tmagscot" >> $ERROR_OUT
	fi
	exit 1
fi

mkdir -p $logout

#################################################################################
#### HMM Files ####
###################
mkdir -p $hmm_dir
if ! test -f $hmm_dir/input.hmm || [ `head $hmm_dir/input.hmm | wc -l` -eq 0 ]
then
	if ! test -f $hmm_dir/prodigal.faa
	then
		echo "####  MAGSCOT WRAPPER ($binning_program): Running PROGIGAL ####" >> $logout/mag_pipe_progress.log
		echo "####  MAGSCOT WRAPPER ($binning_program): Running PRODIGAL ####"
		### ORF detection with prodigal
		singularity exec -B /lustre,/nfs \
			$LOCAL_IMAGES/prodigal.sif \
				prodigal \
					-p meta \
					-a $hmm_dir/prodigal.faa \
					-d $hmm_dir/prodigal.ffn \
					ß-o $hmm_dir/tmpfile
	fi
	if ! test -f $hmm_dir/prodigal.faa || [ `head $hmm_dir/prodigal.faa | wc -l` -eq 0 ]
	then
		echo -e "$prefix\t$binning_program\tmagscot" >> $ERROR_OUT
		echo "#### ERROR --- MAGSCOT WRAPPER ($binning_program): Prodigal produced no output" >> $logout/mag_pipe_progress.log
		echo "#### ERROR --- MAGSCOT WRAPPER ($binning_program): Prodigal produced no output"
		exit 1
	fi
	if ! test -f $hmm_dir/hmm.tigr.out
	then
		echo "####  MAGSCOT WRAPPER ($binning_program): Running HMMSEARCH ####" >> $logout/mag_pipe_progress.log
		echo "####  MAGSCOT WRAPPER ($binning_program): Running HMMSEARCH ####"
		### annotation of protein sequences using HMMer and GTDBtk r207 marker genes
		singularity exec -B /lustre,/nfs \
			$LOCAL_IMAGES/hmmer.sif \
				hmmsearch \
					-o $hmm_dir/hmm.tigr.out \
					--tblout $hmm_dir/hmm.tigr.hit.out \
					--noali --notextw --cut_nc --cpu $threads \
					$GTDBTK_DB/hmm/gtdbtk_rel207_tigrfam.hmm \
					$hmm_dir/prodigal.faa
		singularity exec -B /lustre,/nfs \
			$LOCAL_IMAGES/hmmer.sif \
				hmmsearch \
					-o $hmm_dir/hmm.pfam.out \
					--tblout $hmm_dir/hmm.pfam.hit.out \
					--noali --notextw --cut_nc --cpu $threads \
					$GTDBTK_DB/hmm/gtdbtk_rel207_Pfam-A.hmm \
					$hmm_dir/prodigal.faa
	fi
	cat $hmm_dir/hmm.tigr.hit.out \
		| grep -v "^#" \
		| awk '{print $1"\t"$3"\t"$5}' \
		> $hmm_dir/tigr
	cat $hmm_dir/hmm.pfam.hit.out \
		| grep -v "^#" \
		| awk '{print $1"\t"$4"\t"$5}' \
		> $hmm_dir/pfam
	cat $hmm_dir/pfam $hmm_dir/tigr > $hmm_dir/input.hmm
	if ! test -f $hmm_dir/input.hmm || [ `head $hmm_dir/input.hmm | wc -l` -eq 0 ]
	then
		echo -e "$prefix\t$binning_program\tmagscot" >> $ERROR_OUT
		echo "#### ERROR --- MAGSCOT WRAPPER ($binning_program): HMMSearch produced no output" >> $logout/mag_pipe_progress.log
		echo "#### ERROR --- MAGSCOT WRAPPER ($binning_program): HMMSearch produced no output"
		exit 1
	fi
fi
# cat $hmm_dir/hmm.tigr.hit.out \
# 	| grep -v "^#" \
# 	| awk '{print $1"\t"$3"\t"$5}' \
# 	> $hmm_dir/tigr
# cat $hmm_dir/hmm.pfam.hit.out \
# 	| grep -v "^#" \
# 	| awk '{print $1"\t"$4"\t"$5}' \
# 	> $hmm_dir/pfam
# cat $hmm_dir/pfam $hmm_dir/tigr > $hmm_dir/input.hmm



# Add the x to ensure that unique bin names
# input_maps=`echo $maps | sed "s|,| |g"`
rm -rf $outdir/contig_to_bin.tsv
x=1
for map in $input_maps
do
	program=`dirname $map | awk -F'/' '{print $NF}'`
	if [ "$program" == "dastool" ] || [ "$program" == "magscot" ]
	then
		program=`dirname $map | awk -F'/' '{print $(NF-1)}'`
	fi
	awk -v program=$program -v n=$x '{print $2"."n"\t"$1"\t"program"."n}' $map \
		>> $outdir/contig_to_bin.tsv
	x=`expr $x + 1`
done
if ! test -d $outdir/output_bins || [ `ls $outdir/output_bins/*.fa | wc -l` -eq 0 ]
then
	echo "####  MAGSCOT WRAPPER ($binning_program): Running MAGSCOT ####" >> $logout/mag_pipe_progress.log
	echo "####  MAGSCOT WRAPPER ($binning_program): Running MAGSCOT ####"
	cd $outdir
	singularity run \
		--bind /lustre:/lustre \
		$LOCAL_IMAGES/magscot.sif MagScoT.R \
			-i $outdir/contig_to_bin.tsv \
			--hmm $hmm_dir/input.hmm
	if ! test -f $outdir/MAGScoT.refined.contig_to_bin.out || [ `head $outdir/MAGScoT.refined.contig_to_bin.out | wc -l` -eq 0 ]
	then
		echo -e "$prefix\t$binning_program\tmagscot" >> $ERROR_OUT
		echo "#### ERROR --- MAGSCOT WRAPPER ($binning_program): No output produced by magscot" >> $logout/mag_pipe_progress.log
		echo "#### ERROR --- MAGSCOT WRAPPER ($binning_program): No output produced by magscot"
		exit 1
	fi
	mkdir -p $outdir/output_bins
	contigs2bintofasta.sh \
		<(awk 'BEGIN {OFS="\t"} {print $2,$1}' $outdir/MAGScoT.refined.contig_to_bin.out) \
		$assembly \
		$outdir/output_bins/
fi
contigs2bin.sh \
	$outdir/output_bins \
	$outdir/contigs2bin.tsv \
	fa
touch $outdir/magscot_complete

if [ "$eval" == "TRUE" ]
then
	if ! test -f $outdir/binstats_complete
	then
		echo "#### MAGSCOT WRAPPER ($binning_program): Running Binstat Wrappers for MAGSCOT refined bins ####" >> $logout/mag_pipe_progress.log
		echo "#### MAGSCOT WRAPPER ($binning_program): Running Binstat Wrappers for MAGSCOT refined bins ####"
		# bsub -n$threads -q long -R"span[hosts=1]" \
			# -o $logout/$binning_program.magscot.binstats.%J.o \
			# -e $logout/$binning_program.magscot.binstats.%J.e \
			# -M100000 \
			# -R 'select[mem>100000] rusage[mem=100000]' \
			# 	"binstats.sh -G \
			# 		-a $assembly \
			# 		-o $outdir \
			# 		-i $outdir/output_bins \
			# 		-c $contig_info \
			# 		-d $depth \
			# 		-e fa \
			# 		-b $binning_program \
			# 		-R magscot \
			# 		-p $prefix \
			# 		-t $threads \
			# 		-E $error_out \
			# 		-D $data_out \
			# 		-C $checkm_db \
			# 		-G $gtdbtk_db \
			# 		-l $logout"
		binstats.sh -G \
							-a $assembly \
							-o $outdir \
							-i $outdir/output_bins \
							-c $contig_info \
							-d $depth \
							-e fa \
							-b $binning_program \
							-R magscot \
							-p $prefix \
							-t $threads \
							-E $error_out \
							-D $data_out \
							-C $checkm_db \
							-G $gtdbtk_db \
							-l $logout
	fi
fi
echo "#### MAGSCOT WRAPPER ($binning_program): Complete ####" >> $logout/mag_pipe_progress.log
echo "#### MAGSCOT WRAPPER ($binning_program): Complete ####"
touch $outdir/magscot_complete
