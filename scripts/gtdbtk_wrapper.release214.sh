#!/usr/bin/env bash

GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release214
NCBI_TAXDUMP="$GTDBTK_DB/`ls $GTDBTK_DB | grep 'ncbi_taxdump' | sort -r | head -1`"
IMG_DIR=/software/team311/ng13/local/images
ALTERNATE_KEY=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/ncbi_altkey.csv
bac_prefix="bac120"
ar_prefix="ar53"


usage="SYNOPSIS: gtdbtk_wrapper.sh [OPTIONS]\n\
DESCRIPTION: Classifies bins using GTDB-TK and tranlates lineages to NCBI Taxonomy format.\n\
OPTIONS:\n\
	-d | --directory\tDirectory containing MAG binning output.\n\
	-x | --extension\tBin file extension [fa].\n\
	-g | --gtdb_output_file\tTSV or CSV with GTDB 'classification' column.\n\
	-G | --gtdb_output_dir\tPrevious GTDB output directory (if wishing to convert to NCBI Taxonomy)\n\
	-o | --outdir\t\tOutput directory [./gtdb/output]\n\
	-F | --force\t\tForce overwrite of previous GTDB results.\n\
	-h | --help\t\tPrint help message\n\n\
AUTHOR:\n\
Noah Gettle 2023"


outdir="./gtdb/output"
extension="fa"
force=FALSE
while [ "$1" != "" ]
do
	case $1 in
		-d | --directory)
			if ! `beginswith "-" $2`
			then
				shift
				directory=$1
				if ! test -d $directory
				then
					echo -e "ERROR: Could not find input MAG directory $directory.\n$usage"
					exit 1
				fi
			fi;;
		-g | --gtdb_output_file)
			if ! `beginswith "-" $2`
			then
				shift
				gtdb_output_file=$1
				if ! test -f $gtdb_output_file
				then
					echo -e "ERROR: Could not find suggested GTDB output file $gtdb_output_file.\n$usage"
					exit 1
				fi
			fi;;
		-G | --gtdb_output_dir)
			if ! `beginswith "-" $2`
			then
				shift
				gtdb_output_dir=$1
				if ! test -d $gtdb_output_dir
				then
					echo -e "ERROR: Could not find suggested GTDB output directory $gtdb_output_dir.\n$usage"
					exit 1
				fi
			fi;;
		-x | --extension)
			if ! `beginswith "-" $2`
			then
				shift
				extension=$1
			fi;;
		-o | --outdir)
			if ! `beginswith "-" $2`
			then
				shift
				outdir=$1
			fi;;
		-F | --force)
			force=TRUE;;
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
mkdir -p $outdir
main_dir=`dirname $outdir`

if ([ "$gtdb_output_dir" == "" ] && [ "$gtdb_output_file" == "" ]) || [ "$force" == "TRUE" ]
then
	if [ "$directory" == "" ]
	then
		if test -d $main_dir/genomes
		then
			directory=$main_dir/genomes
		else
			echo -e "ERROR: No input bin directory provided and could not find one.\n$usage"
			exit 1
		fi
	fi
	if [ `ls $directory | grep '.'$extension$'$' | wc -l` -eq 0 ]
	then
		echo -e "ERROR: Could not find any bins with extension .$extension in $directory.\n$usage"
		exit 1
	else
		mkdir -p $main_dir/genomes
		if [ "$directory" != "$main_dir/genomes" ]
		then
			cp $directory/*.$extension $main_dir/genomes/
		fi
	fi
elif [ "$gtdb_output_dir" != "" ] && [ "$gtdb_output_file" == "" ] && ! test -f $gtdb_output_dir/gtdbtk.$bac_prefix.summary.tsv && ! test -f$gtdb_output_dir/gtdbtk.$ar_prefix.summary.tsv 
then
	echo -e "ERROR: Could not find GTDB-TK output files in $gtdb_output_dir.\n$usage"
	exit 1
elif [ "$gtdb_output_file" != "" ] && ! `head -1 $gtdb_output_file | sed "s|[\t,]|\n|g" | grep -q $'^classification$'`
then
	echo -e "ERROR: Could not find 'classification' column in provided GTDB-TK output file $gtdb_output_file.\n$usage"
	exit 1
fi

######################################################################################
# DATABASE CHECK #
##################
# Check for required GTDB data
current_day=`date +'%Y%m%d'`

if ! test -f $GTDBTK_DB/ar53_metadata.tsv
then
	echo "#### DATABASE CHECK: Could not find ${ar_prefix}_metadata in $GTDBTK_DB. Downloading... ####"
	cd $GTDBTK_DB
	rm -rf $GTDBTK_DB/ar53_metadata*.tsv
	rm -rf $GTDBTK_DB/metadata.tsv
	wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz
	gunzip $GTDBTK_DB/ar53_metadata.tsv.gz
	name=`basename $GTDBTK_DB/ar53_metadata.tsv .tsv`
	mv $GTDBTK_DB/ar53_metadata.tsv $GTDBTK_DB/$name.$current_day.tsv
	ln -fs $GTDBTK_DB/$name.$current_day.tsv $GTDBTK_DB/ar53_metadata.tsv
fi
if ! test -f $GTDBTK_DB/bac120_metadata.tsv
then
	echo "#### DATABASE CHECK: Updating ${bac_prefix}_metadata in $GTDBTK_DB. Downloading... ####"
	cd $GTDBTK_DB
	rm -rf $GTDBTK_DB/bac120_metadata*.tsv
	rm -rf $GTDBTK_DB/metadata.tsv
	wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz
	gunzip $GTDBTK_DB/bac120_metadata.tsv.gz
	name=`basename $GTDBTK_DB/bac120_metadata.tsv .tsv`
	mv $GTDBTK_DB/bac120_metadata.tsv $GTDBTK_DB/$name.$current_day.tsv
	ln -fs $GTDBTK_DB/$name.$current_day.tsv $GTDBTK_DB/bac120_metadata.tsv
fi
	

if ! test -f $GTDBTK_DB/metadata.tsv
then
	cat $GTDBTK_DB/bac120_metadata.tsv \
		$GTDBTK_DB/ar53_metadata.tsv \
		| cut -f17,79 > $GTDBTK_DB/metadata.tsv
fi

if [ "$NCBI_TAXDUMP" == "" ] || ! test -d $NCBI_TAXDUMP || [ `ls $GTDBTK_DB | grep 'ncbi_taxdump' | sort -r | head -1 | awk -F'[/.]' '{print $NF}' | cut -c1-6` -ne `date +'%Y%m'` ]
then
	if [ `ls $GTDBTK_DB | grep 'ncbi_taxdump' | sort -r | head -1 | awk -F'[/.]' '{print $NF}' | cut -c1-6` -ne `date +'%Y%m'` ]
	then
		echo "#### DATABASE CHECK: Updating NCBI taxdump... ####"
	else
		echo "#### DATABASE CHECK: Could not find NCBI taxdump. Downloading... ####"
	fi
	NCBI_TAXDUMP=$GTDBTK_DB/ncbi_taxdump.$current_day
	mkdir -p $NCBI_TAXDUMP
	cd $NCBI_TAXDUMP
	wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 
	tar -zxvf $NCBI_TAXDUMP/taxdump.tar.gz 
	rm -rf $NCBI_TAXDUMP/taxdump.tar.gz
fi


###################################################################
##### GTDB-TK ######
####################

if ([ "$gtdb_output_dir" == "" ] && [ "$gtdb_output_file" == "" ]) || [ "$force" == "TRUE" ]
then
	# HMMER chokes on large datasets. Filter out potential hazards
	echo -e "#### GTDB-TK: Checking bin sizes. ####"
	for file in $main_dir/genomes/*.$extension
	do
		if [ `grep -v $'^>' $file | wc -c` -ge 30000000 ]
		then
			mkdir -p $main_dir/genomes.large
			mv $file $main_dir/genomes.large/
		fi
	done
	mkdir -p $main_dir/output
	if [ `ls $main_dir/genomes | wc -l ` -eq 0 ]
	then
		echo "#### ERROR: No output fastas left after filtering out bins > 30Mbp. ####"
		touch $main_dir/output/complete
		exit 1
	fi
	if (! test -f $main_dir/output/gtdbtk.$bac_prefix.summary.tsv && ! test -f $outdir/gtdb/output/gtdbtk.$ar_prefix.summary.tsv) || [ "$force" == "TRUE" ]
	then
		rm -rf $main_dir/output/*
		echo -e "#### GTDB-TK: Starting GTDB-TK. ####"
		singularity exec \
			-B $main_dir:/data \
			-B $GTDBTK_DB:/refdata \
			$IMG_DIR/gtdbtk.sif \
				gtdbtk classify_wf \
					--genome_dir /data/genomes \
					-x $extension \
					--out_dir /data/output \
					--mash_db /refdata/mash_out
	fi
	if ! test -f $main_dir/output/gtdbtk.$bac_prefix.summary.tsv && ! test -f $main_dir/output/gtdbtk.$ar_prefix.summary.tsv
	then
		echo "#### ERROR: GTDB_TK produced no output. ####"
		touch $main_dir/output/complete
		exit 1
	fi
	gtdb_output_dir=$main_dir/output
fi

###################################################################

if [ "$gtdb_output_file" == "" ]
then
	rm -rf $gtdb_output_dir/gtdbtk.summary.tsv
	for file in $gtdb_output_dir/gtdbtk.*.summary.tsv
	do
		if ! test -f $gtdb_output_dir/gtdbtk.summary.tsv
		then
			cp $file $gtdb_output_dir/gtdbtk.summary.tsv
		else
			tail -n +2 $file >> $gtdb_output_dir/gtdbtk.summary.tsv
		fi
	done
	gtdb_output_file=$gtdb_output_dir/gtdbtk.summary.tsv
fi

###################################################################

sed "s|,|;|g" $gtdb_output_file > $outdir/gtdbtk.summary.csv
sed -i "s|\t|,|g" $outdir/gtdbtk.summary.csv
sed -i "s|N\/A|NA|g" $outdir/gtdbtk.summary.csv


class_column=`head -1 $outdir/gtdbtk.summary.csv \
	| sed "s|,|\n|g" \
	| grep -n $'^classification$' \
	| cut -d':' -f1`

cut -d',' -f$class_column $outdir/gtdbtk.summary.csv \
	| tail -n +2 \
	> $outdir/gtdb_taxa.tmp.txt
gtdb2ncbi.py \
	-N $NCBI_TAXDUMP \
	-A $ALTERNATE_KEY \
	$outdir/gtdb_taxa.tmp.txt \
	$outdir/tmp.csv
cut -d',' -f2- $outdir/tmp.csv > $outdir/tmp2.csv
paste -d',' $outdir/gtdbtk.summary.csv $outdir/tmp2.csv \
	> $outdir/tmp.comb.csv
mv $outdir/tmp.comb.csv $outdir/gtdbtk.summary.csv

touch $outdir/complete

rm -rf $outdir/*tmp*

