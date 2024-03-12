#!/usr/bin/env bash

GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release207_v2
NCBI_TAXDUMP="$GTDBTK_DB/`ls $GTDBTK_DB | grep 'ncbi_taxdump' | sort -r | head -1`"
IMG_DIR=/software/team311/ng13/local/images
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

if ! test -f $GTDBTK_DB/ar53_metadata_r*.tsv
then
	echo "#### DATABASE CHECK: Could not find ${ar_prefix}_metadata in $GTDBTK_DB. Downloading... ####"
	cd $GTDBTK_DB
	wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz 
	tar -zxvf $GTDBTK_DB/ar53_metadata.tar.gz
	rm -rf $GTDBTK_DB/ar53_metadata.tar.gz
fi
if ! test -f $GTDBTK_DB/bac120_metadata_r*.tsv
then
	echo "#### DATABASE CHECK: Could not find ${bac_prefix}_metadata in $GTDBTK_DB. Downloading... ####"
	wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz
	tar -zxvf $GTDBTK_DB/bac120_metadata.tar.gz
	rm -rf $GTDBTK_DB/bac120_metadata.tar.gz
fi
	

if ! test -f $GTDBTK_DB/metadata.tsv
then
	cat $GTDBTK_DB/${bac_prefix}_metadata_r*.tsv \
		$GTDBTK_DB/${ar_prefix}_metadata_r*.tsv \
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
		singularity run \
			-B $main_dir:/data \
			-B $GTDBTK_DB:/refdata \
			$IMG_DIR/gtdbtk.2.2.4.sif \
				classify_wf \
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
rm -rf $outdir/tmp.csv
n=0
while read line
do
	n=`expr $n + 1`
	echo $n
	if ! test -f $outdir/tmp.csv
	then
		echo "$line,ncbi_classification,ncbi_taxon,taxon_id" \
			> $outdir/tmp.csv
	else
		gtdb_lineage=`echo $line \
			| cut -d',' -f$class_column \
			| sed "s|p__Proteobacteria|p__Pseudomonadota|" \
			| sed "s|Desulfobacterota |Desulfobacterota_|"`
			# | sed "s|__|%|g" \
			# | sed "s|_| |g" \
			# | sed "s|%|__|g"`
		if [ "$gtdb_lineage" != "Unclassified" ] && [ "$gtdb_lineage" != "" ]
		then
			# Remove trailing clade spots w/o id
			while `echo $gtdb_lineage | grep -q $'__$'`
			do
				gtdb_lineage=`echo $gtdb_lineage | sed "s|;.__$||g"`
			done
			# Replace spaces in species name with _
			gtdb_lineage=`echo $gtdb_lineage | sed "s|Unclassified |d__|g"`
			# Should be unnecessary. Cut the classification until a match is found in $GTDBTK_DB/metadata.tsv.
			# while ! `grep -q $'^'"$gtdb_lineage"$'[;\t]' $GTDBTK_DB/metadata.tsv`
			# do
			# 	gtdb_lineage=`echo $gtdb_lineage | awk 'BEGIN{FS=OFS=";"}{NF--; print}'`
			# done
			# First see if there is direct match in NCBI database
			N=`echo $gtdb_lineage | awk -F';' '{print NF}'`
			ncbi_taxon=""
			ncbi_lineage=""
			if [ `echo "$gtdb_lineage" | cut -d';' -f$N | awk -F'__' '{print $1}'` != "s" ]
			then
				while [ $N -gt 1 ] && [ "$ncbi_taxon" == "" ]
				do
					test_taxon=`echo "$gtdb_lineage" | cut -d';' -f$N | awk -F'__' '{print $NF}' | sed "s|_A | |g"`
					# Look first in manually generated synonym file to see
					if test -f $GTDBTK_DB/ncbi_synonyms.tsv && `grep -q $'^'"$test_taxon"$'\t' $GTDBTK_DB/ncbi_synonyms.tsv`
					then
						ncbi_taxon=`grep $'^'"$test_taxon"$'\t' $GTDBTK_DB/ncbi_synonyms.tsv | cut -f2`
						ncbi_lineage=`grep $'^'"$test_taxon"$'\t' $GTDBTK_DB/ncbi_synonyms.tsv | cut -f3`
					else
						# Current taxon name not in NCBI taxdump move to higher level
						if ! `grep -q $'\t'"$test_taxon"$'\t' $NCBI_TAXDUMP/names.dmp` && ! `grep -q $'\tCandidatus '"$test_taxon"$'\t' $NCBI_TAXDUMP/names.dmp`
						then
							N=`expr $N - 1`
						else
							if `grep -q $'\tCandidatus '"$test_taxon"$'\t' $NCBI_TAXDUMP/names.dmp`
							then
								ncbi_taxon="Candidatus $test_taxon"
							else
								ncbi_taxon=$test_taxon
							fi
							ncbi_lineage=`cut -f2 $GTDBTK_DB/metadata.tsv \
								| grep "__$ncbi_taxon"$'[$;]' \
								| head -1`
							taxon_place=`echo $ncbi_lineage \
								| sed "s|;|\n|g" \
								| grep -n "$ncbi_taxon$" \
								| tail -1 \
								| cut -d':' -f1`
							ncbi_lineage=`echo $ncbi_lineage \
								| cut -d';' -f1-$taxon_place`
						fi
					fi
				done
			fi
			if [ "$ncbi_lineage" == "" ]
			then
				gtdb_taxon=`echo $gtdb_lineage | awk -F';' '{print $NF}' | awk -F'__' '{print $NF}'`
				field_number=`echo $gtdb_lineage | awk -F';' '{print NF}'`
				grep $'^'"$gtdb_lineage"$'[;\t]' $GTDBTK_DB/metadata.tsv \
					| sort -u | cut -f2 \
					| cut -d';' -f1-$field_number \
					> $outdir/tmp.ncbi.txt
				while `grep -q $'__$' $outdir/tmp.ncbi.txt`
				do
					sed -i "s|;.__$||" $outdir/tmp.ncbi.txt
				done
				if [ `awk -F';' -v N=$field_number '{if(NF == N){print}}' $outdir/tmp.ncbi.txt | wc -l` == 1 ]
				then
					ncbi_lineage=`awk -F';' -v N=$field_number '{if(NF == N){OFS=";"; print}}' $outdir/tmp.ncbi.txt`
				elif `grep -q "__$gtdb_taxon"$'$' $outdir/tmp.ncbi.txt`
				then
					ncbi_lineage=`grep "__$gtdb_taxon"$'$' $outdir/tmp.ncbi.txt | head -1`
				else
					sort $outdir/tmp.ncbi.txt | uniq -c | sort -rnk1 > $outdir/tmp.ncbi.count.txt
					top_count=`awk '{print $1}' $outdir/tmp.ncbi.count.txt | head -1`
					N=`expr $field_number - 1`
					while [ `awk -v count=$top_count '{if($1 == count){print $2}}' $outdir/tmp.ncbi.count.txt | wc -l` -gt 1 ] && [ $N -gt 0 ]
					do
						cut -d';' -f1-$N $outdir/tmp.ncbi.txt \
							| sort | uniq -c | sort -rnk1 > $outdir/tmp.ncbi.count.txt
						top_count=`awk '{print $1}' $outdir/tmp.ncbi.count.txt | head -1`
						N=`expr $N - 1`
					done
					ncbi_lineage=`head -1 $outdir/tmp.ncbi.count.txt | sed "s|^.*$top_count ||"`
				fi
			fi
			level=`echo "$ncbi_lineage" | awk -F';' '{print $NF}' | cut -d'_' -f1`
			if [ "$ncbi_taxon" == "" ]
			then
				ncbi_taxon=`echo "$ncbi_lineage" | awk -F'[;_]' '{print $NF}' | sed "s|\[\]||g"`
			fi
			if [ "$level" != "s" ]
			then
				domain=`echo "$ncbi_lineage" | awk -F'[;_]' '{print $3}'`
				if [ "$level" == "g" ]
				then
					ncbi_taxon="$ncbi_taxon sp."
				else
					if [ "$domain" == "Archaea" ]
					then
						ncbi_taxon="$ncbi_taxon archaeon"
					elif [ "$domain" == "Bacteria" ]
					then
						ncbi_taxon="$ncbi_taxon bacterium"
					fi
				fi
			fi
			taxid=`grep $'\t'"uncultured $ncbi_taxon"$'\t' $NCBI_TAXDUMP/names.dmp | head -1 | cut -f1`
			if [ "$taxid" == "" ]
			then
				ncbi_taxon_stripped=`echo "$ncbi_taxon" | sed "s|Candidatus ||g"`
				taxid=`grep $'\t'"uncultured $ncbi_taxon_stripped"$'\t' $NCBI_TAXDUMP/names.dmp | head -1 | cut -f1`
				if [ "$taxid" == "" ] 
				then
					taxid=`grep $'\t'"$ncbi_taxon"$'\t' $NCBI_TAXDUMP/names.dmp | head -1 | cut -f1`
					if [ "$taxid" == "" ] 
					then
						taxid=`grep $'\t'"$ncbi_taxon_stripped"$'\t' $NCBI_TAXDUMP/names.dmp | head -1 | cut -f1`
						if [ "$taxid" == "" ] 
						then
							taxid="NA"
						else
							ncbi_taxon=$ncbi_taxon_stripped
						fi
					fi
				else
					ncbi_taxon=$ncbi_taxon_stripped
				fi
			fi	
		else
			ncbi_taxon="NA"
			taxid="NA"
		fi
		echo "$line,$ncbi_lineage,$ncbi_taxon,$taxid" \
			>> $outdir/tmp.csv
	fi
done < $outdir/gtdbtk.summary.csv
mv $outdir/tmp.csv $outdir/gtdbtk.summary.csv
touch $outdir/complete

rm -rf $outdir/*tmp*

