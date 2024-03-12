#!/usr/bin/env bash

GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release214
NCBI_TAXDUMP="$GTDBTK_DB/`ls $GTDBTK_DB | grep 'ncbi_taxdump' | sort -r | head -1`"
ALTERNATE_KEY=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/ncbi_altkey.csv
IMG_DIR=/software/team311/ng13/local/images
bac_prefix="bac120"
ar_prefix="ar53"


usage="SYNOPSIS: gtdb2ncbi.sh [OPTIONS]\n\
DESCRIPTION: Takes GTDB-format lineage and returns corresponding NCBI lineage, taxon, and taxon ID.\n\
OPTIONS:\n\
	-l | --gtdb_lineage\tGTDB-formatted lineage (should be surrounded by \"\")\n\
	-G | --GTDBTK_DB\tPath to GTDB database [$GTDBTK_DB]\n\
	-N | --NCBI_TAXDUMP\tPath to NCBI database [$NCBI_TAXDUMP] \n\
	-h | --help\t\tPrint help message\n\n\
AUTHOR:\n\
Noah Gettle 2023"


while [ "$1" != "" ]
do
	case $1 in
		-l | --gtdb_lineage)
			if ! `beginswith "-" $2`
			then
				shift
				gtdb_lineage=$1
			fi;;
		-G | --GTDBTK_DB)
			if ! `beginswith "-" $2`
			then
				shift
				GTDBTK_DB=$1
				if ! test -d $GTDBTK_DB
				then
					echo -e "ERROR: Could not find suggested GTDB database directory $GTDBTK_DB.\n"
					exit 1
				fi
			fi;;
		-N | --NCBI_TAXDUMP)
			if ! `beginswith "-" $2`
			then
				shift
				NCBI_TAXDUMP=$1
				if ! test -d $NCBI_TAXDUMP
				then
					echo -e "ERROR: Could not find suggested NCBI database directory $NCBI_TAXDUMP.\n$usage"
					exit 1
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


######################################################################################
# DATABASE CHECK #
##################
# Check for required GTDB data
current_day=`date +'%Y%m%d'`

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


######################################
gtdb_taxa=`echo $gtdb_lineage | sed "s|;|\n|g" | sed "/^.__$/d" | tac | sed "s| |@|g"`
for gtax in $gtdb_taxa
do
	gtax=`echo $gtax | sed "s|@| |g"`
	level=`echo $gtax | cut -c1`
	gtax=`echo $gtax | sed "s|.__||g"`
	if [ "$level" == "s" ] && ! `grep -q $'\t'"$gtax"$'\t' $NCBI_TAXDUMP/names.dmp` && ! `grep -q $'\t'"Candidatus $gtax"$'\t' $NCBI_TAXDUMP/names.dmp`
	then
		gtax=`echo "$gtax" | cut -d' ' -f1`
		level="g"
	fi
	if ! `grep $'\tscientific name\t'  $NCBI_TAXDUMP/names.dmp | grep -q $'\t'"$gtax"$'\t'` && `grep $'\tscientific name\t'  $NCBI_TAXDUMP/names.dmp | grep -q $'\t'"Candidatus $gtax"$'\t'`
	then
		gtax="Candidatus $gtax"
	fi
	if [ "$level" != "s" ] && [ "$level" != "g" ]
	then
		ncbi_taxon="$gtax bacterium"
		if `echo $gtdb_taxa | grep -q 'd__Archaea'`
		then
			ncbi_taxon="$gtax archaeon"
		fi	
	elif [ "$level" == "g" ]
	then
		ncbi_taxon="uncultured $gtax sp."
		if ! `grep $'\tscientific name\t'  $NCBI_TAXDUMP/names.dmp | grep -q $'\t'"$ncbi_taxon"$'\t'` && `grep $'\tscientific name\t'  $NCBI_TAXDUMP/names.dmp | grep -q $'\t'"$gtax sp."$'\t'`
		then
			ncbi_taxon="$gtax sp."
		fi
	else
		ncbi_taxon="$gtax"
	fi
	if `grep -q $'^'"$gtax," $ALTERNATE_KEY`
	then
		alt_line=`grep $'^'"$gtax," $ALTERNATE_KEY`
		ncbi_taxid=`echo $alt_line | cut -d',' -f3`
		ncbi_taxon=`echo $alt_line | cut -d',' -f2`
		alt=true
	else
		ncbi_taxid=`grep $'\tscientific name\t'  $NCBI_TAXDUMP/names.dmp | grep $'\t'"$ncbi_taxon"$'\t' | awk '{print $1}'`
		alt=false
	fi
	if [ "$ncbi_taxid" == "" ]
	then
		if `grep $'\tscientific name\t'  $NCBI_TAXDUMP/names.dmp | grep -q $'\t'"$gtax"$'\t'` && ! `grep -q $'^'"$gtax," $ALTERNATE_KEY`
		then
			ncbi_taxid="NA"
		fi
	fi
	if [ "$ncbi_taxid" != "" ]
	then
		tmp_taxon=$gtax
		if [ "$alt" == "true" ] 
		then
			tmp_taxon=`echo $ncbi_taxon | sed "s| bacterium||g" | sed "s|uncultured ||g" | sed "s| cyanobacterium||g" | sed "s| archaeon||g" | sed "s| sp\.||g"`
		fi
		tmp_taxid=`grep $'\t'"$tmp_taxon"$'\t' $NCBI_TAXDUMP/names.dmp | grep $'\tscientific name\t'| awk '{print $1}'`
		node=`grep $'^'$tmp_taxid$'\t' $NCBI_TAXDUMP/nodes.dmp | head -1 |  awk -F'\t' '{print $3" "$5}'`
		tmp_parent=`echo $node | cut -d' ' -f1`
		tmp_level=`echo $node | cut -d' ' -f2 | sed "s|superkingdom|domain|g" | cut -c1`
		if [ "$tmp_level" != "s" ]
		then
			ncbi_lineage="${tmp_level}__${tmp_taxon}"
		else
			ncbi_lineage=""
		fi
		while [ "$tmp_level" != "d" ]
		do
			tmp_taxid="$tmp_parent"
			tmp_taxon=`grep $'^'$tmp_taxid$'\t' $NCBI_TAXDUMP/names.dmp | grep $'\tscientific name\t' | awk -F'\t' '{print $3}'`
			node=`grep $'^'$tmp_taxid$'\t' $NCBI_TAXDUMP/nodes.dmp | head -1 | awk -F'\t' '{OFS=","; print $3,$5}'`
			tmp_parent=`echo $node | cut -d',' -f1`
			tmp_level_full=`echo $node | cut -d',' -f2 | sed "s|superkingdom|domain|g"`
			tmp_level=`echo $node | cut -d',' -f2 | sed "s|superkingdom|domain|g" | cut -c1`
			if [ "$tmp_level_full" == "domain" ] || [ "$tmp_level_full" == "kingdom" ] || [ "$tmp_level_full" == "phylum" ] || [ "$tmp_level_full" == "class" ] || [ "$tmp_level_full" == "order" ] || [ "$tmp_level_full" == "family" ] || [ "$tmp_level_full" == "genus" ]
			then
				ncbi_lineage="${tmp_level}__${tmp_taxon};$ncbi_lineage"
			fi
		done
		echo -e "$ncbi_lineage,$ncbi_taxon,$ncbi_taxid" | sed "s|,bacterium,|,Bacteria bacterium,|g"
		break
	fi
done
