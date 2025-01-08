#!/usr/bin/env bash


usage="SYNOPSIS: bin_data2csvs.sh <command> [OPTIONS]\n\
DESCRIPTION: Collates bin data information into separate CSVs for use w/ BTK.\n\
OPTIONS:\n\
    -a | --assembly              Assembly (can be gzipped) (Required)\n\
    -d | --directory         Directory containing MAG binning output (Required).\n\
    -b | --bin_data              Bin data file (csv). (Required).\n\
	-R | --bin_dir               Directory containing bin fasta files[\$bin_directory/output_bins]\n\
    -x | --biosample_accessions  CSV of BioSample accessions in the format \"bin_name,tolid,accession\"\n\
	-C | --contig_info           Contig info TSV in format \"contig\tcontig_length\tcontig_depth\tcircular\"\n\
    -o | --outdir                Output directory [./bin_data2csvs_output]\n\
    -h | --help                  Print help message\n\n\
\n\n\
AUTHOR:\n\
Noah Gettle 2023"


outdir="./bin_data2csvs_output"
biosample_accessions=""
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
		-d | --directory)
			if ! `beginswith "-" $2`
			then
				shift
				directory=$1
				if ! test -d $directory
				then
					echo -e "ERROR: Could not find input MAG directory $directory.\n"
					exit 1
				fi
			fi;;
		-b | --bin_data)
			if ! `beginswith "-" $2`
			then
				shift
				bin_data=$1
				if ! test -f $bin_data
				then
					echo -e "ERROR: Could not find bin data file $bin_data.\n\n$usage"
					exit 1
				fi
			fi;;
		-R | --bin_dir)
			if ! `beginswith "-" $2`
			then
				shift
				bin_dir=$1
				if ! test -d $bin_dir
				then
					echo -e "ERROR: Could not find bin directory $bin_dir.\n\n$usage"
					exit 1
				fi
			fi;;
		-x | --biosample_accessions)
			if ! `beginswith "-" $2`
			then
				shift
				biosample_accessions=$1
				if ! test -f $biosample_accessions
				then
					echo -e "#### ERROR: Could not find BioSample Accession file $biosample_accessions. ####\n\n$usage"
				fi
			fi;;
		-C | --contig_info)
			if ! `beginswith "-" $2`
			then
				shift
				contig_info=$1
				if ! test -f $contig_info
				then
					echo -e "#### ERROR: Could not find input contig info file $contig_info. ####\n\n$usage"
					exit 1
				fi
			fi;;
		-o | --outdir)
			if ! `beginswith "-" $2`
			then
				shift
				outdir=$1
				mkdir -p $outdir
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
	echo -e "ERROR: Assemby file required."
	exit 1
fi


if [ "$directory" == "" ]
then
	echo -e "ERROR: Binning output directory required."
	exit 1
fi

if [ "$bin_data" == "" ]
then
	echo -e "ERROR: Bin data file required."
	exit
fi

if [ "$bin_dir" == "" ]
then
	if test -d $directory/output_bins
	then
		bin_dir=$directory/output_bins
	else
		echo -e "ERROR: Directory containing binned MAGs could not be inferred and was not provided."
		exit 1
	fi
fi

if [ `echo $assembly | awk -F'.' '{print $NF}'` == "gz" ]
then
	gunzip -c $assembly > $directory/assembly.tmp.fa
	assembly=$directory/assembly.tmp.fa
fi

# assembly=/lustre/scratch124/tol/projects/darwin/users/ng13/asg/data/protists/Heterometopus_palaeformis/working/piHetPala1.meta-mdbg.20230726
# directory=/lustre/scratch124/tol/projects/darwin/users/ng13/asg/data/protists/Heterometopus_palaeformis/working/piHetPala1.meta-mdbg.20230726/mags/all/magscot_raw
# bin_data=/lustre/scratch124/tol/projects/darwin/users/ng13/asg/data/protists/Heterometopus_palaeformis/working/piHetPala1.meta-mdbg.20230726/map_output/output_data.csv
# bin_dir=$directory/output_bins
# contig_info=/lustre/scratch124/tol/projects/darwin/users/ng13/asg/data/protists/Heterometopus_palaeformis/working/piHetPala1.meta-mdbg.20230726/contig_info.tsv
# outdir=/lustre/scratch124/tol/projects/darwin/users/ng13/asg/data/protists/Heterometopus_palaeformis/working/piHetPala1.meta-mdbg.20230726/map_output/bin_data_csvs

bin_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^bin_id$' | cut -d':' -f1`
qual_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^quality$' | cut -d':' -f1`
gtdb_class_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^classification$' | cut -d':' -f1`
ncbi_class_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^ncbi_classification$' | cut -d':' -f1`
species_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^ncbi_taxon$' | cut -d':' -f1`
drep_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^drep$' | cut -d':' -f1`
completeness_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^Completeness$' | cut -d':' -f1`
contamination_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^Contamination$' | cut -d':' -f1`

rm -rf $outdir/*
mkdir -p $outdir/tmp
bins=`cut -d',' -f$bin_col $bin_data | tail -n +2`
echo '"identifier","checkm_predicted_ssu_count"' > $outdir/checkm_predicted_ssu_count.csv
for bin in $bins
do
	echo $bin
	rm -rf $outdir/tmp/*
	qual=`grep ",$bin," $bin_data | cut -d',' -f$qual_col`
	gtdb_class=`grep ",$bin," $bin_data | cut -d',' -f$gtdb_class_col`
	ncbi_class=`grep ",$bin," $bin_data | cut -d',' -f$ncbi_class_col`
	species=`grep ",$bin," $bin_data | cut -d',' -f$species_col`
	drep=`grep ",$bin," $bin_data | cut -d',' -f$drep_col`
	completeness=`grep ",$bin," $bin_data | cut -d',' -f$completeness_col`
	contamination=`grep ",$bin," $bin_data | cut -d',' -f$contamination_col`
	groups="domain phylum class order family genus species"
	rm -rf $outdir/tmp/gtdb.class.csv
	for group in $groups
	do
		X=${group:0:1}
		G=`echo $gtdb_class | sed "s|;|\n|g" | grep $'^'$X"_" | awk -F'__' '{print $2}'`
		if [ "$G" != "" ]
		then
			echo -e "$group,$G" >> $outdir/tmp/gtdb.class.csv
		else
			echo -e "$group,unknown" >> $outdir/tmp/gtdb.class.csv
		fi
	done
	groups="domain phylum class order family genus"
	rm -rf $outdir/tmp/ncbi.class.csv
	for group in $groups
	do
		X=${group:0:1}
		G=`echo $ncbi_class | sed "s|;|\n|g" | grep $'^'$X"_" | awk -F'__' '{print $2}'`
		if [ "$G" != "" ]
		then
			echo -e "$group,$G" >> $outdir/tmp/ncbi.class.csv
		else
			echo -e "$group,unknown" >> $outdir/tmp/ncbi.class.csv
		fi
	done
	echo -e "species,$species" >> $outdir/tmp/ncbi.class.csv
	if [ "$biosample_accessions" != "" ] && test -f $biosample_accessions
	then
		biosample=`grep $'^'$bin$',' $biosample_accessions | cut -d',' -f3`
		tolid=`grep $'^'$bin$',' $biosample_accessions | cut -d',' -f2`
	else
		biosample=""
		tolid=""
	fi
	contigs=`grep $'^>' $bin_dir/$bin.fa | cut -d'>' -f2 | sort -u`
	for contig in $contigs
	do
		if [ "$contig_info" != "" ]
		then
			if ! test -f $outdir/circular.csv
			then
				echo '"identifier","circular"' > $outdir/circular.csv
			fi
			circular=`grep $'^'$contig$'\t' $contig_info | cut -f4 | sed "s|Y|TRUE|g" | sed "s|N|FALSE|g"`
			echo "\"$contig\",\"$circular\"" >> $outdir/circular.csv
		fi
		while read line
		do
			group=`echo $line | cut -d',' -f1`
			if ! test -f $outdir/gtdb_$group.csv
			then
				echo "\"identifier\",\"gtdb_$group\"" > $outdir/gtdb_$group.csv
			fi
			tax=`echo $line | cut -d',' -f2`
			echo "\"$contig\",\"$tax\"" >> $outdir/gtdb_$group.csv
		done < $outdir/tmp/gtdb.class.csv
		while read line
		do
			group=`echo $line | cut -d',' -f1`
			if ! test -f $outdir/ncbi_$group.csv
			then
				echo "\"identifier\",\"ncbi_$group\"" > $outdir/ncbi_$group.csv
			fi
			tax=`echo $line | cut -d',' -f2`
			echo "\"$contig\",\"$tax\"" >> $outdir/ncbi_$group.csv
		done < $outdir/tmp/ncbi.class.csv
		if [ "$biosample" != "" ]
		then
			if ! test -f $outdir/biosample.csv
			then
				echo '"identifier","biosample_id"' > $outdir/biosample.csv
			fi
			echo "\"$contig\",\"$biosample\"" >> $outdir/biosample.csv
		fi
		if [ "$tolid" != "" ]
		then
			if ! test -f $outdir/tolid.csv
			then
				echo '"identifier","tolid"' > $outdir/tolid.csv
			fi
			echo "\"$contig\",\"$tolid\"" >> $outdir/tolid.csv
		fi
		if test -d $directory/checkm_ssu && test -f $directory/checkm_ssu/ssu_summary.tsv
		then
			if ! test $outdir/checkm_predicted_ssu_count.csv
			then
				echo '"identifier","checkm_predicted_ssu_count"' > $outdir/checkm_predicted_ssu_count.csv
			fi
			checkm_ssus=`grep $'\t'$contig$'[\t-]' $directory/checkm_ssu/ssu_summary.tsv | wc -l`
			echo "\"$contig\",$checkm_ssus" >> $outdir/checkm_predicted_ssu_count.csv
		fi
		if ! test -f $outdir/other_bin_info.csv
		then
			echo '"identifier","bin","passed_dereplication","bin_completeness","bin_contamination"' > $outdir/other_bin_info.csv
		fi
		echo "\"$contig\",\"$bin\",\"$drep\",$completeness,$contamination" >> $outdir/other_bin_info.csv
	done
done
contigs=`grep $'^>' $assembly | cut -d'>' -f2`
for contig in $contigs
do
	if ! `grep -q $'^\"'$contig'\",' $outdir/other_bin_info.csv`
	then
		for file in $outdir/*.csv
		do
			if [ "$file" != "$outdir/other_bin_info.csv" ]
			then
				if [ "$file" == "$outdir/checkm_predicted_ssu_count.csv" ]
				then
					echo "\"$contig\",0" >> $file
				else
					echo "\"$contig\",\"NA\"" >> $file
				fi
			else
				echo "\"$contig\",\"NA\",\"NA\",0,0" >> $file
			fi
		done
	fi
done
		
rm -rf $outdir/tmp

