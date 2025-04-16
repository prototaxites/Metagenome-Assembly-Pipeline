#!/usr/bin/env bash

GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release214
NCBI_TAXDUMP="$GTDBTK_DB/`ls $GTDBTK_DB | grep 'ncbi_taxdump' | sort -r | head -1`"
group_key=/lustre/scratch124/tol/projects/darwin/users/ng13/group_key.tsv
metagenome_taxids=/lustre/scratch124/tol/projects/darwin/users/ng13/organismal_metagenome_taxids.tsv
LOCAL_IMAGES=/software/team311/ng13/local/images

usage="SYNOPSIS: metagenome_submission.sh <command> [OPTIONS]\n\
DESCRIPTION: Takes input information to prep and send metagenome curation requests.\n\
  Commands:\n\
    config_file        --> Generate BioSample submission config file\n\
    check              --> Classify bins with GTDB-TK, if necessary, and create submission form for new taxon ids\n\
    biosample_prep     --> Create directories containing information for biosample creation\n\
    spreadsheet_update --> Create CSV to update the ToL Informatics Sheet\n\
	generate_chromosome_list --> Annotate MAG contigs as plasmids/chromosomes\n\
    curation_request   --> Create YAMLs and send curation requests\n\
OPTIONS:\n\
    -a | --assembly              Assembly (can be gzipped)\n\
    -d | --directory             Directory containing MAG binning output.\n\
    -i | --tol_id                Host ToL id\n\
    -A | --assembler             Assembly program\n\
    -B | --binning_program       Binning program\n\
    -R | --refining_program      Refining program (if using only one binning program). Default = 'dastool'.\n\
    -b | --bin_data              Bin data file (csv). Required for 'check','biosample_prep', and 'curation_request'\n\
	-R | --bin_dir               Directory containing bin fasta files\n\
    -P | --part                  Part of assembly used for binning [primary, all]\n\
    -x | --biosample_accessions  CSV of BioSample accessions in the format \"bin_name,tolid,accession\"\n\
    -p | --project               Sanger project (see /lustre/scratch124/tol/projects)\n\
    -g | --group                 Organismal group\n\
    -T | --primary_taxon         Primary metagenome taxon (e.g. 'sponge metagenome')\n\
    -s | --species               Species\n\
    -c | --biosample_configs     BioSample config file\n\
	-C | --contig_info           Contig info TSV in format \"contig\tcontig_length\tcontig_depth\tcircular\"\n\
	-l | --chromosome_list		 Contig chromosome annotations in format \"bin\tcontig\tchromosome_type\"\n\
	-S | --submit				 Submit curation requests for MAGS rather than just creating submission directories\"\n\
    -o | --outdir                Output directory\n\
	-O | --cr_outdir			 Curation request output directory [ <tol_id>.metagenome.<yyyymmdd>]\n\
	-V | --metagenome_version	 Version of the metagenome\n\
	-K | --btk_dir				 BTK directory\n\
	-F | --force				 Overwrite final outputs for submission\n\
    -h | --help                  Print help message\n\n\
\n\n\
AUTHOR:\n\
Noah Gettle 2023"

if ! `beginswith "-" $1`
then
	command=$1
	shift
	if [ "$command" != "config_file" ] && [ "$command" != "check" ] && [ "$command" != "biosample_prep" ] && [ "$command" != "spreadsheet_update" ] && [ "$command" != "curation_request" ] && [ "$command" != "generate_chromosome_list" ] && [ "$command" != "submit_cr" ] && [ "$command" != "create_btk" ]
	then
		echo -e "#### UNKNOWN COMMAND: $command ####\n\n$usage"
		exit 1
	fi
else
	echo -e $usage
fi

biosample_configs="./metagenome_pipeline/biosample_configs.txt"
primary_taxon=""
threads=1
refining_program="dastool"
metagenome_taxids=/lustre/scratch124/tol/projects/darwin/users/ng13/organismal_metagenome_taxids.tsv
outdir="."
submit="false"
force="false"
metagenome_version=1
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
		-i | --tol_id)
			if ! `beginswith "-" $2`
			then
				shift
				prefix=$1
			fi;;
		-A | --assembler)
			if ! `beginswith "-" $2`
			then
				shift
				assembler=$1
			fi;;
		-B | --binning_program)
			if ! `beginswith "-" $2`
			then
				shift
				binning_program=$1
			fi;;
		-R | --refining_program)
			if ! `beginswith "-" $2`
			then
				shift
				refining_program=$1
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
		-P | --part)
			if ! `beginswith "-" $2`
			then
				shift
				part=$1
				if [ "$part" != "primary" ] && [ "$part" != "all" ] && [ "$part" != "consensus" ]
				then
					echo -e "ERROR: Value for -p/--part must be either 'primary' or 'all' (Value $part not accepted)\n\n$usage"
					exit 1
				fi
			fi;;
		-p | --project)
			if ! `beginswith "-" $2`
			then
				shift
				project=$1
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
		-g | --group)
			if ! `beginswith "-" $2`
			then
				shift
				group=$1
			fi;;
		-V | --metagenome_version)
			if ! `beginswith "-" $2`
			then
				shift
				metagenome_version=$1
			fi;;
		-s | --species)
			if ! `beginswith "-" $2`
			then
				shift
				species=$1
			fi;;
		-c | --biosample_configs)
			if ! `beginswith "-" $2`
			then
				shift
				biosample_configs=$1
				if [ "$command" != "config_file" ] && ! test -f $biosample_configs
				then
					echo -e "#### ERROR: Could not find input BioSample config file $biosample_configs. ####\n\n$usage"
					exit 1
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
		-K | --btk_dir)
			if ! `beginswith "-" $2`
			then
				shift
				btk_dir=$1
				if ! test -d $btk_dir
				then
					echo -e "#### ERROR: Could not find btk directory $btk_dir. ####\n\n$usage"
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
		-O | --cr_outdir)
			if ! `beginswith "-" $2`
			then
				shift
				cr_outdir=$1
			fi;;
		-l | --chromosome_list)
			if ! `beginswith "-" $2`
			then
				shift
				chromosome_list=$1
				if ! test -f $chromosome_list
				then
					echo -e "#### ERROR: Could not find input chromosome list file $chromosome_list. ####\n\n$usage"
					exit 1
				fi
			fi;;
		-S | --submit)
			shift
			submit="true";;
		-F | --force)
			shift
			force="true";;
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

if [ "$binning_program" == "" ] && [ "$command" != "spreadsheet_update" ] && [ "$command" != "curation_request" ]
then
	echo -e "ERROR: Binning program must be provided.\n\n$usage"
	exit 1
fi

if ([ "$command" != "config_file" ] || [ "$command" != "biosample_prep" ]) && ! test -f $bin_data
then
	echo -e "ERROR: Bin data file required to run '$command'.\n\n$usage"
	exit -1
fi

if [ "$command" == "biosample_prep" ] && ([ "$biosample_configs" == "" ] || ! test -f $biosample_configs)
then
	echo -e "ERROR: BioSample configuration file is required to run '$command'.\n\n$usage"
	exit 1
fi

if [ "$directory" != "" ]
then
	if [ "$project" == "" ]
	then
		project=`echo $directory | awk -F'/' '{print $(NF - 7)}'`
	fi
	if [ "$group" == "" ]
	then
		group=`echo $directory | awk -F'/' '{print $(NF - 5)}'`
	fi
	if [ "$species" == "" ]
	then 
		species=`echo $directory | sed "s|_sp._|_sp_|g" | awk -F'[/]' '{print $(NF - 4)}' | sed "s|_sp_|_sp._|g"`
	fi
	if [ "$tol_id" == "" ]
	then
		tol_id=`echo $directory | sed "s|_sp._|_sp_|g" | awk -F'[/.]' '{print $(NF - 4)}'`
	fi
	if [ "$assembler" == "" ]
	then
		assembler=`echo $directory | sed "s|_sp._|_sp_|g" | awk -F'[/.]' '{print $(NF - 3)}'`
		if [ "$assembler" == "metaflye" ]
		then
			assembler="flye"
		fi
	fi
	if [ "$part" == "" ]
	then
		part=`echo $directory | sed "s|_sp._|_sp_|g" | awk -F'[/.]' '{print $NF}'`
	fi
	if [ "$assembly" == "" ]
	then
		assembly_dir=`dirname $directory`
		assembly_dir=`dirname $assembly_dir`
		if [ "$assembler" == "metaflye" ] || [ "$assembler" == "flye" ]
		then
			assembly=$assembly_dir/assembly.fasta
		elif [ "$assembler" == "hifiasm" ] || [ "$assembler" == "hifiasm-meta" ]
		then
			if [ "$part" == "all" ] || [ "$part" == "alternate_primary" ]
			then
				assembly=$assembly_dir/$tol_id.all.fa
			elif [ "$part" == "primary" ]
			then
				assembly=$assembly_dir/$tol_id.p_ctg.fa
			elif  [ "$part" == "alternate" ]
			then
				assembly=$assembly_dir/$tol_id.a_ctg.fa
			fi
		elif [ "$assembler" == "hifiasm-hic" ]
		then
			if [ "$part" == "all" ] || [ "$part" == "alternate_primary" ]
			then
				assembly=$assembly_dir/$tol_id.all.fa
			elif [ "$part" == "primary" ]
			then
				assembly=$assembly_dir/$tol_id.hic.hap1.p_ctg.fa
			elif [ "$part" == "alternate" ]
			then
				assembly=$assembly_dir/$tol_id.hic.hap2.p_ctg.fa
			fi
		elif [ "$assembler" == "hicanu" ]
		then
			assembly=$assembly_dir/$tol_id.contigs.fasta
		elif [ "$assembler" == "meta-mdbg" ]
		then
			assembly=$assembly_dir/contigs.fasta
		else
			echo -e "ERROR: Could not determine appropriate assembly based on inferred assembler ($assembler) information.\n\n$usage"
			# exit 1
		fi
		if ! test -f $assembly && test -f $assembly.gz
		then
			assembly=$assembly.gz
		fi
		if ! test $assembly
		then
			echo -e "ERROR: No assembly given and none could be inferred based on information given.\n\n$usage"
			# exit 1
		fi
	fi
fi


if [ "$assembler" == "" ]
then
	echo -e "Assembly program required.\n\n$usage"
	exit 1
fi

if [ "$part" == "" ] && ! `echo $assembler | grep -q "hifiasm"`
then
	part="all"
fi

if [ "$part" == "" ]
then
	echo -e "ERROR: Assembly part not provided or could not be inferred.\n\n$usage"
	exit 1
fi

if [ "$species" == "" ]
then
	echo -e "ERROR: No species provided.\n\n$usage"
	exit 1
else
	species_proper=`echo $species | sed "s|_| |g"`
fi

if [ "$project" == "" ]
then
	echo -e "ERROR: No project provided.\n$usage"
	exit 1
fi

if [ "$group" == "" ]
then
	echo -e "ERROR: No organismal group provided.\n\n$usage"
	exit 1
fi

	
if [ "$directory" == "" ]
then
	if ! test -d /lustre/scratch124/tol/projects/$project/data/$group/$species/working
	then
		echo -e "ERROR: Could not find directory /lustre/scratch124/tol/projects/$project/data/$group/$species/working.\n\n$usage"
		exit 1
	fi			
	directory="ls /lustre/scratch124/tol/projects/$project/data/$group/$species/working \
		| grep $'^'$assembler$'\.' | sort -r | head -1"
	directory=$directory/mags/$part
	if ! test -d $directory
	then
		echo -e "ERROR: Could not find directory $directory.\n\n$usage"
		exit 1
	fi
fi
		
if [ "$tol_id" == "" ]
then
	tol_id=`echo $directory | awk -F'[/.]' '{print $11}'`
fi

# if ! test -f $directory/bin_stats.csv
# then
# 	echo -e "ERROR: Could not find bin statistics file $directory/bin_stats.csv.\n\n$usage"
# 	exit 1
# fi
if ! test -d $directory/$binning_program
then
	echo -e "ERROR: Could not find directory $directory/$binning_program.\n\n$usage"
	exit 1
fi

if [ "$GTDBTK_DB" == "" ]
then
	export GTDBTK_DB=/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/release214
fi

if [ "$NCBI_TAXDUMP" == "" ] || ! test -d $NCBI_TAXDUMP || [ `ls $GTDBTK_DB | grep 'ncbi_taxdump' | sort -r | head -1 | awk -F'[/.]' '{print $NF}' | cut -c1-6` -ne `date +'%Y%m'` ]
then
	if [ `ls $GTDBTK_DB | grep 'ncbi_taxdump' | sort -r | head -1 | awk -F'[/.]' '{print $NF}' | cut -c1-6` -ne `date +'%Y%m'` ]
	then
		echo "#### DATABASE CHECK: Updating NCBI taxdump... ####"
	else
		echo "#### DATABASE CHECK: Could not find NCBI taxdump. Downloading... ####"
	fi
	current_day=`date +%Y%m%d`
	NCBI_TAXDUMP=$GTDBTK_DB/ncbi_taxdump.$current_day
	mkdir -p $NCBI_TAXDUMP
	cd $NCBI_TAXDUMP
	wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 
	tar -zxvf $NCBI_TAXDUMP/taxdump.tar.gz 
	rm -rf $NCBI_TAXDUMP/taxdump.tar.gz
fi



if [ "$primary_taxon" == "" ]
then
	tol_group=`echo ${tol_id:0:2}`
	primary_taxon=`grep $'[\t,]'$tol_group $metagenome_taxids | cut -f1 | sed "s| \[species\]||g"`
	if [ "$primary_taxon" == "" ]
	then
		echo -e "ERROR: Could not find appropriate primary metagenome taxon ID for $tol_id,$tol_group"
		exit 1
	fi
	informal_group=`grep $'[\t,]'$tol_group $metagenome_taxids | cut -f4`
fi
primary_taxid=`grep "$primary_taxon" $NCBI_TAXDUMP/names.dmp | cut -f1`
if [ "$primary_taxid" == "" ]
then
	echo -e "ERROR: Could not find a TaxID associated with \"$primary_taxon\" in $NCBI_TAXDUMP/names.dmp"
	exit 1
fi

if `echo $binning_program | grep -q 'magscot'` || `echo $binning_program | grep -q 'dastool'`
then
	bindir=$directory/$binning_program
else
	bindir=$directory/$binning_program/$refining_program
fi
if [ "$command" != "config_file" ] && ! test -d $bindir/output_bins
then
	echo -e "ERROR: Could not find expected bin directory ($bindir/output_bins).\n\n$usage"
	exit 1
fi

if [ "$command" == "check" ]
then
	# Subset bin data and remove empty classification columns
	head -1 $bin_data > $outdir/output_data.csv
	if `echo $binning_program | grep -q 'magscot'` || `echo $binning_program | grep -q 'dastool'`
	then
		grep ",$assembler,$part,$binning_program," $bin_data >> $outdir/output_data.csv
	else
		grep ",$assembler,$part,$binning_program,$refining_program," $bin_data >> $outdir/output_data.csv
	fi
	classification_col=`head -1 $outdir/output_data.csv | sed "s|,|\n|g" | grep -n $'^classification$' | cut -d':' -f1`
	taxon_col=`head -1 $outdir/output_data.csv | sed "s|,|\n|g" | grep -n $'^taxon_id$' | cut -d':' -f1`
	if [ `tail -n +2 $outdir/output_data.csv | awk -F',' -v N=$classification_col '{print $N}' | sort -u | wc -l` -eq 1 ] && [ `tail -n +2 $outdir/output_data.csv | awk -F',' -v N=$classification_col '{print $N}' | sort -u` == "NA" ]
	then
		n=`expr $classification_col - 1`
		m=`expr $taxon_col + 1`
		cut -d',' -f1-$n $outdir/output_data.csv > $outdir/a.tmp.csv
		cut -d',' -f$m- $outdir/output_data.csv > $outdir/b.tmp.csv
		paste -d','  $outdir/a.tmp.csv $outdir/b.tmp.csv > $outdir/output_data.csv
	fi	
fi

if test -d $assembly_dir && [ "$contig_info" == "" ]
then
	if test -f $assembly_dir/contig_info.tsv
	then
		contig_info=$assembly_dir/contig_info.tsv
	fi
fi

######################################################################################
#### CREATE BIOSAMPLE CONFIG FILE ####
######################################
hifiasm_version="Hifiasm (version 0.19.8)"
hifiasm_meta_version="Hifiasm-Meta (version 0.0.2); Hifiasm (vversion 0.19.8)"
flye_version="Flye (version 2.9.3)"
canu_version="Canu (version 2.2)"
meta_mdbg_version="metaMDBG (version 0.3)"
metabat2_version="MetaBat (version 2.15)"
bin3c_version="bin3C (version 0.3.3)"
maxbin2_version="MaxBin (version 2.2.7)"
metator_version="MetaTOR"
dastool_version="DASTOOL (version 1.1.5)"
magscot_version="MAGScoT (version 1.0.0)"
checkm_version="checkM (version 1.2.2); checkM_DB (release 2015-01-16)"
prokka_version="PROKKA (version 1.14.6)"
gtdbtk_version="GTDB-TK (version 2.4.0); GTDB (release 214)"
minimap_version="Minimap2 (version 2.24-r1122)"

if [ "$command" == "config_file" ]
then
	isolation_source="${informal_group^}: $species_proper"
	as=`echo "${assembler}_version" | sed "s|-|_|g" | sed "s|metaflye|flye|g"`
	asm_soft=${!as}
	binner=`echo $binning_program | sed "s|_raw||g" | sed "s|_filtered||g"`
	
	bs="${binner}_version"
	binning_soft=${!bs}
	if `echo $binning_program | grep -q 'dastool'` || `echo $binning_program | grep -q 'magscot'`
	then
		binning_soft="$binning_soft; $metabat2_version; $bin3c_version; $maxbin2_version; $metator_version"
	else
		binning_soft="$binning_soft; $dastool_version"
	fi
	coverage_soft=$minimap_version
	if [ "$binning_program" == "bin3c" ] || [ "$binning_program" == "metator" ]
	then
		binning_params="hic-mapping"
	elif [ "$binning_program" == "metabat2" ]
	then
		binning_params="coverage; graph"
	elif [ "$binning_program" == "maxbin2" ]
	then
		binning_params="coverage"
	else
		binning_params="coverage; graph; hic-mapping"
	fi
	if [ "$project" == "asg" ]
	then
		broadscale_environmental_context="aquatic biome"
	else
		broadscale_environmental_context=""
	fi
	if [ "$group" == "protists" ]
	then
		environmental_medium="$informal_group culture"
	else
		environmental_medium="adult $informal_group tissue"
	fi
	sts_line=`grep $'\t'"$tol_id"$'\t' /lustre/scratch123/tol/tolqc/track/tol_sts.tsv | sed "s|\t|;|g"`
	host_biospecimen=`echo $sts_line | awk -F';' '{print $8}'`
	host_taxid=`echo $sts_line | awk -F';' '{print $15}'`
	metagenomic_source="$primary_taxon"
	echo -e "\
# SOFTWARE #\n\
asm_soft=\"$asm_soft\"\n\
binning_soft=\"$binning_soft\"\n\
binning_params=\"$binning_params\"\n\
coverage_soft=\"$minimap_version\"\n\
completeness_soft=\"$checkm_version\"\n\
s16_soft=\"$prokka_version\"\n\
tRNA_soft=\"$prokka_version\"\n\
taxonomic_classification=\"$gtdbtk_version\"\n\n\
# SAMPLE INFO #\n\
host_biospecimen=\"$host_biospecimen\"\n\
host_taxname=\"$species_proper\"\n\
host_taxid=\"$host_taxid\"\n\
primary_taxon=\"$primary_taxon\"\n\
primary_taxid=\"$primary_taxid\"     # Inferred from file \"$metagenome_taxids\". See also https://www.ebi.ac.uk/ena/browser/view/410656?dataType=TAXON&show=tax-tree\n\
sequencing_method=\"Pacbio Sequel II\"\n\
investigation_type=\"metagenome-assembled genome\"\n\
broadscale_environmental_context=\"$broadscale_environmental_context\"   # Use EnvO \"biome\" term (http://purl.obolibrary.org/obo/ENVO_00000428)\n\
local_environmental_context=\"\"    # Subclass of \"broadscale_environmental_context\". Anatomical sites can work as well (https://www.ebi.ac.uk/ols/ontologies)\n\
isolation_source=\"$isolation_source\"\n\
environmental_medium=\"$environmental_medium\"\n\
metagenomic_source=\"$metagenomic_source\"\n\
" > $biosample_configs
fi


######################################################################################
#### CHECK BINS ####
####################
if [ "$command" == "check" ]
then
	# Expecting to find 'taxon_id' field in bin data. If not, run GTDB-TK classification
	taxon_id_col=`head -1 $outdir/output_data.csv | sed "s|,|\n|g" | grep -n $'^taxon_id$' | cut -d':' -f1`
	if ! `head -1 $outdir/output_data.csv | sed "s|,|\n|g" | grep -q $'^taxon_id$'`
	then
		mkdir -p $bindir/gtdb/output
		if ! test -f $bindir/gtdb/output/complete
		then
			rm -rf $outdir/gtdb-final_outputs.*
			if ! `ls $bindir/gtdb/output | grep -q $'gtdbtk.bac*.summary.tsv'` && ! `ls $bindir/gtdb/output | grep -q $'gtdbtk.ar*.summary.tsv'`
			then
				test_bin=`tail -n +2 $outdir/output_data.csv | cut -d',' -f4`
				extension=`ls $bindir/output_bins | grep $test_bin$'.' | awk -F'.' '{print $NF}'`
				echo -e "Classifying bins with GTDB-TK. This may take some time...."
				bsub -n$threads -q long -R"span[hosts=1]" \
					-o $outdir/gtdb-final_outputs.o \
					-e $outdir/gtdb-final_outputs.e \
					-M100000 \
					-R 'select[mem>100000] rusage[mem=100000]' \
						"gtdbtk_wrapper.sh \
							-d $bindir/output_bins \
							-x $extension \
							-o $bindir/gtdb/output"
				n=0
				while ! test -f $bindir/gtdb/output/complete && (! test -f $outdir/gtdb-final_outputs.o || ! `grep -q 'ERROR' $outdir/gtdb-final_outputs.o`) && [ $n -le 288 ]
				do
					sleep 5m
					n=`expr $n + 1`
				done
				if ! test -f $bindir/gtdb/output/complete && ! test -f $bindir/gtdb/output/complete && ! `grep -q 'ERROR' $outdir/gtdb-final_outputs.o`
				then
					echo "#### ERROR --- GTDB-TK did not finish within 24 hrs. Exiting ####"
					exit 1
				fi
			else
				gtdbtk_wrapper.sh \
					-d $input_bins \
					-x $extension \
					-o $bindir/gtdb/output
			fi			
		fi
		gtdbtk_file=$bindir/gtdb/output/gtdbtk.summary.csv
		gtdbtk_file_header="classification,fastani_reference,fastani_reference_radius,fastani_taxonomy,fastani_ani,fastani_af,closest_placement_reference,closest_placement_radius,closest_placement_taxonomy,closest_placement_ani,closest_placement_af,pplacer_taxonomy,classification_method,note,other_related_references(genome_id;species_name;radius;ANI;AF),msa_percent,translation_table,red_value,warnings,ncbi_classification,ncbi_taxon,taxon_id"
		main_file_header=`head -1 $outdir/output_data.csv`
		echo "$main_file_header,$gtdbtk_file_header" > $outdir/tmp.csv
		bins=`cut -d',' -f6 $outdir/output_data.csv | tail -n +2 | sort -u`
		for bin in $bins
		do
			part1=`grep $','$bin$',' $outdir/output_data.csv`
			part2=`grep $','$bin$',' $gtdbtk_file | cut -d',' -f2-`
			echo "$part1,$part2" >> $outdir/tmp.csv
		done
		mv $outdir/tmp.csv $outdir/output_data.csv
	fi
	taxon_field=`head -1 $outdir/output_data.csv | sed "s|,|\n|g" | grep -n $'^taxon_id$' | cut -d':' -f1`
	awk -F',' -v N=$taxon_field '{if($N == "NA"){print $(N-2)","$(N-1)","$N}}' $outdir/output_data.csv \
		| sort -u > $outdir/missing_taxon_ids.csv
	host=`echo $species | sed "s|_| |g"`
	name_type="Unidentified Species"
	project_id="NA"
	echo -e "proposed_name\tname_type\thost\tproject_id\tdescription" > $outdir/prokaryotic_taxonomy_registration.$tol_id.tsv
	while read line
	do
		proposed_name=`echo $line | cut -d',' -f2`
		description=`echo $line | cut -d',' -f1`
		if ! `grep -q $'\t'"$proposed_name"$'\t' $NCBI_TAXDUMP/names.dmp`
		then
			echo -e "$proposed_name\t$name_type\t$host\t$project_id\t$description" >>  $outdir/prokaryotic_taxonomy_registration.$tol_id.tsv
		else
			echo "Found $proposed_name"
			taxon_id=`grep $'\t'"$proposed_name"$'\t' $NCBI_TAXDUMP/names.dmp | head -1 | cut -f1`
			sed -i "s|,$proposed_name,NA,|,$proposed_name,$taxon_id,|g" $outdir/output_data.csv
		fi
	done < $outdir/missing_taxon_ids.csv
fi

######################################################################################
#### PUT TOGETHER DATA FOR BIOSAMPLE GENERATION ####
####################################################

if [ "$command" == "biosample_prep" ]
then
	. $biosample_configs
	# Primary metagenome
	echo -e "host_biospecimen,host_taxname,host_taxid,metagenome_taxname,metagenome_taxid,metagenome_tolid,broad-scale environmental context,local environmental context,environmental medium,binned_path,mag_path" \
		> $outdir/primary_biosample.csv
	echo -e "$host_biospecimen,$host_taxname,$host_taxid,$primary_taxon,$primary_taxid,$tol_id.metagenome,$broadscale_environmental_context,$local_environmental_context,$environmental_medium,$outdir/binned_biosample_metadata.csv,$outdir/mag_biosample_metadata.csv" \
		>> $outdir/primary_biosample.csv
	# echo -e "tol_id,taxon,taxon_id,broad-scale environmental context,local environmental context,environmental medium" \
	# 	> $outdir/primary_biosample.csv
	# echo -e "$tol_id.metagenome,$primary_taxon,$primary_taxid,$broadscale_environmental_context,$local_environmental_context,$environmental_medium" \
	# 	>> $outdir/primary_biosample.csv
	# MAGs
	header="bin_name,tol_id,taxon,taxon_id,number of standard tRNAs extracted,assembly software,16S recovered,16S recovery software,tRNA extraction software,completeness score,completeness software,contamination score,binning software,MAG coverage software,binning parameters,taxonomic identity marker,taxonomic classification,assembly quality,sequencing method,investigation type,isolation_source,broad-scale environmental context,local environmental context,environmental medium,metagenomic source"
	echo -e "$header" \
		> $outdir/binned_biosample_metadata.csv
	echo -e "$header" \
		> $outdir/mag_biosample_metadata.csv
	HIGH="Single contiguous sequence without gaps or ambiguities with a consensus error rate equivalent to Q50 or better"
	LOW="Many fragments with little to no review of assembly other than reporting of standard assembly statistics."
	bins=`tail -n +2 $bin_data | cut -d',' -f6`
	drep_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^drep$' | cut -d':' -f1`
	taxon_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^taxon_id$' | cut -d':' -f1`
	trna_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^total_trnas$' | cut -d':' -f1`
	completeness_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^Completeness$' | cut -d':' -f1`
	contamination_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^Contamination$' | cut -d':' -f1`
	s16_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^rrna_16s$' | cut -d':' -f1`
	marker_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^classification_method$' | cut -d':' -f1`
	quality_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^quality$' | cut -d':' -f1`
	rm -rf $outdir/taxon_tally.tmp.txt
	for bin in $bins
	do
		# echo $bin
		line=`grep ","$bin"," $bin_data`
		tolid_taxon=`echo $line \
			| awk -F',' -v N=$taxon_field '{print $(N-1)}' \
			| sed "s|[uU]ncultured ||g" \
			| sed "s| bacterium||g" \
			| sed "s|[- ]|_|g" \
			| sed "s| archaeon||g" \
			| sed "s|[Cc]andidatus_||g" \
			| sed "s|[().:]||g" | sed "s|[[]]||g" | sed "s|\"||g"`
		if [ `echo "$tolid_taxon" | wc -c` -gt 32 ]
		then
			tolid_taxon=${tolid_taxon:0:32}
		fi
		if ! test -f $outdir/taxon_tally.tmp.txt || ! `grep -q $tolid_taxon$'\t' $outdir/taxon_tally.tmp.txt`
		then
			n=1
			echo -e "$tolid_taxon\t$n" >> $outdir/taxon_tally.tmp.txt
		else
			m=`grep $tolid_taxon$'\t' $outdir/taxon_tally.tmp.txt | cut -f2`
			n=`expr $m + 1`
			sed -i "s|^$tolid_taxon\t$m$|$tolid_taxon\t$n|g" $outdir/taxon_tally.tmp.txt
		fi
		bin_tolid=$tol_id.${tolid_taxon}_${n}
		outtaxon=`echo $line | awk -F',' -v N=$taxon_field '{print $(N-1)}'`
		outtaxon_id=`echo $line | awk -F',' -v N=$taxon_field '{print $N}'`
		n_trna=`echo $line | awk -F',' -v N=$trna_field '{print $N}'`
		s16='No'
		if [ `echo $line | awk -F',' -v N=$s16_field '{print $N}'` == "Y" ]
		then
			s16='Yes'
		fi
		completeness=`echo $line | awk -F',' -v N=$completeness_field '{print $N}'`
		contamination=`echo $line | awk -F',' -v N=$contamination_field '{print $N}'`
		taxon_id_marker=`echo $line |awk -F',' -v N=$marker_field '{print $N}'`
		if [ `echo $line | awk -F',' -v N=$quality_field '{print $N}'` == "HIGH" ] && [ `echo $line | awk -F',' -v N=$drep_field '{print $N}'` == "PASSED" ]
		then
			asm_qual=$HIGH
			outline="$bin,$bin_tolid,$outtaxon,$outtaxon_id,$n_trna,$asm_soft,$s16,$s16_soft,$tRNA_soft,$completeness,$completeness_soft,$contamination,$binning_soft,$coverage_soft,$binning_params,$taxon_id_marker,$taxonomic_classification,$asm_qual,$sequencing_method,$investigation_type,$isolation_source,$broadscale_environmental_context,$local_environmental_context,$environmental_medium,$metagenomic_source"
			echo "$outline" \
				>>  $outdir/mag_biosample_metadata.csv
		else
			asm_qual=$LOW
			outline="$bin,$bin_tolid,$outtaxon,$outtaxon_id,$n_trna,$asm_soft,$s16,$s16_soft,$tRNA_soft,$completeness,$completeness_soft,$contamination,$binning_soft,$coverage_soft,$binning_params,$taxon_id_marker,$taxonomic_classification,$asm_qual,$sequencing_method,$investigation_type,$isolation_source,$broadscale_environmental_context,$local_environmental_context,$environmental_medium,$metagenomic_source"
			echo "$outline" \
				>> $outdir/binned_biosample_metadata.csv
		fi
	done
fi

######################################################################################
#### SPREADSHEET UPDATE ####
############################
# host_tolid,primary_tolid,primary_biosample,tolid,taxname,taxid,bin_type,length,contigs,circular_contigs,completeness,contamination,mean_coverage,SSUs,total_trnas,unique_trnas,23s,16s,5s,biosample,bioproject,assembly_accession
if [ "$command" == "spreadsheet_update" ]
then
	. $biosample_configs
	quality_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^quality$' | cut -d':' -f1`
	drep_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^drep$' | cut -d':' -f1`
	primary_biosample=`grep $'^'primary"," $biosample_accessions | cut -d',' -f3`
	primary_tolid=`grep $'^'primary"," $biosample_accessions | cut -d',' -f2`
	bins=`tail -n +2 $biosample_accessions | grep -v 'primary' | cut -d',' -f1`
	header_line=`head -1 $bin_data | sed "s|$|,primary_tolid,primary_biosample,tolid,bin_type,biosample,group,project|" \
		| sed "s|^host,|host,host_species,|g"`
	host_tolid=`tail -n +2 $bin_data | head -1 | cut -d',' -f1`
	first=`echo $host_tolid | cut -c1-1`
	group=`grep $'^'$first',' /lustre/scratch124/tol/projects/darwin/users/ng13/tolid_groups.csv \
		| cut -d',' -f2`
	project=`grep $'\t'$host_tolid$'\t' /lustre/scratch123/tol/tolqc/track/tol_sts.tsv \
		| head -1 | cut -f3`
	project=${project,,}
	echo $header_line > $outdir/bin_data.full.csv
	for bin in $bins
	do
		biosample=`grep $'^'$bin"," $biosample_accessions | cut -d',' -f3`
		tolid=`grep $'^'$bin"," $biosample_accessions | cut -d',' -f2`
		bin_line=`grep ",$bin," $bin_data`
		bin_type="binned metagenome"
		if [ `echo $bin_line | awk -F',' -v N=$quality_field '{print $N}'` == "HIGH" ] && [ `echo $bin_line | awk -F',' -v N=$drep_field '{print $N}'` == "PASSED" ]
		then
			bin_type="MAG"
		fi
		bin_line=`echo $bin_line | sed "s|^$host_tolid,|$host_tolid,$host_taxname,|"`
		echo "$bin_line,$primary_tolid,$primary_biosample,$tolid,$bin_type,$biosample,$group,$project" >> $outdir/bin_data.full.csv
	done
	format_bindata.py -p $outdir/bin_data.spreadsheet_update $outdir/bin_data.full.csv
	mags=`grep ',MAG,' $outdir/bin_data.spreadsheet_update.csv | wc -l`
	binned_metagenomes=`grep ',binned metagenome,' $outdir/bin_data.spreadsheet_update.csv | wc -l`
	echo -e "host_tolid,host_taxname,host_taxid,host_biosample,tolid,taxname,taxid,mags,binned_metagenomes,biosample,bioproject,assembly_accession,status" \
		> $outdir/primary.spreadsheet_update.csv
	echo -e "$tol_id,$host_taxname,$host_taxid,$host_biospecimen,$tol_id.metagenome,$primary_taxon,$primary_taxid,$mags,$binned_metagenomes,$primary_biosample,,,," \
		>> $outdir/primary.spreadsheet_update.csv
fi

######################################################################################
#### CHROMOSOME LIST ####
##########################
#
if [ "$command" == "generate_chromosome_list" ]
then
	. $biosample_configs
	bin_field=`head -1  $bin_data \
		| sed "s|,|\n|g" \
		| grep -n $'^bin_id$' \
		| cut -d':' -f1`
	quality_field=`head -1  $bin_data \
		| sed "s|,|\n|g" \
		| grep -n $'^quality$' \
		| cut -d':' -f1`
	drep_field=`head -1  $bin_data \
		| sed "s|,|\n|g" \
		| grep -n $'^drep$' \
		| cut -d':' -f1`
	mags=`tail -n +2 $bin_data \
		| awk -F',' -v B=$bin_field -v Q=$quality_field -v D=$drep_field '{if($Q == "HIGH" && $D == "PASSED") {print $B}}'`
	rm -rf $outdir/circ.tmp.list
	rm -rf $outdir/single.tmp.list
	for mag in $mags
	do
		count=`grep $'\t'$mag$'$' $directory/$binning_program/contigs2bin.tsv | wc -l`
		# grep $'\t'$mag$'$' $directory/$binning_program/contigs2bin.tsv \
		# 	| grep $'c\t' \
		# 	| cut -f1 \
		# 	>> $outdir/circ.tmp.list
		contigs=`grep $'\t'$mag$'$' $directory/$binning_program/contigs2bin.tsv \
			| cut -f1`
		for contig in $contigs
		do
			if ! `grep -q $contig$'\t' $contig_info`
			then
				echo "ERROR: $contig not found in $contig_info"
				exit 1
			fi
			if [ `grep $contig$'\t' $contig_info | cut -f4` == "Y" ]
			then
				echo $contig >> $outdir/circ.tmp.list
			fi
		done
		if [ $count -eq 1 ]
		then
			grep $'\t'$mag$'$' $directory/$binning_program/contigs2bin.tsv \
				| cut -f1 \
					>> $outdir/single.tmp.list
		fi
	done
	if test -f $outdir/circ.tmp.list && [ `cat $outdir/circ.tmp.list | wc -l` -gt 0 ]
	then
		if ! test -d $outdir/genomad_output || ! test -d $outdir/genomad_output/circular_contigs_summary || ! test -f $outdir/genomad_output/circular_contigs_summary/circular_contigs_plasmid_summary.tsv || [ "$force" == "true" ]
		then
			assembly_dir=`dirname $assembly`
			if ! test -f $assembly_dir/contigs.fasta
			then
				gunzip -c $assembly \
					> $assembly_dir/contigs.fasta
			fi
			select_fasta_by_list.pl \
				-i $assembly_dir/contigs.fasta \
				-l $outdir/circ.tmp.list \
				-o $outdir/circular_contigs.fa
			fastalength $outdir/circular_contigs.fa \
				> $outdir/circular_contigs.len
			singularity run -B /lustre,/nfs \
				$LOCAL_IMAGES/genomad.sif end-to-end \
					$outdir/circular_contigs.fa \
					$outdir/genomad_output \
					$GENOMAD_DB
		fi
	fi
	rm -rf $outdir/chromosome_list.tsv
	if test -f $outdir/circ.tmp.list
	then
		while read contig
		do
			len=`grep $' '$contig$'$' $outdir/circular_contigs.len | cut -d' ' -f1`
			mag=`grep $'^'$contig$'\t' $directory/$binning_program/contigs2bin.tsv | cut -f2`
			plasmid=`grep $contig$'\t' $outdir/genomad_output/circular_contigs_summary/circular_contigs_plasmid_summary.tsv \
				| awk -F'\t' '{if($6 >= .95){print "true"}else{print "false"}}'`
			if [ "$plasmid" == "true" ]
			then
				echo -e "$mag\t$contig\tCircular-Plasmid" >> $outdir/chromosome_list.tsv
			elif [ $len -ge 500000 ]
			then
				echo -e "$mag\t$contig\tCircular-Chromosome" >> $outdir/chromosome_list.tsv
			fi
		done < $outdir/circ.tmp.list
	fi
	if test -f $outdir/single.tmp.list
	then
		while read contig
		do
			if ! `grep -q $'\t'$contig$'\t' $outdir/chromosome_list.tsv`
			then
				mag=`grep $'^'$contig$'\t' $directory/$binning_program/contigs2bin.tsv | cut -f2`
				if `echo $contig | grep -q $'c$'`
				then
					echo -e "$mag\t$contig\tCircular-Chromosome" >> $outdir/chromosome_list.tsv
				else
					echo -e "$mag\t$contig\tLinear-Chromosome" >> $outdir/chromosome_list.tsv
				fi
			fi
		done < $outdir/single.tmp.list
	fi
fi

######################################################################################
#### CURATION REQUEST ####
##########################

if [ "$command" == "curation_request" ]
then
	. $biosample_configs
	if `echo $project | grep -iq 'darwin'`
	then
		outproject='Darwin'
	else
		outproject=`echo $project | tr '[a-z]' '[A-Z]'`
	fi
	# Primary metagenome
	if [ "$cr_outdir" == "" ]
	then
		current_day=`date +%Y%m%d`
		cr_outdir=/lustre/scratch124/tol/projects/$project/data/$group/$species/assembly/draft
	fi
	# if [ "$force" == "true" ]
	# then
	# 	rm -rf $cr_outdir
	# fi
	mainout=/lustre/scratch124/tol/projects/$project/data/$group/$species/assembly/draft
	pb_read_dir=/lustre/scratch124/tol/projects/$project/data/$group/$species/genomic_data/$tol_id/pacbio
	filesout=$cr_outdir/$tol_id.metagenome
	mkdir -p $filesout
	echo "Storing final output in $cr_outdir"
	if ! test -f $filesout/$tol_id.metagenome.fa.gz \
		|| ! test -f $filesout/$tol_id.metagenome.yaml \
		|| ! `grep -q "largest" $filesout/$tol_id.metagenome.yaml` \
		|| [ "$force" == "true" ]
	then
		echo "Prepping primary metagenome."
		if [ `echo $assembly | awk -F[.] '{print $NF}'` != "gz" ]
		then
			gzip -c $assembly \
				> $filesout/$tol_id.metagenome.fa.gz
			ln -fs $assembly $outdir/tmp.fa
		else
			cp $assembly $filesout/$tol_id.metagenome.fa.gz
			gunzip -c $filesout/$tol_id.metagenome.fa.gz \
				> $outdir/tmp.fa
		fi
	#################################################
	# YAML for Curation Request: Primary Metagenome #
		echo -e \
"---
species: \"$primary_taxon\"
specimen: $tol_id.metagenome
projects:
  - $outproject
data_location: Sanger RW
cobiont_status: cobiont
primary: $filesout/$tol_id.metagenome.fa.gz
assembly_type: primary metagenome
pacbio_read_dir: $pb_read_dir
jira_queue: DS
pipeline:
  - $asm_soft
stats: |" > $filesout/$tol_id.metagenome.yaml
		singularity exec -B /lustre:/lustre \
			/software/tola/images/asmstats.sif \
				asmstats $outdir/tmp.fa \
					| sed "s|^|  |g" | tail -n +2 \
					>> $filesout/$tol_id.metagenome.yaml
	fi
	# echo "Submitting curation request for \"$tol_id.metagenome\"."
	# metagenome_curation_request.sh $filesout/$tol_id.metagenome.yaml
	# ################################################
	# YAML for Curation Request: Binned Assemblies
	bin_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^bin_id$' | cut -d':' -f1`
	species_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^ncbi_taxon$' | cut -d':' -f1`
	quality_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^quality$' | cut -d':' -f1`
	drep_field=`head -1  $bin_data | sed "s|,|\n|g" | grep -n $'^drep$' | cut -d':' -f1`
	bins=`tail -n +2 $bin_data | awk -F',' -v N=$bin_field '{print $N}'`
	for bin in $bins
	do
		bin_file=`ls $bindir/output_bins | grep $bin.$'fa.*'`
		if [ "$bin_file" == "" ]
		then
			echo "ERROR: Could not find bin file in \'$bindir/output_bins\'."
			exit 1
		elif [ `head $bindir/output_bins/$bin_file | wc -l` -eq 0 ]
		then
			echo "ERROR: $binfile is empty."
			exit 1
		fi
		echo "Prepping $bin"
		rm -rf $outdir/bin.tmp.fa
		bin_tolid=`grep $'^'$bin"," $biosample_accessions | cut -d',' -f2`
		if [ `echo $bin_tolid | wc -c` -gt 47 ]
		then
			echo "ERROR: Sequence name length $bin_tolid exceeds limit (47 characters)"
			exit 1
		fi
		bin_species=`grep ",$bin," $bin_data | awk -F',' -v N=$species_field '{print $N}'`
		bin_quality=`grep ",$bin," $bin_data | awk -F',' -v N=$quality_field '{print $N}'`
		bin_drep=`grep ",$bin," $bin_data | awk -F',' -v N=$drep_field '{print $N}'`
		assembly_type="binned metagenome"
		bin_tolid_base=`echo $bin_tolid | cut -d'.' -f2`
		submission_dir="$cr_outdir/${bin_tolid}"
		mkdir -p $submission_dir
		if ! test -f $submission_dir/$bin_tolid.fa.gz \
			|| ! test -f $submission_dir/$bin_tolid.yaml \
			|| ! `grep -q "largest" $submission_dir/$bin_tolid.yaml` \
			|| [ "$force" == "true" ]
		then
			if [ "$bin_quality" == "HIGH" ] && [ "$bin_drep" == "PASSED" ]
			then
				assembly_type="Metagenome-Assembled Genome (MAG)"
			fi
			if [ `echo $bin_file | awk -F'.' '{print $NF}'` == "gz" ]
			then
				# cp $bindir/output_bins/$bin_file $submission_dir/$bin_tolid.fa.gz
				gunzip -c $bindir/output_bins/$bin_file \
					> $outdir/bin.tmp.fa
			else
				# gzip -c $bindir/output_bins/$bin_file > $submission_dir/$bin_tolid.fa.gz
				cp $bindir/output_bins/$bin_file $outdir/bin.tmp.fa
			fi
			contigs=`grep $'^>' $outdir/bin.tmp.fa | cut -d'>' -f2`
			if [ "$assembly_type" == "Metagenome-Assembled Genome (MAG)" ] && [ `echo $contigs | wc -w` -eq 1 ] && ( [ "$chromosome_list" == "" ] || ! `grep -q "$contigs" $chromosome_list` )
			then
				echo "ERROR: $bin is classified as MAG and has only one contig. A chromosome label is required for submission."
			fi
			contig_num=1
			chrom_num=1
			rm -rf $submission_dir/chromosome_list.tsv
			for contig in $contigs
			do
				if [ "$chromosome_list" != "" ] && `grep -q $'\t'$contig$'\t' $chromosome_list`
				then
					sed -i "s|>$contig$|>chromosome_$chrom_num contig=$contig|g" $outdir/bin.tmp.fa
					chrom_type=`grep $'\t'$contig$'\t' $chromosome_list | cut -f3`
					echo -e  "chromosome_$chrom_num\t$chrom_num\t$chrom_type" >> $submission_dir/chromosome_list.tsv
					chrom_num=`expr $chrom_num + 1`
				else
					sed -i "s|>$contig$|>contig_$contig_num contig=$contig|g" $outdir/bin.tmp.fa
					contig_num=`expr $contig_num + 1`
				fi
			done
			gzip -c $outdir/bin.tmp.fa > $submission_dir/$bin_tolid.fa.gz
			echo -e \
"---
species: \"$bin_species\"
specimen: $bin_tolid
projects:
  - $outproject
data_location: Sanger RW
cobiont_status: cobiont
primary: $submission_dir/$bin_tolid.fa.gz
assembly_type: $assembly_type
pacbio_read_dir: $pb_read_dir
jira_queue: DS
pipeline:
  - $asm_soft
stats: |" > $submission_dir/$bin_tolid.yaml
			singularity exec -B /lustre:/lustre \
				/software/tola/images/asmstats.sif \
					asmstats $outdir/bin.tmp.fa \
						| sed "s|^|  |g" | tail -n +2 \
						>> $submission_dir/$bin_tolid.yaml
			if test -f $submission_dir/chromosome_list.tsv
			then
				sed -i "s|Sanger RW$|Sanger RW\nchromosome_list: $submission_dir/chromosome_list.tsv|g" $submission_dir/$bin_tolid.yaml
			fi
		fi
	done
	echo "Output directory: $cr_outdir"
	if [ "$submit" == "true" ]
	then
		for magdir in $cr_outdir/*
		do
			bin_tolid=`basename $magdir`
			echo "Submitting curation request for \"$bin_tolid\"."
			if test -f $magdir/$bin_tolid.yaml
			then
				metagenome_curation_request.sh $magdir/$bin_tolid.yaml
			else
				echo "ERROR: Could not find yaml $magdir/$bin_tolid.yaml"
			fi
		done
	fi		
fi

if [ "$command" == "submit_cr" ]
then
	if [ "$cr_outdir" == "" ] || ! test -d $cr_outdir
	then
		echo "Could not find curation request directory: $cr_outdir"
		exit 1
	else
		for magdir in $cr_outdir/*
		do
			bin_tolid=`basename $magdir`
			if [ "$bin_tolid" != "blocklist.txt" ]
			then
				if ! test -f $cr_outdir/blocklist.txt || ( test -f $cr_outdir/blocklist.txt && ! `grep -q $'^'$bin_tolid$'$' $cr_outdir/blocklist.txt` )
				then
					if test -f $magdir/$bin_tolid.yaml
					then
						echo "Submitting curation request for \"$bin_tolid\"."
						sleep_time=`shuf -i20-100 -n1`
						sleep_time=300
						sleep $sleep_time
						metagenome_curation_request.sh $magdir/$bin_tolid.yaml
					else
						echo "ERROR: Could not find yaml $magdir/$bin_tolid.yaml"
					fi
				else
					echo "Skipping $bin_tolid"
				fi
			fi
		done
	fi
fi
# rm -rf $outdir/*tmp*

if [ "$command" == "create_btk" ]
then
	host_tolid=`tail -n +2 $bin_data | head -1 | cut -d',' -f1`
	curl -L "https://docs.google.com/spreadsheets/d/1RKubj10g13INd4W7alHkwcSVX_0CRvNq0-SRe21m-GM/export?gid=1641921323&format=csv" \
		> ~/primary_accession_data.csv
	grep $host_tolid ~/primary_accession_data.csv \
		> $outdir/primary_accession_data.csv
	bioproject_col=`head -1 ~/primary_accession_data.csv | sed "s|,|\n|g" | grep -n $'^bioproject$' | cut -d':' -f1`
	bioproject=`head -1 $outdir/primary_accession_data.csv | cut -d',' -f$bioproject_col`
	curl -X 'GET' \
	  "https://www.ebi.ac.uk/ena/portal/api/filereport?result=analysis&accession=$bioproject&fields=accession%2Canalysis_title&format=tsv" \
	  -H 'accept: */*' \
		  | tail -n +2 | awk '{OFS=","; print $4,$1}' \
		  > $outdir/accessions.csv
	curl -X 'GET' \
	  "https://www.ebi.ac.uk/ena/portal/api/filereport?result=assembly&accession=$bioproject&fields=accession%2Cassembly_title&format=tsv" \
	  -H 'accept: */*' \
		  | tail -n +2 | awk '{OFS=","; print $2,$1}' \
			  >> $outdir/accessions.csv
	if ! test -f $outdir/bin_data/other_bin_info.csv || [ `cat $outdir/bin_data/other_bin_info.csv | wc -l` -eq 0 ] || [ "$force" == "true" ]
	then
		tolid_col=`head -1 $bin_data | sed "s|,|\n|g" | grep -n $'^tolid$' | cut -d':' -f1`
		header=`head -1 $bin_data`
		echo "$header,assembly_accession" > $outdir/bindata.updated.csv
		while read line
		do
			if ! `echo $line | grep -q $'^host,'`
			then
				bin_tolid=`echo $line | cut -d',' -f$tolid_col`
				echo $bin_tolid
				asm_accession=`grep $'^'$bin_tolid$'\.' $outdir/accessions.csv \
					| cut -d',' -f2`
				if [ "$asm_accession" == "" ]
				then
					echo "ERROR: Could not find accession for $bin_tolid"
					# exit -1
					asm_accession="NOT_FOUND"
				fi
				echo "$line,$asm_accession" >> $outdir/bindata.updated.csv
			fi
		done < $bin_data				
		bin_data2csvs.sh \
			-a $assembly \
			-d $bindir \
			-R $bindir/output_bins \
			-b $outdir/bindata.updated.csv \
			-x $biosample_accessions \
			-C $contig_info \
			-o $outdir/bin_data
	fi
	if [ "$btk_dir" != "" ]
	then		
		mkdir -p $outdir/btk_dataset
		cp -R $btk_dir/* $outdir/btk_dataset
		singularity exec -B /lustre,/nfs \
			$LOCAL_IMAGES/blobtoolkit.sif \
				blobtools add \
					--text $outdir/bin_data/other_bin_info.csv \
					--text-delimiter ',' \
					--text-header \
					--replace \
					$outdir/btk_dataset \
					--text-no-array
		files=`ls $outdir/bin_data/*.csv \
			| grep -v 'ncbi' \
			| grep -v 'gtdb' \
			| grep -v "other_bin_info"`
		for file in $files
		do
			singularity exec -B /lustre,/nfs \
				$LOCAL_IMAGES/blobtoolkit.sif \
					blobtools add \
						--text $file \
						--text-delimiter ',' \
						--text-header \
						--replace \
						$outdir/btk_dataset \
						--text-no-array
		done
		groups="domain phylum class order family genus species"
		dbs="gtdb ncbi"
		for db in $dbs
		do
			for group in $groups
			do
				singularity exec -B /lustre,/nfs \
					$LOCAL_IMAGES/blobtoolkit.sif \
						blobtools add \
							--text $outdir/bin_data/${db}_$group.csv \
							--text-delimiter ',' \
							--text-header \
							--replace \
							$outdir/btk_dataset \
							--text-no-array
			done
		done
		singularity exec -B /lustre,/nfs \
			$LOCAL_IMAGES/blobtoolkit.sif \
				blobtools create \
					$outdir/btk_dataset
		primary_accession=`grep "$host_tolid.metagenome" $outdir/accessions.csv | cut -d',' -f2`
		primary_biosample_col=`head -1 ~/primary_accession_data.csv | sed "s|,|\n|g" | grep -n $'^biosample$' | cut -d':' -f1`
		primary_biosample=`grep ','$bioproject"," ~/primary_accession_data.csv | cut -d',' -f$primary_biosample_col`
		singularity exec -B /lustre,/nfs \
			$LOCAL_IMAGES/blobtoolkit.sif \
				blobtools replace \
					--key assembly.accession=$primary_accession \
					--key assembly.bioproject=$bioproject \
					--key assembly.biosample=$primary_biosample \
					--key assembly.alias="$host_tolid.metagenome.1" \
					--key assembly.file="$dir/assembly/draft/$host_tolid.metagenome.$date/$host_tolid.metagenome.fa.gz" \
					--key plot.cat="ncbi_family" \
					$outdir/btk_dataset
		cp -R $outdir/btk_dataset \
			/lustre/scratch123/tol/share/mg-btk/blobplots/$host_tolid.metagenome.$metagenome_version
		curl -s 'https://metagenomes-api.genomehubs.org/api/v1/search/reload/testkey%20npm%20start'
	fi
fi

	
