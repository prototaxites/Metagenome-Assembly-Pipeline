# Metagenome Assembly Pipeline
## Overview
Running the pipeline follow the following path...

### 1. Setting up the environment
Currently, all requisite executables should be located in my `bin` directory or can be run \
via Singularity images stored locally. 

```bash
export PATH=$PATH:/software/team311/ng13/local/bin
export LOCAL_IMAGES=/software/team311/ng13/local/images

module load samtools/1.20--h50ea8bc_0
```

We'll also set some basic variables (largely ToL-related) including our output directory...
```bash
project=<PROJECT> # e.g ASG
tolid_group=<TOLID_GROUP> # e.g. chordates
species=<SPECIES> # e.g. Trididemnum_miniatum
tolid=<TOLID> # kaTriMin1
assembler=<ASSEMBLER> # e.g. meta-mdbg
date=<DATE>
dir=/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species
outdir=/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date
```

NOTE: The pipeline is primarily designed to work with ToL Farm directory structures \
(`/lustre/scratch124/tol/projects/<PROJECT>/data/<TOLID_GROUP>/<SPECIES>/working/<TOLID>.<ASSEMBLER>.<DATE>`) \
For example, `/lustre/scratch124/tol/projects/asg/data/sponges/Sycon_ciliatum/working/ocSycCili1.meta-mdbg.20231108`.


Create output directory and combine read files
```bash
mkdir -p $outdir

for file in $dir/genomic_data/$tolid/pacbio/fasta/*.fasta.gz
do
	echo $file >> $outdir/$tolid.pacbio.fa.fofn
	gunzip -c $file >> $outdir/$tolid.pacbio.fa.gz
done
```

If HiC exists and is from the same sample, we need to convert from cram to fastq
```bash
for cram in $dir/genomic_data/$tolid/hic-arima2/*.cram
do
	bsub -n14 -R"span[hosts=1]" \
		-o samtools_fastq.%J.o \
		-e samtools_fastq.%J.e \
		-M10000 \
		-R 'select[mem>10000] rusage[mem=10000]' \
			"samtools fastq -@14 -n -f0x40 -F0xB00 $cram \
			| gzip -c - >> $outdir/$tolid.hic.1.fastq.gz"
	bsub -n14 -q long -R"span[hosts=1]" \
		-o samtools_fastq.%J.o \
		-e samtools_fastq.%J.e \
		-M10000 \
		-R 'select[mem>10000] rusage[mem=10000]' \
			"samtools fastq -@14 -n -f0x80 -F0xB00 $cram \
			| gzip -c - >> $outdir/$tolid.hic.2.fastq.gz"
done
```

### 2. Metagenome assembly

The pipeline should be compatible with MetaFlye, Hifiasm-Meta, and metaMDBG. \
In some basic benchmarking, metaMDBG consistently outperformed the other tools so we'll \
stick with that for now...


Run metaMDBG. 
```bash
outdir=/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date
cd $outdir
metaMDBG asm \
	--out-dir $outdir \
	--in-hifi $outdir/$tolid.pacbio.fa
```
Or as a job (metaMDBG is not super memory intensive relative to other assemblers like Hifiasm. \
Adjust accordingly)
```bash
bsub -n12 -q long -R"span[hosts=1]" \
	-o meta-mdbg.%J.o \
	-e meta-mdbg.%J.e \
	-M30000 \
	-R 'select[mem>30000] rusage[mem=30000]' \
		"metaMDBG asm \
			--threads 12 \
			--out-dir $outdir \
			--in-hifi $outdir/$specimen.pacbio.fa"
```

Output...
```
$outdir/contigs.fasta.gz
```


### 3. Run the pipeline
Run binning pipeline via the runner script. The runner script takes as input a file of \
assembly directories (one per line) in the format...
`<ASSEMBLER>	<ASSEMBLY DIRECTORY>`


Set up output directory and text file with assembly directories
```bash
outdir=/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/map_output
mkdir -p $outdir
echo -e "$assembler\t/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date" > $outdir/asm_dirs.txt
reads=/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/$tolid.pacbio.fa
hic1=/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/$tolid.hic.1.fastq.gz
hic2=/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/$tolid.hic.2.fastq.gz
```

Run runner script.
```bash
map_runner.sh \
	-D $outdir/asm_dirs.txt \
	-p $tolid \
	-A all \
	-r  \
	-R dastool,magscot \
	-O $outdir \
	-c $hic1,$hic2 \
	-B metabat2,maxbin2,bin3c
```
NOTE: `bin3c` requires HiC data.


Output..
```bash
/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/mags/all/bin_stats.csv
```

To check if all of the outputs were generated...
```bash
cut -d',' -f1-3 /lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/mags/all/bin_stats.csv \
	| sort -u
```


### 3. Score the output
```bash
score_mag_outputs.py \
	-p $outdir/final \
	/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/mags/all/bin_stats.csv
```

This will produce two files:
- final.scored_bindata.csv 		<- The output bin data with scores and score-breakdowns for each bin
- final_counts.csv 				<- A file containing a summary of each pipeline that was run including overall score

The `final_counts.csv` output gives counts of the number of bins of varying quality, the
number of bins in each category that represent a unique taxon as determined by the
program `drep`,the number of total MAGs (of all quality) and the number of SSUs (should 
roughly match the MAG count) found by each pipeline.


### 4. Post-binning processing
Once you have chosen the best pipeline, we need to check the bins and associate the chosen 
output with the metadata required for submission to ENA.

First create some variables to point to the output of the pipeline you've chosen.

```bash
mag_dir=/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/mags/all
binning_program=magscot_raw
```

#### A. BioSample Config Setup
We can now set up the metadata associated with our sample that will be attached to our
BioSample.

```bash
metagenome_submission.sh config_file \
	-d $mag_dir \
	-B $binning_program \
	-c $outdir/biosample_configs.txt
```

This script will do its best to generate the correct information for each required 
BioSample field in the file `$outdir/biosample_configs.txt`, however, editing the 
information is often required so make sure to **check the output file before proceeding**.

The BioSample will require a TaxID for the metagenome as a whole. These can be found 
[here](https://www.ebi.ac.uk/ena/browser/view/410656?dataType=TAXON&show=tax-tree). NCBI 
suggests that they will only add new categories if absolutely necessary so pick the closest 
thing and know that, given the current structure of stored metagenomic data, finding these 
data will require special insider knowledge of its existence, regardless.

For whatever reason, the BioSample associated with the metagenomic data also requires 
recording two levels of the biome specificity with which it the sample is associated. The
list of these can be founnd [here](http://purl.obolibrary.org/obo/ENVO_00000428). As with 
the TaxIDs for metagenomes, these are poorly designed so I'd recommend just finding two 
obvious expansive/general biome matches (e.g. aquatic biome -- marine biome) and spending 
your time on something that actually matters).

#### B. Check Bin Taxon Info
Next, we pull out just the bin info associated with our pipeline of choice...

```bash
metagenome_submission.sh check \
	-d $mag_dir \
	-B $binning_program \
	-c $outdir/biosample_configs.txt \
	-b $outdir/final.scored_bindata.csv
```

The output of this command is `output_data.csv`, which is just the relevant subset of
$outdir/final.scored_bindata.csv.

The taxonomic assignments associated with each bin are based on those assigned by GTDB-TK 
according to information in the GTDB[https://gtdb.ecogenomic.org/]. This database follows 
what is arguably the most up-to-date taxonomic classifications available in the world of 
prokaryotic research. NCBI is much slower to adopt change, however, but as we are depositing 
our data to ENA it must follow NCBI's taxonomic system. The pipeline will attempt to 
translate GTDB taxonomic classifications into the most detailed NCBI classification. When 
clade names don't align the pipeline will pick the most specific NCBI lineage details to 
assign to a bin. Occasionally a more specific NCBI classification can be made despite a lack 
of cross database correspondence, often due to disagreements in naming conventions.

Compare the GTDB and NCBI lineage details in the file `output_data.csv` and update the NCBI 
lineage where necessary. It is also recommended that the these details be added to the file 
`/lustre/scratch124/tol/projects/darwin/users/ng13/gtdb-tk/ncbi_altkey.csv` which contains 
alternate keys for the GTDB to NCBI translation script to look at when unable to find clear 
correspondance. Columns for this file are as follows:

`<GTDB TAXON>,<CORRESPONDING NCBI TAXON ASSIGNMENT>,<CORRESPONDING NCBI TAXONID>`

Note: The translation script will only refer to this as a last resort.

---
**A note on metagenome bin taxonomy assignments**
GTDB-TK will not be able to classify most bins to the species level. Ambiguous classifications 
should accord to the following rules...

- Genus level assignments: <Genus> sp. (e.g. _Nitrosopumilus sp._)
- Non-genus level assignments: <Taxon> [bacterium,cyanobacterium,archaeon]

To check if a taxon assignment exists and if you can submit under that name, you can
check by going to 
https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/<TAXON NAME>

For example,
[https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/Candidatus Woesearchaeota archaeon](https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/Candidatus Woesearchaeota archaeon)

_NOTE: In the past, ENA recommended that non-specific taxa from metagenomes_ 
_be assigned to an 'uncultured' taxon (e.g. uncultured _Nitrosopumilus sp._). This no_ 
_longer appears to be the recommendation although such designations are still okay._

---


#### C. Create BioSamples
The final step before submitting our metagenome assembly is to create BioSample Accessions 
for our _primary metagenome_, our _binned metagenome assemblies_, and our 
_metagenome assembled genomes (MAGS)_.

---
**A note on metagenomic assembly/bin categorization**
There are 3 types of objects that comprise a metagenomic assembly submission.
([_see ENA guidelines_](https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome.html)). I've adjusted some submission criteria to reflect the higher 
quality nature of PacBio HiFi-based assemblies:

[_Primary metagenome_](https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome/primary.html): The raw assembly of all the reads, unbinned and unfiltered. The 
bins/mags are derived directly from this assembly

[_Binned metagenome assemblies_](https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome/binned.html): Medium quality and duplicate bins (bins of an OTU such as differen strains with a higher quality representative bin present)
	Critera
		- Contamination ≤ 10%
		- Completeness ≥ 50%

[_Metagenome-assembled genomes (MAGs)_](https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome/mag.html): High quality bins that represent a unique set of OTUs from the sample. 
These are meant to represent full genomes.
	Critera
		- 5s, 16s, & 23s rRNAs present
		- ≥ 18 unique tRNAs present
		- Contamination ≤ 5%
		- Completeness ≥ 90%
				OR
		  Completeness ≥ 50% and all contigs circular
	
---

We start by setting up all the necessary data to create the BioSamples.

```bash
metagenome_submission.sh biosample_prep \
	-d $mag_dir \
	-B $binning_program \
	-c $outdir/biosample_configs.txt \
	-b $outdir/output_data.csv
```

This will create three files:
 - `primary_biosample.csv`
 - `binned_biosample_metadata.csv`
 - `mag_biosample_metadata.csv`

We then do a test BioSample creation run to make sure everything is in order.

```
generate_metagenome_biosampleId.py \
	-a /nfs/users/nfs_n/ng13/mag_submission_credentials.dev.json \
	-p asg \
	-d $outdir/primary_biosample.csv \
	-o $outdir/${tolid}_biosamples.csv
```

If everything runs okay, you should now have a csv, `${tolid}_biosamples.csv`, that lists 
each metagenome object and a dummy BioSample accession.

We can now create real BioSamples...

```bash
generate_metagenome_biosampleId.py \
	-a /nfs/users/nfs_n/ng13/mag_submission_credentials.json \
	-p asg \
	-d $outdir/primary_biosample.csv \
	-o $outdir/${tolid}_biosamples.csv
```

#### D. Circular contig classification
A recent addition to this pipeline that has not been integrated well yet is the 
classification of circular contigs as plasmids/chromosomes.

```bash
metagenome_submission.sh generate_chromosome_list \
	-B $binning_program \
	-d $mag_dir \
	-c $outdir/biosample_configs.txt \
	-b $outdir/output_data.csv \
	-x $outdir/${specimen}_biosamples.csv
```

You may need to queue this command as it will run [geNomad](https://github.com/apcamargo/genomad) on all circular contigs.

The output from this command with be a file of circular contig annotations, `chromosome_list.tsv`.


#### E. GRIT Ticket creation
While metagenomic assembly objects do not undergo any of the standard curation steps a 
typical eukaryotic assembly would, we still submit to 'GRIT' (really one of James's pipelines) 
in a similar manner i.e. each object will get its own YAML and be assigned a ticket.

The pipeline for processing these objects requires that all relevant information be present 
and up-to-date in the [Tree of Life assembly informatics spreadsheet](https://docs.google.com/spreadsheets/d/1RKubj10g13INd4W7alHkwcSVX_0CRvNq0-SRe21m-GM/edit?gid=1442224132#gid=1442224132).

We can create the required entry lines by running the following...

```bash
metagenome_submission.sh spreadsheet_update \
	-c $outdir/biosample_configs.txt \
	-B $binning_program \
	-d $mag_dir \
	-x $outdir/${tolid}_biosamples.csv \
	-b $outdir/output_data.csv
```

This will create two csv's.
 - `primary.spreadsheet_update.csv` : The info from this file can be copy and pasted into the **Primary Metagenome Submission** tab.
 - `bin_data.spreadsheet_update.csv` : The info from this file can be copy and pasted into the **Binned Metagenome Submission** tab.
 
Before the data can be successfully submitted to ENA, a metagenome BioProject will 
need to be created. Currently, Shane is the only one who can do this so ask him to create 
this before putting in a ticket request.

When ready, we can create 'Curation Request' tickets for all the samples by running the 
following...

```
metagenome_submission.sh curation_request \
	-B $binning_program \
	-d $mag_dir \
	-c $outdir/biosample_configs.txt \
	-b $outdir/output_data.csv \
	-x $outdir/${specimen}_biosamples.csv \
	-l $outdir/chromosome_list.tsv
```

Each metagenome object will get its own subdirectory in the assembly draft directory for 
that sample (`$dir/assembly/draft`).

