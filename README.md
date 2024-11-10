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
```
score_mag_outputs.py \
	-p $outdir/final \
	/lustre/scratch124/tol/projects/$project/data/$tolid_group/$species/working/$tolid.$assembler.$date/mags/all/bin_stats.csv
```





