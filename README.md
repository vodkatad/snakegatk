# Installation

All references and annotation files can be obtained in a separate task from the repo https://github.com/vodkatad/snaketry, in `annotation/dataset/gnomad`.
Required packages can be obtained via dockerfiles (look in local/share/envs/dockerfiles), the branch occam has all snakefile rules with a new directive, 'docker' - a patched
snakemake is available at https://github.com/vodkatad/snakemake_docker which is able to call rules inside docker containers.
Otherwise you can use a standard snakemake and run all rules in a conda environment obtained with:

`conda config --add channels defaults && conda config --add channels bioconda && conda config --add channels conda-forge`

`conda create -n CONDA_ENV_NAME bcftools=1.9 bedtools=2.27 picard=2.18.15 samtools=1.9 trimmomatic=0.38  ucsc-liftover=357 bwa=0.7.17 fastqc=0.11.7 mosdepth=0.2.3 ensembl-vep=94.5 python=3.7.3 multiqc=1.7`

For GATK rules use --use-singularity, the singularity image of GATK version 4.1.4.0 is referenced with the singularity directive and can be obtained on your local
systems via:

`singularity pull --name gatk.img docker://broadinstitute/gatk:4.1.4.0`

# Setup

Create a new directory inside dataset and add a symbolic link to the appropriate Snakefile, eg:

`ln -s ../../local/share/snakerule/Snakefile_WES_base Snakefile`

then copy one of the example conf file, adapt it to your needs and link it there as conf.sk 
(look in dataset/example_* for minimal examples of setup).

### Major variables that needs to be adapted to your needs:
  - ROOT path where you put your task/annotation/dataset/gnomad
  - TYPE= "WES" | "WGS"
  - PAIRED=True | False
  - CORES=12 # how many cores you want to use

### Variables with docker images or path to singularity images:
 - GATK_SING="/home/egrassi/gatk4140/gatk.img"
 - GATK_ODOCKER="egrassi/occamsnakes/gatk:1"
 - XENOME_ODOCKER="egrassi/occamsnakes/xenome:1"
 - ALIGN_ODOCKER="egrassi/occamsnakes/align_recalibrate:1"

### Specific variables with info on your samples:
 - DATA=PRJ_ROOT+"/local/share/data/example" # where the fastq are stored
 - SAMPLES_ORIG=["example","examplen"] # names of the fastq files in the fastq_dir
 - SAMPLES=["pippo","pluto"] # samples name, same order as SAMPLES_ORIC
 - FASTQ_SUFFIX=".{pair}.fastq.gz" # structure of your fastq names
 - FASTQ_PRE="."
 - FASTQ_POST=".fastq.gz"
 - PAIRS=['1','2']

### Only for WES, the targeted regions:
 - EXONS=DATA+'/targed.bed'

### If you want to process some samples via xenome, their list of names:
 - XENOMED_SAMPLES=[]

**TODO** Xenome indices path are still hardcoded in rules and not configured via `conf.sk`!

*WARNING* due to a bug xenome does not end but hangs in multithreaded mode, here some steps need to be carried out manually changing rules targets (first xenome, kill
after some time manually, change targets to call the checkxenome rule).

### Only for paired setups sample names for each normal-tumor pair, in the same order
 - NORMAL=["pluto"]
 - TUMOR=["pippo"]

# Run
The all target will generate all vcf for your samples and some basic qc/multiqc files.