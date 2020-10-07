# SEquence Quality Control (SEQC -- /sek-si:/)

## Overview:

SEQC is a python package that processes single-cell sequencing data in the cloud and analyzes it interactively on your local machine.

To faciliate easy installation and use, we have made available Amazon Machine Images (AMIs) that come with all of SEQC's dependencies pre-installed. In addition, we have uploaded common genome indices (`-i/--index parameter`) and barcode data (`--barcode-files`) to public Amazon S3 repositories. These links can be provided to SEQC and it will automatically fetch them prior to initiating an analysis run. Finally, it can fetch input data directly from BaseSpace or amazon s3 for analysis.

For users with access to in-house compute clusters, SEQC can be installed on your systems and run using the `--local` parameter.

## Dependencies:

### Python 3

Python3 must be installed on your local machine to run SEQC. We recommend installing Python3 through Miniconda (https://docs.conda.io/en/latest/miniconda.html).

### Python 3 Libraries

 We recommend creating a virtual environment before installing anything:

```bash
conda create -n seqc python=3.7.7 pip
conda activate seqc
```

```bash
pip install Cython
pip install numpy
pip install bhtsne
```

### STAR, Samtools, and HDF5

To process data locally using SEQC, you must install the <a href=https://github.com/alexdobin/STAR>STAR Aligner</a>, <a href=http://www.htslib.org/>Samtools</a>, and <a href=https://support.hdfgroup.org/HDF5/>hdf5</a>. If you only intend to use SEQC to trigger remote processing on AWS, these dependencies are optional. We recommend installing samtools and hdf5 through your package manager, if possible.

## SEQC Installation

Once all dependencies have been installed, SEQC can be installed on any machine by running:

```bash
export SEQC_VERSION="0.2.6"
wget https://github.com/hisplan/seqc/archive/v${SEQC_VERSION}.tar.gz
tar xvzf v${SEQC_VERSION}.tar.gz
cd seqc-${SEQC_VERSION}
pip install .
```

## Hardware Requirements:

For processing a single lane (~200M reads) against human- and mouse-scale genomes, SEQC requires 30GB RAM, approximately 200GB free hard drive space, and scales linearly with additional compute cores. If running on AWS (see below), jobs are automatically scaled up or down according to the size of the input. There are no hardware requirements for the computer used to launch remote instances.

## Running SEQC on Local Machine:

Download an example dataset (1k PBMCs from a healthy donor; freely available at 10x Genomics https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3):

```bash
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar xvf pbmc_1k_v3_fastqs.tar
```

Move R1 FASTQ files to the `barcode` folder and R2 FASTQ files to the `genomic` folder:

```bash
mkdir barcode
mkdir genomic
mv ./pbmc_1k_v3_fastqs/*R1*.fastq.gz barcode
mv ./pbmc_1k_v3_fastqs/*R2*.fastq.gz genomic/
```

Download the 10x barcode whitelist file:

```bash
mkdir whitelist
wget https://seqc-public.s3.amazonaws.com/barcodes/ten_x_v3/flat/3M-february-2018.txt
mv 3M-february-2018.txt ./whitelist
```

The resulting directory structure should look something like this:

```
.
├── barcode
│   ├── pbmc_1k_v3_S1_L001_R1_001.fastq.gz
│   └── pbmc_1k_v3_S1_L002_R1_001.fastq.gz
├── genomic
│   ├── pbmc_1k_v3_S1_L001_R2_001.fastq.gz
│   └── pbmc_1k_v3_S1_L002_R2_001.fastq.gz
├── pbmc_1k_v3_fastqs
│   ├── pbmc_1k_v3_S1_L001_I1_001.fastq.gz
│   └── pbmc_1k_v3_S1_L002_I1_001.fastq.gz
├── pbmc_1k_v3_fastqs.tar
└── whitelist
    └── 3M-february-2018.txt
```

Create a reference package (STAR index + gene annotation):

```bash
SEQC index \
  --organism homo_sapiens \
  --ensemble-release 93 \
  --valid-biotypes protein_coding lincRNA antisense IG_V_gene IG_D_gene IG_J_gene IG_C_gene TR_V_gene TR_D_gene TR_J_gene TR_C_gene \
  --read-length 101 \
  --folder index \
  --local
```

Run SEQC:

```bash
export AWS_DEFAULT_REGION=us-east-1
export SEQC_MAX_WORKERS=7

SEQC run ten_x_v3 \
  --index ./index/ \
  --barcode-files ./whitelist/ \
  --barcode-fastq ./barcode/ \
  --genomic-fastq ./genomic/ \
  --upload-prefix ./seqc-results/ \
  --output-prefix PBMC \
  --no-filter-low-coverage \
  --min-poly-t 0 \
  --star-args runRNGseed=0 \
  --local
```

## Running SEQC on Amazon Web Services:

SEQC can be run on any unix-based operating system, however it also features the ability to automatically spawn Amazon Web Services instances to process your data.

1. <a href=http://aws.amazon.com>Set up an AWS account</a>
2. <a href=https://aws.amazon.com/cli/>Install and configure AWS CLI</a>
3. <a href=http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html>Create and upload an rsa-key for AWS</a>

Run SEQC:

```bash
SEQC run ten_x_v2 \
  --ami-id ami-08652ee2477761403 \
  --user-tags Job:Test,Project:PBMC-Test,Sample:pbmc_1k_v3 \
  --index s3://seqc-public/genomes/hg38_long_polya/ \
  --barcode-files s3://seqc-public/barcodes/ten_x_v2/flat/ \
  --genomic-fastq s3://.../genomic/ \
  --barcode-fastq s3://.../barcode/ \
  --upload-prefix s3://.../seqc-results/ \
  --output-prefix PBMC \
  --no-filter-low-coverage \
  --min-poly-t 0 \
  --star-args runRNGseed=0
```
