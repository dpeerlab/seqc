# Running Test

## Setup

Set the following environment variables:

```bash
export SEQC_TEST_RSA_KEY=/Users/chunj/dpeerlab-chunj.pem
export SEQC_TEST_EMAIL=jaeyoung.chun@gmail.com
export SEQC_TEST_AMI_ID=ami-07d15d7b3bb78506d
```

For local test, download test data in S3 to your test machine:

```
aws s3 sync s3://seqc-public/test/ten_x_v2/ ./test-data/datasets/ten_x_v2/
aws s3 sync s3://seqc-public/barcodes/ten_x_v2/ ./test-data/datasets/barcodes/ten_x_v2/
aws s3 sync s3://seqc-public/genomes/hg38_chr19/ ./test-data/datasets/genomes/hg38_chr19/
```

## SEQC index

```bash
nose2 -s src/seqc test_index
```

Besides the nose2 test results, actual SEQC output files can be found here, for example:

```
s3://dp-lab-cicd/seqc/index-ciona_intestinalis-0d19e818-7623-4a1d-bac3-a8c9e3be1e3e/
```

## SEQC run

### Local

SEQC will run with `--local`.

```bash
nose2 -s src/seqc test_run_e2e_local
```

### Remote

SEQC will run on AWS.

The following will generate a package that can be uploaded to AWS EC2 for testing:

```bash
python repackage.py
```

```bash
nose2 -s src/seqc test_run_e2e_remote
```

Besides the nose2 test results, actual SEQC output files can be found here, for example:

```
s3://dp-lab-cicd/seqc/run-in_drop_v2-a997b408-f883-4ba2-9941-8b541e319850/
```

### Clean Up

```bash
aws s3 rm s3://dp-lab-cicd/seqc/ --recursive
```