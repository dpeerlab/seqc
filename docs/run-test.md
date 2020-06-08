# Running Test

## Setup

Set the following environment variables:

```bash
export SEQC_TEST_RSA_KEY=/Users/chunj/dpeerlab-chunj.pem
export SEQC_TEST_EMAIL=jaeyoung.chun@gmail.com
export SEQC_TEST_AMI_ID=ami-07d15d7b3bb78506d
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

This will generate a package that can be uploaded to EC2 for testing:

```bash
python repackage.py
```

```bash
nose2 -s src/seqc test_run_e2e.TestRunRemote
```

Besides the nose2 test results, actual SEQC output files can be found here, for example:

```
s3://dp-lab-cicd/seqc/run-in_drop_v2-a997b408-f883-4ba2-9941-8b541e319850/
```
