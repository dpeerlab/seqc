# Running Test

## Setup

Set the following environment variables:

```bash
export SEQC_TEST_RSA_KEY=/Users/chunj/dpeerlab-chunj.pem
export SEQC_TEST_EMAIL=jaeyoung.chun@gmail.com
export SEQC_TEST_AMI_ID=ami-0a4d2955fe21dee72
```

## SEQC index

```bash
nose2 -s src/seqc test.TestIndex
```
