# docs

## Developers

- [Environment setup for development](./install-dev.md)
- [Running test](./run-test.md)


## Generating Reference Packages

This generates a reference package (STAR index and GTF) using SEQC v0.2.8.

- Ensembl 86
- Gene annotation file that contains only the reference chromosomes (no scaffolds, no patches)
- Only these biotypes: 'protein_coding', 'lincRNA', 'IG_V_gene', 'IG_C_gene', 'IG_J_gene', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_D_gene', 'IG_D_gene'
- Not passing anything to `--additional-id-types`
- Setting the read length to 101 (internally, this becomes 100)

### Local

```bash
SEQC index \
    -o homo_sapiens \
    -f homo_sapiens \
    --ensemble-release 93 \
    --valid-biotypes protein_coding lincRNA antisense IG_V_gene IG_D_gene IG_J_gene IG_C_gene TR_V_gene TR_D_gene TR_J_gene TR_C_gene \
    --read-length 101 \
    --folder ./test-data/index/ \
    --local
```

### AWS

```bash
SEQC index \
    -o homo_sapiens \
    -f homo_sapiens \
    --ensemble-release 93 \
    --valid-biotypes protein_coding lincRNA antisense IG_V_gene IG_D_gene IG_J_gene IG_C_gene TR_V_gene TR_D_gene TR_J_gene TR_C_gene \
    --read-length 101 \
    --upload-prefix s3://dp-lab-test/seqc/index/86/ \
    --rsa-key ~/dpeerlab-chunj.pem \
    --ami-id ami-037cc8c1417e197c1
```
