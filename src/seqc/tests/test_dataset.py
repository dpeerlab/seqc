from collections import namedtuple

TestDataset = namedtuple(
    "datasets",
    ["barcode_fastq", "genomic_fastq", "merged_fastq", "bam", "index", "barcodes",],
)

dataset_s3 = TestDataset(
    barcode_fastq="s3://seqc-public/test/%s/barcode/",  # platform
    genomic_fastq="s3://seqc-public/test/%s/genomic/",  # platform
    merged_fastq="s3://seqc-public/test/%s/%s_merged.fastq.gz",  # platform, platform
    bam="s3://seqc-public/test/%s/Aligned.out.bam",  # platform
    index="s3://seqc-public/genomes/hg38_chr19/",
    barcodes="s3://seqc-public/barcodes/%s/flat/",  # platform
)

dataset_local = TestDataset(
    barcode_fastq="test-data/datasets/%s/barcode/",  # platform
    genomic_fastq="test-data/datasets/%s/genomic/",  # platform
    merged_fastq=None,
    bam="test-data/datasets/%s/Aligned.out.bam",  # platform
    index="test-data/datasets/genomes/hg38_chr19/",
    barcodes="test-data/datasets/barcodes/%s/flat/",  # platform
)
