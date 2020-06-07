import os
import unittest
import gzip
import pandas as pd
import numpy as np
import ftplib
import nose2
from nose2.tools import params
from seqc.sequence import index, gtf
from seqc.core import main
from seqc.read_array import ReadArray
import seqc
import logging
import shutil


logging.basicConfig()

seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'

# fill and uncomment these variables to avoid having to provide input to tests
TEST_BUCKET = "dp-lab-test/seqc"  # None
EMAIL = os.environ["SEQC_TEST_EMAIL"]
RSA_KEY = os.environ["SEQC_TEST_RSA_KEY"]
AMI_ID = os.environ["SEQC_TEST_AMI_ID"]

# define some constants for testing
BARCODE_FASTQ = 's3://seqc-public/test/%s/barcode/'  # platform
GENOMIC_FASTQ = 's3://seqc-public/test/%s/genomic/'  # platform
MERGED = 's3://seqc-public/test/%s/%s_merged.fastq.gz'  # platform, platform
SAMFILE = 's3://seqc-public/test/%s/Aligned.out.bam'  # platform
INDEX = 's3://seqc-public/genomes/hg38_chr19/'
LOCAL_OUTPUT = os.environ['TMPDIR'] + 'seqc/%s/test'  # test_name
UPLOAD = 's3://%s/%s/'  # bucket_name, test_folder
PLATFORM_BARCODES = 's3://seqc-public/barcodes/%s/flat/'  # platform


def makedirs(output):
    """make directories based on OUTPUT"""
    stem, _ = os.path.split(output)
    os.makedirs(stem, exist_ok=True)


class TestSEQCLocal(unittest.TestCase):

    email = globals()['EMAIL'] if 'EMAIL' in globals() else None
    bucket = globals()['TEST_BUCKET'] if 'TEST_BUCKET' in globals() else None
    rsa_key = globals()['RSA_KEY'] if 'RSA_KEY' in globals() else None

    # @params('in_drop', 'in_drop_v2', 'drop_seq', 'ten_x', 'mars_seq')
    def test_local(self, platform='in_drop_v2'):
        """test seqc after pre-downloading all files"""
        with open('seqc_log.txt', 'w') as f:
            f.write('Dummy log. nose2 captures input, so no log is produced. This causes pipeline errors.\n')
        test_name = 'test_no_aws_%s' % platform
        makedirs(LOCAL_OUTPUT % test_name)
        if self.email is None:
            self.email = input('please provide an email address for SEQC to mail results: ')
        argv = [
            'run',
            platform,
            '-o', test_name,
            '-i', INDEX,
            '-b', BARCODE_FASTQ % platform,
            '-g', GENOMIC_FASTQ % platform,
            '--barcode-files', PLATFORM_BARCODES % platform,
            '-e', self.email,
            '--local']
        main.main(argv)
        # os.remove('./seqc_log.txt')  # clean up the dummy log we made.


class TestSEQCRemote(unittest.TestCase):

    email = globals()['EMAIL'] if 'EMAIL' in globals() else None
    bucket = globals()['TEST_BUCKET'] if 'TEST_BUCKET' in globals() else None
    rsa_key = globals()['RSA_KEY'] if 'RSA_KEY' in globals() else None
    ami_id = globals()['AMI_ID'] if 'AMI_ID' in globals() else None

    # @params('in_drop', 'in_drop_v2', 'drop_seq', 'ten_x', 'mars_seq')
    def test_remote_from_raw_fastq(self, platform='ten_x_v2'):
        test_name = 'test_remote_%s' % platform
        argv = [
            'run',
            platform,
            '-o', 'from-raw-fastq',
            '-u', UPLOAD % (self.bucket, test_name),
            '-i', INDEX,
            '-e', self.email,
            '-b', BARCODE_FASTQ % platform,
            '-g', GENOMIC_FASTQ % platform,
            '--instance-type', 'c4.large',
            '--spot-bid', '1.0',
            '-k', self.rsa_key,
            '--debug',
            '--ami-id', self.ami_id
        ]
        if platform != 'drop_seq':
            argv += ['--barcode-files', PLATFORM_BARCODES % platform]
        main.main(argv)

    # @params('in_drop', 'in_drop_v2', 'drop_seq', 'ten_x', 'mars_seq')
    def test_remote_from_merged(self, platform='in_drop_v2'):
        test_name = 'test_remote_%s' % platform
        argv = [
            'run',
            platform,
            '-o', 'from-merged',
            '-u', UPLOAD % (self.bucket, test_name),
            '-i', INDEX,
            '-e', self.email,
            '-m', MERGED % (platform, platform),
            '-k', self.rsa_key,
            '--instance-type', 'c4.large',
            '--ami-id', self.ami_id
            # '--spot-bid', '1.0'
        ]
        if platform != 'drop_seq':
            argv += ['--barcode-files', PLATFORM_BARCODES % platform]
        main.main(argv)

    # @params('in_drop', 'in_drop_v2', 'drop_seq', 'ten_x', 'mars_seq')
    def test_remote_from_samfile(self, platform='in_drop_v2'):
        test_name = 'test_remote_%s' % platform
        argv = [
            'run',
            platform,
            '-o', 'from-samfile',
            '-u', UPLOAD % (self.bucket, test_name),
            '-i', INDEX,
            '-e', self.email,
            '-a', SAMFILE % platform,
            '-k', self.rsa_key,
            '--instance-type', 'r5.2xlarge',
            '--debug',
            '--ami-id', self.ami_id
            # '--spot-bid', '1.0'
        ]
        if platform != 'drop_seq':
            argv += ['--barcode-files', PLATFORM_BARCODES % platform]
        main.main(argv)

class TestReadArrayCreation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # cls.bamfile = LOCAL_OUTPUT % platform + '_bamfile.bam'
        # cls.annotation = LOCAL_OUTPUT % platform + '_annotations.gtf'
        # S3.download(SAMFILE % platform, cls.bamfile, recursive=False)
        # S3.download(INDEX + 'annotations.gtf', cls.annotation, recursive=False)

        cls.bamfile = os.path.expanduser('~/Downloads/mm_test_short.bam')
        cls.annotation = os.path.expanduser('~/Downloads/annotations.gtf')
        cls.summary = os.path.expanduser('~/Downloads/mm_test_summary.txt')
        cls.total_input_reads = 12242659
        cls.translator = gtf.GeneIntervals(cls.annotation, 10000)

    def test_read_array_creation(self):
        ra = ReadArray.from_alignment_file(self.bamfile, self.translator,
                                           required_poly_t=0)
        print(repr(ra.data.shape[0]))
        print(repr(ra.genes))
        ra.save(os.path.expanduser('~/Downloads/test.ra'))


class TestTranslator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.annotation = os.path.expanduser('~/Downloads/annotations.gtf')

    def test_construct_translator(self):
        translator = gtf.GeneIntervals(self.annotation)
        print(len(translator._chromosomes_to_genes))

    def get_length_of_gtf(self):
        rd = gtf.Reader(self.annotation)
        # print(len(rd))
        print(sum(1 for _ in rd.iter_transcripts()))


#########################################################################################

if __name__ == "__main__":
    nose2.main()

#########################################################################################
