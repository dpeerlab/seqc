import unittest
import os
import uuid
import shutil
import subprocess
import re
from nose2.tools import params
from seqc.core import main
from test_dataset import dataset_local, dataset_s3


def get_output_file_list(test_id, test_folder):

    proc = subprocess.Popen(
        ["find", test_folder, "-type", "f"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, _ = proc.communicate()
    files = stdout.decode().splitlines()

    # extract only filenames (i.e. remove directory hierarchy)
    # convert to a set for easy comparison
    files = set(map(lambda filename: filename.replace(test_folder + "/", ""), files))

    return files


def expected_output_files(file_prefix):

    files = set(
        [
            f"{file_prefix}.h5",
            f"{file_prefix}_alignment_summary.txt",
            f"{file_prefix}_cell_filters.png",
            f"{file_prefix}_de_gene_list.txt",
            f"{file_prefix}_dense.csv",
            f"{file_prefix}_merged.fastq.gz",
            f"{file_prefix}_mini_summary.json",
            f"{file_prefix}_mini_summary.pdf",
            f"{file_prefix}_seqc_log.txt",
            f"{file_prefix}_sparse_counts_barcodes.csv",
            f"{file_prefix}_sparse_counts_genes.csv",
            f"{file_prefix}_sparse_molecule_counts.mtx",
            f"{file_prefix}_sparse_read_counts.mtx",
            f"{file_prefix}_summary.tar.gz",
            f"{file_prefix}_Aligned.out.bam",
        ]
    )

    return files


class TestRunLocal(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.test_id = str(uuid.uuid4())
        cls.path_temp = os.path.join(
            os.environ["TMPDIR"], "seqc-test", str(uuid.uuid4())
        )
        os.makedirs(cls.path_temp, exist_ok=True)
        with open("seqc_log.txt", "wt") as f:
            f.write("Dummy log.\n")
            f.write("nose2 captures input, so no log is produced.\n")
            f.write("This causes pipeline errors.\n")

    @classmethod
    def tearDown(self):
        if os.path.isdir(self.path_temp):
            shutil.rmtree(self.path_temp, ignore_errors=True)

    def test_using_dataset_in_s3(self, platform="ten_x_v2"):
        # must NOT end with a slash
        file_prefix = "test"
        output_prefix = os.path.join(self.path_temp, file_prefix)

        params = [
            ("run", platform),
            ("--local",),
            ("--output-prefix", output_prefix),
            ("--index", dataset_s3.index),
            ("--barcode-files", dataset_s3.barcodes % platform),
            ("--barcode-fastq", dataset_s3.barcode_fastq % platform),
            ("--genomic-fastq", dataset_s3.genomic_fastq % platform),
            ("--star-args", "runRNGseed=0"),
        ]

        argv = [element for tupl in params for element in tupl]

        if platform != "drop_seq":
            argv += ["--barcode-files", dataset_s3.barcodes % platform]

        main.main(argv)

        # get output file list
        files = get_output_file_list(self.test_id, self.path_temp)

        # check if each expected file is found in the list of files generated
        for file in expected_output_files(file_prefix):
            self.assertIn(file, files)

    def test_using_local_dataset(self, platform="ten_x_v2"):
        # must NOT end with a slash
        file_prefix = "test"
        output_prefix = os.path.join(self.path_temp, file_prefix)

        params = [
            ("run", platform),
            ("--local",),
            ("--output-prefix", output_prefix),
            ("--index", dataset_local.index),
            ("--barcode-files", dataset_local.barcodes % platform),
            ("--barcode-fastq", dataset_local.barcode_fastq % platform),
            ("--genomic-fastq", dataset_local.genomic_fastq % platform),
            ("--star-args", "runRNGseed=0"),
        ]

        argv = [element for tupl in params for element in tupl]

        if platform != "drop_seq":
            argv += ["--barcode-files", dataset_local.barcodes % platform]

        main.main(argv)

        # get output file list
        files = get_output_file_list(self.test_id, self.path_temp)

        # check if each expected file is found in the list of files generated
        for file in expected_output_files(file_prefix):
            self.assertIn(file, files)
