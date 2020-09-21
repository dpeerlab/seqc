import unittest
import os
import uuid
import shutil
import re
from seqc.core import main
from seqc import io
import boto3
from nose2.tools import params
from test_dataset import dataset_s3


def get_instance_by_test_id(test_id):

    ec2 = boto3.resource("ec2")
    instances = ec2.instances.filter(
        Filters=[{"Name": "tag:TestID", "Values": [test_id]}]
    )
    instances = list(instances)

    if len(instances) != 1:
        raise Exception("Test ID is not found or not unique!")

    return instances[0]


def expected_output_files(output_prefix):

    files = set(
        [
            f"{output_prefix}.h5",
            f"{output_prefix}_Aligned.out.bam",
            f"{output_prefix}_alignment_summary.txt",
            f"{output_prefix}_cell_filters.png",
            f"{output_prefix}_de_gene_list.txt",
            f"{output_prefix}_dense.csv",
            f"{output_prefix}_merged.fastq.gz",
            f"{output_prefix}_mini_summary.json",
            f"{output_prefix}_mini_summary.pdf",
            f"{output_prefix}_seqc_log.txt",
            f"{output_prefix}_sparse_counts_barcodes.csv",
            f"{output_prefix}_sparse_counts_genes.csv",
            f"{output_prefix}_sparse_molecule_counts.mtx",
            f"{output_prefix}_sparse_read_counts.mtx",
            f"{output_prefix}_summary.tar.gz",
            f"seqc_log.txt",
        ]
    )

    return files


def expected_output_files_run_from_merged(output_prefix):

    files = expected_output_files(output_prefix)

    excludes = set([f"{output_prefix}_merged.fastq.gz"])

    return files - excludes


def expected_output_files_run_from_bam(output_prefix):

    files = expected_output_files(output_prefix)

    excludes = set(
        [
            f"{output_prefix}_Aligned.out.bam",
            f"{output_prefix}_alignment_summary.txt",
            f"{output_prefix}_merged.fastq.gz",
        ]
    )

    return files - excludes


def get_output_file_list(test_id, s3_bucket, test_folder):

    # get instance and wait until terminated
    instance = get_instance_by_test_id(test_id)
    instance.wait_until_terminated()

    # check files generated in S3
    files = io.S3.listdir(s3_bucket, test_folder)

    # extract only filenames (i.e. remove directory hierarchy)
    # convert to a set for easy comparison
    files = set(map(lambda filename: filename.replace(test_folder, ""), files))

    return files


def check_for_success_msg(s3_seqc_log_uri, path_temp):

    # download seqc_log.txt
    io.S3.download(
        link=s3_seqc_log_uri, prefix=path_temp, overwrite=True, recursive=False
    )

    # check if seqc_log.txt has a successful message
    with open(os.path.join(path_temp, "seqc_log.txt"), "rt") as fin:
        logs = fin.read()
        match = re.search(r"Execution completed successfully", logs, re.MULTILINE)

        return True if match else False


class TestRunRemote(unittest.TestCase):

    email = os.environ["SEQC_TEST_EMAIL"]
    rsa_key = os.environ["SEQC_TEST_RSA_KEY"]
    ami_id = os.environ["SEQC_TEST_AMI_ID"]

    s3_bucket = "dp-lab-cicd"

    @classmethod
    def setUp(cls):
        cls.test_id = str(uuid.uuid4())
        cls.path_temp = os.path.join(
            os.environ["TMPDIR"], "seqc-test", str(uuid.uuid4())
        )
        os.makedirs(cls.path_temp, exist_ok=True)

    @classmethod
    def tearDown(self):
        if os.path.isdir(self.path_temp):
            shutil.rmtree(self.path_temp, ignore_errors=True)

    @params("in_drop_v2", "ten_x_v2")
    def test_remote_from_raw_fastq(self, platform="ten_x_v2"):
        output_prefix = "from-raw-fastq"
        # must end with a slash
        test_folder = f"seqc/run-{platform}-{self.test_id}/"

        params = [
            ("run", platform),
            ("--output-prefix", "from-raw-fastq"),
            ("--upload-prefix", f"s3://{self.s3_bucket}/{test_folder}"),
            ("--index", dataset_s3.index),
            ("--email", self.email),
            ("--barcode-fastq", dataset_s3.barcode_fastq % platform),
            ("--genomic-fastq", dataset_s3.genomic_fastq % platform),
            ("--instance-type", "r5.2xlarge"),
            ("--spot-bid", "1.0"),
            ("--rsa-key", self.rsa_key),
            ("--debug",),
            ("--remote-update",),
            ("--ami-id", self.ami_id),
            ("--user-tags", f"TestID:{self.test_id}"),
        ]

        argv = [element for tupl in params for element in tupl]

        if platform != "drop_seq":
            argv += ["--barcode-files", dataset_s3.barcodes % platform]

        main.main(argv)

        # wait until terminated
        # get output file list
        files = get_output_file_list(self.test_id, self.s3_bucket, test_folder)

        # check for the exact same filenames
        self.assertSetEqual(files, expected_output_files(output_prefix))

        # check for success message in seqc_log.txt
        has_success_msg = check_for_success_msg(
            s3_seqc_log_uri="s3://{}/{}".format(
                self.s3_bucket, os.path.join(test_folder, "seqc_log.txt")
            ),
            path_temp=self.path_temp,
        )

        self.assertTrue(
            has_success_msg, msg="Unable to find the success message in the log"
        )

    def test_remote_from_merged(self, platform="in_drop_v2"):
        output_prefix = "from-merged"
        # must end with a slash
        test_folder = f"seqc/run-{platform}-{self.test_id}/"

        params = [
            ("run", platform),
            ("--output-prefix", output_prefix),
            ("--upload-prefix", f"s3://{self.s3_bucket}/{test_folder}"),
            ("--index", dataset_s3.index),
            ("--email", self.email),
            ("--merged-fastq", dataset_s3.merged_fastq % (platform, platform)),
            ("--rsa-key", self.rsa_key),
            ("--instance-type", "r5.2xlarge"),
            ("--ami-id", self.ami_id),
            ("--remote-update",),
            ("--user-tags", f"TestID:{self.test_id}")
            # ('--spot-bid', '1.0')
        ]

        argv = [element for tupl in params for element in tupl]

        if platform != "drop_seq":
            argv += ["--barcode-files", dataset_s3.barcodes % platform]

        main.main(argv)

        # wait until terminated
        # get output file list
        files = get_output_file_list(self.test_id, self.s3_bucket, test_folder)

        # check for the exact same filenames
        self.assertSetEqual(files, expected_output_files_run_from_merged(output_prefix))

        # check for success message in seqc_log.txt
        has_success_msg = check_for_success_msg(
            s3_seqc_log_uri="s3://{}/{}".format(
                self.s3_bucket, os.path.join(test_folder, "seqc_log.txt")
            ),
            path_temp=self.path_temp,
        )

        self.assertTrue(
            has_success_msg, msg="Unable to find the success message in the log"
        )

    def test_remote_from_bamfile(self, platform="in_drop_v2"):
        output_prefix = "from-bamfile"
        # must end with a slash
        test_folder = f"seqc/run-{platform}-{self.test_id}/"

        params = [
            ("run", platform),
            ("--output-prefix", output_prefix),
            ("--upload-prefix", f"s3://{self.s3_bucket}/{test_folder}"),
            ("--index", dataset_s3.index),
            ("--email", self.email),
            ("--alignment-file", dataset_s3.bam % platform),
            ("--rsa-key", self.rsa_key),
            ("--instance-type", "r5.2xlarge"),
            ("--debug",),
            ("--ami-id", self.ami_id),
            ("--remote-update",),
            ("--user-tags", f"TestID:{self.test_id}")
            # ('--spot-bid', '1.0')
        ]

        argv = [element for tupl in params for element in tupl]

        if platform != "drop_seq":
            argv += ["--barcode-files", dataset_s3.barcodes % platform]

        main.main(argv)

        # wait until terminated
        # get output file list
        files = get_output_file_list(self.test_id, self.s3_bucket, test_folder)

        # check for the exact same filenames
        self.assertSetEqual(files, expected_output_files_run_from_bam(output_prefix))

        # check for success message in seqc_log.txt
        has_success_msg = check_for_success_msg(
            s3_seqc_log_uri="s3://{}/{}".format(
                self.s3_bucket, os.path.join(test_folder, "seqc_log.txt")
            ),
            path_temp=self.path_temp,
        )

        self.assertTrue(
            has_success_msg, msg="Unable to find the success message in the log"
        )
