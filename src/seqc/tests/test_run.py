from unittest import TestCase, mock
import os
import uuid
import shutil
import nose2
import numpy as np
from nose2.tools import params
from seqc.sequence import index, gtf
from seqc.sequence.encodings import DNA3Bit
from seqc.read_array import ReadArray
from seqc import rmt_correction
from test_dataset import dataset_local


class TestReadArray(TestCase):
    @classmethod
    def setUp(cls):
        cls.test_id = str(uuid.uuid4())
        cls.path_temp = os.path.join(
            os.environ["TMPDIR"], "seqc-test", str(uuid.uuid4())
        )
        cls.annotation = os.path.join(dataset_local.index, "annotations.gtf")
        cls.translator = gtf.GeneIntervals(cls.annotation, 10000)

    @classmethod
    def tearDown(self):
        if os.path.isdir(self.path_temp):
            shutil.rmtree(self.path_temp, ignore_errors=True)

    def test_read_array_creation(self, platform="ten_x_v2"):
        ra = ReadArray.from_alignment_file(
            dataset_local.bam % platform, self.translator, required_poly_t=0
        )
        self.assertIsNotNone(ra)

    def test_read_array_rmt_decode_10x_v2(self):
        platform = "ten_x_v2"

        # create a readarray
        ra = ReadArray.from_alignment_file(
            dataset_local.bam % platform, self.translator, required_poly_t=0
        )

        # see if we can decode numeric UMI back to nucleotide sequence
        dna3bit = DNA3Bit()
        for rmt in ra.data["rmt"]:
            decoded = dna3bit.decode(rmt).decode()
            # ten_x_v2 UMI length = 10 nt
            self.assertEqual(len(decoded), 10)

    def test_read_array_rmt_decode_10x_v3(self):
        platform = "ten_x_v3"

        # create a readarray
        ra = ReadArray.from_alignment_file(
            dataset_local.bam % platform, self.translator, required_poly_t=0
        )

        # see if we can decode numeric UMI back to nucleotide sequence
        dna3bit = DNA3Bit()
        for rmt in ra.data["rmt"]:
            decoded = dna3bit.decode(rmt).decode()
            # ten_x_v3 UMI length = 12 nt
            self.assertEqual(len(decoded), 12)


class TestTranslator(TestCase):
    @classmethod
    def setUp(cls):
        cls.test_id = str(uuid.uuid4())
        cls.path_temp = os.path.join(
            os.environ["TMPDIR"], "seqc-test", str(uuid.uuid4())
        )
        cls.annotation = os.path.join(dataset_local.index, "annotations.gtf")

    @classmethod
    def tearDown(self):
        if os.path.isdir(self.path_temp):
            shutil.rmtree(self.path_temp, ignore_errors=True)

    def test_construct_translator(self):
        translator = gtf.GeneIntervals(self.annotation)
        self.assertIsNotNone(translator)

    def test_num_of_transcripts(self):
        rd = gtf.Reader(self.annotation)
        num_transcripts = sum(1 for _ in rd.iter_transcripts())
        # awk -F'\t' '$3=="transcript" { print $0 }' annotations.gtf | wc -l
        self.assertEqual(num_transcripts, 12747)

    def test_iter_transcripts(self):
        rd = gtf.Reader(self.annotation)
        (transcript_chromosome, transcript_strand, transcript_gene_id), exons = next(
            rd.iter_transcripts()
        )

        # this should give us 3 exons of the first transcript of the first gene found in inverse order:
        #
        # chr19  HAVANA  gene        60951  71626  .  -  .  gene_id  "ENSG00000282458.1";  gene_type      "transcribed_processed_pseudogene";  gene_status  "KNOWN";                             gene_name    "WASH5P";  level      2;         havana_gene      "OTTHUMG00000180466.8";
        # chr19  HAVANA  transcript  60951  70976  .  -  .  gene_id  "ENSG00000282458.1";  transcript_id  "ENST00000632506.1";                 gene_type    "transcribed_processed_pseudogene";  gene_status  "KNOWN";   gene_name  "WASH5P";  transcript_type  "processed_transcript";              transcript_status  "KNOWN";  transcript_name  "WASH5P-008";  level        2;  tag                       "basic";              transcript_support_level  "1";                     havana_gene               "OTTHUMG00000180466.8";  havana_transcript         "OTTHUMT00000471217.2";
        # chr19  HAVANA  exon        70928  70976  .  -  .  gene_id  "ENSG00000282458.1";  transcript_id  "ENST00000632506.1";                 gene_type    "transcribed_processed_pseudogene";  gene_status  "KNOWN";   gene_name  "WASH5P";  transcript_type  "processed_transcript";              transcript_status  "KNOWN";  transcript_name  "WASH5P-008";  exon_number  1;  exon_id                   "ENSE00003781173.1";  level                     2;                       tag                       "basic";                 transcript_support_level  "1";                     havana_gene        "OTTHUMG00000180466.8";  havana_transcript         "OTTHUMT00000471217.2";
        # chr19  HAVANA  exon        66346  66499  .  -  .  gene_id  "ENSG00000282458.1";  transcript_id  "ENST00000632506.1";                 gene_type    "transcribed_processed_pseudogene";  gene_status  "KNOWN";   gene_name  "WASH5P";  transcript_type  "processed_transcript";              transcript_status  "KNOWN";  transcript_name  "WASH5P-008";  exon_number  2;  exon_id                   "ENSE00003783498.1";  level                     2;                       tag                       "basic";                 transcript_support_level  "1";                     havana_gene        "OTTHUMG00000180466.8";  havana_transcript         "OTTHUMT00000471217.2";
        # chr19  HAVANA  exon        60951  61894  .  -  .  gene_id  "ENSG00000282458.1";  transcript_id  "ENST00000632506.1";                 gene_type    "transcribed_processed_pseudogene";  gene_status  "KNOWN";   gene_name  "WASH5P";  transcript_type  "processed_transcript";              transcript_status  "KNOWN";  transcript_name  "WASH5P-008";  exon_number  3;  exon_id                   "ENSE00003783010.1";  level                     2;                       tag                       "basic";                 transcript_support_level  "1";                     havana_gene        "OTTHUMG00000180466.8";  havana_transcript         "OTTHUMT00000471217.2";

        self.assertEqual(transcript_chromosome, "chr19")
        self.assertEqual(transcript_strand, "-")
        self.assertEqual(transcript_gene_id, 282458)
        self.assertEqual(len(exons), 3)

        # 8th column has exon ID
        self.assertIn("ENSE00003783010.1", exons[0][8])  # exon number 3
        self.assertIn("ENSE00003783498.1", exons[1][8])  # exon number 2
        self.assertIn("ENSE00003781173.1", exons[2][8])  # exon number 1

    def test_translate(self):
        translator = gtf.GeneIntervals(self.annotation)
        # chr19	HAVANA	gene	60951	71626	.	-	.	gene_id "ENSG00000282458.1"; gene_type "transcribed_processed_pseudogene"; gene_status "KNOWN"; gene_name "WASH5P"; level 2; havana_gene "OTTHUMG00000180466.8";
        gene_id = translator.translate("chr19", "-", 60951)
        self.assertEqual(gene_id, 282458)


class TestRmtCorrection(TestCase):
    @classmethod
    def setUp(self):
        # pre-allocate arrays
        n_barcodes = 183416337
        data = np.recarray((n_barcodes,), ReadArray._dtype)
        genes = np.zeros(n_barcodes, dtype=np.int32)
        positions = np.zeros(n_barcodes, dtype=np.int32)
        self.ra = ReadArray(data, genes, positions)

    @classmethod
    def tearDown(self):
        pass

    def test_should_return_correct_ra_size(self):

        ra_size = self.ra.data.nbytes + self.ra.genes.nbytes + self.ra.positions.nbytes

        self.assertEqual(4768824762, ra_size)

    # 64GB
    @mock.patch("seqc.rmt_correction._get_total_memory", return_value=64 * 1024 ** 3)
    def test_should_return_correct_max_workers(self, mock_mem):

        n_workers = rmt_correction._calc_max_workers(self.ra)

        self.assertEqual(n_workers, 9)

    # 1TB
    @mock.patch("seqc.rmt_correction._get_total_memory", return_value=1079354630144)
    def test_should_return_correct_max_workers2(self, mock_mem):

        n_workers = rmt_correction._calc_max_workers(self.ra)

        self.assertEqual(n_workers, 156)

    # having less memory than ra size
    @mock.patch("seqc.rmt_correction._get_total_memory")
    def test_should_return_one_if_ra_larger_than_mem(self, mock_mem):

        ra_size = self.ra.data.nbytes + self.ra.genes.nbytes + self.ra.positions.nbytes

        # assume the memory size is a half of ra
        mock_mem.return_value = int(ra_size) / 2

        n_workers = rmt_correction._calc_max_workers(self.ra)

        self.assertEqual(n_workers, 1)


if __name__ == "__main__":
    nose2.main()
