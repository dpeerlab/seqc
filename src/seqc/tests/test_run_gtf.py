from unittest import TestCase, mock
import os
import uuid
import shutil
import nose2
from seqc.sequence import gtf
from test_dataset import dataset_local


class TestGtf(TestCase):
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


if __name__ == "__main__":
    nose2.main()
