import nose2
import unittest

import seqc
from seqc.core import main


# class TestSEQC(unittest.TestCase):
#     def setUp(self):
#         pass

#     def tearDown(self):
#         pass

#     def test_args(self):

#         argv = ["start", "-k", "/Users/dchun/dpeerlab-chunj.pem", "-t", "t2.micro"]

#         self.assertRaises(ValueError, lambda: main.main(argv))

# class MyUnitTest(unittest.TestCase):
#     def setUp(self):
#         pass

#     def tearDown(self):
#         pass

#     def test_args(self):

#         # argv = [
#         #     "run", "ten_x_v2", "--local",
#         #     "--index", "s3://seqc-public/genomes/hg38_chr19/",
#         #     "--barcode-files", "s3://seqc-public/barcodes/ten_x_v2/flat/",
#         #     "--genomic-fastq", "./test-data/genomic/",
#         #     "--barcode-fastq", "./test-data/barcode/",
#         #     "--output-prefix", "./test-data/seqc-results/",
#         #     "--email", "jaeyoung.chun@gmail.com",
#         #     "--star-args", "\"runRNGseed=0\""
#         # ]

#         argv = [
#             "run"
#         ]        

#         try:
#             main.main(argv)
#             # self.assertRaises(BaseException, lambda: main.main(argv))
#         except:
#             pass
#         # self.assertRaises(ValueError, lambda: main.main(argv))


# class TestSEQC(unittest.TestCase):
#     def setUp(self):
#         pass

#     def tearDown(self):
#         pass

#     def test_args(self):

#         from seqc.sequence import gtf

#         # remove any invalid ids from the annotation file
#         gr = gtf.Reader("./test-data/homo_sapiens.gtf.gz")

#         for line_fields in gr:
#             record = gtf.Record(line_fields)
#             print(record)
#             biotype = record.attribute("gene_biotype")
#             print(biotype)

#         # self.assertRaises(ValueError, lambda: main.main(argv))


if __name__ == "__main__":

    unittest.main()
