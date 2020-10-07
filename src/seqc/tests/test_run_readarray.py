from unittest import TestCase, mock
import os
import uuid
import shutil
import nose2
from test_dataset import dataset_local
from seqc.sequence.encodings import DNA3Bit
from seqc.read_array import ReadArray
from seqc.sequence import gtf


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
        ra, _ = ReadArray.from_alignment_file(
            dataset_local.bam % platform, self.translator, required_poly_t=0
        )
        self.assertIsNotNone(ra)

    def test_read_array_rmt_decode_10x_v2(self):
        platform = "ten_x_v2"

        # create a readarray
        ra, _ = ReadArray.from_alignment_file(
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
        ra, _ = ReadArray.from_alignment_file(
            dataset_local.bam % platform, self.translator, required_poly_t=0
        )

        # see if we can decode numeric UMI back to nucleotide sequence
        dna3bit = DNA3Bit()
        for rmt in ra.data["rmt"]:
            decoded = dna3bit.decode(rmt).decode()
            # ten_x_v3 UMI length = 12 nt
            self.assertEqual(len(decoded), 12)


if __name__ == "__main__":
    nose2.main()
