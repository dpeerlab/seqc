from unittest import TestCase, mock
import nose2
import os
import numpy as np
from seqc.read_array import ReadArray
from seqc import rmt_correction


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
    @mock.patch(
        "seqc.rmt_correction._get_available_memory", return_value=50 * 1024 ** 3
    )
    def test_should_return_correct_max_workers(self, mock_mem):

        n_workers = rmt_correction._calc_max_workers(self.ra)

        self.assertEqual(n_workers, 7)

    # 1TB
    @mock.patch("seqc.rmt_correction._get_available_memory", return_value=1079354630144)
    def test_should_return_correct_max_workers2(self, mock_mem):

        n_workers = rmt_correction._calc_max_workers(self.ra)

        self.assertEqual(n_workers, 156)

    # having less memory than ra size
    @mock.patch("seqc.rmt_correction._get_available_memory")
    def test_should_return_one_if_ra_larger_than_mem(self, mock_mem):

        ra_size = self.ra.data.nbytes + self.ra.genes.nbytes + self.ra.positions.nbytes

        # assume the available memory is a half of ra
        mock_mem.return_value = int(ra_size) / 2

        n_workers = rmt_correction._calc_max_workers(self.ra)

        self.assertEqual(n_workers, 1)


class TestRmtCorrection2(TestCase):
    @classmethod
    def setUp(self):
        # pre-allocate arrays
        n_barcodes = 183416337
        data = np.recarray((n_barcodes,), ReadArray._dtype)
        genes = np.zeros(n_barcodes, dtype=np.int32)
        positions = np.zeros(n_barcodes, dtype=np.int32)
        self.ra = ReadArray(data, genes, positions)

        import pickle

        with open("pre-correction-ra.pickle", "wb") as fout:
            pickle.dump(self.ra, fout)

    @classmethod
    def tearDown(self):
        import os

        try:
            os.remove("pre-correction-ra.pickle")
        except:
            pass

    @mock.patch("seqc.rmt_correction._correct_errors_by_cell_group", return_value=0)
    def test_correct_errors_by_chunks(self, mock_correct):
        cell_group = [1, 2, 3]
        x = rmt_correction._correct_errors_by_cell_group_chunks(
            self.ra, cell_group, 0.02, 0.05
        )
        mock_correct.assert_called()
        self.assertEquals(len(cell_group), mock_correct.call_count)
        self.assertEquals([0, 0, 0], x)


if __name__ == "__main__":
    nose2.main()
