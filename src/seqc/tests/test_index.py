import os
import shutil
import uuid
import unittest
import gzip
import pandas as pd
import numpy as np
import ftplib
import nose2
from nose2.tools import params
import seqc
from seqc.sequence import index, gtf
from seqc import io


def expected_output_files():

    files = set(
        [
            "Genome",
            "SA",
            "SAindex",
            "annotations.gtf",
            "chrLength.txt",
            "chrName.txt",
            "chrNameLength.txt",
            "chrStart.txt",
            "exonGeTrInfo.tab",
            "exonInfo.tab",
            "geneInfo.tab",
            "genomeParameters.txt",
            "sjdbInfo.txt",
            "sjdbList.fromGTF.out.tab",
            "sjdbList.out.tab",
            "transcriptInfo.tab",
        ]
    )

    return files


class TestIndexRemote(unittest.TestCase):

    s3_bucket = "dp-lab-cicd"

    @classmethod
    def setUp(cls):
        cls.test_id = str(uuid.uuid4())
        cls.outdir = os.path.join(os.environ["TMPDIR"], "seqc-test", cls.test_id)
        os.makedirs(cls.outdir, exist_ok=True)

    @classmethod
    def tearDown(self):
        if os.path.isdir(self.outdir):
            shutil.rmtree(self.outdir, ignore_errors=True)

    def test_upload_star_index_correctly_places_index_on_s3(self):
        organism = "ciona_intestinalis"
        # must end with a slash
        test_folder = f"seqc/index-{organism}-{self.test_id}/"

        idx = index.Index(
            organism, ["external_gene_name"], index_folder_name=self.outdir
        )
        index_directory = os.path.join(self.outdir, organism) + "/"
        idx._download_ensembl_files(ensemble_release=None)
        idx._subset_genes()
        idx._create_star_index()
        idx._upload_index(index_directory, f"s3://{self.s3_bucket}/{test_folder}")

        # check files generated in S3
        files = io.S3.listdir(self.s3_bucket, test_folder)

        # extract only filenames (i.e. remove directory hierarchy)
        # convert to a set for easy comparison
        files = set(map(lambda filename: filename.replace(test_folder, ""), files))

        # check for the exact same filenames
        self.assertSetEqual(files, expected_output_files())

    def test_create_index_produces_and_uploads_an_index(self):
        organism = "ciona_intestinalis"
        # must end with a slash
        test_folder = f"seqc/index-{organism}-{self.test_id}/"

        idx = index.Index(
            organism, ["external_gene_name"], index_folder_name=self.outdir
        )
        idx.create_index(
            s3_location=f"s3://{self.s3_bucket}/{test_folder}",
            ensemble_release=None,
            read_length=101,
        )

        # check files generated in S3
        files = io.S3.listdir(self.s3_bucket, test_folder)

        # extract only filenames (i.e. remove directory hierarchy)
        # convert to a set for easy comparison
        files = set(map(lambda filename: filename.replace(test_folder, ""), files))

        # check for the exact same filenames
        self.assertSetEqual(files, expected_output_files())


class MyUnitTest(unittest.TestCase):

    s3_bucket = "dp-lab-cicd"

    @classmethod
    def setUp(cls):
        cls.test_id = str(uuid.uuid4())
        cls.outdir = os.path.join(os.environ["TMPDIR"], "seqc-test", cls.test_id)
        os.makedirs(cls.outdir, exist_ok=True)

    @classmethod
    def tearDown(self):
        if os.path.isdir(self.outdir):
            shutil.rmtree(self.outdir, ignore_errors=True)

    def test_Index_raises_ValueError_when_organism_is_not_provided(self):
        self.assertRaises(ValueError, index.Index, organism="", additional_id_types=[])

    def test_Index_raises_ValueError_when_organism_isnt_lower_case(self):
        self.assertRaises(
            ValueError, index.Index, organism="Homo_sapiens", additional_id_types=[]
        )
        self.assertRaises(
            ValueError, index.Index, organism="Homo_Sapiens", additional_id_types=[]
        )
        self.assertRaises(
            ValueError, index.Index, organism="hoMO_Sapiens", additional_id_types=[]
        )

    def test_Index_raises_ValueError_when_organism_has_no_underscore(self):
        self.assertRaises(
            ValueError, index.Index, organism="homosapiens", additional_id_types=[]
        )

    def test_Index_raises_TypeError_when_additional_id_fields_is_not_correct_type(self):
        self.assertRaises(
            TypeError,
            index.Index,
            organism="homo_sapiens",
            additional_id_types="not_an_array_tuple_or_list",
        )
        self.assertRaises(
            TypeError, index.Index, organism="homo_sapiens", additional_id_types=""
        )

    def test_False_evaluating_additional_id_fields_are_accepted_but_set_empty_list(
        self,
    ):
        idx = index.Index("homo_sapiens", [])
        self.assertEqual(idx.additional_id_types, [])
        idx = index.Index("homo_sapiens", tuple())
        self.assertEqual(idx.additional_id_types, [])
        idx = index.Index("homo_sapiens", np.array([]))
        self.assertEqual(idx.additional_id_types, [])

    def test_converter_xml_contains_one_attribute_line_per_gene_list(self):
        idx = index.Index("homo_sapiens", ["hgnc_symbol", "mgi_symbol"])
        self.assertEqual(idx._converter_xml.count("Attribute name"), 3)
        idx = index.Index("homo_sapiens", [])
        self.assertEqual(idx._converter_xml.count("Attribute name"), 1)

    def test_converter_xml_formats_genome_as_first_initial_plus_species(self):
        idx = index.Index("homo_sapiens", ["hgnc_symbol", "mgi_symbol"])
        self.assertIn("hsapiens", idx._converter_xml)
        idx = index.Index("mus_musculus")
        self.assertIn("mmusculus", idx._converter_xml)

    def test_identify_gtf_file_should_return_correct_file(self):

        files = [
            "CHECKSUMS",
            "Homo_sapiens.GRCh38.86.abinitio.gtf.gz",
            "Homo_sapiens.GRCh38.86.chr.gtf.gz",
            "Homo_sapiens.GRCh38.85.chr.gtf.gz",
            "Homo_sapiens.GRCh38.86.chr_patch_hapl_scaff.gtf.gz",
            "Homo_sapiens.GRCh38.86.gtf.gz",
            "README",
        ]
        release_num = 86

        filename = index.Index._identify_gtf_file(files, release_num)

        self.assertEqual(filename, "Homo_sapiens.GRCh38.86.chr.gtf.gz")

    def test_identify_gtf_file_should_throw_exception(self):

        files = [
            "CHECKSUMS",
            "Homo_sapiens.GRCh38.86.abinitio.gtf.gz",
            "Homo_sapiens.GRCh38.86.chr_patch_hapl_scaff.gtf.gz",
            "Homo_sapiens.GRCh38.86.gtf.gz",
            "README",
        ]
        release_num = 86

        self.assertRaises(
            FileNotFoundError,
            index.Index._identify_gtf_file,
            files=files,
            release_num=release_num,
        )

    def test_can_login_to_ftp_ensembl(self):
        with ftplib.FTP(host="ftp.ensembl.org") as ftp:
            ftp.login()

    def test_download_converter_gets_output_and_is_pandas_loadable(self):
        idx = index.Index("ciona_intestinalis", ["external_gene_name"])
        filename = os.path.join(self.outdir, "ci.csv")
        idx._download_conversion_file(filename)
        converter = pd.read_csv(filename, index_col=0)
        self.assertGreaterEqual(len(converter), 10)
        self.assertEqual(converter.shape[1], 1)

    def test_identify_newest_release_finds_a_release_which_is_gt_eq_85(self):
        idx = index.Index("ciona_intestinalis", ["external_gene_name"])
        with ftplib.FTP(host="ftp.ensembl.org") as ftp:
            ftp.login()
            newest = idx._identify_newest_release(ftp)
        self.assertGreaterEqual(int(newest), 85)  # current=85, they only get bigger

    def test_identify_genome_file_finds_primary_assembly_when_present(self):
        idx = index.Index("homo_sapiens", ["entrezgene"])
        with ftplib.FTP(host="ftp.ensembl.org") as ftp:
            ftp.login()
            newest = idx._identify_newest_release(ftp)
            ftp.cwd("/pub/release-%d/fasta/%s/dna" % (newest, idx.organism))
            filename = idx._identify_genome_file(ftp.nlst())
        self.assertIn("primary_assembly", filename)

    def test_identify_genome_file_defaults_to_toplevel_when_no_primary_assembly(self):
        idx = index.Index("ciona_intestinalis", ["entrezgene"])
        with ftplib.FTP(host="ftp.ensembl.org") as ftp:
            ftp.login()
            newest = idx._identify_newest_release(ftp)
            ftp.cwd("/pub/release-%d/fasta/%s/dna" % (newest, idx.organism))
            filename = idx._identify_genome_file(ftp.nlst())
        self.assertIn("toplevel", filename)

    def test_download_fasta_file_gets_a_properly_formatted_file(self):
        idx = index.Index(
            "ciona_intestinalis", ["external_gene_name"], index_folder_name=self.outdir
        )
        with ftplib.FTP(host="ftp.ensembl.org") as ftp:
            ftp.login()
            filename = os.path.join(self.outdir, "ci.fa.gz")
            idx._download_fasta_file(ftp, filename, ensemble_release=None)
        with gzip.open(filename, "rt") as f:
            self.assertIs(
                f.readline()[0], ">"
            )  # starting character for genome fa record

    def test_identify_annotation_file_finds_a_gtf_file(self):
        idx = index.Index("ciona_intestinalis", ["external_gene_name"])
        with ftplib.FTP(host="ftp.ensembl.org") as ftp:
            ftp.login()
            newest = idx._identify_newest_release(ftp)
            ftp.cwd("/pub/release-%d/gtf/%s/" % (newest, idx.organism))
            filename = idx._identify_gtf_file(ftp.nlst(), newest)
        self.assertIsNotNone(filename)

    def test_download_gtf_file_gets_a_file_readable_by_seqc_gtf_reader(self):

        idx = index.Index("ciona_intestinalis", ["entrezgene"])

        with ftplib.FTP(host="ftp.ensembl.org") as ftp:
            ftp.login()
            filename = self.outdir + "ci.gtf.gz"
            idx._download_gtf_file(ftp, filename, ensemble_release=99)

        rd = gtf.Reader(filename)
        (transcript_chromosome, transcript_strand, transcript_gene_id), exons = next(
            rd.iter_transcripts()
        )

        # (('1', '+', 17842), [['1', 'ensembl', 'exon', '1636', '1902', '.', '+', '.', 'gene_id "ENSCING00000017842"; gene_version "1"; transcript_id "ENSCINT00000030147"; transcript_version "1"; exon_number "1"; gene_name "RNaseP_nuc"; gene_source "ensembl"; gene_biotype "misc_RNA"; transcript_name "RNaseP_nuc-201"; transcript_source "ensembl"; transcript_biotype "misc_RNA"; exon_id "ENSCINE00000207263"; exon_version "1";\n']])
        self.assertEqual(transcript_chromosome, "1")
        self.assertEqual(transcript_strand, "+")
        self.assertEqual(transcript_gene_id, 17842)
        self.assertEqual(len(exons), 1)

    def test_subset_genes_should_returns_original_if_no_additional_fields_or_valid_biotypes(
        self,
    ):

        fasta_name = os.path.join(self.outdir, "ci.fa.gz")
        gtf_name = os.path.join(self.outdir, "ci.gtf.gz")
        conversion_name = os.path.join(self.outdir, "ci_ids.csv")

        idx = index.Index("ciona_intestinalis", index_folder_name=self.outdir)

        idx._download_ensembl_files(
            ensemble_release=None,
            fasta_name=fasta_name,
            gtf_name=gtf_name,
            conversion_name=conversion_name,
        )
        truncated_gtf = os.path.join(self.outdir, "test.gtf")
        idx._subset_genes(conversion_name, gtf_name, truncated_gtf, valid_biotypes=None)

        # expect the same file as the original file
        self.assertTrue(os.path.isfile(truncated_gtf))

        # the current implementation of GTF Reader doesn't allow this:
        # for gr1, gr2 in zip(gtf.Reader(gtf_name), gtf.Reader(truncated_gtf)):

        records = []
        for gr in gtf.Reader(gtf_name):
            records.append(gtf.Record(gr))

        for i, gr in enumerate(gtf.Reader(truncated_gtf)):
            rec1 = records[i]
            rec2 = gtf.Record(gr)
            self.assertEqual(rec1, rec2)

    def test_subset_genes_produces_a_reduced_annotation_file_when_passed_fields(self):
        organism = "ciona_intestinalis"
        idx = index.Index(
            organism, ["external_gene_name"], index_folder_name=self.outdir
        )
        idx._download_ensembl_files(ensemble_release=None)
        self.assertTrue(
            os.path.isfile(os.path.join(self.outdir, "%s.fa.gz" % organism)),
            "fasta file not found",
        )
        self.assertTrue(
            os.path.isfile(os.path.join(self.outdir, "%s.gtf.gz" % organism)),
            "gtf file not found",
        )
        self.assertTrue(
            os.path.isfile(os.path.join(self.outdir, "%s_ids.csv" % organism)),
            "id file not found",
        )

        valid_biotypes = {"protein_coding", "lincRNA"}
        idx._subset_genes(valid_biotypes=valid_biotypes)

        self.assertTrue(
            os.path.isfile(os.path.join(self.outdir, organism, "annotations.gtf"))
        )
        gr_subset = gtf.Reader(os.path.join(self.outdir, organism, "annotations.gtf"))
        gr_complete = gtf.Reader(os.path.join(self.outdir, "%s.gtf.gz" % organism))
        self.assertLess(
            len(gr_subset),
            len(gr_complete),
            "Subset annotation was not smaller than the complete annotation",
        )

        # make sure only valid biotypes are returned
        complete_invalid = False

        for r in gr_complete:
            record = gtf.Record(r)
            if record.attribute("gene_biotype") not in valid_biotypes:
                complete_invalid = True
                break
        self.assertTrue(complete_invalid)

        subset_invalid = False
        for r in gr_subset:
            record = gtf.Record(r)
            if record.attribute("gene_biotype") not in valid_biotypes:
                subset_invalid = True
                break
        self.assertFalse(subset_invalid)
        self.assertGreater(len(gr_subset), 0)

    def test_create_star_index_produces_an_index(self):
        organism = "ciona_intestinalis"
        idx = index.Index(
            organism, ["external_gene_name"], index_folder_name=self.outdir
        )
        idx._download_ensembl_files(ensemble_release=None)
        idx._subset_genes()
        idx._create_star_index()
        expected_file = os.path.join(
            self.outdir,
            "{outdir}/{organism}/SAindex".format(outdir=self.outdir, organism=organism),
        )
        self.assertTrue(os.path.isfile(expected_file))


#########################################################################################

if __name__ == "__main__":
    nose2.main()

#########################################################################################
