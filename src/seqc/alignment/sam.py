from collections import namedtuple
from subprocess import Popen, PIPE
import shutil
import gzip


def get_version():

    proc = Popen(["samtools", "--version"], stderr=PIPE, stdout=PIPE)
    out, err = proc.communicate()
    if err:
        raise ChildProcessError(err)

    # e.g.
    # samtools 1.9
    # Using htslib 1.9
    # Copyright (C) 2018 Genome Research Ltd.
    # --> 'samtools 1.9'
    version = out.decode().strip().split("\n")[0]

    # --> '1.9'
    version = version.split(" ")[1]

    return version


class SamRecord:
    """Simple record object allowing access to Sam record properties"""

    __slots__ = ["_record", "_parsed_name_field"]

    NameField = namedtuple("NameField", ["pool", "cell", "rmt", "poly_t", "name"])

    def __init__(self, record):
        self._record = record
        self._parsed_name_field = None

    def __repr__(self):
        return "<SamRecord{!s}>".format("\t".join(self._record))

    def __bytes__(self):
        return "\t".join(self._record) + "\n"

    @property
    def qname(self) -> str:
        return self._record[0]

    @property
    def flag(self) -> int:
        return int(self._record[1])

    @property
    def rname(self) -> str:
        return self._record[2]

    @property
    def pos(self) -> int:
        return int(self._record[3])

    @property
    def mapq(self) -> int:
        return int(self._record[4])

    @property
    def cigar(self) -> str:
        return self._record[5]

    @property
    def rnext(self) -> str:
        return self._record[6]

    @property
    def pnext(self) -> int:
        return int(self._record[7])

    @property
    def tlen(self) -> int:
        return int(self._record[8])

    @property
    def seq(self) -> str:
        return self._record[9]

    @property
    def qual(self) -> str:
        return self._record[10]

    @property
    def optional_fields(self):
        flags_ = {}
        for f in self._record[11:]:
            k, _, v = f.split(":")
            flags_[k] = int(v)
        return flags_

    def _parse_name_field(self):
        fields, name = self.qname.split(";")
        processed_fields = fields.split(":")
        processed_fields.append(name)
        self._parsed_name_field = self.NameField(*processed_fields)

    @property
    def pool(self) -> str:
        try:
            return self._parsed_name_field.pool
        except AttributeError:
            self._parse_name_field()
            return self._parsed_name_field.pool

    @property
    def rmt(self) -> str:
        try:
            return self._parsed_name_field.rmt
        except AttributeError:
            self._parse_name_field()
            return self._parsed_name_field.rmt

    @property
    def cell(self) -> str:
        try:
            return self._parsed_name_field.cell
        except AttributeError:
            self._parse_name_field()
            return self._parsed_name_field.cell

    @property
    def poly_t(self) -> str:
        try:
            return self._parsed_name_field.poly_t
        except AttributeError:
            self._parse_name_field()
            return self._parsed_name_field.poly_t

    @property
    def name(self):
        try:
            return self._parsed_name_field.name
        except AttributeError:
            self._parse_name_field()
            return self._parsed_name_field.name

    @property
    def is_mapped(self):
        return False if (int(self.flag) & 4) else True

    @property
    def is_unmapped(self):
        return not self.is_mapped

    @property
    def is_multimapped(self):
        return True if self.optional_fields["NH"] > 1 else False

    @property
    def is_uniquely_mapped(self):
        return True if self.optional_fields["NH"] == 1 else False

    @property
    def strand(self):
        minus_strand = int(self.flag) & 16
        return "-" if minus_strand else "+"

    # # todo this takes up 66% of the processing time for parsing the sam record
    # @property
    # def dust_low_complexity_score(self) -> int:
    #
    #     # Counts of 3-mers in the sequence
    #     counts = {}
    #     for i in range(len(self.seq) - 2):
    #         kmer = self.seq[i:i + 3]
    #         counts[kmer] = counts.get(kmer, 0) + 1
    #
    #     # Calculate dust score  # todo this is 30% faster when vectorized
    #     score = sum([i * (i - 1) / 2 for i in counts.values()]) / (len(self.seq) - 3)
    #
    #     # Scale score (Max score possible is no. of 3mers/2)
    #     score = int(score / ((len(self.seq) - 2) / 2) * 100)
    #
    #     return score


class Reader:
    """Simple sam reader, optimized for utility rather than speed"""

    def __init__(self, samfile: str):
        """
        :param samfile: str, location of a .sam file

        usage:
        if rd = Reader(samfile)
        :method __iter__: iterate over the .sam file's records (also usable in for loop)
        :method __len__: return the number of alignments in the file
        :method itermultialignments: return tuples of multiple alignments, all from the
           same fastq record
        """

        self._samfile = samfile
        try:
            samfile_iterator = iter(self)
            next(samfile_iterator)
        except RuntimeError as ex:
            raise ex
        except:
            raise ValueError(
                "%s is an invalid samfile. Please check file formatting." % samfile
            )

    @property
    def samfile(self):
        return self._samfile

    def _open(self):
        """
        seamlessly open self._samfile, whether gzipped or uncompressed
        :returns: open file object
        """
        if self.samfile.endswith(".gz"):
            fobj = gzip.open(self.samfile, "rb")
        elif self.samfile.endswith(".bam"):
            if not shutil.which("samtools"):
                raise RuntimeError("samtools utility must be installed to run bamfiles")
            p = Popen(["samtools", "view", self.samfile], stdout=PIPE)
            fobj = p.stdout
        else:
            fobj = open(self.samfile, "rb")
        return fobj

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """return an iterator over all non-header records in samfile"""
        fobj = self._open()
        try:
            for line in fobj:
                line = line.decode()
                # todo move this if statement to execute only until header is exhausted
                if line.startswith("@"):
                    continue
                yield SamRecord(line.strip().split("\t"))
        finally:
            fobj.close()

    def iter_multialignments(self):
        """yields tuples of all alignments for each fastq record"""
        sam_iter = iter(self)
        fq = [next(sam_iter)]
        for record in sam_iter:
            if record.qname == fq[0].qname:
                fq.append(record)
            else:
                yield tuple(fq)
                fq = [record]
        yield tuple(fq)
