import sys
import resource
import logging

from whatshap.bam import (
    AlignmentFileNotIndexedError,
    EmptyAlignmentFileError,
    SampleNotFoundError,
    ReferenceNotFoundError,
)
from whatshap.variants import ReadSetReader, ReadSetError
from whatshap.utils import IndexedFasta, FastaNotIndexedError, detect_file_format
from whatshap.core import ReadSet
from whatshap.vcf import VcfReader

logger = logging.getLogger(__name__)


class CommandLineError(Exception):
    """An anticipated command-line error occurred. This ends up as a user-visible error message"""

class PhasedInputReader:
    def __init__(
        self,
        vcf_paths,
        numeric_sample_ids,
        indels,
        **kwargs,  # passed to ReadSetReader constructor
    ):
        self._vcf_paths = self._split_input_file_list(vcf_paths)

        # TODO exit stack!
        self._numeric_sample_ids = numeric_sample_ids

        vcf_readers = [VcfReader(f, indels=indels, phases=True) for f in self._vcf_paths]

        self._vcf_readers = vcf_readers

        if not self._vcf_readers:
            self._vcfs = []
        else:
            self._vcfs = None  # None means uninitialized, call .read_vcf() first

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self._fasta is not None:
            self._fasta.close()

    @property
    def has_vcfs(self):
        return bool(self._vcf_paths)


    @staticmethod
    def _split_input_file_list(paths):
        vcfs = []
        for path in paths:
            try:
                file_format = detect_file_format(path)
            except OSError as e:
                raise CommandLineError(e)
            if file_format == "VCF":
                vcfs.append(path)
            else:
                raise CommandLineError("Unable to determine type of input file {!r}".format(path))
        return vcfs

    def read_vcfs(self):
        # Read phase information provided as VCF files, if provided.
        # TODO: do this chromosome- and/or sample-wise on demand to save memory.
        self._vcfs = []
        for reader in self._vcf_readers:
            # create dict mapping chromosome names to VariantTables
            m = dict()
            logger.info("Reading phased blocks from %r", reader.path)
            for variant_table in reader:
                m[variant_table.chromosome] = variant_table
            self._vcfs.append(m)

    def read(self, chromosome, variants, sample, *, read_vcf=True, regions=None):
        """
        Return a pair (readset, vcf_source_ids) where readset is a sorted ReadSet.

        Set read_vcf to False to not read phased blocks from the VCFs
        """
        readset = ReadSet()
        vcf_source_ids = set()
        if read_vcf:
            # TODO this is a bit clumsy
            if self._vcfs is None:
                raise ValueError("call PhasedInputReader.read_vcfs() first")
            # Add phasing information from VCF files, if present
            sample_id = self._numeric_sample_ids[sample]
            for i, vcf in enumerate(self._vcfs):
                if chromosome in vcf:
                    variant_table = vcf[chromosome]
                    source_id = i
                    vcf_source_ids.add(source_id)
                    for read in variant_table.phased_blocks_as_reads(
                        sample, variants, source_id, sample_id
                    ):
                        readset.add(read)

        # TODO is this necessary?
        for read in readset:
            read.sort()
        readset.sort()

        logger.info(
            "Found %d reads covering %d variants", len(readset), len(readset.get_positions())
        )
        return readset, vcf_source_ids


def log_memory_usage(include_children=False):
    if sys.platform == "linux":
        if include_children:
            memory_kb = (
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            )
        else:
            memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum memory usage: %.3f GB", memory_kb / 1e6)
