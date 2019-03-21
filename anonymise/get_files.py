'''
File manipulation code.

The anonymisation program deals with a number of file types
and naming conventions.

This module tries to abstract over the file finding and manipulation
routines.
'''

import os.path
import logging
from metadata import BATCHES_DIR_NAME
from error import print_error, ERROR_BAD_FILENAME

BAM_SUFFIX = "merge.dedup.realign.recal.bam"
BAI_SUFFIX = "merge.dedup.realign.recal.bai"
FASTQ_SUFFIX = "fastq.gz"
# XXX Assuming unfiltered VCF file
VCF_SUFFIX = "merge.dedup.realign.recal.vcf"
FASTQ_DIR_NAME = "data"
ANALYSIS_DIR_NAME = "analysis"
ALIGN_DIR_NAME = "align"
VCF_DIR_NAME = "variants"


def get_files(data_dir, file_types, metadata):
    fastqs = []
    bams = []
    bais = []
    vcfs = []
    if "fastq" in file_types:
        fastqs = get_files_by_type(data_dir, metadata, FASTQ_filename) 
    if "bam" in file_types:
        bams = get_files_by_type(data_dir, metadata, BAM_filename) 
        bais = get_files_by_type(data_dir, metadata, BAI_filename)
    if "vcf" in file_types:
        vcfs = get_files_by_type(data_dir, metadata, VCF_filename) 
    return fastqs, bams, bais, vcfs


def get_files_by_type(datadir, metadata, file_type):
    results = []
    for batch in metadata.batches:
        directory = file_type.make_batch_dir(datadir, batch)
        all_filenames = os.listdir(directory)
        logging.info("Searching for files in: {}".format(directory))
        for filename in all_filenames:
            full_path = os.path.join(directory, filename)
            try:
                file_handler = file_type(full_path)
            except FileTypeException:
                # ignore this file because it does not match what
                # we are looking for
                pass
            else:
                filename_sample_id = file_handler.get_sample_id()
                if filename_sample_id in metadata.sample_ids:
                    results.append(full_path)
    return results


class FileTypeException(Exception):
    pass


class Data_filename(object):
    def __init__(self, absolute_path, suffix):
        if absolute_path.endswith(suffix):
            self.absolute_path = absolute_path
        else:
            raise FileTypeException

    def get_directory(self):
        directory, filename = os.path.split(self.absolute_path)
        return directory

    def get_filename(self):
        directory, filename = os.path.split(self.absolute_path)
        return filename

    def get_sample_id(self):
        sample_id, _rest = self.split_sample_id()
        return sample_id

    # Most file names use dot for fields, can be overridden for special cases
    @staticmethod
    def get_fields(filename):
        return filename.split(".")

    def split_sample_id(self):
        filename = self.get_filename()
        fields = self.get_fields(filename) 
        if len(fields) > 0:
            prefix_len = len(fields[0])
            return fields[0], filename[prefix_len:] 
        else:
            print_error("Cannot find sample ID in filename: {}".format(self.absolute_path))
            exit(BAD_FILENAME)

    def replace_sample_id(self, new_id):
        '''Modifies self.absolute_path'''
        sample_id, rest = self.split_sample_id()
        directory = self.get_directory()
        new_path = os.path.join(directory, new_id + rest)
        self.absolute_path = new_path


class FASTQ_filename(Data_filename):
    def __init__(self, absolute_path):
        Data_filename.__init__(self, absolute_path, FASTQ_SUFFIX)

    @staticmethod
    def get_fields(filename):
        """
        A field is a substring in filename separated by '_'

        """
        return filename.split("_")

    @staticmethod
    def make_batch_dir(data_dir, batch):
        return os.path.join(data_dir, BATCHES_DIR_NAME, batch, FASTQ_DIR_NAME)

    def split_sample_id(self):
        filename = self.get_filename()
        fields = self.get_fields(filename)
        if len(fields) > 0:
            # SAMPLEID_AGRF_111_HHMN7BCXX_TAAGGCGA_L001_R1.fastq.gz
            # get rid of _AGRF_111
            # prefix_len = len(fields[0]) + 1 + len(fields[1]) + 1 + len(fields[2])

            prefix_len = len(fields[0])
            return fields[0], filename[prefix_len:]
        else:
            print_error("Cannot find sample ID in filename: {}".format(self.absolute_path))
            exit(BAD_FILENAME)

    def get_sample_id(self):
        filename = self.get_filename()
        fields = self.get_fields(filename)
        return fields[0]

    def replace_sample_id(self, new_id):
        filename = self.get_filename()
        fields = self.get_fields(filename)
        directory = self.get_directory()
        prefix_len = len(fields[0])
        new_path = os.path.join(directory, new_id + filename[prefix_len:])
        self.absolute_path = new_path

    def replace_field(self, new_id, field_index_start, field_index_end=None):
        """
        A field is a substring in filename separated by '_'

        """
        filename = self.get_filename()
        fields = self.get_fields(filename)
        directory = self.get_directory()

        if field_index_end == None:
            field_index_end = field_index_start

        # batch id is in this format AGRF_001
        prefix = ''
        for index in range(field_index_start):
            # len of the field + len of '_'
            prefix += fields[index] + '_'

        replaced = ''
        for index in range(field_index_start, field_index_end+1):
            replaced += fields[index] + '_'

        suffix_start = len(prefix + replaced) - 1

        new_path = os.path.join(directory, prefix + new_id + filename[suffix_start:])
        self.absolute_path = new_path


class BAM_filename(Data_filename):
    def __init__(self, absolute_path):
        Data_filename.__init__(self, absolute_path, BAM_SUFFIX)

    @staticmethod
    def make_batch_dir(data_dir, batch):
        return os.path.join(data_dir, BATCHES_DIR_NAME, batch, ANALYSIS_DIR_NAME, ALIGN_DIR_NAME)


class BAI_filename(Data_filename):
    def __init__(self, absolute_path):
        Data_filename.__init__(self, absolute_path, BAI_SUFFIX)

    @staticmethod
    def make_batch_dir(data_dir, batch):
        return os.path.join(data_dir, BATCHES_DIR_NAME, batch, ANALYSIS_DIR_NAME, ALIGN_DIR_NAME)


class VCF_filename(Data_filename):
    def __init__(self, absolute_path):
        Data_filename.__init__(self, absolute_path, VCF_SUFFIX)

    @staticmethod
    def make_batch_dir(data_dir, batch):
        return os.path.join(data_dir, BATCHES_DIR_NAME, batch, ANALYSIS_DIR_NAME, VCF_DIR_NAME)
