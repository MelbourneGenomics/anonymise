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
        return filename.split("_")

    @staticmethod
    def make_batch_dir(data_dir, batch):
        return os.path.join(data_dir, BATCHES_DIR_NAME, batch, FASTQ_DIR_NAME)


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
