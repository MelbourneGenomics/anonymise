import os
import csv
from constants import BATCHES_DIR_NAME

METADATA_FILENAME = "samples.txt"
DEFAULT_METADATA_OUT_FILENAME = "samples.out.txt"
METADATA_HEADINGS = ['Batch', 'Sample_ID', 'DNA_Tube_ID', 'Sex',
    'DNA_Concentration', 'DNA_Volume', 'DNA_Quantity', 'DNA_Quality',
    'DNA_Date', 'Cohort', 'Sample_Type', 'Fastq_Files',
    'Prioritised_Genes', 'Consanguinity', 'Variants_File',
    'Pedigree_File', 'Ethnicity', 'VariantCall_Group',
    'Capture_Date', 'Sequencing_Date', 'Mean_Coverage',
    'Duplicate_Percentage', 'Machine_ID', 'DNA_Extraction_Lab',
    'Sequencing_Lab', 'Exome_Capture', 'Library_Preparation',
    'Barcode_Pool_Size', 'Read_Type', 'Machine_Type',
    'Sequencing_Chemistry', 'Sequencing_Software',
    'Demultiplex_Software', 'Hospital_Centre',
    'Sequencing_Contact', 'Pipeline_Contact', 'Notes']

# Metadata desribing samples is in:
#    $datadir/batches/$batchNum/samples.txt
# batch numbers are directory names with three digits in their name
# (seems limiting, so we will assume any directory with all digits
# in its name is a batch)

class Metadata(object):

    def __init__(self, datadir, cohorts):
        '''Return a dictionary mapping batch number to a list of sample
        metadata, for all samples in the desired cohort'''
        self.samples = []
        batches_dir = os.path.join(datadir, BATCHES_DIR_NAME)
        batches_dir_contents = os.listdir(batches_dir)
        for file in batches_dir_contents:
            if is_batch_dir(file):
                batch_number = file
                metadata_path = os.path.join(batches_dir, batch_number, METADATA_FILENAME)
                samples_infos = get_batch_metadata(cohorts, metadata_path)
                self.samples.extend(samples_infos)
        # set of all batches used by all samples in all cohorts
        self.batches = { sample['Batch'] for sample in self.samples }
        # set of all sample IDs used by all samples in all cohorts
        self.sample_ids = { sample['Sample_ID'] for sample in self.samples }

    def filter_consent(self, consent_file, allowed_data_types):
        # XXX fixme, don't forget to filter the self.sample_ids as well
        pass


def get_batch_metadata(cohorts, metadata_filename):
    result = []
    with open(metadata_filename) as metadata_file:
        reader = csv.DictReader(metadata_file, delimiter='\t')
        for row in reader:
            if row["Cohort"] in cohorts:
                result.append(row)
        return result


# We assume batch filenames are all digits and nothing else
def is_batch_dir(filename):
    return filename.isdigit()

def write_metadata(output_filename, cohort_sample_infos):
    with open(output_filename, "w") as output_file:
        writer = csv.DictWriter(output_file, METADATA_HEADINGS, delimiter='\t')
        writer.writeheader()
        for batch in sorted(cohort_sample_infos):
            for info in cohort_sample_infos[batch]:
                writer.writerow(info)
