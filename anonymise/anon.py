#!/usr/bin/env python

'''
Orchestrate the anonymisation process for Melbourne Genomics

1) For each request, create a directory Application ID_request ID e.g. APPID3_REQID6
2) Identify what data is available to researcher
3) Identify all patients for conditions requested, or a particular patient in a "condition"
4) Create metadata file for patients in "condition", and edit columns accordingly. 
5) Symbolic link to re-identifiable data, or anonymise data depending on request combination.
6) md5.txt
7) Create upload files, and send links to requestor and PI of Application ID

Usage:

Authors: Bernie Pope, Gayle Philip
'''

from __future__ import print_function
import json
import os
import sys
import csv
import random
from argparse import ArgumentParser
from jsonschema import validate, ValidationError
from pkg_resources import resource_filename
from collections import namedtuple
import sqlite3


FILE_TYPES = ["fastq", "bam", "vcf"]
COHORTS = ["AML", "EPIL", "CS", "CRC", "CMT"]
METADATA_FILENAME = "samples.txt"
BATCHES_DIR_NAME = "batches"
FASTQ_DIR_NAME = "data"
ANALYSIS_DIR_NAME = "analysis"
ALIGN_DIR_NAME = "align"
VCF_DIR_NAME = "variants"
JSON_SCHEMA = "application_json_schema.txt"
PROGRAM_NAME = "anonymise"
ERROR_MAKE_DIR = 1
ERROR_JSON_SCHEMA_DEFINE = 2
ERROR_JSON_SCHEMA_OPEN = 3
ERROR_JSON_SCHEMA_INVALID = 4
ERROR_INVALID_APPLICATION = 5
BAM_SUFFIX = "merge.dedup.realign.recal.bam"
BAI_SUFFIX = "merge.dedup.realign.recal.bai" 
DEFAULT_METADATA_OUT_FILENAME = "samples.out.txt"
DEFAULT_USED_IDS_DATABASE = "used_random_sample_ids.db"
# Start the minimum random sample ID at some non-low number.
# Users might think there is something special about low number IDs
MIN_RANDOM_ID = 1000
# XXX This may not be portable: we might want to pick a high
# maximum value ourselves
MAX_RANDOM_ID = sys.maxint


def print_error(message):
    print("{}: ERROR: {}".format(PROGRAM_NAME, message), file=sys.stderr)


def parse_args():
    """Orchestrate the anonymisation process for Melbourne Genomics"""
    parser = ArgumentParser(description="Orchestrate the anonymisation process for Melbourne Genomics")
    parser.add_argument("--app", required=True, type=str,
        help="name of input application JSON file")
    parser.add_argument("--data", required=True,
        type=str, help="directory containing production data")
    parser.add_argument("--metaout",
        required=False, default=DEFAULT_METADATA_OUT_FILENAME, type=str,
        help="name of output metadatafile, defaults to {}".format(DEFAULT_METADATA_OUT_FILENAME))
    parser.add_argument("--usedids", required=False,
        default=DEFAULT_USED_IDS_DATABASE, type=str,
        help="sqlite3 databsae of previously used randomised sample ids defaults to {}".format(DEFAULT_USED_IDS_DATABASE))
    return parser.parse_args() 


def validate_json(application_filename, application):
    '''Validate input JSON file against schema. This function exits the program if
    the validation fails, otherwise it returns, but does not have a return result.
    '''
    try:
        json_schema_filename = resource_filename(PROGRAM_NAME, os.path.join('data', JSON_SCHEMA))
    except Exception as e:
        print_error("JSON schema file not defined, program not installed correctly")
        exit(ERROR_JSON_SCHEMA_DEFINE)
    try:
        json_schema_file = open(json_schema_filename)
        json_schema = json.load(json_schema_file)
    except OSError as e:
        print_error("Cannot open JSON schema file: {}".format(json_schema_filename))
        print_error(str(e), file=sys.stderr)
        exit(ERROR_JSON_SCHEMA_OPEN)
    finally:
        json_schema_file.close()
    try:
        validate(application, json_schema)
    except ValidationError as e:
        print_error("JSON input file is not valid: {}".format(application_filename))
        print(e, file=sys.stderr)
        exit(ERROR_JSON_SCHEMA_INVALID)


# Assumes application JSON is validated
def create_app_dir(application):
    path = os.path.join(application['application id'], application['request id'])
    try:
        os.makedirs(path)
    except OSError as e:
        print_error("failed to make directory {}".format(PROGRAM_NAME, dir))
        print(e, file=sys.stderr)
        exit(ERROR_MAKE_DIR)
    return path


def get_data_available(application):
    '''Based on the input application, decide what data is available
    to the requestor.'''
    request = application_to_request(application)
    try:
        result_combinations = REQUEST_COMBINATIONS[request]
    except KeyError:
        print_error("Application does not have a valid interpretation")
        print(format(json.dumps(application, indent=4)), file=sys.stderr)
        exit(ERROR_INVALID_APPLICATION)
    return result_combinations 
         

def application_to_request(application):
    return Request(ethics=application['ethics'],
               research_related=application['research_related'],
               filter_results=application['filter_results'],
               method_dev=application['method_dev'],
               return_results=application['return_results'],
               genes_approved=application['genes_approved'],
               reconsent_patient=application['reconsent_patient'])


Request = namedtuple("Request",
    [ "ethics", "research_related", "filter_results"
    , "method_dev", "return_results", "genes_approved"
    , "reconsent_patient" ])


# We define this table of combinations to be explicit
# The equivalent conditional statement is hard to read and
# easy to get wrong.
REQUEST_COMBINATIONS = {

    Request(ethics='MGHA',
        research_related='TRUE',
        filter_results='TRUE',
        method_dev='FALSE',
        return_results='FALSE',
        genes_approved='FALSE',
        reconsent_patient='FALSE'): ['Re-identifiable'],

    Request(ethics='MGHA',
        research_related='TRUE',
        filter_results='FALSE',
        method_dev='FALSE',
        return_results='FALSE',
        genes_approved='TRUE',
        reconsent_patient='FALSE'): ['Re-identifiable'],

    Request(ethics='MGHA',
        research_related='TRUE',
        filter_results='FALSE',
        method_dev='FALSE',
        return_results='FALSE',
        genes_approved='FALSE',
        reconsent_patient='FALSE'): [],

    Request(ethics='MGHA',
        research_related='FALSE',
        filter_results='FALSE',
        method_dev='TRUE',
        return_results='FALSE',
        genes_approved='FALSE',
        reconsent_patient='FALSE'): ['Anonymised'],

    Request(ethics='MGHA',
        research_related='FALSE',
        filter_results='FALSE',
        method_dev='FALSE',
        return_results='FALSE',
        genes_approved='FALSE',
        reconsent_patient='FALSE'): [],

    Request(ethics='HREC',
        research_related='FALSE',
        filter_results='FALSE',
        method_dev='FALSE',
        return_results='FALSE',
        genes_approved='FALSE',
        reconsent_patient='FALSE'): ['Anonymised', 'Future'],

    Request(ethics='HREC',
        research_related='TRUE',
        filter_results='FALSE',
        method_dev='FALSE',
        return_results='FALSE',
        genes_approved='FALSE',
        reconsent_patient='FALSE'): ['Anonymised', 'Future'],

    Request(ethics='HREC',
        research_related='TRUE',
        filter_results='TRUE',
        method_dev='FALSE',
        return_results='FALSE',
        genes_approved='FALSE',
        reconsent_patient='TRUE'): ['Re-identifiable'],

    Request(ethics='HREC',
        research_related='TRUE',
        filter_results='TRUE',
        method_dev='FALSE',
        return_results='FALSE',
        genes_approved='FALSE',
        reconsent_patient='FALSE'): ['Re-identifiable', 'Future'],

    Request(ethics='HREC',
        research_related='TRUE',
        filter_results='TRUE',
        method_dev='FALSE',
        return_results='TRUE',
        genes_approved='FALSE',
        reconsent_patient='FALSE'): ['Re-identifiable', 'Future', 'Return'],
}

def get_samples_metadata(cohorts, metadata_filename):
    result = []
    with open(metadata_filename) as metadata_file:
        reader = csv.DictReader(metadata_file, delimiter='\t')
        for row in reader:
            if row["Cohort"] in cohorts:
                result.append(row)
        return reader.fieldnames, result


# We assume batch filenames are all digits and nothing else
def is_batch_dir(filename):
    return filename.isdigit()

# Metadata desribing samples is in:
#    $datadir/batches/$batchNum/samples.txt
# batch numbers are directory names with three digits in their name
# (seems limiting, so we will assume any directory with all digits
# in its name is a batch)
def get_samples_metadata_for_cohorts(datadir, cohorts):
    '''Return a dictionary mapping batch number to a list of sample
    metadata, for all samples in the desired cohort'''
    sample_ids = set([])
    batch_sample_infos = {}
    batches_dir = os.path.join(datadir, BATCHES_DIR_NAME)
    batches_dir_contents = os.listdir(batches_dir)
    for file in batches_dir_contents:
        if is_batch_dir(file):
            batch_number = file 
            samples_metadata_path = os.path.join(batches_dir, batch_number, METADATA_FILENAME)
            fieldnames, samples_infos = get_samples_metadata(cohorts, samples_metadata_path)
            batch_sample_infos[batch_number] = samples_infos
            for info in samples_infos:
                sample_ids.add(info['Sample_ID'])
    return fieldnames, batch_sample_infos, sample_ids

def get_requested_cohorts(application):
    '''Return a list of all the cohorts requested in an application'''
    return [cohort for cohort in COHORTS if application["condition"][cohort] == "TRUE"]
            
def get_requested_file_types(application):
    '''Return a list of all the file types requested in an application'''
    return [file_type for file_type in FILE_TYPES if application["file types"][file_type] == "TRUE"]

def get_fastq_files(data_dir, batch_sample_infos):
    results = [] 
    for batch in batch_sample_infos:
        fastq_dir = os.path.join(data_dir, BATCHES_DIR_NAME, batch, FASTQ_DIR_NAME)
        requested_sample_ids = set() 
        for sample in batch_sample_infos[batch]:
            requested_sample_ids.add(sample["Sample_ID"])
        all_filenames = os.listdir(fastq_dir)
        for filename in all_filenames:
            fields = filename.split("_")
            if len(fields) > 0:
                filename_sample_id = fields[0]
                if filename_sample_id in requested_sample_ids:
                    full_fastq_path = os.path.join(fastq_dir, filename)
                    results.append(full_fastq_path)
    return results

def get_bam_files(data_dir, batch_sample_infos):
    results = [] 
    for batch in batch_sample_infos:
        bam_dir = os.path.join(data_dir, BATCHES_DIR_NAME, batch, ANALYSIS_DIR_NAME, ALIGN_DIR_NAME)
        requested_sample_ids = set() 
        for sample in batch_sample_infos[batch]:
            requested_sample_ids.add(sample["Sample_ID"])
        all_filenames = os.listdir(bam_dir)
        for filename in all_filenames:
            head, tail = os.path.split(filename)
            if tail.endswith(BAM_SUFFIX) or tail.endswith(BAI_SUFFIX):
                fields = filename.split(".")
                if len(fields) > 0:
                    filename_sample_id = fields[0]
                    if filename_sample_id in requested_sample_ids:
                        full_bam_path = os.path.join(bam_dir, filename)
                        results.append(full_bam_path)
    return results

# XXX fixme
def get_vcf_files(data_dir, batch_sample_infos):
    return []

def get_files(data_dir, batch_sample_infos, file_types):
    fastqs = bams = vcfs = []
    if "fastq" in file_types:
        fastqs = get_fastq_files(data_dir, batch_sample_infos)
    if "bam" in file_types:
        bams = get_bam_files(data_dir, batch_sample_infos)
    if "vcf" in file_types:
        vcfs = get_vcf_files(data_dir, batch_sample_infos)
    return fastqs, bams, vcfs

def link_files(application_dir, filepaths):
    for path in filepaths:
        _, filename = os.path.split(path)
        link_name = os.path.join(application_dir, filename)
        os.symlink(path, link_name)

def output_metadata(output_filename, fieldnames, batch_sample_infos):
    with open(output_filename, "w") as output_file:
        writer = csv.DictWriter(output_file, fieldnames, delimiter='\t')
        writer.writeheader()
        for batch in sorted(batch_sample_infos):
            for info in batch_sample_infos[batch]:
                writer.writerow(info)

def make_random_id():
    return random.randint(MIN_RANDOM_ID, MAX_RANDOM_ID)

# We keep an sqlite database of previously used random IDs. This ensures that
# we never re-use the same random ID. It may seem like overkill to use a
# database, but the benefit is that sqlite will handle the locking of
# database access, so that we do not have race conditions if multiple instances
# of this program run at the same time.
# If the database does not exist we will create a new empty one.
def make_random_ids(used_ids_database, sample_ids):
    conn = sqlite3.connect(used_ids_database)
    cursor = conn.cursor()
    cursor.execute('CREATE TABLE IF NOT EXISTS unique_ids (id integer)')
    # The set of IDs that have already been used
    used_ids = set([])
    # Grab all the IDs from the database
    for (next_used_id,) in cursor.execute('SELECT * from unique_ids'):
        used_ids.add(next_used_id)
    # The newly generated IDs by the call to this function
    new_ids = []
    # A list of pairs containing the original ID and its new randomised ID
    result = []
    for old_sample in sample_ids:
        # Make sure the newly generated ID has not been seen before
        new_id = make_random_id()
        while new_id in used_ids:
            new_id = make_random_id()
        result.append((old_sample, new_id))
        # Record this new ID in the set of previously used IDs so we don't
        # use it again
        used_ids.add(new_id)
        new_ids.append(new_id)
    # Write the newly created IDs out to the database
    # XXX should be able to do this as a single INSERT statement
    for new_id in new_ids:
        cursor.execute('INSERT into unique_ids (id) VALUES ({})'.format(new_id))
    conn.commit()
    return result


def main():
    args = parse_args()
    with open(args.app) as application_filename:
        application = json.load(application_filename) 
        validate_json(application_filename, application)
        data_available = get_data_available(application) 
        if len(data_available) > 0:
            cohorts = get_requested_cohorts(application)
            fieldnames, batch_sample_infos, sample_ids = get_samples_metadata_for_cohorts(args.data, cohorts)
            file_types = get_requested_file_types(application)
            vcfs, bams, fastqs = get_files(args.data, batch_sample_infos, file_types)
            #application_dir = create_app_dir(application)
            #link_files(application_dir, vcfs + bams + fastqs)
            output_metadata(args.metaout, fieldnames, batch_sample_infos)
            randomised_ids = make_random_ids(args.usedids, sample_ids)
            print(randomised_ids)
        else:
            print("No data available for this application")
        

if __name__ == '__main__':
    main()

