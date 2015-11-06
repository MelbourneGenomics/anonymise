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
from argparse import ArgumentParser
from jsonschema import validate, ValidationError
from pkg_resources import resource_filename
from collections import namedtuple


COHORTS = ["AML", "EPIL", "CS", "CRC", "CMT"]
METADATA_FILENAME = "samples.txt"
BATCHES_DIR_NAME = "batches"
JSON_SCHEMA = "application_json_schema.txt"
PROGRAM_NAME = "anonymise"
ERROR_MAKE_DIR = 1
ERROR_JSON_SCHEMA_DEFINE = 2
ERROR_JSON_SCHEMA_OPEN = 3
ERROR_JSON_SCHEMA_INVALID = 4
ERROR_INVALID_APPLICATION = 5


def print_error(message):
    print("{}: ERROR: {}".format(PROGRAM_NAME, message), file=sys.stderr)


def parse_args():
    """Orchestrate the anonymisation process for Melbourne Genomics"""
    parser = ArgumentParser(description="Orchestrate the anonymisation process for Melbourne Genomics")
    parser.add_argument("--app", required=True, type=str, help="name of input application JSON file")
    parser.add_argument("--data", required=True, type=str, help="directory containing production data")
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
    return result


# We assume batch filenames are all digits and nothing else
def is_batch_dir(filename):
    return filename.isdigit()

# Metadata desribing samples is in:
#    $datadir/batches/$batchNum/samples.txt
# batch numbers are directory names with three digits in their name
# (seems limiting, so we will assume any directory with all digits
# in its name is a batch)
def get_sample_metadata_for_cohorts(datadir, cohorts):
    '''Return a dictionary mapping batch number to a list of sample
    metadata, for all samples in the desired cohort'''
    batch_sample_infos = {}
    batches_dir = os.path.join(datadir, BATCHES_DIR_NAME)
    batches_dir_contents = os.listdir(batches_dir)
    for file in batches_dir_contents:
        if is_batch_dir(file):
            batch_number = file 
            samples_metadata_path = os.path.join(batches_dir, batch_number, METADATA_FILENAME)
            samples_infos = get_samples_metadata(cohorts, samples_metadata_path)
            batch_sample_infos[batch_number] = samples_infos
    return batch_sample_infos

def get_requested_cohorts(application):
    '''Return a list of all the cohorts requested in an application'''
    return [cohort for cohort in COHORTS if application["condition"][cohort] == "TRUE"]
            

def main():
    args = parse_args()
    with open(args.app) as application_filename:
        application = json.load(application_filename) 
        validate_json(application_filename, application)
        data_available = get_data_available(application) 
        if len(data_available) > 0:
            cohorts = get_requested_cohorts(application)
            batch_sample_infos = get_sample_metadata_for_cohorts(args.data, cohorts)
            for batch in batch_sample_infos:
                print(batch)
                for info in batch_sample_infos[batch]:
                    print(info)
        else:
            print("No data available for this application")
        


if __name__ == '__main__':
    main()

