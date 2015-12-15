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

TODO:
    - decide if we are using filtered or unfiltered VCF files
'''

from __future__ import print_function
import os
import sys
import random
import logging
from argparse import ArgumentParser
import sqlite3
from error import print_error, ERROR_MAKE_DIR, ERROR_BAD_ALLOWED_DATA
from application import Application
from random_id import make_random_ids, DEFAULT_USED_IDS_DATABASE
from constants import BATCHES_DIR_NAME
from metadata import Metadata, DEFAULT_METADATA_OUT_FILENAME
from get_files import get_files, FileTypeException, VCF_filename, BAM_filename
from vcf_edit import vcf_edit
from bam_edit import bam_edit


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
    parser.add_argument("--consent", required=True, type=str,
        help="file path of consent metadata{}")
    parser.add_argument('--log', metavar='FILE', type=str, \
        help='Log progress in FILENAME, defaults to stdout.')
    return parser.parse_args() 


# Assumes application JSON is validated
def create_app_dir(application):
    path = os.path.join(application.fields['application id'], application.fields['request id'])
    try:
        os.makedirs(path)
    except OSError as e:
        print_error("failed to make directory {}".format(path))
        print(e, file=sys.stderr)
        exit(ERROR_MAKE_DIR)
    return path


def link_files(application_dir, filepaths):
    for path in filepaths:
        _, filename = os.path.split(path)
        link_name = os.path.join(application_dir, filename)
        os.symlink(path, link_name)


def anonymise_files(filenames, randomised_ids, application_dir, filename_type, file_editor):
    for file_path in filenames:
        try:
            file_handler = filename_type(file_path)
        except FileTypeException:
            # skip this file 
            continue 
        else:
            old_id = file_handler.get_sample_id()
            if old_id is None or old_id not in randomised_ids:
                print_error("Cannot randomise this file: {}".format(file_path))
                exit(ERROR_RANDOMISE_ID)
            else:
                new_id = str(randomised_ids[old_id])
                file_handler.replace_sample_id(new_id)
                new_filename = file_handler.get_filename()
                new_path = os.path.join(application_dir, new_filename)
                file_editor(old_id, new_id, file_path, new_path)
                logging.info("Anonymised {} to {}".format(file_path, new_path))


def init_log(log_file):
    '''Set up log output, if log_file is None, output does to stderr'''
    logging.basicConfig(
        filename=log_file,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('Program started')
    logging.info('Command line: {0}'.format(' '.join(sys.argv)))


def main():
    args = parse_args()
    init_log(args.log)
    with open(args.app) as app_file:
        # parse and validate the requested data application JSON file
        application = Application(app_file) 
        logging.info("Input data application parsed: {}".format(args.app))
        # check what data types are allowed for this application
        allowed_data_types = application.allowed_data_types()
        logging.info("Allowed data types: {}".format(' '.join(allowed_data_types)))
        if len(allowed_data_types) > 0:
            # Get all the sample metadata for all requested cohorts
            requested_cohorts = application.cohorts()
            metadata = Metadata(args.data, requested_cohorts)
            logging.info("Metadata collected for requested cohorts: {}".format(' '.join(requested_cohorts)))
            metadata_sample_ids = sorted(metadata.get_sample_ids())
            logging.info("Metadata for sample IDs: {}".format(' '.join(metadata_sample_ids)))
            # Filter the sample metadata based on patient consent
            metadata.filter_consent(args.consent, allowed_data_types)
            logging.info("Metadata filtered by consent")
            # Find all the file paths for requested file types for each
            # consented sample
            requested_file_types = application.file_types()
            logging.info("Requested file types: {}".format(' '.join(requested_file_types)))
            fastqs, bams, bais, vcfs = get_files(args.data, requested_file_types, metadata)
            logging.info("VCF files:\n{}".format('\n'.join(vcfs)))
            logging.info("BAM files:\n{}".format('\n'.join(bams)))
            logging.info("BAI files:\n{}".format('\n'.join(bais)))
            logging.info("FASTQ files:\n{}".format('\n'.join(fastqs)))
            # Create output directory for the results
            application_dir = create_app_dir(application)
            if 'Anonymised' in allowed_data_types:
                randomised_ids = make_random_ids(args.usedids, metadata.sample_ids)
                # XXX randomise sample IDs in metadata
                metadata.anonymise(randomised_ids)
                metadata.write(args.metaout)
                logging.info("Anonymised metadata written to: {}".format(args.metaout))
                anonymise_files(vcfs, randomised_ids, application_dir, VCF_filename, vcf_edit)
                anonymise_files(bams, randomised_ids, application_dir, BAM_filename, bam_edit)
            elif 'Re-identifiable' in allowed_data_types:
                link_files(application_dir, vcfs + bams_bais + fastqs)
                logging.info("Files linked in directory: {}".format(application_dir))
                metadata.write(args.metaout)
            else:
                print_error("Allowed data is neither anonymised nor re-identifiable")
                exit(ERROR_BAD_ALLOWED_DATA)
        else:
            print("No data available for this application")
        

if __name__ == '__main__':
    main()

