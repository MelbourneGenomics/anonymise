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
    - check if BAI files need anonymising
    - check consent
    - file uploading
'''

from __future__ import print_function
import os
import sys
import random
import logging
import string
from argparse import ArgumentParser
import sqlite3
from error import print_error, ERROR_MAKE_DIR, ERROR_BAD_ALLOWED_DATA, ERROR_MD5, ERROR_RANDOMISE_ID
from application import Application
from random_id import make_random_ids, DEFAULT_USED_IDS_DATABASE
from constants import BATCHES_DIR_NAME
from metadata import Metadata, DEFAULT_METADATA_OUT_FILENAME
from get_files import get_files, FileTypeException, VCF_filename, BAM_filename, BAI_filename, FASTQ_filename
from vcf_edit import vcf_edit
from bam_edit import bam_edit
from version import program_version
from subprocess import call


DEFAULT_MD5_COMMAND = "openssl md5"

def parse_args():
    """Orchestrate the anonymisation process for Melbourne Genomics"""
    parser = ArgumentParser(description="Orchestrate the anonymisation process for Melbourne Genomics, version {}".format(program_version))
    parser.add_argument('--version', action='version', version='%(prog)s ' + program_version)
    parser.add_argument("--app", required=True, type=str,
        help="Name of input application JSON file")
    parser.add_argument("--data", required=True,
        type=str, help="Directory containing production data")
    parser.add_argument("--metaout",
        required=False, default=DEFAULT_METADATA_OUT_FILENAME, type=str,
        help="Name of output metadatafile, defaults to {}".format(DEFAULT_METADATA_OUT_FILENAME))
    parser.add_argument("--usedids", required=False,
        default=DEFAULT_USED_IDS_DATABASE, type=str,
        help="Sqlite3 databsae of previously used randomised sample ids defaults to {}".format(DEFAULT_USED_IDS_DATABASE))
    parser.add_argument("--consent", required=True, type=str,
        help="File path of consent metadata")
    parser.add_argument("--md5", required=False, type=str, default=DEFAULT_MD5_COMMAND,
        help="MD5 checksum command, defaults to {}".format(DEFAULT_MD5_COMMAND))
    parser.add_argument('--log', metavar='FILE', type=str,
        help='Log progress in FILENAME, defaults to stdout')
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
    output_files = []
    for path in filepaths:
        _, filename = os.path.split(path)
        link_name = os.path.join(application_dir, filename)
        output_files.append(link_name)
        os.symlink(path, link_name)
    return output_files


def anonymise_files(filenames, randomised_ids, application_dir, filename_type, file_editor=None):
    output_files = []

    randomised_batch_ids = {}
    # randomised_flowcell_ids = {}

    for file_path in filenames:
        try:
            file_handler = filename_type(file_path)
        except FileTypeException:
            # skip this file 
            continue 
        else:
            # sample id
            old_id = file_handler.get_sample_id()
            if old_id is None or old_id not in randomised_ids:
                print_error("Cannot randomise this file: {}".format(file_path))
                exit(ERROR_RANDOMISE_ID)
            else:
                new_id = str(randomised_ids[old_id])
                file_handler.replace_sample_id(new_id)

            # replace batch id with XXXXX: AGRF_024
            # 010108101_AGRF_024_HG3JKBCXX_CGTACTAG_L001_R1.fastq.gz
            fields = file_handler.get_fields(file_handler.get_filename())
            old_batch_id = fields[1] + '_' + fields[2]
            if old_batch_id not in randomised_batch_ids:
                new_batch_id = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(5))
                randomised_batch_ids[old_batch_id] = new_batch_id

            file_handler.replace_field(randomised_batch_ids[old_batch_id], 1, 2)

            # # replace flow cell id: HG3JKBCXX
            # # 010108101_AGRF_024_HG3JKBCXX_CGTACTAG_L001_R1.fastq.gz
            # # after batch renamed from AGRF_027 to XXXXX
            # # 010108101_xxxxx_HG3JKBCXX_CGTACTAG_L001_R1.fastq.gz
            # fields = file_handler.get_fields(file_handler.get_filename())
            # old_flowcell_id = fields[2]
            # if old_flowcell_id not in randomised_flowcell_ids:
            #     new_flowcell_id = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(9))
            #     randomised_flowcell_ids[old_flowcell_id] = new_flowcell_id
            #
            # file_handler.replace_field(randomised_flowcell_ids[old_flowcell_id], 2)

            # file_handler has updated file name at this point
            new_filename = file_handler.get_filename()
            new_path = os.path.join(application_dir, new_filename)
            output_files.append(new_path)
            if file_editor is not None:
                file_editor(old_id, new_id, file_path, new_path)
                logging.info("Anonymised {} to {}".format(file_path, new_path))
            else:
                os.symlink(file_path, new_path)
                logging.info("Linked {} to {}".format(new_path, file_path))

    return output_files


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


def md5_files(md5_command, filenames):
    for filename in filenames:
        output_filename = filename + ".md5"
        logging.info("{} {} > {}".format(md5_command, filename, output_filename))
        with open(output_filename, "w") as out_file:
            try:
                command = md5_command.split() + [filename]
                call(command, stdout=out_file)
            except OSError as e:
                print_error(e)
                exit(ERROR_MD5)


def main():
    args = parse_args()
    init_log(args.log)
    with open(args.app) as app_file:
        # parse and validate the requested data application JSON file
        application = Application(app_file) 
        logging.info("Input data application parsed: {}".format(args.app))
        # Create output directory for the results
        application_dir = create_app_dir(application)
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
            logging.warning("Consent not handled yet. FIXME")
            # Find all the file paths for requested file types for each
            # consented sample
            requested_file_types = application.file_types()
            logging.info("Requested file types: {}".format(' '.join(requested_file_types)))
            fastqs, bams, bais, vcfs = get_files(args.data, requested_file_types, metadata)
            logging.info("VCF files selected:\n{}".format('\n'.join(vcfs)))
            logging.info("BAM files selected:\n{}".format('\n'.join(bams)))
            logging.info("BAI files selected:\n{}".format('\n'.join(bais)))
            logging.info("FASTQ files selected:\n{}".format('\n'.join(fastqs)))
            output_files = []
            if 'Anonymised' in allowed_data_types:
                # generate random IDs for all output samples
                randomised_ids = make_random_ids(args.usedids, metadata.sample_ids)
                metadata.anonymise(randomised_ids)
                metadata.write(args.metaout)
                logging.info("Anonymised metadata written to: {}".format(args.metaout))
                new_vcfs = anonymise_files(vcfs, randomised_ids, application_dir, VCF_filename, vcf_edit)
                new_bams = anonymise_files(bams, randomised_ids, application_dir, BAM_filename, bam_edit)
                # BAIs and FASTQs are just sym-linked to output with randomised name
                new_bais = anonymise_files(bais, randomised_ids, application_dir, BAI_filename)
                new_fastqs = anonymise_files(fastqs, randomised_ids, application_dir, FASTQ_filename)
                output_files.extend(new_vcfs + new_bams + new_bais + new_fastqs)
                logging.info("Output files are anonymised")
            elif 'Re-identifiable' in allowed_data_types:
                new_links = link_files(application_dir, vcfs + bams + bais + fastqs)
                output_files.extend(new_links)
                logging.info("Files linked in directory: {}".format(application_dir))
                metadata.write(args.metaout)
                logging.info("Output files are re-identifiable")
            else:
                print_error("Allowed data is neither anonymised nor re-identifiable")
                exit(ERROR_BAD_ALLOWED_DATA)
            logging.info("Generating MD5 checksums on output files")
            md5_files(args.md5, output_files)
        else:
            logging.warning("No data available for this application")
        

if __name__ == '__main__':
    main()
