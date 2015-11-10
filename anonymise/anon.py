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
import os
import sys
import random
from argparse import ArgumentParser
import sqlite3
from error import print_error, ERROR_MAKE_DIR
from application import Application
from random_id import make_random_ids, DEFAULT_USED_IDS_DATABASE
from constants import BATCHES_DIR_NAME
from metadata import Metadata, DEFAULT_METADATA_OUT_FILENAME

FASTQ_DIR_NAME = "data"
ANALYSIS_DIR_NAME = "analysis"
ALIGN_DIR_NAME = "align"
VCF_DIR_NAME = "variants"
BAM_SUFFIX = "merge.dedup.realign.recal.bam"
BAI_SUFFIX = "merge.dedup.realign.recal.bai" 
FASTQ_SUFFIX = "fastq.gz"

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
    return parser.parse_args() 


# Assumes application JSON is validated
def create_app_dir(application):
    path = os.path.join(application.fields['application id'], application.fields['request id'])
    try:
        os.makedirs(path)
    except OSError as e:
        print_error("failed to make directory {}".format(PROGRAM_NAME, dir))
        print(e, file=sys.stderr)
        exit(ERROR_MAKE_DIR)
    return path


def get_files_by_type(metadata, make_dir_name, get_sample_from_filename):
    results = [] 
    for batch in metadata.batches:
        dir = make_dir_name(batch)
        all_filenames = os.listdir(dir)
        for filename in all_filenames:
            filename_sample_id = get_sample_from_filename(filename)
            if (filename_sample_id is not None) and \
               (filename_sample_id in metadata.sample_ids):
                full_path = os.path.join(dir, filename)
                results.append(full_path)
    return results


def get_fastq_files(data_dir, metadata):
    def make_dir_name(batch):
        return os.path.join(data_dir, BATCHES_DIR_NAME, batch, FASTQ_DIR_NAME)
    def get_sample_from_filename(filename):
        if filename.endswith(FASTQ_SUFFIX):
            fields = filename.split("_")
            if len(fields) > 0:
                return fields[0]
    return get_files_by_type(metadata, make_dir_name, get_sample_from_filename)


def get_bam_files(data_dir, metadata):
    def make_dir_name(batch):
        return os.path.join(data_dir, BATCHES_DIR_NAME, batch, ANALYSIS_DIR_NAME, ALIGN_DIR_NAME)
    def get_sample_from_filename(filename):
        if filename.endswith(BAM_SUFFIX) or filename.endswith(BAI_SUFFIX):
            fields = filename.split(".")
            if len(fields) > 0:
                filename_sample_id = fields[0]
                return filename_sample_id
    return get_files_by_type(metadata, make_dir_name, get_sample_from_filename)

# XXX fixme
def get_vcf_files(data_dir, metadata):
    return []


def get_files(data_dir, file_types, metadata):
    fastqs = bams = vcfs = []
    if "fastq" in file_types:
        fastqs = get_fastq_files(data_dir, metadata)
    if "bam" in file_types:
        bams = get_bam_files(data_dir, metadata)
    if "vcf" in file_types:
        vcfs = get_vcf_files(data_dir, metadata)
    return fastqs, bams, vcfs


def link_files(application_dir, filepaths):
    for path in filepaths:
        _, filename = os.path.split(path)
        link_name = os.path.join(application_dir, filename)
        os.symlink(path, link_name)


def main():
    args = parse_args()
    with open(args.app) as app_filename:
        # parse and validate the requested data application JSON file
        application = Application(app_filename) 
        # check what data types are allowed for this application
        allowed_data_types = application.allowed_data_types()
        if len(allowed_data_types) > 0:
            # Get all the sample metadata for all requested cohorts
            metadata = Metadata(args.data, application.cohorts())
            # Filter the sample metadata based on patient consent
            metadata.filter_consent(args.consent, allowed_data_types)
            # Find all the file paths for requested file types for each
            # consented sample
            vcfs, bams, fastqs = get_files(args.data, application.file_types(), metadata)
            # Create output directory for the results
            #application_dir = create_app_dir(application)
            if 'Anonymised' in allowed_data_types:
                randomised_ids = make_random_ids(args.usedids, metadata.sample_ids)
                print(randomised_ids)
                # XXX randomise sample IDs in metadata
                meta_data.write(args.metaout)
            #elif 'Re-identifiable' in allowed_data_types:
            #    link_files(application_dir, vcfs + bams + fastqs)
        else:
            print("No data available for this application")
        

if __name__ == '__main__':
    main()

