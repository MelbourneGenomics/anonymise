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
from argparse import ArgumentParser
from jsonschema import validate, ValidationError
from pkg_resources import resource_filename
import sys

JSON_SCHEMA = "application_json_schema.txt"
PROGRAM_NAME = "anonymise"
ERROR_MAKE_DIR = 1
ERROR_JSON_SCHEMA_DEFINE = 2
ERROR_JSON_SCHEMA_OPEN = 3
ERROR_JSON_SCHEMA_INVALID = 4

def print_error(message):
    print("{}: ERROR: {}".format(PROGRAM_NAME, message), file=sys.stderr)


def parse_args():
    """Orchestrate the anonymisation process for Melbourne Genomics"""
    parser = ArgumentParser(description="Orchestrate the anonymisation process for Melbourne Genomics")
    parser.add_argument("--app", required=True, type=str, help="name of input application JSON file")
    return parser.parse_args() 


def validate_json(application_filename, application):
    '''Validate input JSON file against schema. This function exits the program if
    the validation fails, otherwise it returns, but does not have a return result.
    '''
    try:
        json_schema_filename = resource_filename(PROGRAM_NAME, os.path.join('data', JSON_SCHEMA))
    except Exception as e:
        print(e)
        print_error("JSON schema file not defined, program not installed correctly")
        exit(ERROR_JSON_SCHEMA_DEFINE)
    try:
        json_schema_file = open(json_schema_filename)
        json_schema = json.load(json_schema_file)
    except OSError as e:
        print_error("Cannot open JSON schema file: {}".format(json_schema_filename))
        print_error(str(e))
        exit(ERROR_JSON_SCHEMA_OPEN)
    finally:
        json_schema_file.close()
    try:
        validate(application, json_schema)
    except ValidationError as e:
        print_error("JSON input file is not valid: {}".format(application_filename))
        print(e)
        exit(ERROR_JSON_SCHEMA_INVALID)

# Assumes application JSON is validated
def create_app_dir(application):
    path = os.path.join(application['application id'], application['request id'])
    try:
        os.makedirs(path)
    except OSError as e:
        print("{}: ERROR: failed to make directory {}".format(PROGRAM_NAME, dir))
        print(e)
        exit(ERROR_MAKE_DIR)


def main():
    args = parse_args()
    with open(args.app) as application_filename:
        application = json.load(application_filename) 
        validate_json(application_filename, application)
        create_app_dir(application)

if __name__ == '__main__':
    main()

