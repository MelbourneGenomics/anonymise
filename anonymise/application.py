'''
Manipulation of data request application

Usage:

Authors: Bernie Pope, Gayle Philip
'''

from __future__ import print_function
import json
import os
import sys
from jsonschema import validate, ValidationError
from pkg_resources import resource_filename
from collections import namedtuple
from program_name import PROGRAM_NAME
from error import print_error, ERROR_JSON_SCHEMA_DEFINE, \
    ERROR_INCOMPATIBLE_REQUEST, ERROR_INVALID_APPLICATION, \
    ERROR_JSON_SCHEMA_INVALID, ERROR_JSON_SCHEMA_OPEN

JSON_SCHEMA = "application_json_schema.txt"
FILE_TYPES = ["fastq", "bam", "vcf"]
COHORTS = ["AML", "EPIL", "CS", "CRC", "CMT"]

class Application(object):

    def __init__ (self, filename):
        fields = json.load(filename)
 
        # Validate input JSON file against schema. This function exits the program if
        # the validation fails, otherwise we return a dictionary representing
        # the JSON file
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
            validate(fields, json_schema)
        except ValidationError as e:
            print_error("JSON input file is not valid: {}".format(filename))
            print(e, file=sys.stderr)
            exit(ERROR_JSON_SCHEMA_INVALID)

        # We only reach here if everything succeeded
        self.fields = fields 

    def allowed_data_types(self):
        '''Based on the input application, decide what data is available
        to the requestor.'''
        request = application_to_request(self.fields)
        try:
            allowed_combinations = REQUEST_COMBINATIONS[request]
        except KeyError:
            print_error("Application does not have a valid interpretation")
            print(format(json.dumps(self.fields, indent=4)), file=sys.stderr)
            exit(ERROR_INVALID_APPLICATION)
        # check that what was requested is compatible with what is available
        requested_identifiability = self.fields['identifiability']
        if requested_identifiability in allowed_combinations:
            return allowed_combinations 
        else:
            print_error("Requested identifiability {} is not compatible with allowed results {}".format(requested_identifiability, allowed_combinations))
            exit(ERROR_INCOMPATIBLE_REQUEST)

    def cohorts(self):
        '''Return a list of all the cohorts requested in an application'''
        return [cohort for cohort in COHORTS if self.fields["condition"][cohort] == "TRUE"]

    def file_types(self):
        '''Return a list of all the file types requested in an application'''
        return [file_type for file_type in FILE_TYPES if self.fields["file types"][file_type] == "TRUE"]


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
