#!/usr/bin/env python

'''
Replace text in a BAM file in certain fields for the purposes
of anonymisation.

Replaces --old with --new in the following places:

    - BAM header read group (RG) field
    - Each alignment record in the body:
        - the query name
        - the RG field in tags

Usage:

    bam_edit.py --old oldtext --new newtext --input example_input.bam --output example_output.bam 

Authors: Bernie Pope, Gayle Philip

'''

import pysam

from argparse import ArgumentParser

def parse_args():
    """Replace old text in a BAM file"""
    parser = ArgumentParser(description="Edit reads in a BAM file")
    parser.add_argument("--old", required=True, type=str, help="old string (to be replaced)")
    parser.add_argument("--new", required=True, type=str, help="new string (to replace old)")
    parser.add_argument("--output", required=True, type=str, help="output BAM file path")
    parser.add_argument("--input", required=True, type=str, help="input BAM file path")
    return parser.parse_args() 

def bam_edit(old, new, input_filename, output_filename):
    print((old, new, input_filename, output_filename))
    with pysam.AlignmentFile(input_filename, "r") as bam_input:
        input_header = bam_input.header
        # replace old with new in the ID field of RG in the header
        if 'RG' in input_header:
            for read_group in input_header['RG']:
                if 'ID' in read_group:
                    new_id = read_group['ID'].replace(old, new)
                    read_group['ID'] = new_id
                if 'SM' in read_group:
                    new_sm = read_group['SM'].replace(old, new)
                    read_group['SM'] = new_sm
        with pysam.AlignmentFile(output_file, "wb", header=input_header) as bam_output:
            # replace old with new in the query name for each read
            for read in bam_input:
                read.query_name = read.query_name.replace(old, new)
                if read.has_tag('RG'):
                    new_tag = read.get_tag('RG').replace(old, new)
                    read.set_tag('RG', new_tag)
                bam_output.write(read)

def main():
    args = parse_args()
    bam_edit(args.old, args.new, args.input, args.output)


if __name__ == '__main__':
    main()
