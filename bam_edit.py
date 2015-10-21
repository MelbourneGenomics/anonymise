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
    parser.add_argument("--input", type=str, help="input BAM file path")
    return parser.parse_args() 

def main():
    args = parse_args()
    with pysam.AlignmentFile(args.input, "r") as bam_input:
        input_header = bam_input.header
        # replace old with new in the ID field of RG in the header
        if 'RG' in input_header:
            for read_group in input_header['RG']:
                if 'ID' in read_group:
                    new_id = read_group['ID'].replace(args.old, args.new)
                    read_group['ID'] = new_id
        with pysam.AlignmentFile(args.output, "wb", header=input_header) as bam_output:
            # replace old with new in the query name for each read
            for read in bam_input:
                read.query_name = read.query_name.replace(args.old, args.new)
                if read.has_tag('RG'):
                    new_tag = read.get_tag('RG').replace(args.old, args.new)
                    read.set_tag('RG', new_tag)
                bam_output.write(read)
    bam_input.close()
    bam_output.close()

if __name__ == '__main__':
    main()
