#!/usr/bin/env python

'''
Replace text in the header of a VCF file in certain fields for the purposes
of anonymisation.

Replaces --old with --new in the following places:

    - VCF header up until and including the column header line 

We assume that only header lines start with '#'.

Usage:

    vcf_edit.py --old oldtext --new newtext --input example_input.vcf --output example_output.vcf

Authors: Bernie Pope, Gayle Philip

'''

from argparse import ArgumentParser

def parse_args():
    """Replace old text with new text in the header of a VCF file"""
    parser = ArgumentParser(description="Replace old text with new text in the header of a VCF file")
    parser.add_argument("--old", required=True, type=str, help="old string (to be replaced)")
    parser.add_argument("--new", required=True, type=str, help="new string (to replace old)")
    parser.add_argument("--output", required=True, type=str, help="output VCF file path")
    parser.add_argument("--input", required=True, type=str, help="input VCF file path")
    return parser.parse_args() 

def vcf_edit(old, new, input_filename, output_filename):
    with open(input_filename) as input_file, \
         open(output_filename, "w") as output_file:
        for line in input_file:
            if line.startswith('#'):
                # this replaces all occurrences of old with new
                # on the input line
                # We assume header lines start with at least one '#'
                # character.
                new_line = line.replace(old, new)
            else:
                new_line = line
            output_file.write(new_line)

def main():
    args = parse_args()
    vcf_edit(args.old, args.new, args.input, args.output)

if __name__ == '__main__':
    main()
