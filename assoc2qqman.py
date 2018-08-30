#!/usr/bin/python

import sys
import os
import getopt
import re

'''
        AUTHOR: Johan Zicola

	USAGE: assoc2qqman.py <prefix.assoc.txt>
    Read the output of gemma 'prefix.assoc.txt' derived 
    from the association analysis command and performs 
    readjustment of the file to be compatible for Read
    R analysis with the package 'qqman'.
'''


def main():

    # Open the file
    input_file = sys.argv[1]
    
    # Check if input file exists
    if os.path.isfile(input_file):
        input_file = open(input_file, "r")
        lines_file = input_file.read().splitlines()
    else:
        sys.exit("Error: input file does not exist")


    for line in lines_file:
        if line[:3] == "chr":  #First line should be header
            header = "SNP\tCHR\tBP\tP\tzscore"
            print header
        else:
            line = line.strip().split("\t")  # Get rid of EOL and Create a list based on \t separation
            SNP = line[2]  # third column
            chr_name = line[1].split(":")
            #Keep only lines matching chromosomes (exclude Mt, Pt, ...)
            if chr_name[0][0:3] == "Chr":
                CHR = chr_name[0].replace("Chr", "")
                BP = line[2]
                P = line[8]
                zscore = line[7]
                new_line = SNP, CHR, BP, P, zscore
                print "\t".join(new_line)
    input_file.close()

if __name__ == "__main__":
    sys.exit(main())

