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
    R analysis with the package 'qqman'. Note that chromosomes
    should named as digits.
    Compatible Python27 and Python37
'''


# Function to test whether the string following "Chr" 
# in the chromosome name is a digit or not. This allow
# to discard lines with ChrM or ChrC for instance
def is_digit(x):
    try:
        int(x)
        return True
    except ValueError:
        return False

def main():

    # Open the file
    input_file = sys.argv[1]
    
    # Check if input file exists
    if os.path.isfile(input_file):
        input_file = open(input_file, "r")
        lines_file = input_file.read().splitlines()
    else:
        sys.exit("Error: input file does not exist")

    # Get through the file and parse the fiels to reorganize them as expected
    # by manahattan R package
    for line in lines_file:
        if line[:3] == "chr":  #First line should be header (chr per default in gemma output)
            header = "SNP\tCHR\tBP\tP\tzscore"
            print(header)
        else:
            line = line.strip().split("\t")  # Get rid of EOL and Create a list based on \t separation
            
            # Get chromosome name (should be a digit) and the name of the SNP
            # on the chromosome
            CHR = line[0]
            
            SNP = line[1]
            
            # If CHR=0, it means the chromosome names in the VCF file are not pure digit
            # then, retrieve chromosome name from the SNP field and remove all non digits
            if int(CHR) == 0:
                CHR = re.split(':', SNP)[0] 
                CHR = re.sub("[^0-9]", "", CHR)
            

            # Skip the SNP if the CHR is not a digit)
            if not is_digit(CHR):
                continue

            # Get the information for the other fields 
            # (note that it works when GEMMA option -lmm 2 is set)
            BP = line[2]
            P = line[8]
            zscore = line[7]
            new_line = SNP, CHR, BP, P, zscore
            print("\t".join(new_line))
    input_file.close()

if __name__ == "__main__":
    sys.exit(main())

