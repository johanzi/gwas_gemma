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
    should either be named as digits or digits with prefix "Chr"
    as e.g Chr3. All other nomenclatures will be removed from the
    output file.
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
            SNP = line[1].replace(":","_")  # Second column but replace column by underscore
            CHR = line[1].split(":")[0] # Get the chromosome ID
            # Note that the CHR field should be numeric for plotting R 
            # manhattan package. Therefore, prefix "Chr" should be removed and organelles
            # with string names should be skipped
            
            # Nothing to do if the chromosome is already a digit
            if is_digit(CHR):
                CHR = CHR            
            # Remove "Chr or chr" prefix if present and check if suffix is a digit
            elif CHR[0:3].lower == "chr" and is_digit(CHR[3::]):
                CHR = CHR.lower().replace("chr","")
            # Skip the line if CHR is something else than e.g. 
            # a digit or "Chr" + digit
            else:
               continue
        
            # Get the information for the other fields
            BP = line[2]
            P = line[8]
            zscore = line[7]
            new_line = SNP, CHR, BP, P, zscore
            print("\t".join(new_line))
    input_file.close()

if __name__ == "__main__":
    sys.exit(main())

