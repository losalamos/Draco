#======================================================================
# Module: milagro.input
#
# Contains functions for working with Milagro input files,
#======================================================================

import re, string, array, os

#----------------------------------------------------------------------

def extract_field(file_name, field):

    """Find the first appearance of the field in a milagro imput file

    Field is a string field indicator found in a Milagro input file
    (no colon required). This function looks for the first appearence
    of the field in the input file and returns the pattern match
    object"""

    # Build a regex object
    regexField     = re.compile(field + ':[\s]+(.+)')
    regexFirstChar = re.compile('^[ ]*(\w)')
    returnString = None

    # Read the lines of the file and look for a match
    file_in = open(file_name, 'r')
    contents = file_in.readlines()
    for line in contents:

        fieldMatch = regexField.search(line)

        if fieldMatch:
            firstCharMatch = regexFirstChar.search(line)
            firstChar = firstCharMatch.start(1)

            # Disallow lines that begin with 'c', unless it's the first character in the tag:
            if (line[firstChar]!="c" or fieldMatch.start()==firstChar): 
                returnString = fieldMatch.group(1)
                break

    file_in.close()
    return returnString

