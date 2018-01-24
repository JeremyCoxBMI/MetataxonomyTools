__author__ = 'COX1KB'



#############
#
#   Takes a file with columns featuring taxon ID and possibly counts associated with that ID.
#
#   filename: name of file to import
#   column_taxon:  0-indexed column number of the taxon ID
#   delimiter:     defaults to splitting on whitespace, however sometimes spaces are in values (NOT regex)
#                  (input to the string.split() function   https://docs.python.org/2/library/string.html
#   column_counts: 0-indexed column number of the associated counts
#############

def readFileToDictionary( filename, column_taxon, delimiter=None, column_counts=None):
    result = {}

    if column_counts == None:
        counts = False
    else:
        counts = True

    if delimiter==None and not counts:
        for line in open(filename):
            splits=line.split()
            result[ int(splits[column_taxon]) ]=1
    elif delimiter==None and counts:
        for line in open(filename):
            splits=line.split()
            result[ int(splits[column_taxon]) ]= int(splits[column_counts])
    #else delimiter exists
    elif not counts:
        for line in open(filename, delimiter):
            splits=line.split()
            result[ int(splits[column_taxon]) ]=1
    elif counts:
        for line in open(filename, delimiter):
            splits=line.split()
            result[ int(splits[column_taxon]) ]= int(splits[column_counts])

