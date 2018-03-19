# -*- coding: utf-8 -*-
__author__ = 'COX1KB'

import sys


def initializeDoubleDict(aList):
    result = {}
    for e in aList:
        result[e] = {}
    return result


colnames = []
clades =[]

firstTaxonToLength = {}

# maps clade name to
#   dictionary of taxonIDS at that clade to LIST of firstTaxonID
firstTaxonListAtCladeAndTaxonID = {}



first=True
for line in open(sys.argv[1]):
    splits = line.split('\t')
    if first:
        colnames=splits

#typo previous program
	#clades=splits[1:-1] #skip first (accession) and last (seq length)
        clades=splits[2:-1] #skip first (accession) and last (seq length)

	#print splits
	#print '\n'

        firstTaxonListAtCladeAndTaxonID = initializeDoubleDict(clades)
        first=False
    else:
        ft = splits[1]
        if not ft in firstTaxonToLength:
            firstTaxonToLength[ft] = int(splits[-1])
        else:
            firstTaxonToLength[ft] += int(splits[-1])

        for x in range(len(clades)):
#            print x, "\t", len(splits), "\t", len(clades), "\t", line

	    taxon = splits[x+1]
            clade = clades[x]
            if taxon in firstTaxonListAtCladeAndTaxonID[clade]:
                firstTaxonListAtCladeAndTaxonID[clade][taxon].append(ft)
            else:
                firstTaxonListAtCladeAndTaxonID[clade][taxon] = [ ft ]

outFileByClade = []
total_size=0
for clade in clades:
    outF = open(sys.argv[1]+"."+clade,'w')

    outF.write("taxon_within_clade_"+clade+"\tfirstTaxonLongestGenome\tgenome_length\n")

    for taxonID in firstTaxonListAtCladeAndTaxonID[clade]:
        list_find_max = firstTaxonListAtCladeAndTaxonID[clade][taxonID]
        for k in range(len(list_find_max)):
            ft = list_find_max[k]
            list_find_max[k] = (firstTaxonToLength[ft],ft)
        list_find_max.sort()
        (l3ntgh, ft) = list_find_max[-1] #maximum
        total_size += l3ntgh
        ###output taxonID, ft, length
        outLine = taxonID+"\t"+ft+"\t"+str(l3ntgh)+"\n"
        outF.write(outLine)
    print >> sys.stderr, "Longest genomes only at clade\t"+clade+"\thas total length\t"+str(total_size)
    total_size = 0

