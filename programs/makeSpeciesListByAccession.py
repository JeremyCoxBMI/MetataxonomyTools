__author__ = 'COX1KB'

import taxonomyLib as tl
import sys

if __name__ == "__main__":

    print >> sys.stderr, "01 Monitoring Start Up"
    files_list_file=sys.argv[1]
    # print str(len(files))
    # print files

    print >> sys.stderr, "10 Building databases"
    (nodes, levels) = tl.buildNodes()
    names = tl.buildNames()


    accessionToFind = {}

    print >> sys.stderr, "20 Build list accession numbers"
    for l in open(files_list_file):
        for line in open(l.strip()):
            if line[0]==">":
                accessionToFind[ line.split()[0][1:] ] = 1

    print >> sys.stderr, "30 Find Taxa from Accession Number"
    taxa = tl.findTaxaAccessionNumbers(accessionToFind)

    species = {}

    print >> sys.stderr, "40 Find species and Kingdom"
    for taxon in taxa:
        speciesID = tl.getSpeciesID(taxon, nodes)
        if not speciesID in species:
            species[speciesID] = tl.findKingdom(speciesID, names, nodes)

    print >> sys.stderr, "50 write file"

    outF = open("species.db.list.txt","w")
    for sID in species:
        outF.write(str(sID)+"\t"+species[sID]+"\n")
    outF.close()

