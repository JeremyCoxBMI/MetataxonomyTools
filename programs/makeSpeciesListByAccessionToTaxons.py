__author__ = 'COX1KB'

import taxonomyLib as tl
import sys

if __name__ == "__main__":

    print >> sys.stderr, "01 Monitoring Start Up"
    # files_list_file=sys.argv[1]
    # print str(len(files))
    # print files

    print >> sys.stderr, "10 Building databases"
    (nodes, levels) = tl.buildNodes()
    names = tl.buildNames()


    accessionNumbersToFind = {}

    print >> sys.stderr, "20 Build list accession numbers"
    # for l in open(files_list_file):
    for line in open(sys.argv[1]):
        accessionNumbersToFind[ line.split()[0].split('.')[0] ] = 1

    print >> sys.stderr, "Number of Accession to Find", str(len(accessionNumbersToFind))

    print >> sys.stderr, "30 Find Taxa from Accession Number"
    accession2taxa = tl.findTaxaAccessionNumbersMap(accessionNumbersToFind)

    print >> sys.stderr, "TAXA found", str(len(accession2taxa.keys()))

    accession2firstTaxonSpeciesTaxon = {}

    print >> sys.stderr, "40 Find species and Kingdom, write file real time"


    outF = open(sys.argv[1]+"accession2firstTaxonSpeciesTaxon.db.list.txt","w")
    outF.write("Accession\tfirstTaxon\tspeciesTaxon\tGenusTaxon\tKingdom")
    for accession in accession2taxa:
        taxon = int(accession2taxa[accession])
        (speciesID,genusID) = tl.getSpeciesIDGenusID(taxon, nodes)
        if not accession in accession2firstTaxonSpeciesTaxon:
            kingID = tl.findKingdom(speciesID, names, nodes)
            accession2firstTaxonSpeciesTaxon[accession] = (taxon, speciesID,genusID,kingID)
            outF.write(str(accession)+"\t"+str(taxon)+"\t"+str(speciesID)+"\t"+str(genusID)+"\n"+str(kingID))

    # print >> sys.stderr, "50 write file"
    #
    # outF = open("species.db.list.txt","w")
    # for sID in species:
    #     outF.write(str(sID)+"\t"+species[sID]+"\n")
    # outF.close()

