# -*- coding: utf-8 -*-
__author__ = 'COX1KB'

import taxonomyLib as tL
import sys
import os
import cPickle as pickle


# takes a list of (integer, text) and returns pair wwith largest integer

def findLargestKey(adict):
    temp = []
    for key in adict:
        temp.append(( adict[key], key ))
    temp.sort()

    return temp[-1]  # (length, taxon)


# takes accession numbers to length
# maps accession to length,

def findAccessionNumbers2TaxaLength(adict):
    result = {}
    # iterate over file
    for line in open(tL.AccessionDBpath + "master.accession2taxid"):
        splits = line.split()
        acc = splits[0].split('.')[0]
        if acc == "NZ_GG666849":
            DEBUG = "LINE"

        # if len(acc) > 0 and acc[0] == 'N':
        #     cheese = acc in adict
        #     DEBUG = "LINE"

        if acc in adict:
            #print splits[0], " found ", splits[2]
            result[acc] = (splits[2], adict[acc])  #accession to (taxon, length)
    return result


if __name__ == "__main__":

    print >> sys.stderr, "01 Monitoring Start Up"
    clades = ["firstTaxon"] + tL.LINNAEUS_TAXONOMY_REVERSE
    taxaList = {}
    for clade in tL.LINNAEUS_TAXONOMY_REVERSE:
        taxaList[clade] = {}


    files_list_file = sys.argv[1:]
    # print str(len(files))
    # print files




    print >> sys.stderr, "10 Building databases"
    (nodes, levels) = tL.buildNodes()
    names = tL.buildNames()



    print >> sys.stderr, "BEGIN LOOP"
    print >> sys.stderr, "20 Build list accession numbers"
    for f in files_list_file:

        accessionNumbersToFind = {}

        print >> sys.stderr, "\t\tFILE:\t" + f

        for line in open(f.strip()):
            (acc, l) = line.split()
            acc = acc.split('.')[0]
            accessionNumbersToFind[acc] = l

        print >> sys.stderr, "Number of Accession to Find", str(len(accessionNumbersToFind))

        print >> sys.stderr, "30 Find Taxa from Accession Number"

        taxa = findAccessionNumbers2TaxaLength(accessionNumbersToFind)

        print >> sys.stderr, "TAXA found", str(len(taxa.keys()))

        species = {}

        print >> sys.stderr, "40 Find species and Kingdom, write file real time"

        # outF = open(f + ".species.db.list.txt", "w")
        # outFerr = open(f + ".species.db.list.txt.Accession.errors", "w")
        # outF.write("ACCESSION\t|\tFIRST_TAXON\t|\tSPECIES_ID\t|\tKINGDOM\t|\tLENGTH\n")
        # outFerr.write("ACCESSION\t|\tFIRST_TAXON\t|\tSPECIES_ID\t|\tKINGDOM\t|\tLENGTH\n")

        k = 0
        kingdom = "Bacteria"


        header = "ACCESSION\t" + "\t".join(clades) + "\tSEQ_LENGTH\n"
        outF2 = open(f + ".taxonomy_as_text.txt", "w")
        outF2.write(header)


        # print taxa
        for acc in taxa:
            temp = taxa[acc]
            (taxon, L) = temp
            taxon = int(taxon)
            taxonomy_txt = tL.buildTaxaLevels2(int(taxon), names, nodes, tL.LINNAEUS_TAXONOMY_REVERSE)
            taxa = taxonomy_txt.split()
            for clade in tL.LINNAEUS_TAXONOMY_REVERSE:
                k = tL.LINNAEUS_TAXONOMY_REVERSE.index(clade)
                if taxa[k] in taxaList[clade]:
                    taxaList[clade][taxa[k]][taxon] += int(L)
                else:
                    taxaList[clade][taxa[k]] = {}
                    taxaList[clade][taxa[k]][taxon] = int(L)

            outLine2 = acc + "\t" + str(taxon) + "\t" + taxonomy_txt + "\t" + str(L) + "\n"

            outF2.write(outLine2)

            if k % (100 * 1000) == 0:
                print >> sys.stderr, "Processing " + str(k / (1000 * 1000.0)) + " millionth Accession"
            k += 1
            #end for l in file_list
        outF2.close()

        print >> sys.stderr, "Outputting one taxon per clade files"

        for clade in tL.LINNAEUS_TAXONOMY_REVERSE:

            filename = f + "." + str(clade) + '.longest.only.txt'
            outF = open(filename, "w")
            total = 0
            k = tL.LINNAEUS_TAXONOMY_REVERSE.index(clade)

            print >> sys.stderr, "\tOutputting one taxon per clade files"
            print >> sys.stderr, "\t\t", str(clade), "\t" + filename

            for taxon in taxaList[clade]:
                (L, taxon) = findLargestKey(taxaList[clade][taxon])
                outF.write(str(taxon) + "\t" + str(L) + "\n")

            outF.close()
            print >> sys.stderr, "\t\t", str(clade), "\tcount\t" + str(total)