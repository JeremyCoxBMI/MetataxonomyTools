# -*- coding: utf-8 -*-
__author__ = 'COX1KB'

import taxonomyLib as tL
import sys
import os
import cPickle as pickle


# takes accession numbers to length
# maps accession to length,

def findLargestKey(adict):
    temp = []
    for key in adict:
        temp.append(  ( adict[key], key )   )
    temp.sort()
    
    return temp[-1]     # (length, taxon)

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
    files_list_file = sys.argv[1:]
    # print str(len(files))
    # print files

    clades = ["firstTaxon"]+tL.LINNAEUS_TAXONOMY_REVERSE
    taxaList ={}
    for clade in tL.LINNAEUS_TAXONOMY_REVERSE:
        taxaList[clade] = {}


    # indexToLevels = [-1 for x in range(len(tL.height.keys()) + 1)]
    # indexToLevels[0] = "firstTaxon"
    # for level in tL.height:
    #     indexToLevels[tL.height[level]] = level
    #
    # for key in range(1, 9):
    #     clades.append(indexToLevels[key])

    print >> sys.stderr, "10 Building databases"
    (nodes, levels) = tL.buildNodes()
    names = tL.buildNames()

    print >> sys.stderr, "BEGIN LOOP"
    print >> sys.stderr, "20 Build list accession numbers"
    for f in files_list_file:

        pickleName = f + ".accession2taxa.pickle"
        if os.path.exists(pickleName):
            print >> sys.stderr, "30 Loading TAXA object from pickle"
            taxa = pickle.load(open(pickleName, 'rb'))
        else:

            accessionNumbersToFind = {}

            print >> sys.stderr, "\t\tFILE:\t" + f

            for line in open(f.strip()):
                spl = line.split()
                acc = spl[0]
                l = int(spl[1])
                acc = acc.split('.')[0]
                accessionNumbersToFind[acc] = l

            print >> sys.stderr, "Number of Accession to Find", str(len(accessionNumbersToFind))

            print >> sys.stderr, "30 Find Taxa from Accession Number"

            taxa = findAccessionNumbers2TaxaLength(accessionNumbersToFind)

            pickle.dump(taxa, open(pickleName, "wb"))

        print >> sys.stderr, "TAXA found", str(len(taxa.keys()))

        accessionToTaxonomy = {}

        print >> sys.stderr, "40 Find species and Kingdom, write file real time"

        header = "ACCESSION\t" + "\t".join(clades) + "\tSEQ_LENGTH\n"
        outF2 = open(f + ".taxonomy_as_text.txt", "w")
        outF2.write(header)

        k = 0
        kingdom = "Bacteria"

        for acc in taxa:
            (taxon, L) = taxa[acc]
            taxon = int(taxon)
            taxonomy_txt = tL.buildTaxaLevels2(int(taxon), names, nodes, tL.LINNAEUS_TAXONOMY_REVERSE)
            taxa = taxonomy_txt.split()
            for clade in tL.LINNAEUS_TAXONOMY_REVERSE:
                k=tL.LINNAEUS_TAXONOMY_REVERSE.index(clade)
                if taxa[k] in taxaList[clade]:
                    taxaList[clade][taxa[k]][ taxon ] += int(L)
                else:
                    taxaList[clade][taxa[k]] = {}
                    taxaList[clade][taxa[k]][ taxon ] = int(L)
                
            
            outLine2 = acc + "\t" + str(taxon)+ "\t" + taxonomy_txt + "\t" + str(L) + "\n"

            outF2.write(outLine2)

            if k % (100 * 1000) == 0:
                print >> sys.stderr, "Processing " + str(k / (1000 * 1000.0)) + " millionth Accession"
            k += 1
    #end for l in file_list
        outF2.close()

    print >> sys.stderr, "Outputting one taxon per clade files"

    for clade in tL.LINNAEUS_TAXONOMY_REVERSE:

        filename=f + "."+str(clade)+'.longest.only.txt'
        outF = open(filename, "w")
        total=0
        k=tL.LINNAEUS_TAXONOMY_REVERSE.index(clade)

        print >> sys.stderr, "\tOutputting one taxon per clade files"
        print >> sys.stderr, "\t\t",str(clade),"\t"+filename

        for taxon in taxaList[clade]:
            (L, taxon) = findLargestKey(taxaList[clade][taxon])
            outF.write(str(taxon)+"\t"+str(L)+"\n")

        outF.close()
        print >> sys.stderr, "\t\t",str(clade),"\tcount\t"+str(total)