__author__ = 'COX1KB'

import taxonomyLib as tl
import sys



# takes accession numbers to length
# maps accession to length,

def findAccessionNumbers2TaxaLength( adict ):
    result = {}
    #iterate over file
    for line in open(tl.AccessionDBpath+"master.accession2taxid"):
        splits = line.split()
        acc = splits[0].split('.')[0]
        if acc == "NZ_GG666849":
            DEBUG = "LINE"

        # if len(acc) > 0 and acc[0] == 'N':
        #     cheese = acc in adict
        #     DEBUG = "LINE"

        if acc in adict:
            #print splits[0], " found ", splits[2]
            result[acc]=(splits[2], adict[acc])  #accession to (taxon, length)
    return result


if __name__ == "__main__":

    print >> sys.stderr, "01 Monitoring Start Up"
    files_list_file=sys.argv[1:]
    # print str(len(files))
    # print files




    print >> sys.stderr, "10 Building databases"
    (nodes, levels) = tl.buildNodes()
    names = tl.buildNames()

    print >> sys.stderr, "BEGIN LOOP"
    print >> sys.stderr, "20 Build list accession numbers"
    for f in files_list_file:

        accessionNumbersToFind = {}

        print >> sys.stderr, "\t\tFILE:\t"+f

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

        outF = open(f+"species.db.list.txt","w")
        outF.write("ACCESSION\tFIRST_TAXON\tSPECIES_ID\tKINGDOM\tLENGTH\n")

        k=0
        for acc in taxa:
            (taxon, l) = taxa[acc]
            speciesID = tl.getSpeciesID(taxon, nodes)

            kingdom = tl.findKingdom(speciesID, names, nodes)
            outF.write(str(acc)+"\t"+taxon+"\t"+str(speciesID)+"\t"+str(kingdom)+"\t"+str(l)+"\n")
            if k % (100*1000) == 0:
                print >> sys.stderr, "Processing "+str(k/(1000*1000.0))+" millionth Accession"
            k+=1
    #end for l in file_list