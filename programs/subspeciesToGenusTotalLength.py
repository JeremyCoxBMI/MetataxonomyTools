# -*- coding: utf-8 -*-
__author__ = 'COX1KB'

import sys

#data structure
# maps speciesID to, DICT firstTaxon to length)

#program modified from species ==> genus  VARIABLE NAMES NOT CHANGED

speciesToFirstTaxonLength = {}

first = True
for line in open(sys.argv[1]):
    splits = line.split('\t')
    if first:
      first=False
    elif len(splits) >= 7:
      acc = splits[0]
      firstTaxon = splits[1]
      speciesTaxon = splits[2]
      genusTaxon = splits[3]
      length = int(splits[4])
    
      if not speciesTaxon in speciesToFirstTaxonLength:
        speciesToFirstTaxonLength[speciesTaxon] = {}
        specDict = speciesToFirstTaxonLength[genusTaxon]
      if firstTaxon in specDict:
        specDict[firstTaxon] += length
      else:
        specDict[firstTaxon] = length



outF = open(sys.argv[1]+".genus.to.length", "w")
outF.write("GENUS_ID\tFIRST_TAXON\tLENGTH\n")
outF2 = open(sys.argv[1]+".one_species.to.length", "w")
outF2.write("GENUS_ID\tFIRST_TAXON\tLENGTH\n")

all_total=0
total=0
for specID in speciesToFirstTaxonLength:
    firstDict = speciesToFirstTaxonLength[specID]
    maxLength = (0,"fakeName")
    for f in firstDict:
      if maxLength[0] < firstDict[f]:
        maxLength = (firstDict[f],f)
      outF.write(specID+"\t"+f+"\t"+str(firstDict[f])+"\n")
      all_total += firstDict[f]
    outF2.write(specID+"\t"+maxLength[1]+"\t"+str(maxLength[0])+"\n")
    total += maxLength[0]

print "total length of Bacteria Metagenome per species is "+str(all_total)
print "total length using longest genome per species only is "+str(total)