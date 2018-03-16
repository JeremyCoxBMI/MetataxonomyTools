# -*- coding: utf-8 -*-
__author__ = 'COX1KB'

import sys

#data structure
# maps speciesID to, DICT firstTaxon to length)

speciesToFirstTaxonLength = {}

first = True
for line in open(sys.argv[1]):
    splits = line.split('\t')
    if first:
      first=False
    elif len(splits) >= 7:
      acc = splits[0]
      firstTaxon = splits[2]
      speciesTaxon = splits[4]
      length = int(splits[8])

      if not speciesTaxon in speciesToFirstTaxonLength:
	  speciesToFirstTaxonLength[speciesTaxon] = {}
      specDict = speciesToFirstTaxonLength[speciesTaxon]
      if firstTaxon in specDict:
	  specDict[firstTaxon] += length
      else:
	  specDict[firstTaxon] = length



outF = open(sys.argv[1]+".species.to.length", "w")
outF.write("SPECIES_IT\tFIRST_TAXON\tLENGTH\n")

for specID in speciesToFirstTaxonLength:
    firstDict = speciesToFirstTaxonLength[specID]
    for f in firstDict:
        outF.write(specID+"\t"+f+"\t"+str(firstDict[f])+"\n")