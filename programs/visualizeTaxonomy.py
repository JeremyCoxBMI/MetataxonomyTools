# -*- coding: utf-8 -*-
__author__ = 'COX1KB'

import sys
import taxonomyLib as tl

def makeTabs( k):
    result=""
    for x in range(k):
        result+='\t'
    return result

def recurseOutput( aString, atree, numTabs):
    if isinstance(atree, dict):
        first=True
        for key in atree:
            if first:
                first=False
                string1 = aString+'\t'+str(key)
            else:
                string1 = makeTabs(numTabs)+str(key)
            recurseOutput(string1, atree[key],numTabs+1)

    elif isinstance(atree, list):
        #species level, next level is firstTaxon level
        #sort from greatest to least for output
        atree.sort()
        atree.reverse()

        first=True
        for value in atree:
            (l3ngth, name) = value
            if first:
                first=False
                string1 = aString+'\t'+str(name)+"\t"+str(l3ngth)
                splits = string1.split("\t")
                for k in range(len(splits)):
                    if len(splits[k]) > 0:
                        cladeBaseCounts[k] += l3ngth

            else:
                string1 = makeTabs(numTabs)+str(name)+"\t"+str(l3ngth)
            print string1
    # else:
    #     ##at leaves
    #     print aString+'\t'+str(atree)


def initializeDoubleDict(aList):
    result = {}
    for e in aList:
        result[e] = {}
    return result



taxonomyTree = {}


first=True
clades_to_proces = range(1,6)
clades_to_proces.reverse()

cladeBaseCounts = [0 for x in range(len(clades_to_proces)+4)]

cladeIndex = [tl.LINNAEUS_TAXONOMY[x] for x in clades_to_proces]

for line in open(sys.argv[1]):
    splits = line.split('\t')
    if first:
        colnames=splits
        clades=colnames[1:-1]
        first=False
    else:
        curr_dict = taxonomyTree
        l3ngth = int(splits[-1])
        if len(splits) >= 9:
            for x in clades_to_proces:
                v=splits[x]
                if x == clades_to_proces[-1]: #firstTaxon level; is a list
                    found=False
                    for entry in curr_dict:
                        if v == entry[1]:
                            entry = (entry[0]+l3ngth,entry[1])
                            found = True
                            break
                    if not found:
                        curr_dict.append( (l3ngth, v) )
                elif v in curr_dict:
                    curr_dict = curr_dict[v]
                else:
                    if x == clades_to_proces[-2]: #species level
                        curr_dict[v] = []
                    else:
                        curr_dict[v] = {}
                    curr_dict = curr_dict[v]


for key in taxonomyTree:
    output = str(key)
    recurseOutput(output, taxonomyTree[key],1)


for x in range(len(clades_to_proces)):
    clade = clades[ clades_to_proces[x] ]
    print >> sys.stderr, "Clade by longest only\t"+clade+"\ttotal length\t"+str(cladeBaseCounts[x])# -*- coding: utf-8 -*-
__author__ = 'COX1KB'

import sys

def makeTabs( k):
    result=""
    for x in range(k):
        result+='\t'
    return result

def recurseOutput( aString, atree, numTabs):
    if isinstance(atree, dict):
        first=True
        for key in atree:
            if first:
                first=False
                string1 = aString+'\t'+str(key)
            else:
                string1 = makeTabs(numTabs)+str(key)
            recurseOutput(string1, atree[key],numTabs+1)

    elif isinstance(atree, list):
        #species level, next level is firstTaxon level
        #sort from greatest to least for output
        atree.sort()
        atree.reverse()

        first=True
        for value in atree:
            (l3ngth, name) = value
            if first:
                first=False
                string1 = aString+'\t'+str(name)+"\t"+str(l3ngth)
                splits = string1.split("\t")

                for k in range(len(cladeBaseCounts)):
                    if len(splits[k+1]) > 0:
                        cladeBaseCounts[k] += l3ngth

            else:
                string1 = makeTabs(numTabs)+str(name)+"\t"+str(l3ngth)
            print string1
    # else:
    #     ##at leaves
    #     print aString+'\t'+str(atree)


def initializeDoubleDict(aList):
    result = {}
    for e in aList:
        result[e] = {}
    return result



taxonomyTree = {}


first=True
clades_to_proces = range(1,6)
clades_to_proces.reverse()

cladeBaseCounts = [0 for x in clades_to_proces]


for line in open(sys.argv[1]):
    splits = line.split('\t')
    if first:
        colnames=splits
        clades=colnames[1:-1]
        first=False
    else:
        curr_dict = taxonomyTree
        l3ngth = int(splits[-1])
        if len(splits) >= 9:
            for x in clades_to_proces:
                v=splits[x]
                if x == clades_to_proces[-1]: #firstTaxon level; is a list
                    found=False
                    for entry in curr_dict:
                        if v == entry[1]:
                            entry = (entry[0]+l3ngth,entry[1])
                            found = True
                            break
                    if not found:
                        curr_dict.append( (l3ngth, v) )
                elif v in curr_dict:
                    curr_dict = curr_dict[v]
                else:
                    if x == clades_to_proces[-2]: #species level
                        curr_dict[v] = []
                    else:
                        curr_dict[v] = {}
                    curr_dict = curr_dict[v]


for key in taxonomyTree:
    output = str(key)
    recurseOutput(output, taxonomyTree[key],1)


for x in range(len(clades_to_proces)):
    clade = clades[ clades_to_proces[x] ]
    print >> sys.stderr, "Clade by longest only\t"+clade+"\ttotal length\t"+str(cladeBaseCounts[x])