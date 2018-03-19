# -*- coding: utf-8 -*-
from config import *
import sys

# Thanks to NCBI, the structure of the taxonomy table is constantly changing
# 1 is no longer the root node
# 2018-02-02 noticed tree no longer rooted to Taxon=1 (ROOT)
            # root (old school)     #bacteria#cellular organisms    #viruses           #eukaroyta
BREAK_IDS = {1 : True,      #old root
             2 : True,      #bacteria
             131567 : True, #celular organisms
             10239 : True,  #viruses
            2759 : True}    #eukaryota
# There seems to be some broken edges in the nodes table
# Making max iteration to stop any problems
MAX_ITERATIONS=20

def buildNames():
    return buildNamesLocal(DBpath+"names.dmp")

def buildNamesLocal(file):
    inFile = open(file)

    result = dict()

    for line in inFile:
        line = line.replace("\t","")
        line = line.replace("|\n","")
        splits = line.split("|")
        if splits[-1] == 'scientific name':
            key = int( line.split("|")[0] )
            name = line.split("|")[1]
            result[key] = name

    return result

def buildNodes():
    return buildNodesLocal(DBpath+"nodes.dmp")

def buildNodesLocal(file):
    inFile = open(file)

    result = dict()
    levels = dict()

    for line in inFile:
        #print line
        line = line.replace("\t","")
        line = line.replace("|\n","")
        splits = line.split("|")
        #print splits
        taxaID = int(splits[0])
        parentID = int(splits[1])
        level = splits[2]
        result[taxaID] = (parentID, level)
        levels[level]=1

    return (result,levels)

def buildTaxaLevels( x , names, nodes, filter = None ):
    open = ""
    result = ""
    k = x
    result += names[k]
    (next, level) = nodes[k]
    k = next
    #EXAMPLE
    #  (((((Gammaretrovirus)Orthoretrovirinae)Retroviridae)Retro-transcribing_viruses)Viruses,
    z=0
    while True:
        (next, level) = nodes[k]
        if filter == None or level in filter:
            result = result + (")_%s_" % level) + names[k]
            open += "("
        elif z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::buildTaxaLevels could not complete for taxon\t"+str(x)
            break
        if k in BREAK_IDS:
            break
        k = next
        z+=1


    return open+result

height = dict()
height["superkingdom"]=8
height["kingdom"]=7
height["phylum"]=6
height["class"]=5
height["order"]=4
height["family"]=3
height["genus"]=2
height["species"] = 1
#height["firstTaxon"] = 0

#TODO - update for MAX_ITERATIONS
def buildTaxaLevelList2( x , nodes, filter):
    #Returns index key list from highest level to most specific
    #reverse order is needed to build tree
    result = []
    k = x
    last = "species"
    while True:
        if k in nodes:
            (next, level) = nodes[k]
        else:
            (next, level) = 1, "no rank"
        if level in filter:
            distance = height[level]-height[last]
            result.append( (k, distance) )
            #last updates only when a level is used
            last = level
        if k in BREAK_IDS:
            break
        #k always updates to move next
        k = next

    #add prokaryota level
    if len(result) > 0 and result[-1] == 2:
        result.append( (1712345,1) )

    #add archaebacteria
    if len(result) > 0 and result[-1] == 2157:
        result.append( (18181818,1) )

    result.reverse()

    return result

def buildTaxaLevelList3( x , nodes, filter=height):
    #Returns index key list from x to species to kingdom
    # -1 if taxon is not found
    #reverse order is needed to build tree
    result = [-1 for x in range(len(filter)+1)]
    k = int(x)
    result[0] = k
    last = "species"
    z=0
    while True:
        if k in nodes:
            (next, level) = nodes[k]
        else:
            (next, level) = 1, "no rank"
        if level in filter:
            # distance = height[level]-height[last]
            index = filter[level]
            result[index] = k
            #last updates only when a level is used
            last = level
        if k in BREAK_IDS:
            break
        #k always updates to move next
        k = next
        z+=1
        if z >= MAX_ITERATIONS:
            break

    #add prokaryota level
    if len(result) > 0 and result[-1] == 2:
        result.append( (1712345,1) )

    #add archaebacteria
    if len(result) > 0 and result[-1] == 2157:
        result.append( (18181818,1) )

    result.reverse()

    return result



def buildTaxaLevelList( x , nodes):
    #Returns index key list from highest level to most specific
    #reverse order is needed to build tree
    result = []
    k = x

    z=0
    while True:
        result.append(k)
        if k in nodes:
            (next, level) = nodes[k]
        else:
            (next, level) = 1, "no rank"
        if k in BREAK_IDS:
            break
        k = next
        z+=1
        if z >= MAX_ITERATIONS:
            break

    #add prokaryota level
    if result[-1] == 2:
        result.append(1712345)

    #add archaebacteria
    if result[-1] == 2157:
        result.append(18181818)
    result.reverse()

    return result

def levelToText( level ):
    #note "u" used for 'super-" as in the German 'ueber'
    if level=="superkingdom":
        return "uk"
    elif level=="tribe":
        return "t"
    elif level=="subgenus":
        return "sg"
    elif level=="family":
        return "f"
    elif level=="species subgroup":
        return "ssg"
    elif level=="species group":
        return "sg"
    elif level=="phylum":
        return "p"
    elif level=="superclass":
        return "uc"
    elif level=="subphylum":
        return "sp"
    elif level=="subspecies":
        return "ss"
    elif level=="no rank":
        return "NA"
    elif level=="superorder":
        return "uo"
    elif level=="infraorder":
        return "io"
    elif level=="subclass":
        return "sc"
    elif level=="species":
        return "s"
    elif level=="superphylum":
        return "up"
    elif level=="kingdom":
        return "k"
    elif level=="subtribe":
        return "st"
    elif level=="subkingdom":
        return "sk"
    elif level=="forma":
        return "fr"
    elif level=="infraclass":
        return "ic"
    elif level=="varietas":
        return "v"
    elif level=="subfamily":
        return "sf"
    elif level=="class":
        return "c"
    elif level=="superfamily":
        return "uf"
    elif level=="parvorder":
        return "po"
    elif level=="suborder":
        return "so"
    elif level=="genus":
        return "g"
    elif level=="order":
        return "o"
    return '-1'


# return an integer
def getSpeciesID(taxonID, nodes):
    result = -1
    currID = int(taxonID)
    z=0
    while (currID in nodes):
        if nodes[currID][1] == "no rank":
            result = currID

        if nodes[currID][1] == "species":
            result = currID
            break

            #top nodes
            #2018-02-02 noticed tree no longer rooted to Taxon=1 (ROOT)
            # root (old school)     #bacteria#cellular organisms    #viruses           #eukaroyta
        if currID in BREAK_IDS:
            break
        currID = nodes[currID][0]

        z+=1
        if z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::getSpeciesID could not find taxonID\t"+str(taxonID)
            break
    return result

#####
#  Supposed to get Species and Genus taxon ID.  However, if there is no genus, it will run amok
#
####
def getSpeciesIDGenusID(taxonID, nodes):

    specID=-1
    genusID=-1

    z=0
    currID = int(taxonID)
    while (currID in nodes):
        if specID == -1 and nodes[currID][1] == "no rank":  #capture FIRST TAXON, hopefully to be overwritten
             result = currID

        if nodes[currID][1] == "species":
            specID = currID

        if nodes[currID][1] == "genus":
            genusID = currID
            break

            #top nodes
            #2018-02-02 noticed tree no longer rooted to Taxon=1 (ROOT)
            # root (old school)     #bacteria#cellular organisms    #viruses           #eukaroyta
        if currID in BREAK_IDS:
            break
        z+=1
        if z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::findSpeciesIDGenusID coudl not find taxon\t"+str(taxonID)+"\tresults spec\t"+str(specID)+"\tgenus\t"+str(genusID)
            break
        currID = nodes[currID][0]
    return (specID, genusID)



# returns dictionary mapping accession number to speciesID (int)
def extractSpeciesID(accessionNos, nodes):
    result = dict()
    for line in open(AccessionDBpath+"nucl_gb.accession2taxid"):
        splits = line.split()
        # 0 accession       1 accession.version     2 taxid     3 gid
        if splits[0] in accessionNos:
            result[splits[0]] = getSpeciesID( int(splits[2]),nodes)

    return result

def findTaxaAccessionNumbers( adict ):
    result = {}
    #iterate over file
    for line in open(AccessionDBpath+"master.accession2taxid"):
        splits = line.split()

        if splits[0] in adict:
            #print splits[0], " found ", splits[2]
            result[splits[2]]=1
    return result


def findTaxaAccessionNumbersMap( adict ):
    result = {}
    #iterate over file
    for line in open(AccessionDBpath+"master.accession2taxid"):
        splits = line.split()

        if splits[0] in adict:
            #print splits[0], " found ", splits[2]
            result[splits[0]]=splits[2]
    return result


def findKingdom(speciesID, names, nodes):
    result = "" #text name
    next=speciesID

    z=0
    while (True):
        if next in nodes:
            (next, level) = nodes[next]
            if level == "kingdom":
                result=names[next]
            if level == "superkingdom":
                result+="\t"+names[next]
                break
        if next in BREAK_IDS:
            break
        z+= 1
        if z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::findKingdom could not find speciesID\t"+str(speciesID)
            break
    return result
# Thanks to NCBI, the structure of the taxonomy table is constantly changing
# 1 is no longer the root node
# 2018-02-02 noticed tree no longer rooted to Taxon=1 (ROOT)
            # root (old school)     #bacteria#cellular organisms    #viruses           #eukaroyta
BREAK_IDS = {1 : True,      #old root
             2 : True,      #bacteria
             131567 : True, #celular organisms
             10239 : True,  #viruses
            2759 : True}    #eukaryota
# There seems to be some broken edges in the nodes table
# Making max iteration to stop any problems
MAX_ITERATIONS=20

def buildNames():
    return buildNamesLocal(DBpath+"names.dmp")

def buildNamesLocal(file):
    inFile = open(file)

    result = dict()

    for line in inFile:
        line = line.replace("\t","")
        line = line.replace("|\n","")
        splits = line.split("|")
        if splits[-1] == 'scientific name':
            key = int( line.split("|")[0] )
            name = line.split("|")[1]
            result[key] = name

    return result

def buildNodes():
    return buildNodesLocal(DBpath+"nodes.dmp")

def buildNodesLocal(file):
    inFile = open(file)

    result = dict()
    levels = dict()

    for line in inFile:
        #print line
        line = line.replace("\t","")
        line = line.replace("|\n","")
        splits = line.split("|")
        #print splits
        taxaID = int(splits[0])
        parentID = int(splits[1])
        level = splits[2]
        result[taxaID] = (parentID, level)
        levels[level]=1

    return (result,levels)

def buildTaxaLevels( x , names, nodes, filter = None ):
    open = ""
    result = ""
    k = x
    result += names[k]
    (next, level) = nodes[k]
    k = next
    #EXAMPLE
    #  (((((Gammaretrovirus)Orthoretrovirinae)Retroviridae)Retro-transcribing_viruses)Viruses,
    z=0
    while True:
        (next, level) = nodes[k]
        if filter == None or level in filter:
            result = result + (")_%s_" % level) + names[k]
            open += "("
        elif z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::buildTaxaLevels could not complete for taxon\t"+str(x)
            break
        if k in BREAK_IDS:
            break
        k = next
        z+=1


    return open+result

#TODO restore to previous version
def buildTaxaLevels( x , names, nodes, filter = None ):
    open = ""
    result = ""
    k = x
    result += names[k]
    (next, level) = nodes[k]
    k = next
    #EXAMPLE
    #  (((((Gammaretrovirus)Orthoretrovirinae)Retroviridae)Retro-transcribing_viruses)Viruses,
    z=0
    while True:
        (next, level) = nodes[k]
        if filter == None or level in filter:
            result = result + ("\t_%s_" % level) + names[k]

        elif z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::buildTaxaLevels could not complete for taxon\t"+str(x)
            break
        if k in BREAK_IDS:
            break
        k = next
        z+=1


    return result

#returns a tab in front
#TODO can set if puts (x) first; remove the \t first
def buildTaxaLevels2(x , names, nodes, filter = None, xIsFirst=False ):
    open = ""
    result = ""
    #resultInt = []
    k = int(x)

    if xIsFirst:
      if k in names:
        result += "\t"+str(x)+"_"+names[k]+"\t"
    else:
      result += "\t-1_unknown\t"

    z=0
    while True:
	if k in nodes:
	  (next, level) = nodes[k]
	  if filter == None or level in filter:
	      result = result + ("\t%d_%s_" % (k, level)) + names[k]
	  elif z >= MAX_ITERATIONS:
	      print >> sys.stderr,"taxonomyLib::buildTaxaLevels could not complete for taxon\t"+str(x)
	      break
	  if k in BREAK_IDS:
	      break
	  k = next
	  z+=1
	else:
	  result+= str(k)+"_has_no_parent_node"
	  break
    return result


height = dict()
height["superkingdom"]=8
height["kingdom"]=7
height["phylum"]=6
height["class"]=5
height["order"]=4
height["family"]=3
height["genus"]=2
height["species"] = 1
#height["firstTaxon"] = 0

#TODO - update for MAX_ITERATIONS
def buildTaxaLevelList2( x , nodes, filter):
    #Returns index key list from highest level to most specific
    #reverse order is needed to build tree
    result = []
    k = x
    last = "species"
    while True:
        if k in nodes:
            (next, level) = nodes[k]
        else:
            (next, level) = 1, "no rank"
        if level in filter:
            distance = height[level]-height[last]
            result.append( (k, distance) )
            #last updates only when a level is used
            last = level
        if k in BREAK_IDS:
            break
        #k always updates to move next
        k = next

    #add prokaryota level
    if len(result) > 0 and result[-1] == 2:
        result.append( (1712345,1) )

    #add archaebacteria
    if len(result) > 0 and result[-1] == 2157:
        result.append( (18181818,1) )

    result.reverse()

    return result

def buildTaxaLevelList3( x , nodes, filter=height):
    #Returns index key list from x to species to kingdom
    # -1 if taxon is not found
    #reverse order is needed to build tree
    result = [-1 for z in range(len(filter)+1)]

    k = int(x)
    result[0] = k
    last = "species"
    z=0
    if k in nodes:
      (next, level) = nodes[k]
      if (level == "species"):
        index = filter[level]
        result[index] = k

    while True:
        if k in nodes:
            (next, level) = nodes[k]
        else:
            (next, level) = (1, "no rank")
        if level in filter:
            # distance = height[level]-height[last]
            index = filter[level]
            result[index] = k
            #last updates only when a level is used
            last = level
        if k in BREAK_IDS:
            break
        #k always updates to move next
        k = next
        z+=1
        if z >= MAX_ITERATIONS:
            break

    #add prokaryota level
    if len(result) > 0 and result[-1] == 2:
        result.append( (1712345,1) )

    #add archaebacteria
    if len(result) > 0 and result[-1] == 2157:
        result.append( (18181818,1) )

    result.reverse()

    return result



def buildTaxaLevelList( x , nodes):
    #Returns index key list from highest level to most specific
    #reverse order is needed to build tree
    result = []
    k = x

    z=0
    while True:
        result.append(k)
        if k in nodes:
            (next, level) = nodes[k]
        else:
            (next, level) = 1, "no rank"
        if k in BREAK_IDS:
            break
        k = next
        z+=1
        if z >= MAX_ITERATIONS:
            break

    #add prokaryota level
    if result[-1] == 2:
        result.append(1712345)

    #add archaebacteria
    if result[-1] == 2157:
        result.append(18181818)
    result.reverse()

    return result

def levelToText( level ):
    #note "u" used for 'super-" as in the German 'ueber'
    if level=="superkingdom":
        return "uk"
    elif level=="tribe":
        return "t"
    elif level=="subgenus":
        return "sg"
    elif level=="family":
        return "f"
    elif level=="species subgroup":
        return "ssg"
    elif level=="species group":
        return "sg"
    elif level=="phylum":
        return "p"
    elif level=="superclass":
        return "uc"
    elif level=="subphylum":
        return "sp"
    elif level=="subspecies":
        return "ss"
    elif level=="no rank":
        return "NA"
    elif level=="superorder":
        return "uo"
    elif level=="infraorder":
        return "io"
    elif level=="subclass":
        return "sc"
    elif level=="species":
        return "s"
    elif level=="superphylum":
        return "up"
    elif level=="kingdom":
        return "k"
    elif level=="subtribe":
        return "st"
    elif level=="subkingdom":
        return "sk"
    elif level=="forma":
        return "fr"
    elif level=="infraclass":
        return "ic"
    elif level=="varietas":
        return "v"
    elif level=="subfamily":
        return "sf"
    elif level=="class":
        return "c"
    elif level=="superfamily":
        return "uf"
    elif level=="parvorder":
        return "po"
    elif level=="suborder":
        return "so"
    elif level=="genus":
        return "g"
    elif level=="order":
        return "o"
    return '-1'


# return an integer
def getSpeciesID(taxonID, nodes):
    result = -1
    currID = int(taxonID)
    z=0
    while (currID in nodes):
        if nodes[currID][1] == "no rank":
            result = currID

        if nodes[currID][1] == "species":
            result = currID
            break

            #top nodes
            #2018-02-02 noticed tree no longer rooted to Taxon=1 (ROOT)
            # root (old school)     #bacteria#cellular organisms    #viruses           #eukaroyta
        if currID in BREAK_IDS:
            break
        currID = nodes[currID][0]

        z+=1
        if z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::getSpeciesID could not find taxonID\t"+str(taxonID)
            break
    return result

#####
#  Supposed to get Species and Genus taxon ID.  However, if there is no genus, it will run amok
#
####
def getSpeciesIDGenusID(taxonID, nodes):

    specID=-1
    genusID=-1

    z=0
    currID = int(taxonID)
    while (currID in nodes):
        if specID == -1 and nodes[currID][1] == "no rank":  #capture FIRST TAXON, hopefully to be overwritten
             result = currID

        if nodes[currID][1] == "species":
            specID = currID

        if nodes[currID][1] == "genus":
            genusID = currID
            break

            #top nodes
            #2018-02-02 noticed tree no longer rooted to Taxon=1 (ROOT)
            # root (old school)     #bacteria#cellular organisms    #viruses           #eukaroyta
        if currID in BREAK_IDS:
            break
        z+=1
        if z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::findSpeciesIDGenusID coudl not find taxon\t"+str(taxonID)+"\tresults spec\t"+str(specID)+"\tgenus\t"+str(genusID)
            break
        currID = nodes[currID][0]
    return (specID, genusID)



# returns dictionary mapping accession number to speciesID (int)
def extractSpeciesID(accessionNos, nodes):
    result = dict()
    for line in open(AccessionDBpath+"nucl_gb.accession2taxid"):
        splits = line.split()
        # 0 accession       1 accession.version     2 taxid     3 gid
        if splits[0] in accessionNos:
            result[splits[0]] = getSpeciesID( int(splits[2]),nodes)

    return result

def findTaxaAccessionNumbers( adict ):
    result = {}
    #iterate over file
    for line in open(AccessionDBpath+"master.accession2taxid"):
        splits = line.split()

        if splits[0] in adict:
            #print splits[0], " found ", splits[2]
            result[splits[2]]=1
    return result


def findTaxaAccessionNumbersMap( adict ):
    result = {}
    #iterate over file
    for line in open(AccessionDBpath+"master.accession2taxid"):
        splits = line.split()

        if splits[0] in adict:
            #print splits[0], " found ", splits[2]
            result[splits[0]]=splits[2]
    return result


def findKingdom(speciesID, names, nodes):
    result = "" #text name
    next=speciesID

    z=0
    while (True):
        if next in nodes:
            (next, level) = nodes[next]
            if level == "kingdom":
                result=names[next]
            if level == "superkingdom":
                result+="\t"+names[next]
                break
        if next in BREAK_IDS:
            break
        z+= 1
        if z >= MAX_ITERATIONS:
            print >> sys.stderr,"taxonomyLib::findKingdom could not find speciesID\t"+str(speciesID)
            break
    return result