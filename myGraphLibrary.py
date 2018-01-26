__author__ = 'COX1KB'


# A lot of code borrowed heavily from manuscript_pIMSA_A_OCT2015.py
# copied from 2017-03-manuscript

from ete2 import Tree, TreeStyle, NodeStyle, faces, CircleFace, TextFace
import math
import copy
import re
import cPickle as pickle

LINEWIDTH = 3   #thickness of lines in tree
POINTSIZE = 6   #size of points; obviously, POINTSIZE > LINEWIDTH to appear
FONTSIZE = 18   #fontsize
LINNAEUS_FILTER = ["superkingdom","kingdom","phylum","class","order","family","genus","species"]

# Cladogram rendering inputs
WIDTH = 2400    # width in pixels (4 inches)
DPI = 600


#DBpath = "Z:/microbiome_genome_only/ncbi_gi_taxid/"
DBpath = "Z:/cox1kb/imsa2/data_files/gi_taxid/"
# DBpath = "Z:/pathoscope2/pathoscope/ncbiDB/"


DEBUGtaxa = False

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

def buildNamesPickle():
    outF = open(DBpath+"names.dmp.pickle",'wb')
    names = buildNames()
    pickle.dump(names,outF)


def buildNodesPickle():
    outF = open(DBpath+"nodes.dmp.pickle",'wb')
    (nodes ,levels) = buildNodes()
    pickle.dump(nodes,outF)
    outF.close()
    outF = open(DBpath+"levels.dmp.pickle",'wb')
    pickle.dump(levels, outF)

def loadNamesPickle():
    return pickle.load(open(DBpath+"names.dmp.pickle",'rb'))

def loadNodesPickle():
    nodes = pickle.load(open(DBpath+"nodes.dmp.pickle",'rb'))
    levels = pickle.load(open(DBpath+"levels.dmp.pickle",'rb'))
    return (nodes, levels)

def buildNodes():
    return buildNodes(DBpath+"nodes.dmp")

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
    while k != 1 and k != 131567:
        (next, level) = nodes[k]
        if filter == None or level in filter:
            result = result + (")_%s_" % level) + names[k]
            open += "("
        k = next
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


def buildTaxaLevelList2( x , nodes, filter):
    #Returns index key list from highest level to most specific
    #reverse order is needed to build tree
    result = []
    k = x
    last = "species"
    while k != 1 and k != 131567:
        if k in nodes:
            (next, level) = nodes[k]
        else:
            (next, level) = 1, "no rank"
        if level in filter:
            distance = height[level]-height[last]
            result.append( (k, distance) )
            #last updates only when a level is used
            last = level
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






def buildTaxaLevelList( x , nodes):
    #Returns index key list from highest level to most specific
    #reverse order is needed to build tree
    result = []
    k = x
    while k != 1 and k != 131567:
        result.append(k)
        if k in nodes:
            (next, level) = nodes[k]
        else:
            (next, level) = 1, "no rank"
        k = next

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

def buildTree( filename, names, nodes, filter = None ):
    result = Tree()
    result.name = "root"

    inFile = open(filename,'r')

    convertDist={}
    convertDist[1] = 1
    convertDist[2] = 2.22
    convertDist[3] = 3.43
    convertDist[4] = 4.63
    convertDist[5] = 5.85
    convertDist[6] = 7.05
    convertDist[7] = 8.25
    convertDist[8] = 9.46

    for line in inFile:
        #k = int(line.split()[0])    #2015-08-19
        k = int(line.split(',')[0])    #2016-04-05
        reverseList = buildTaxaLevelList2(k, nodes, filter)
        currentNode = result
        prevDistance=1
        for pair in reverseList:
            (k, distance) = pair
            (junk, level) = nodes[k]
            txt = levelToText(level)
            if DEBUGtaxa:
                name = txt+"_"+names[k]
            else:
                name = " "+names[k]

            #2016-04 change for readability
            #name = names[k]

            if filter == None or level in filter:
                kids = currentNode.get_children() # this is list
                found = False
                #look for name in children
                for m in range(len(kids)):
                    if kids[m].name == name:
                        found = True
                        currentNode = kids[m]
                        break
                if found == False:  # make a new node
                    #print "'%s' not found, adding" % name

                    #add child and returns child node
                    currentNode = currentNode.add_child(name=name, dist=convertDist[prevDistance])
                #Because moving up and down tree to add leaves, need to store previousDistance value as jump around
                prevDistance=distance
            #else skip and go on to next traversal in list

    return result


#takes a standard IMSA+A output counts file
# WARNING: clades lower than the lowest clade in filter will be counted as sum for that lowest clade in filter
# However, this does not count removing ties, so one should always use the report given at the level desired
def buildTreeWithCounts( filename, names, nodes, cutoff = 2, filter = None ):
    result = Tree()
    result.name = "root"

    inFile = open(filename,'r')

    convertDist={}
    convertDist[1] = 1
    convertDist[2] = 2.22
    convertDist[3] = 3.43
    convertDist[4] = 4.63
    convertDist[5] = 5.85
    convertDist[6] = 7.05
    convertDist[7] = 8.25
    convertDist[8] = 9.46

    # skip first line (header)
    inFile.readline()

    totalTaxon=0
    usedTaxon=0

    for line in inFile:
        #k = int(line.split()[0])    #2015-08-19
        k = int(line.split('\t')[0])    #2016-04-05
        count = int(line.split('\t')[2])

        if not k in {9606:"Homo Sapien", 40674:"Mammalia", -1:"Error", 9605:"Homo", 9604:"Homonidae"}:
            reverseList = buildTaxaLevelList2(k, nodes, filter)
            currentNode = result
            prevDistance=1
            #skip zero counts obviously
            #print count
            totalTaxon+=1
            if count < cutoff:
                #print "\ttaxID\t"+str(k)+"\tcutoff\t"+str(cutoff)+"\t>\t"+str(count)
                debug="ON"
            if count >= cutoff:
                usedTaxon+=1
                for pair in reverseList:
                    (k, distance) = pair
                    (junk, level) = nodes[k]
                    txt = levelToText(level)
                    if DEBUGtaxa:
                        name = txt+"_"+names[k]
                    else:
                        name = " "+names[k]

                    #2016-04 change for readability
                    #name = names[k]

                    if filter == None or level in filter:
                        kids = currentNode.get_children() # this is list
                        found = False
                        #look for name in children
                        for m in range(len(kids)):
                            if kids[m].name == name:
                                found = True
                                currentNode = kids[m]
                                if level == filter[-1]:
                                    currentNode.name += "-"+str(count)
                                break
                                #code for error handling, but summation calculation might be misleading
                            elif kids[m].name.split("-")[0] == name:
                                found = True
                                currentNode = kids[m]
                                s=kids[m].name.split("-")
                                n=s[0]
                                if level == filter[-1]:
                                    currentNode.name = n+"-"+str(int(s[1])+count)
                                else:
                                    currentNode.name = n
                                break
                        if found == False:  # make a new node
                            #print "'%s' not found, adding" % name

                            #add child and returns child node
                            currentNode = currentNode.add_child(name=name, dist=convertDist[prevDistance])
                            if level == filter[-1]:
                                currentNode.name += "-"+str(count)

                        #Because moving up and down tree to add leaves, need to store previousDistance value as jump around
                        prevDistance=distance
                    #else skip and go on to next traversal in list

    print filename+"\tused\t"+str(usedTaxon)+"\tof\t"+str(totalTaxon)+"\ttaxons"
    return result


def printNames( filename ):
    inFile = open(filename,'r')

    for line in inFile:
        k = int(line)
        print k, " ", names[k]

def applyNodeStyle(t, rootname="root",linewidth=1,bgcolor="White",dot_size=3,dot_color="Black"):
    #applies node style to sub-tree beginning with text name of root
    style = NodeStyle()
    style["bgcolor"] = bgcolor
    style["hz_line_width"] = linewidth
    style["vt_line_width"] = linewidth
    style["size"] = dot_size
    search = t.search_nodes(name=rootname)
    if len(search) >= 1:    #make sure if not found, there is no error
        search = search[0]
        for n in search.traverse():
            n.set_style(style)
            n.img_style["fgcolor"]=dot_color

def makeGraphsToFile(t, filenameStem, outputpath, count ):
    t_back = copy.deepcopy(t)

    #All Nodes
    applyNodeStyle(t,"root",LINEWIDTH,"White",POINTSIZE,"Black")
    #Other Kingdoms
    applyNodeStyle(t,"uk_Prokaryota",LINEWIDTH,"lightgreen",POINTSIZE,"Black")
    applyNodeStyle(t,"k_Fungi",LINEWIDTH,"wheat",POINTSIZE,"Black")
    #Stramenopiles
    applyNodeStyle(t,"o_Peronosporales",LINEWIDTH,"goldenrod",POINTSIZE,"Black")
    applyNodeStyle(t,"o_Saprolegniales",LINEWIDTH,"goldenrod",POINTSIZE,"Black")

    t2 = copy.deepcopy(t)

    for n in t.iter_leaves():
    #this creates text labels
        (control, infect) = count[n.name]
        #addition of spaces code: +' '
        # helps readability in a few cases, but overall stretches the graph
        T = TextFace(str(control)+' ', fsize=FONTSIZE, fgcolor='MediumBlue')
        n.add_face( T, 0, position="aligned" )
        T = TextFace(str(infect)+' ', fsize=FONTSIZE, fgcolor='FireBrick')
        n.add_face( T, 1, position="aligned" )
        #T = TextFace(str(infect+control)+' ', fsize=10, fgcolor='black')
        #n.add_face( T, 2, position="aligned" )
        T = TextFace(" "+n.name+" ",fsize=(FONTSIZE+2),fgcolor='Black') #add a space so not too crowded
        n.add_face( T, 2, position="aligned" )


    circular_style = TreeStyle()
    circular_style.mode = "c" # draw tree in circular mode
    circular_style.scale = 20
    circular_style.scale = 31   #length of 1 level transition in tree
    circular_style.show_scale = False
    circular_style.show_leaf_name = False
    #circular_style.allow_face_overlap = True
    #t.show(tree_style=circular_style)
    t.render(outputpath+filenameStem+"_color_v1.png", tree_style = circular_style, w=WIDTH, dpi=DPI)

    ### COLOR -- Alternate ordering of text labels
    ### using copied t2

    '''    CAUSING BUG -- don't need to reload from wrong file
    count = dict()
    input = open(filename)
    for line in input:
        line = line.split()
        name = line[1]
        control = (int( line[2]))
        infect = (int(line[3]))
        count[name] = (control, infect)
    '''

    for n in t2.iter_leaves():
    #this creates text labels
        (control, infect) = count[n.name]
        #addition of spaces code: +' '
        # helps readability in a few cases, but overall stretches the graph
        T = TextFace(str(control)+' ', fsize=FONTSIZE, fgcolor='MediumBlue')
        n.add_face( T, 1, position="aligned" )
        T = TextFace(str(infect)+' ', fsize=FONTSIZE, fgcolor='FireBrick')
        n.add_face( T, 2, position="aligned" )
        #T = TextFace(str(infect+control)+' ', fsize=10, fgcolor='black')
        #n.add_face( T, 2, position="aligned" )
        T = TextFace(" "+n.name+" ",fsize=(FONTSIZE+2),fgcolor='Black') #add a space so not too crowded
        n.add_face( T, 0, position="aligned" )
    t2.render(outputpath+filenameStem+"_color_v2.png", tree_style = circular_style, w=WIDTH, dpi=DPI)


    t = copy.deepcopy(t_back)

    ###GREYSCALE
    #t = buildTree( filename, names, nodes, filter = taxa_accepted )
    t = copy.deepcopy(t_back)
    #All Nodes
    applyNodeStyle(t,"root",LINEWIDTH,"White",POINTSIZE,"Black")
    #Other Kingdoms
    #applyNodeStyle(t,"uk_Prokaryota",LINEWIDTH,"White",POINTSIZE,"Black") #already defined by all nodes
    applyNodeStyle(t,"k_Fungi",LINEWIDTH,"Silver",POINTSIZE,"Black")
    applyNodeStyle(t,"f_Retroviridae",LINEWIDTH,"GainsBoro",POINTSIZE,"Black")
    #Stramenopiles
    applyNodeStyle(t,"o_Peronosporales",LINEWIDTH,"DarkGrey",POINTSIZE,"Black")
    applyNodeStyle(t,"o_Saprolegniales",LINEWIDTH,"DarkGrey",POINTSIZE,"Black")

    t2 = copy.deepcopy(t)

    for n in t.iter_leaves():
    #this creates text labels
        #BUG OCCURS HERE: can't find citrobacter, but it isn't in the previous version of t? (above)
        (control, infect) = count[n.name]
        #addition of spaces code: +' '
        # helps readability in a few cases, but overall stretches the graph
        T = TextFace(str(control)+' ', fsize=FONTSIZE, fgcolor='Black')
        n.add_face( T, 1, position="aligned" )
        T = TextFace(str(infect)+' ', fsize=FONTSIZE, fgcolor='DimGray')
        n.add_face( T, 2, position="aligned" )
        #T = TextFace(str(infect+control)+' ', fsize=10, fgcolor='black')
        #n.add_face( T, 2, position="aligned" )
        T = TextFace(" "+n.name+" ",fsize=(FONTSIZE+2),fgcolor='Black') #add a space so not too crowded
        n.add_face( T, 0, position="aligned" )

    #t.show(tree_style=circular_style)
    t.render(outputpath+filenameStem+"_grey_v1.png", tree_style = circular_style, w=WIDTH, dpi=DPI)

    ### GREY -- ALTERNATE ordering of text labels
    ### using copied t2
    for n in t2.iter_leaves():
    #this creates text labels
        (control, infect) = count[n.name]
        #addition of spaces code: +' '
        # helps readability in a few cases, but overall stretches the graph
        T = TextFace(str(control)+' ', fsize=FONTSIZE, fgcolor='Black')
        n.add_face( T, 0, position="aligned" )
        T = TextFace(str(infect)+' ', fsize=FONTSIZE, fgcolor='DimGray')
        n.add_face( T, 1, position="aligned" )
        #T = TextFace(str(infect+control)+' ', fsize=10, fgcolor='black')
        #n.add_face( T, 2, position="aligned" )
        T = TextFace(" "+n.name+" ",fsize=(FONTSIZE+2),fgcolor='Black') #add a space so not too crowded
        n.add_face( T, 2, position="aligned" )

    #t.show(tree_style=circular_style)
    t2.render(outputpath+filenameStem+"_grey_v2.png", tree_style = circular_style, w=WIDTH, dpi=DPI)

def sumListByIndex( aList, indexes ):
    result = 0
    for k in indexes:
        if len(aList[k])>0:   #blank is same as zero
            result += float(aList[k])
    return result


OUTPUTPATH="C:/Users/cox1kb/PycharmProjects/ETE_Tutorial/"
INPUTPATH="C:/Users/cox1kb/PycharmProjects/ETE_Tutorial/"


def recurseAlphabetizeNodes( node ):

    kids = node.get_children()
    if kids == None or len(kids) < 1:
        return

    names = []
    idx=0
    for kid in kids:
        names.append(   (kid.name.replace(" ",""), idx)  )
        idx+=1

    names.sort()
    new_kids = [0 for n in names]

    nidx = 0
    for n in names:
        (name, idx) = n
        new_kids[nidx] = kids[idx]
        recurseAlphabetizeNodes(new_kids[nidx])
        nidx+=1

    debug = -1
    return

# countsByTaxon: dictionary taxonID ==> count (int)
# filter: list of clades (text) to include in descending order. Thus last one is the outer level of cladogram.
def buildTreeGeneral(countsByTaxon, names, nodes, filter):  ##need sort nodes by
    result = Tree()
    result.name = "root"

    convertDist={}
    convertDist[1] = 1
    convertDist[2] = 2.22
    convertDist[3] = 3.43
    convertDist[4] = 4.63
    convertDist[5] = 5.85
    convertDist[6] = 7.05
    convertDist[7] = 8.25
    convertDist[8] = 9.46

    cladeToGraph = filter[-1]
    countsByTaxonBaseClade = {}
    for k in countsByTaxon:
        myCount = int(countsByTaxon[k])
        reverseList = buildTaxaLevelList2(k, nodes, filter)
        for pair in reverseList:
            (k, distance) = pair
            (junk, level) = nodes[k]
            if level == cladeToGraph:
                if k in countsByTaxonBaseClade:
                    countsByTaxonBaseClade[k] += myCount
                else:
                    countsByTaxonBaseClade[k] = myCount
                break

    for k in countsByTaxonBaseClade:
        #k = int(line.split()[0])    #2015-08-19
        #k = int(line.split(',')[0])    #2016-04-05
        myCount = int(countsByTaxonBaseClade[k])
        reverseList = buildTaxaLevelList2(k, nodes, filter)
        currentNode = result
        prevDistance=1
        for pair in reverseList:
            (k, distance) = pair
            (junk, level) = nodes[k]
            txt = levelToText(level)
            if DEBUGtaxa:
                name = " "+txt+"_"+names[k]
            else:
                if int(myCount) >= 0:
                    if level == cladeToGraph:
                        if k in names:
                            name = " "+names[k]+" "+str(myCount)
                        else:
                            name = "+"+str(k)+" "+str(myCount)
                    else:
                        if k in names:
                            name = names[k]
                        else:
                            name = str(k)

            #2016-04 change for readability
            #name = names[k]

            if filter == None or level in filter:
                kids = currentNode.get_children() # this is list
                #kids = kids.sort()
                found = False
                #look for name in children
                for m in range(len(kids)):
                    if kids[m].name == name:
                        found = True
                        currentNode = kids[m]
                        break
                if found == False:  # make a new node
                    #print "'%s' not found, adding" % name

                    #add child and returns child node
                    currentNode = currentNode.add_child(name=name, dist=convertDist[prevDistance])
                #Because moving up and down tree to add leaves, need to store previousDistance value as jump around
                prevDistance=distance
            #else skip and go on to next traversal in list

    return result


# countsByTaxon: dictionary taxonID ==> count (int)
# filter: list of clades (text) to include in descending order. Thus last one is the outer level of cladogram.
# subgroup must be a child of the root node
def buildTreeSubgroup(countsByTaxon, names, nodes, filter, subgroup, includeZeroes = False):  ##need sort nodes by
    result = Tree()
    result.name = "root"

    convertDist={}
    convertDist[1] = 1
    convertDist[2] = 2.22
    convertDist[3] = 3.43
    convertDist[4] = 4.63
    convertDist[5] = 5.85
    convertDist[6] = 7.05
    convertDist[7] = 8.25
    convertDist[8] = 9.46

    cladeToGraph = filter[-1]
    countsByTaxonBaseClade = {}
    for k in countsByTaxon:
        myCount = int(countsByTaxon[k])
        reverseList = buildTaxaLevelList2(k, nodes, filter)
        for pair in reverseList:
            (k, distance) = pair
            (junk, level) = nodes[k]
            if level == cladeToGraph:
                if k in countsByTaxonBaseClade:
                    countsByTaxonBaseClade[k] += myCount
                else:
                    countsByTaxonBaseClade[k] = myCount
                break

    for k in countsByTaxonBaseClade:
        #k = int(line.split()[0])    #2015-08-19
        #k = int(line.split(',')[0])    #2016-04-05
        myCount = int(countsByTaxonBaseClade[k])
        reverseList = buildTaxaLevelList2(k, nodes, filter)
        currentNode = result
        prevDistance=1
        for pair in reverseList:
            (k, distance) = pair
            (junk, level) = nodes[k]
            txt = levelToText(level)
            if DEBUGtaxa:
                name = " "+txt+"_"+names[k]
            else:
                #default for failing logic below
                name = " "+names[k]+" "+str(myCount)

                if int(myCount) > 0 or (includeZeroes and int(myCount) == 0):
                    if level == cladeToGraph:
                        if k in names:
                            name = " "+names[k]+" "+str(myCount)
                        else:
                            name = "+"+str(k)+" "+str(myCount)
                    else:
                        if k in names:
                            name = names[k]
                        else:
                            name = str(k)

            #2016-04 change for readability
            #name = names[k]

            if filter == None or level in filter:
                kids = currentNode.get_children() # this is list
                #kids = kids.sort()
                found = False
                #look for name in children
                for m in range(len(kids)):
                    if kids[m].name == name:
                        found = True
                        currentNode = kids[m]
                        break
                if found == False:  # make a new node
                    #print "'%s' not found, adding" % name

                    #add child and returns child node
                    if name == None or name =="":
                        debug = "test"
                    currentNode = currentNode.add_child(name=name, dist=convertDist[prevDistance])
                #Because moving up and down tree to add leaves, need to store previousDistance value as jump around
                prevDistance=distance
            #else skip and go on to next traversal in list

    #very wasteful of computations, filter afterwards
    for x in result.get_children():
        if x.name == subgroup:
            result = x
            break

    return result




def findTaxonAtClade(taxon, clade, nodes ):
    k = taxon
    while k != 1 and k != 131567:
        if k in nodes:
            (next, level) = nodes[k]
            k = next
        else:
            return None
        if k == 1:
            return None
        if level == clade:
            return k

def findPhylum( taxon, names,nodes ):

    k = taxon
    last = "species"
    while k != 1 and k != 131567:
        if k in nodes:
            (next, level) = nodes[k]
        else:
            print "ERROR 36: \t" + str(k) + "\t taxon \t"+str(taxon)
            print buildTaxaLevels( k , names, nodes, None)
            return None

        if level == "phylum" or k == 2 or k == 1:
            if names[k] == "None" or names[k] == None:
                print "ERROR 26: " + str(k)
            return k
        k = next


def findGenus( taxon, names,nodes ):

    k = taxon
    last = "species"
    while k != 1 and k != 131567:
        if k in nodes:
            (next, level) = nodes[k]
        else:
            print "ERROR 36: \t" + str(k) + "\t taxon \t"+str(taxon)
            print buildTaxaLevels( k , names, nodes, None)
            return None

        if level == "genus" or k == 2 or k == 1:
            if names[k] == "None" or names[k] == None:
                print "ERROR 26: " + str(k)
            return k
        k = next


# apply style to node and all descendants based on substring search
def applyNodeStyle2(t, rootname="root",linewidth=1,bgcolor="White",dot_size=3,dot_color="Black"):
    #applies node style to sub-tree beginning with text name of root
    style = NodeStyle()
    style["bgcolor"] = bgcolor
    style["hz_line_width"] = linewidth
    style["vt_line_width"] = linewidth
    style["size"] = dot_size

    for n in t.traverse():
        if rootname in n.name:
            # n.set_style(style)
            # n.img_style["fgcolor"]=dot_color
            for x in n.traverse():
                x.set_style(style)
                x.img_style["fgcolor"]=dot_color
            break

# apply style to node and all descendants based on substring search
def applyNodeStylePhyla(t, rootname="root",linewidth=1,bgcolor="White",dot_size=3,dot_color="Black"):
    #applies node style to sub-tree beginning with text name of root
    style = NodeStyle()
    style["bgcolor"] = bgcolor
    style["hz_line_width"] = linewidth
    style["vt_line_width"] = linewidth
    style["size"] = dot_size

    even = True

    for n in t.traverse():
        if rootname in n.name:
            # n.set_style(style)
            # n.img_style["fgcolor"]=dot_color

            for x in n.get_children():
            #for x in n.traverse():
                #color every other phylum differently
                if even:
                    x.set_style(style)
                    x.img_style["fgcolor"]=dot_color
                    for y in x.traverse():
                        y.set_style(style)
                        y.img_style["fgcolor"]=dot_color
                even = not even
            break

leafNameColor = 'MediumBlue'
leafCountsColor = 'FireBrick'

def labelsOnTree( tree, splitChar=" " ):
    if tree == None:
        return

    for n in tree.iter_leaves():
        #print n.name
        # if re.search("-", n.name) != None:
        #     splits = n.name.split("-")
        # else:
        #     splits = n.name.split(" ")
        splits = n.name.split(splitChar)
        if len(splits) == 1:
            T1 = TextFace(splits[0], fsize=FONTSIZE, fgcolor=leafNameColor)
            T2 = TextFace("0", fsize=FONTSIZE, fgcolor=leafCountsColor)
        else:
            T1 = TextFace(" ".join(splits[:-1]), fsize=FONTSIZE, fgcolor=leafNameColor)
            T2 = TextFace(splits[-1], fsize=FONTSIZE, fgcolor=leafCountsColor)
        n.add_face( T1, 0, position="aligned" )

        n.add_face( T2, 1, position="aligned" )

def importFileData(file, countsByTreatmentCladeTaxonTOcounts, orderByTreatmentCladeTOtaxon, filterClade, levels, names, nodes, filterCutOff = 20):

    first = True
    labels = []


    for line in open(file):  # raw species data file
        splits = line.strip().split('\t')
        if first:
            labels = splits[1:]
            first = False
        else:
            taxonID = int(splits[0])

            k = taxonID

            while k != 1 and k != 131567 and k != 2:
                if k in nodes:
                    (next, level) = nodes[k]
                else:
                    level = ""
                    next = -1
                    print "ERROR 36: \t" + str(k) + "\t taxon \t"+str(taxonID)
                    print buildTaxaLevels( k , names, nodes, None)
                    return None

                # counted = []


                if level in filterClade:
                    for i in range(len(countsByTreatmentCladeTaxonTOcounts)): #number columns
                        c = int(splits[i+1])
                        idx = levels.index(level)
                        orderByTreatmentCladeTOtaxon[i][idx].append(taxonID)
                        # if level == "phylum":
                        #     debug = 1
                        # if taxonID in countsByTreatmentCladeTaxonTOcounts[i][idx]:
                        #     debug = 1
                        # if taxonID == 735:


                        if c > filterCutOff:
                            #list of treatment indexes to dictionaries of taxon to list of counts
                            if not taxonID in countsByTreatmentCladeTaxonTOcounts[i][idx]:
                                countsByTreatmentCladeTaxonTOcounts[i][idx][taxonID] = 0
                            countsByTreatmentCladeTaxonTOcounts[i][idx][taxonID] += c
                        # else:
                        #     countsByTreatmentCladeTaxonTOcounts[i][idx][taxonID] += 0

                    # idx = levels.index(level)
                    # countsByTreatmentCladeTaxonTOcounts[idx][k] = counted

                k = next
                if k == -1:
                    break

    return labels






if __name__== "__main__":

    print "Step 100"
    #BEGIN NEW CODE
    file1 = "Z:/cox1kb/2016.01.simulations/2016.04.specieslists/B50.species.list.csv"
    file2 = "Z:/cox1kb/2016.01.simulations/2016.04.specieslists/E5A.species.list.csv"
    file3 = "Z:/cox1kb/2016.01.simulations/2016.04.specieslists/E5B.species.list.csv"
    file4 = "Z:/cox1kb/2016.01.simulations/2016.04.specieslists/E5C.species.list.csv"

    names = buildNames()
    (nodes, levels) = buildNodes()

    print "Step 200"

    #Apparently, the NCBI taxa tree is missing "Prokaryota"
    #made up my own number
    names[1712345] = "Prokaryota"
    nodes[1712345] = (1,"kingdom")
    names[18181818] = "Archaebacteria"
    nodes[18181818] = (2157,"kingdom")



    names[1649298] = "Siccibacter"
    nodes[1649298] = (543,"genus")
    names[1616789] = "Salinispira"
    nodes[1616789] = (137,"genus")
    names[1621534] = "Paraglaciecola"
    nodes[1621534] = (72275, "genus")
    names[1637257] = "Mageeibacillus"
    nodes[1637257] = (541000,"genus")
    names[1652133] = "Halobacteriovorax"
    nodes[1652133] = (1652133,"genus")

    taxa_accepted = ["superkingdom","kingdom","phylum","class","order","family","genus","species"]


    #TESTS
    print buildTaxaLevels( 2426, names, nodes, filter = taxa_accepted )
    print buildTaxaLevels( 1541743, names, nodes, filter = taxa_accepted )
    print buildTaxaLevels( 301375, names, nodes, filter = taxa_accepted )
    print buildTaxaLevels( 4792, names, nodes, filter = taxa_accepted )
    print buildTaxaLevels( 5180, names, nodes, filter = taxa_accepted )
    print buildTaxaLevels( 36427, names, nodes, filter = taxa_accepted )

    print buildTaxaLevelList2( 1541743 ,nodes, filter = taxa_accepted )

    print "Step 300"
    (tB50) = buildTree( file1, names, nodes, filter = taxa_accepted )


    print "Step 400"
    (tE5A) = buildTree( file2, names, nodes, filter = taxa_accepted )

    print "Step 420"
    (tE5B) = buildTree( file3, names, nodes, filter = taxa_accepted )

    print "Step 440"
    (tE5C) = buildTree( file4, names, nodes, filter = taxa_accepted )


    print "Step 500"
    #tB50.show()
    #tE5A.show()


    square_style = TreeStyle()
    square_style.scale = 31   #length of 1 level transition in tree
    square_style.show_scale = False

    tB50square = copy.deepcopy(tB50)
    applyNodeStyle(tB50square,"root",LINEWIDTH,"White",POINTSIZE,"Black")
    applyNodeStyle(tB50square," Bacteria",LINEWIDTH,"Silver",POINTSIZE,"Black")
    #applyNodeStyle(tB50square,"Archaebacteria",LINEWIDTH,"LightGrey",POINTSIZE,"Black")
    #Archaebacteria is white
    applyNodeStyle(tB50square," Eukaryota",LINEWIDTH,"Gainsboro",POINTSIZE,"Black")
    applyNodeStyle(tB50square," Viruses",LINEWIDTH,"DarkGrey",POINTSIZE,"Black")

    applyNodeStyle(tE5A,"root",LINEWIDTH,"White",POINTSIZE,"Black")
    applyNodeStyle(tE5A," Bacteria",LINEWIDTH,"Silver",POINTSIZE,"Black")
    #applyNodeStyle(tE5A,"Archaebacteria",LINEWIDTH,"LightGrey",POINTSIZE,"Black")
    #Archaebacteria is white
    applyNodeStyle(tE5A," Eukaryota",LINEWIDTH,"Gainsboro",POINTSIZE,"Black")
    applyNodeStyle(tE5A," Viruses",LINEWIDTH,"DarkGrey",POINTSIZE,"Black")

    applyNodeStyle(tE5B,"root",LINEWIDTH,"White",POINTSIZE,"Black")
    applyNodeStyle(tE5B," Bacteria",LINEWIDTH,"Silver",POINTSIZE,"Black")
    applyNodeStyle(tE5B," Eukaryota",LINEWIDTH,"Gainsboro",POINTSIZE,"Black")
    applyNodeStyle(tE5B," Viruses",LINEWIDTH,"DarkGrey",POINTSIZE,"Black")

    applyNodeStyle(tE5C,"root",LINEWIDTH,"White",POINTSIZE,"Black")
    applyNodeStyle(tE5C," Bacteria",LINEWIDTH,"Silver",POINTSIZE,"Black")
    applyNodeStyle(tE5C," Eukaryota",LINEWIDTH,"Gainsboro",POINTSIZE,"Black")
    applyNodeStyle(tE5C," Viruses",LINEWIDTH,"DarkGrey",POINTSIZE,"Black")




    FONTSIZE = 24
    for leaf in tB50square.iter_leaves():
        T = TextFace(' '+leaf.name, fsize=(FONTSIZE), fgcolor='Black')
        leaf.add_face( T, 0, position="aligned" )

    #square_style.show_branch_length = True;
    #tB50square.show(tree_style = square_style)
    tB50 = copy.deepcopy(tB50square)

    for leaf in tE5A.iter_leaves():
        T = TextFace(' '+leaf.name, fsize=(FONTSIZE), fgcolor='Black')
        leaf.add_face( T, 0, position="aligned" )

    for leaf in tE5B.iter_leaves():
        T = TextFace(' '+leaf.name, fsize=(FONTSIZE), fgcolor='Black')
        leaf.add_face( T, 0, position="aligned" )

    for leaf in tE5C.iter_leaves():
        T = TextFace(' '+leaf.name, fsize=(FONTSIZE), fgcolor='Black')
        leaf.add_face( T, 0, position="aligned" )

    square_style.show_leaf_name = False

    tB50square.render(OUTPUTPATH+"B50_square_BW.png", tree_style = square_style, w=1200, dpi=DPI)
    #tB50square.render(OUTPUTPATH+"B50_circle_BW.png", tree_style = circular_style, w=2400, dpi=DPI)
    tE5A.render(OUTPUTPATH+"E5A_square_BW.png", tree_style = square_style, w=1200, dpi=DPI)

    tE5B.render(OUTPUTPATH+"E5B_square_BW.png", tree_style = square_style, w=1200, dpi=DPI)
    tE5C.render(OUTPUTPATH+"E5C_square_BW.png", tree_style = square_style, w=1200, dpi=DPI)


    circular_style = TreeStyle()
    circular_style.mode = "c" # draw tree in circular mode
    circular_style.scale = 20
    circular_style.scale = 31   #length of 1 level transition in tree
    circular_style.show_scale = False
    circular_style.show_leaf_name = False
    #t.render(outputpath+filenameStem+"_color_v1.png", tree_style = circular_style, w=WIDTH, dpi=DPI)

    FONTSIZE = 24
    for leaf in tB50.iter_leaves():
        T = TextFace(' '+leaf.name, fsize=(FONTSIZE), fgcolor='Black')
        leaf.add_face( T, 0, position="aligned" )

    #tB50.render(OUTPUTPATH+"B50_circle_BW.png", tree_style = circular_style, w=2400, dpi=DPI)