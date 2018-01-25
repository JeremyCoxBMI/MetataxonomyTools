DBpath="Z:/cox1kb/imsa2/data_files/gi_taxid/"


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