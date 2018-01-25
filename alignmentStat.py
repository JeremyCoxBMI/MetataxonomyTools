import math
import re
import cigar as cg





# FRIEND functions



## flag is integer
## idx is bit index
def getBit(flag, idx):
    p = math.pow(2,idx)
    q = p * 2
    return (    (flag / q) >= p    )

#place is 1's place, 2's place, 4's place... etc
def getBitPlaceValue(flag, place):
    q = place * 2
    return (    (flag % q) >= place    )


class alignmentStat:

    cigarMErrorPrinted = False

    queryMatchingBases={}
    line = ""
    type = ""
    queryName = ""
    refName = ""

    cigar=""
    mdz=""
    #TODO: difference versus queryIdent?
    pident = 0.0

    # TODO verify
    # intended to be reference (target) length
    length = 0
    mismatch = 0
    match = 0

    #positions are 1-index
    rstart = 0
    rend = 0        #InCLUSIVE
    qstart = 0
    qend = 0        #INCLUSIBE
    queryIdent = 0.0

    queryLength=0

    # line is input line
    # type is note about the kind of line:  "blast6", "sam", "psl"  This is human-readable label only at this point, so whatever you want works
    def __init__(self, line, type):
        self.line = line
        self.type = type
        self.queryName = ""

        self.refName = ""
        self.pident = 0.0

        # the meaning of length is somewhat dubious - length of query or reference alignment region? or length of matches+mismatches?
        self.length = 0
        self.mismatch = 0
        self.match = 0

        #positions are 1-index
        self.rstart = 0
        self.rend = 0        #InCLUSIVE
        self.qstart = 0
        self.qend = 0        #INCLUSIBE
        self.queryIdent = 0.0

        self.queryLength=0
        self.cigar=""
        self.mdz=""
    def blast6ProcessLine(self):
        splits = self.line.strip().split("\t")

        if (len(splits) >= 10):
          self.queryName = splits[0]
          self.refName = splits[1]
          self.pident = float(splits[2])
          self.length = int(splits[3])
          self.mismatch = int(splits[4])

          #self.match = self.length - self.mismatch

          self.qstart = int(splits[6])
          self.qend = int(splits[7])
          self.rstart = int(splits[8])
          self.rend = int(splits[9])
          # qAlignLen =
          self.match = self.qstart - self.qend +1 - self.mismatch
          self.queryIdent = float(self.match)/self.length

    def isHuman(self):
        return self.refName[0] == "c" and self.refName[1] == "h" and self.refName[2] == "r"

    def header(self):
        return "queryName\trefName\tpident\tmatches\tpident\tquery_ident"

    def __str__(self):
        return self.queryName+"\t"+self.refName+"\t"+str(self.pident)+"\t"+str(self.match)+"\t"+str(self.pident) \
                        +"\t"+str(self.queryIdent)

    def isValid(self):
        return self.queryName != "" #not initialized

    #TODO broken to some degree, gotta fix it
    def samProcessLine(self):
        splits = self.line.strip().split()
        if (len(splits) >= 11 and self.line[0] != "@"):

          flag=int(splits[1])
          #based on flag, may be null result
          if not getBitPlaceValue(flag, 4):

              self.queryName = splits[0]
              flag=int(splits[1])
              self.refName = splits[2]
              self.rstart = int(splits[3])


              #parse cigar
              cigar = cg.cigar()
              self.cigar=splits[5]

              #MDZ string
              # mdz="31T6T19C3^TGAT4C3T2C26"
              # numbers = re.search(r'[0-9]+',mdz)
              for s in splits:
                  #print s
                  if s[0:5]=="MD:Z:":
                    self.mdz = s[5:]
                    break
              #numbers = re.findall(r'[0-9]+',mdz)

              # for n in numbers:
              #     self.match += int(n)

              self.samMatchesCoor()

              ## does not work
              # self.queryLength = len(splits[9])

              d = (self.qend - self.qstart)
              if d==0: d=1
              self.pident = float(self.match) / d
              self.length = int(splits[8])
              self.rend = self.rstart + self.length - 1

    def pslProcessLine(self):
        #TODO
        splits=self.line.strip().split()
        self.match=int(splits[0])
        self.mismatch=int(splits[1])
        self.queryName=splits[9]
        self.queryLength=int(splits[10])
        self.refName=splits[13]

        self.qstart=int(splits[11])
        self.qend=int(splits[12])
        self.restart=int(splits[15])
        self.rend=int(splits[16])
        self.length=self.rend-self.rstart+1

        self.queryIdent = float(self.match)/self.queryLength

        #todo correct?
        self.pident = self.queryIdent

    # cutoff is number matches in a "sub alignment" required to be included
    # positions : a dictionary of keys as coordinates
    def samMatchesCoor(self, cutoff=1):
        numbers = re.findall(r'[0-9]+',self.mdz)
        letters = re.findall(r'[\^,A-Z]+',self.mdz)

        pos=1
        l = 0

        #if not (self.mdz[0] >= "0" and self.mdz[0] <= "9"):
        if len(self.mdz)==0:
            print self.line
            print self.mdz
            debuggingLINE="HERE"
        elif not self.mdz[0].isdigit():
            pos += len(letters[l])
            if letters[l][0] == "^":
                pos-=1
            l += 1


        positions = {}

        self.match=0
        for k in range(len(numbers)):
            matches=int(numbers[k])

            if k == 0:
                self.qstart = 0
                #also refStart?  unclear
            if matches >= cutoff:
                for z in range(matches):
                    positions[ pos+z ] = 1
            pos+=matches
            self.match+=matches

            if k == len(numbers)-1:
                self.qend = pos - 1

            if l < len(letters):
                pos += len(letters[l])
                if letters[l][0] == "^":
                    pos-=1
                l+=1
        self.queryMatchingBases = positions

              # setQstart = False
              # setRstart = False
              # rlen = 0
              # prev=1
              # splitCig = re.split("(\d+)",splits[5])
              # for idx in range(1, len(splitCig), 2):
              #     code = splitCig[idx+1]
              #     k = int(splitCig[idx])
              #     # if code == "M":
              #     #     if not alignmentStat.cigarMErrorPrinted:
              #     #         print "\t",str(flag),"\tWARNING: (old-style) code M is ambiguous : match or mismatch.  Assumign match.  Analysis possibly invalid"
              #     #         #print self.line
              #     #         alignmentStat.cigarMErrorPrinted = True
              #     #     self.match += k
              #     # elif code == "=":
              #     #     self.match += k
              #     # elif code == "X":
              #     #     self.mismatch += k
              #
              #
              #     if cigar.consumesQuery(code):
              #         if not setQstart:
              #             self.qstart = prev
              #             setQstart = True
              #         self.queryLength += k
              #         self.qend += prev + k
              #     if cigar.consumesReference(code):
              #         if not setRstart:
              #             #TODO is this right?  or does the pos report correctly?
              #             self.rstart += k
              #             if code == "M":
              #                 setRstart=True
              #         rlen += k
              #     prev+=k
              # self.rend = self.rstart + rlen - 1