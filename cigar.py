
#GLOBAL CONSTANTS

class cigar:

    textToConsumeTuple = {}

    def __init__(self):
        self.textToConsumeTuple["M"] = (True,True)  #match or mismatch (DEPRECATED)
        self.textToConsumeTuple["I"] = (True,False) #insertion to reference
        self.textToConsumeTuple["D"] = (False,True) #deletion from referense
        self.textToConsumeTuple["N"] = (False,True) #skipped region reference
        self.textToConsumeTuple["S"] = (True,False) #soft clipping
        self.textToConsumeTuple["H"] = (False,False) # hard clipping
        self.textToConsumeTuple["P"] = (False,False) # padding
        self.textToConsumeTuple["="] = (True,True)   # sequence match
        self.textToConsumeTuple["X"] = (True,True)   # sequence mismatch

    def consumesQuery(self,char):
        return self.textToConsumeTuple[char][0]

    def consumesReference(self,char):
        return self.textToConsumeTuple[char][1]


