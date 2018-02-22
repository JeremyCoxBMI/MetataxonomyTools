__author__ = 'COX1KB'

import sys

outF = open(sys.argv[1]+".revised")

keeps={}
for line in open(sys.argv[1]):
    splits = line.split()
    if not splits[1] in keeps:
        keeps[splits[1]] = line

for t in keeps:
    outF.write(keeps[t])
