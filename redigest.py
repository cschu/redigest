#!/usr/bin/env python
import sys
import re

HaeIII = ('GGCC', 'GG', 'CC')
EcoRI = ('GAATTC', 'G', 'AATTC')
BamHI = ('GGATCC', 'G', 'GATCC')
XbaI = ('TCTAGA', 'T', 'CTAGA')
XhoI = ('CTCGAG', 'C', 'TCGAG')
RESTRICTION_ENZYMES = [HaeIII, EcoRI, BamHI, XbaI, XhoI]

def readFasta(fn):
    with open(fn) as fi:
        seq, identifier = [], ''
        for line in fi:
            line = line.strip()
            if not line:
                # ignore empty lines
                continue
            if line.startswith('>'):
                # new record starts with >identifier
                if seq:
                    # if there is data in seq, return id, seq and empty seq
                    yield (identifier, ''.join(seq))
                    seq = []
                identifier = line
            else:
                # current line contains sequence information
                seq.append(line)
        if seq:
            # return last sequence record
            yield (identifier, ''.join(seq))

def findRESites(seq, enzymes):
    sites = set()
    for recseq, end1, end2 in enzymes:
        for site in re.finditer(recseq):
            sites.add(site.start() + len(end1))
    return sites

def main():
    fn = sys.argv[1]
    for id_, seq in readFasta(fn):
        resites = findRESites(seq, RESTRICTION_ENZYMES)
        for site in sorted(resites):
            print site
        break
