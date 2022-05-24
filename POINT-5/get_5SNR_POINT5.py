#!/usr/bin/env python

# Import libraries
import os
import pysam
import argparse
import random
import errno
import time
import string
from itertools import groupby
import re

###########################
##                       ##
##       Functions       ##
##                       ##
###########################

#CreateARandomString
def randomString(stringLength=10):
    return ''.join(random.sample(string.ascii_lowercase,stringLength))

def getReadsFirstPosition(flag, line):
    lineSplit = str(line).strip("\n").split("\t")
    if flag == "99":
        lineSplit[1] = str(0)
        lineSplit[2] = str(line.reference_name)
        lineSplit[3] = str(int(line.reference_start) + 1)
        lineSplit[6] = "*"
        lineSplit[7] = str(0)
        lineSplit[8] = str(0)
        lineSplit[9] = str(line.query_alignment_sequence[0])
        lineSplit[10] = str(chr(line.query_alignment_qualities[0] + 33))
        lineSplit[5] = "1M"
        return "\t".join(lineSplit[:11])

    if flag == "83":
        lineSplit[1] = str(16)
        lineSplit[2] = str(line.reference_name)
        lineSplit[3] = str(line.reference_end)
        lineSplit[6] = "*"
        lineSplit[7] = str(0)
        lineSplit[8] = str(0)
        lineSplit[9] = str(line.query_alignment_sequence[-1])
        lineSplit[10] = str(chr(line.query_alignment_qualities[-1] + 33))
        lineSplit[5] = "1M"
        return "\t".join(lineSplit[:11])

# Remove temporary files
def removeTemp(Id):
    filesInDir=os.listdir(os.getcwd())
    for f in filesInDir:
        if Id in f:
            os.remove(f)


###########################
##                       ##
##   Arguments Parsing   ##
##                       ##
###########################


# Input arguments
parser = argparse.ArgumentParser(description="""Creates a BAM file with the left most coordinate
    of Read 1. It expectes a second-stranded protocol, in which the retrived nucleotide correspond to the
    5' end of POINT-5 sequenced fragment.\nOutput file's name is 5pSNR_{POINT-5_BAM}.bam""")
    
parser.add_argument('BAM_POINT5', help='Path to POINT-5 bam file.')
args = parser.parse_args()

BAM_POINT5_basename = os.path.splitext(os.path.basename(args.BAM_POINT5))[0]
randomString=randomString()

print("Processing {} file!".format(BAM_POINT5_basename))

##########################
##                      ##
##    Main Body Code    ##
##                      ##
##########################

#bam to sam
infileHeader = pysam.AlignmentFile(args.BAM_POINT5, mode='rb').header
infile = pysam.AlignmentFile(args.BAM_POINT5, mode='rb', header=infileHeader)
outfile = pysam.AlignmentFile("{}.sam".format(randomString), "wh", header=infileHeader, template=infile)
for s in infile:
    outfile.write(s)
infile.close()
outfile.close()

#get read 1 5'
result=open("snrTemp_{}.sam".format(randomString),'w')
result.writelines([line for line in open("{}.sam".format(randomString), 'r').readlines() if line.startswith("@")])
flags = ["99", "83"]

workingSAMFile = pysam.AlignmentFile("{}.sam".format(randomString), mode='rb', header=infileHeader)

for line in workingSAMFile:
    lineSplit = str(line).strip("\n").split("\t")
    if "I" not in lineSplit[5] and "D" not in lineSplit[5] and lineSplit[1] in flags:
        BAMrowToWrite = getReadsFirstPosition(lineSplit[1], line)
        result.write(BAMrowToWrite)
        result.write("\n")

result.close()

#SAM to BAM
#new sam to bam, then sort and index
infile = pysam.AlignmentFile("snrTemp_{}.sam".format(randomString), "r")
outfile = pysam.AlignmentFile("5pSNR_{}.bam".format(randomString), "wb", template=infile)
for s in infile:
    outfile.write(s)
infile.close()
outfile.close()

pysam.sort('-@', '4', '-o', "Sorted_5pSNR_{}.bam".format(randomString) , "5pSNR_{}.bam".format(randomString))
os.rename("Sorted_5pSNR_{}.bam".format(randomString), "5pSNR_{}.bam".format(BAM_POINT5_basename))

removeTemp(randomString)



