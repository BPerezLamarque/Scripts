#!/usr/bin/python3
# Make stats from OTU fasta file


import sys

UCFileName = sys.argv[1]
File = open(UCFileName)

OTUToReads = {}
LineNr = 0
separator = ";size="
while 1:
	Line = File.readline()
	if len(Line) == 0:
		break
	Line = Line.strip()
	LineNr += 1
	if not Line.startswith(">"):
        	continue
    	
	amplicon, abundance = Line.strip(">;\n").split(separator)

	sys.stdout.write(str(int(((float(LineNr)-1)/2)+1)) + "\t" + abundance + "\t" + amplicon + "\t" + abundance)

	sys.stdout.write("\n")
