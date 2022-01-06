#!/usr/bin/python3
# Converts readmap.uc file from otupipe.bash into QIIME OTU map.


#OTU_PREFIX="OTU"

import sys
import numpy as np

UCFileName = sys.argv[1]
File = open(UCFileName)

#n = len(OTU_PREFIX)
OTUToReads = {}
LineNr = 0
while 1:
	Line = File.readline()
	if len(Line) == 0:
		break
	Line = Line.strip()
	LineNr += 1
	if Line[0] == '#':
		continue
	Fields = Line.split("\t")
	if len(Fields) < 10:
		print >> sys.stderr, "Line %d in .uc file has < 10 fields" % LineNr
		sys.exit(1)

#	if Fields[0] != 'H':
#		continue
	if Fields[9] == '*':
		Fields[9] = Fields[8]

	ReadLabel = Fields[8]
#	TargetLabel = Fields[9]
	OTU = Fields[9]
	
#	if not TargetLabel.startswith(OTU_PREFIX):
#		continue
	
#	nome = "OTU"
#	OTU = nome + (TargetLabel)

	try:
		OTUToReads[OTU].append(ReadLabel)
	except:
		OTUToReads[OTU] = [ ReadLabel ] # [ "" ]

for OTU in OTUToReads.keys():
	sys.stdout.write(str(OTU))
	all_reads=np.unique(OTUToReads[OTU])
#	all_reads=all_reads.remove('')
#    if len(all_reads)>1:
#		print(all_reads)
	for ReadLabel in all_reads:
		if ReadLabel!=str(OTU):
			sys.stdout.write(" " + ReadLabel)

	sys.stdout.write("\n")
