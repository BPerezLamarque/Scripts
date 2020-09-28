#!/usr/bin/env python3

import os


contaminant=[]
for i in open("contamination.txt", "r").readlines():
    contaminant=contaminant+[i.replace("\n","")]



for i in [41, 43, 46, 50, 54, 56, 76, 110, 116, 68, 73, 86, 105, 106, 112, 114, 51, 55, 84, 99, 117, 120, 45, 62, 66, 71, 77, 78, 88, 98, 100, 102, 104, 107, 118, 48, 49, 53, 58, 64, 67, 72, 81, 109, 113, 119, 42, 44, 69, 79, 87, 97, 108, 111, 115, 47, 52, 59, 60, 61, 65, 70, 74, 80, 82, 83, 85, 101, 57, 63, 103]:
    print(i)
    with open("OTU_S"+str(i)+"_16S_table.txt") as oldfile, open("../OTU_table/OTU_S"+str(i)+"_16S_table.txt", 'w') as newfile:
        for line in oldfile:
            if not any(bad_word in line for bad_word in contaminant):
                newfile.write(line)
