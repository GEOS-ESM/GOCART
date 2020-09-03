#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import re
import sys
import numpy
import os

checker = numpy.zeros((5), dtype=bool)
bundle = [None] * 5
#print checker

# Check to see if TRADVI diagnostics are asked for in HISOTRY.rc
# --------------------------------------------------------------
f = open('HISTORY.rc', 'rt')
for line in f:
    line = line.strip()
    if not line.startswith('#'):
#    print line
        if "'TRADVI%" in line:
            checker[0] = True
            bundle[0] = "'TRADVI"
        if "'CHEMTRI%" in line:
            checker[1] = True
            bundle[1] = "'CHEMTRI"
        if "'H2ORTRI%" in line:
            checker[2] = True
            bundle[2] = "'H2ORTRI"
        if "'TRI%" in line:
            checker[3] = True
            bundle[3] = "'TRI"
        if "'MTRI%" in line:
            checker[4] = True
            bundle[4] = "'MTRI"

f.close()


# If no increment species are asked for in HISTORY.rc then exit
#--------------------------------------------------------------
if not any(checker):
    sys.exit()


# Get species names from HISTORY.rc
#----------------------------------
for i in range(len(checker)):
    if checker[i]:
        tag = []
        f = open('HISTORY.rc', 'rt')
        for line in f:
            line = line.strip()
#            print line
            if not line.startswith('#'):
                specName = re.findall(bundle[i]+'%(.*?)\'', line)
                if specName:
                    if specName not in tag:
                        tag.append(specName)
        f.close()

# Write table of species increment names to AGCM.rc or GEOSCTM.rc
#----------------------------------------------------------------  
        exists = os.path.isfile('AGCM.rc')
        if exists:
            with open('AGCM.rc','a') as myfile:
                myfile.write('\n' '\n' '# ' + bundle[i][1:] + ' increment tracers' '\n' 
                             '#-------------------------' '\n' + bundle[i][1:] + '_increments::' 
                             '\n')
                for item in tag:
                    myfile.write(str(item)[2:-4] + '\n')
                myfile.write('::' '\n')

        exists = os.path.isfile('GEOSCTM.rc')
        if exists:
            with open('GEOSCTM.rc','a') as myfile:
                myfile.write('\n' '\n' '# ' + bundle[i][1:] + ' increment tracers' '\n'
                             '#-------------------------' '\n' + bundle[i][1:] + '_increments::'
                             '\n')
                for item in tag:
                    myfile.write(str(item)[2:-4] + '\n')
                myfile.write('::' '\n')


