"""Usage: variant_calling.py (-h --help)  ...
       variant_calling.py (-v <pileup_path>)  ...
        
Options:
  -h --help        | help
  -v pileup_path   | execute variant calling on specified pileup file
"""
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
from enum import Enum
from operator import itemgetter
import math
import sys
import msvcrt
from docopt import docopt

################# CONSTANTS #################
class Constants(object):
    kVCF_PATH = '\\output\VCF_impl.vcf'
    kVCF_HEADER = '##fileformat=VCFv4.2\n'
    kVCF_COLUMNS = '#CHROM\tPOS\tREF\tALT\tFORMAT\tSAMPLE_ID\n'
    #kPILEUP_FILE_PATH = '\\resources\pileup.txt'

############# GLOBAL VARIABLES ##############

class PileupLine(Enum): 
    SEQID = 0
    POSINSEQ = 1
    REFNUCLEOTID = 2
    NUMOFREADS = 3
    BASES = 4
    QUALITYOFBASES = 5
    
def addVarinatToIntVcArray(incVcArrray, insertion):
    #print("addVarinatToIntVcArray")
    insertion = insertion.upper()
    for x in incVcArrray:
        if x[1] == insertion:
            x[0] += 1
            return
    inerList = [1, insertion]
    incVcArrray.append(inerList)
    
def variantCallingCore(x):
    #print("variantCallingCore")
    flag = 0
    insertion = ""
    deletion = ""
    p = ['.', 'T', 'G', 'A', 'C']
    incVcArray = [[0, p[i]] for i in range(len(p))]
    for char in x[PileupLine.BASES.value]:
        if flag == 0:
            if char == '.' or char == ',':
                incVcArray[0][0] += 1
            if char == 'T' or char == 't':
                incVcArray[1][0] += 1
            if char == 'G' or char == 'g':
                incVcArray[2][0] += 1
            if char == 'A' or char == 'a':
                incVcArray[3][0] += 1
            if char == 'C' or char == 'c':
                incVcArray[4][0] += 1
            if char == '+':
                flag = 2
                insertion = insertion + char
            if char == '-':
                flag = 2
                insertion = insertion + char
            if char == '$':
                flag = 0
            if char == '^':
                flag = 1
        elif flag == 1:
            flag = 0
        elif flag == 2:
            count = int(char)
            insertion = insertion + char
            flag = 3
        elif flag == 3:
            count = count - 1
            insertion = insertion + char
            if count == 0:
                flag = 0
                addVarinatToIntVcArray(incVcArray, insertion)
                insertion = ""
    return incVcArray
    
def writeVcfLine(chrom,	pos, ref, alt, format, sample_id):
    #print("writeVcfLine")
    vcfLine = chrom + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + format + "\t" + sample_id + "\n"
    VCFfile.write(vcfLine)
    
def variantCalling():
    #print("variantCalling")
    j = 0
    contenctLen = len(content)
    for x in content:
        sys.stdout.write(str(math.ceil(j*100/contenctLen)) + "% complete         \r")
        j += 1
        incVcArray = variantCallingCore(x)
        incVcArray = sorted(incVcArray, key=itemgetter(0), reverse=True)
        ALT = ""
        GT = "0/1"
        ALTflag = False
        deletionFlag = False
        GT12flag = False
        i = 0
        while i < len(incVcArray) and int(incVcArray[i][0]) > 0:
            if int(incVcArray[i][0]) > int(x[PileupLine.NUMOFREADS.value])*0.5:
                if int(incVcArray[i][0]) > int(x[PileupLine.NUMOFREADS.value])*0.8:
                    GT = "1/1"
                if incVcArray[i][1][0] == '.':
                    i += 1
                    continue
                elif incVcArray[i][1][0] == 'T' or incVcArray[i][1][0] == 'G' or incVcArray[i][1][0] == 'A' or incVcArray[i][1][0] == 'C':
                    if ALTflag:
                        ALT += "," + incVcArray[i][1]
                        GT12flag = True
                    else:
                        ALT += incVcArray[i][1]
                    ALTflag = True
                elif incVcArray[i][1][0] == '+':
                    if ALTflag:
                        ALT += "," + x[PileupLine.REFNUCLEOTID.value] + incVcArray[i][1][2:]
                        GT12flag = True
                    else:
                        ALT += x[PileupLine.REFNUCLEOTID.value] + incVcArray[i][1][2:]
                    ALTflag = True
                elif incVcArray[i][1][0] == '-':
                    if ALTflag:
                        ALT += "," + x[PileupLine.REFNUCLEOTID.value] + incVcArray[i][1][2:]
                        GT12flag = True
                    else:
                        ALT += x[PileupLine.REFNUCLEOTID.value] + incVcArray[i][1][2:]
                    ALTflag = True
                    deletionFlag = True
            i += 1
        if GT12flag:
            GT = "1/2"
        if ALT != "":
            if deletionFlag:
                writeVcfLine(x[PileupLine.SEQID.value], x[PileupLine.POSINSEQ.value], x[PileupLine.REFNUCLEOTID.value] + ALT[1:], ALT[0], "GT", GT)
            else:
                writeVcfLine(x[PileupLine.SEQID.value], x[PileupLine.POSINSEQ.value], x[PileupLine.REFNUCLEOTID.value], ALT, "GT", GT)
        
# Run the program
if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    if arguments['-v']:
        print("Executing variant calling algorithm...")
        print("Initialising, please wait...")
        VCFfile = open(dir_path + Constants.kVCF_PATH, "a+")
        VCFfile.seek(0)
        VCFfile.truncate()
        VCFfile.seek(0)
        VCFfile.write(Constants.kVCF_HEADER)
        VCFfile.write(Constants.kVCF_COLUMNS)
        if not os.path.exists(dir_path + arguments['-v'][0]):
            print("Path of pileup file is not valid: " + dir_path + arguments['-v'][0])
            exit()
        with open(dir_path + arguments['-v'][0]) as f:
            content = f.readlines()
            content = [x.strip() for x in content]
            content = [x.split("\t") for x in content]

            print("Started variant calling...")
            variantCalling()
            
            del content 
            
            
            
    