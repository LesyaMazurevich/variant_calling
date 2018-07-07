import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
from enum import Enum
from operator import itemgetter
import math
import sys
import msvcrt

################# CONSTANTS #################
class Constants(object):
    kVCF_PATH = '\VCF_impl.vcf'
    kVCF_HEADER = '##fileformat=VCFv4.2\n'
    kVCF_COLUMNS = '#CHROM\tPOS\tREF\tALT\tFORMAT\tSAMPLE_ID\n'
    kPILEUP_FILE_PATH = '\pileup.txt'

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
    print("Enter appropriate value to continue:")
    print("1 to execute variant calling algorithm")
    print("2 to execute variantCallingCore function test")
    print("3 to execute variantCalling function test")
    print("4 to execute addVarinatToIntVcArray function test")
    m = msvcrt.getch().decode('utf-8')
    print(m)
    global VCFfile
    if m == '1':
        print("Executing variant calling algorithm...")
        print("Initialising, please wait...")
        VCFfile = open(dir_path + Constants.kVCF_PATH, "a+")
        VCFfile.seek(0)
        VCFfile.truncate()
        VCFfile.seek(0)
        VCFfile.write(Constants.kVCF_HEADER)
        VCFfile.write(Constants.kVCF_COLUMNS)
        with open(dir_path + Constants.kPILEUP_FILE_PATH) as f:
            content = f.readlines()
            content = [x.strip() for x in content]
            content = [x.split("\t") for x in content]

            print("Started variant calling...")
            variantCalling()
            
            del content 
            
    elif m == '2':
        print("Executing variantCallingCore function test...")
        VCFfile = open(dir_path + "\VCF_test.vcf", "a+")
        VCFfile.seek(0)
        VCFfile.truncate()
        VCFfile.seek(0)
        content = [["1", "336", "C", "34", "A$A$A$A$AAAA,.,AA,,..AAAAaaaa...CcTtGg", "test"], ["1", "400", "A", "12", "...,,+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC,.,...,", "test"], ["1", "500", "A", "12", "...,,-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC,.,...,", "test"], ["1", "501", "A", "1", ".", "test"], ["1", "502", "C", "1", ".", "test"]]
        variantCalling()
        VCFfile.seek(0)
        #print(VCFfile.read())
        if VCFfile.read() == "1	336	C	A	GT	0/1\n1	400	A	AAC	GT	1/1\n1	500	AAC	A	GT	1/1\n":
            print("Test run succesufully...")
        else:
            print("Test fail...")
    elif m == '3':
        print("Executing variantCalling function test...")
        incVcArray = variantCallingCore(["1", "336", "C", "34", "A$A$A$A$AAAA,.,.,,,..AAAAaaaa...CcTtGg+2AC", "test"])
        sum = 0
        for x in incVcArray:
            sum += incVcArray[0][0]
        if sum != incVcArray[3][0]:
            print("Test fail...")
        if incVcArray[0][0] == 12 and incVcArray[1][0] == 2 and incVcArray[2][0] == 2 and incVcArray[3][0] == 16 and incVcArray[4][0] == 2 and incVcArray[5][0] == 1:
            print("Test run succesufully...")
        else:
            print("Test fail...")
    elif m == '4':
        print("Executing addVarinatToIntVcArray function test...")
        p = ['.', 'T', 'G', 'A', 'C']
        incVcArray = [[0, p[i]] for i in range(len(p))]
        addVarinatToIntVcArray(incVcArray, "+2AC")
        if incVcArray[5][1] == "+2AC":
            print("Insertion added succesufully...")
            addVarinatToIntVcArray(incVcArray, "-2AT")
            if incVcArray[6][1] == "-2AT":
                print("Deletion added succesufully...")
                print("Function addVarinatToIntVcArray executed succesufully...")
            else:
                print("Test failed due to deletion adding!")
        else:
            print("Test failed due to insertion adding!")
    else:
        print("Entered value is not allowed...")