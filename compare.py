import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
from enum import Enum
import math
import sys

################# CONSTANTS #################
class Constants(object):
    kVCF_PATH = '\VCF_impl.vcf'
    kVCF_SAMTOOLS_PATH = '\VCF_ref.vcf'
    kVCF_RESULT = '\compare_result.txt'
    
############# GLOBAL VARIABLES ##############

class PileupLine(Enum): 
    SEQID = 0
    POSINSEQ = 1
    REFNUCLEOTID = 2
    NUMOFREADS = 3
    BASES = 4
    QUALITYOFBASES = 5
    
# Run the program
if __name__ == '__main__':
    result_file = open(dir_path + Constants.kVCF_RESULT, "a+")
    result_file.seek(0)
    result_file.truncate()
    result_file.seek(0)
    
    with open(dir_path + Constants.kVCF_PATH) as fVCF_mine:
        next(fVCF_mine)
        next(fVCF_mine)
        content_mine = fVCF_mine.readlines()
        content_mine_count = len(content_mine)
        content_mine = [x.strip() for x in content_mine]
        content_mine = [x.split("\t") for x in content_mine]
        
    with open(dir_path + Constants.kVCF_SAMTOOLS_PATH) as fVCF_samtools:
        i = 0
        while i < 117:
            next(fVCF_samtools)
            i += 1
        content_samtools = fVCF_samtools.readlines()
        content_samtools_count = len(content_samtools)
        content_samtools = [x.strip() for x in content_samtools]
        content_samtools = [x.split("\t") for x in content_samtools]
        
    truePositives = 0
    falsePositives = 0
    falseNegatives = 0
    precision = 0
    recall = 0
    f_score = 0
    truePositivesFlag = False
        
    p = 0    
    contentLen = len(content_samtools)
    for x in content_samtools:
        sys.stdout.write(str(math.ceil(p*100/contentLen)) + "% complete         \r")
        p += 1
        truePositivesFlag = False
        for y in content_mine:
            if x[PileupLine.POSINSEQ.value] == y [PileupLine.POSINSEQ.value]:
                truePositives += 1
                truePositivesFlag = True
                break
        if not truePositivesFlag:
            falsePositives += 1
            
    falseNegatives = len(content_mine) - truePositives
    
    precision = truePositives / (truePositives + falsePositives)
    
    recall = truePositives / (truePositives + falseNegatives)
    
    f_score = 2 * precision * recall / (precision + recall)

    result_file.write("Samtools vcf count of variants: " + str(content_samtools_count) + "\n")
    result_file.write("Implemented vcf count of variants: " + str(content_mine_count) + "\n")
    result_file.write("True positives: " + str(truePositives) + "\n")
    result_file.write("False positives: " + str(falsePositives) + "\n")
    result_file.write("False negatives: " + str(falseNegatives) + "\n")
    result_file.write("Precision: " + str(precision) + "\n")
    result_file.write("Recall: " + str(recall) + "\n")
    result_file.write("F-score: " + str(f_score) + "\n")
        
    del content_mine 
    del content_samtools
        