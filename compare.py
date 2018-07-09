"""Usage: test.py (-h --help)  ...
       test.py --compare (<vcf_paths>) ...
       
Options:
  -h --help                            | help
  -compare vcf_impl_path vcf_ref_path  | execute comparison with entered path's 
                                       | to implemented and referent vcf files 
                                       | in that order
"""

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
from enum import Enum
import math
import sys
from docopt import docopt

################# CONSTANTS #################
class Constants(object):
    kVCF_RESULT = '\\output\compare_result.txt'
    
############# GLOBAL VARIABLES ##############

class PileupLine(Enum): 
    SEQID = 0
    POSINSEQ = 1
    REFNUCLEOTID = 2
    NUMOFREADS = 3
    BASES = 4
    QUALITYOFBASES = 5
    
def compareFunction(vcf_path, vcf_samtools_path):
    result_file = open(dir_path + Constants.kVCF_RESULT, "a+")
    result_file.seek(0)
    result_file.truncate()
    result_file.seek(0)
    
    if not os.path.exists(dir_path + vcf_path):
        print("Path of referent vcf file is not valid: " + dir_path + vcf_path)
        exit()
        
    with open(dir_path + vcf_path) as fVCF_mine:
        next(fVCF_mine)
        next(fVCF_mine)
        content_mine = fVCF_mine.readlines()
        content_mine_count = len(content_mine)
        content_mine = [x.strip() for x in content_mine]
        content_mine = [x.split("\t") for x in content_mine]
        
    if not os.path.exists(dir_path + vcf_samtools_path):
        print("Path of implemented vcf file is not valid: " + dir_path + vcf_samtools_path)
        exit()
        
    with open(dir_path + vcf_samtools_path) as fVCF_samtools:
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
            if x[PileupLine.POSINSEQ.value] == y[PileupLine.POSINSEQ.value]:
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
        
# Run the program
if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    vcf_path = arguments['<vcf_paths>'][0]
    vcf_samtools_path = arguments['<vcf_paths>'][1]
    print(vcf_path)
    print(vcf_samtools_path)
    compareFunction(vcf_path, vcf_samtools_path)
    #print(arguments['<ref>'])
    #print(arguments['<impl>'])