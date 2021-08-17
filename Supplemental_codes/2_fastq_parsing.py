###written for python2
###This codes takes the extracted fastq reads from 1_extract_barcode_and_umi_from_fq.py and parses them based on the multiplex tag. Here, the multiplex tags are used to label the envrionmental condition and time point in which the library was generated from.
import os
import sys
from Levenshtein import distance 
import gc

path = '/Volumes/TM/All_segregants/fitness_assay_barcodes/' ##path to the output file from 1_extract_barcode_and_umi_from_fq.py
inf = path + 'crc2_barcode_extraction.csv' ###output file name from 1_extract_barcode_and_umi_from_fq.py

outf1 = open(path + './CuSO4/cuso4_T3.csv', 'w') ###name of the output files
outf2 = open(path + './CuSO4/cuso4_T6.csv', 'w')
outf3 = open(path + './CuSO4/cuso4_T9.csv', 'w')
outf4 = open(path + './CuSO4/cuso4_T12.csv', 'w')

outf5 = open(path + './Rapamycin/rapa_T3.csv', 'w')
outf6 = open(path + './Rapamycin/rapa_T6.csv', 'w')
outf7 = open(path + './Rapamycin/rapa_T9.csv', 'w')
outf8 = open(path + './Rapamycin/rapa_T12.csv', 'w')

outf9 = open(path + './CoCl2/cocl2_T3.csv', 'w')
outf10 = open(path + './CoCl2/cocl2_T6.csv', 'w')
outf11 = open(path + './CoCl2/cocl2_T9.csv', 'w')
outf12 = open(path + './CoCl2/cocl2_T12.csv', 'w')

def process_wrapper(chunkStart, chunkSize):
    with open(inf) as f:
        f.seek(chunkStart) ##read in chunks of umi/barcodes from input file
        chunks = f.read(chunkSize).split('\n') ##split by new line to read in each umi/barcode individually
        for chunk in chunks:
            gc.disable()
            split_read = chunk.split(',') ##split each umi/barcode by comma
            if len(split_read) == 2:
                if (len(split_read[0]) == 26 and len(split_read[1]) == 30):
                    umi = split_read[1]
                    tag = '{0}{1}'.format(umi[2:7], umi[10:15]) ###extract out only the mutiplex tag for each umi/barcode
                    
                    check1 = 'CCCTGCCCTG' ###list of the known multiplex tags
                    check2 = 'GAGAGGAGAG'
                    check3 = 'ATGATATGAT'
                    check4 = 'CTTGGCTTGG'
                    
                    check5 = 'ACACTACACT'
                    check6 = 'AGTAAAGTAA'
                    check7 = 'CGCGCCGCGC'
                    check8 = 'ACCGAACCGA'
                    
                    check9 = 'CATATCATAT'
                    check10 = 'GATGAGATGA'
                    check11 = 'TCTTATCTTA'
                    check12 = 'TAGTCTAGTC'
                    
                    if match(tag, check1):  ###if the multiplex tag from the fastq read has a levenshtein distance of 1 or less with a known multiplex tag, the umi/barcode is written in a specific output file
                        outf1.write(chunk + '\n')
                    elif match(tag, check2):
                        outf2.write(chunk + '\n')
                    elif match(tag, check3):
                        outf3.write(chunk + '\n')
                    elif match(tag, check4):
                        outf4.write(chunk + '\n')
                    
                    elif match(tag, check5):
                        outf5.write(chunk + '\n')    
                    elif match(tag, check6):
                        outf6.write(chunk + '\n')
                    elif match(tag, check7):
                        outf7.write(chunk + '\n')
                    elif match(tag, check8):
                        outf8.write(chunk + '\n')
                    
                    elif match(tag, check9):
                        outf9.write(chunk + '\n')    
                    elif match(tag, check10):
                        outf10.write(chunk + '\n')
                    elif match(tag, check11):
                        outf11.write(chunk + '\n')
                    elif match(tag, check12):
                        outf12.write(chunk + '\n')
                        
                    else: continue
            gc.enable()        
                        

def match(s1, s2):
    return distance(s1, s2) <= 1 ###calculates the levenshtein distance between a known multiplex tag and the multiplex tag from the fastq reads
    
def chunkify(fname, size = 2**25): ###read in chunks of data, instead of reading by line
    fileEnd = os.path.getsize(fname)
    with open(fname, 'r') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break
                
for chunkStart,chunkSize in chunkify(inf):
    process_wrapper(chunkStart, chunkSize)
    
outf1.close()
outf2.close()
outf3.close()
outf4.close()
outf5.close()
outf6.close()
outf7.close()
outf8.close()
outf9.close()
outf10.close()
outf11.close()
outf12.close()
