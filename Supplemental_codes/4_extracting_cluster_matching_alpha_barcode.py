####Takes the clustered barcodes from 3_bartender_cluster and assign them to a parental MATalpha strain
import glob
import os
import sys
from Levenshtein import distance 
import time

path = '/Volumes/TM/All_segregants/fitness_assay_barcodes/Glucose/Barcodes/' ###path to where the parental MATalpha strain barcodes are kept
path2 = '/Volumes/TM/All_segregants/fitness_assay_barcodes/Rapamycin/' ###path to where the output files from 3_bartender_cluster are kept
barcodes = path + 'only_MATalpha_barcodes.txt' ###name of th parental MATalpha barcode file
start_time = time.time()

def match(s1, s2):
    return distance(s1, s2) <= 1 ###calculates levenshtein distance between a barcode and a parental MATalpha strain barcode
    
for filename in glob.iglob(path2 + 'unmatched_cluster/*.csv'):
    filename = os.path.basename(filename)
    inf = path2 + 'unmatched_cluster/' + filename
    outf = open(path2 + 'alpha_matching/' + filename.split('.')[0] + '_alpha_matching_cluster.csv', 'w')  ###output file name
    
    with open(inf) as f:
        next(f)
        with open(barcodes) as g:
            checks = g.readlines()
            bc_list = [i.split('\t')[1].split('\n')[0][7:] for i in checks] ###extract out the last 13bps of each parental MATalpha strain barcode
        
        for line in enumerate(f):
            split = line[1].split(',')
            bc = split[1]
             
            if(bc[13:] in bc_list): ##if a barcode is an exact match with one of the parental MATalpha strain barcode
                for check in checks:
                    bc_check = check.split('\t')
                    if(bc[13:] in bc_check[1]):
                        outf.write('{0},{1},{2}'.format(bc_check[0], bc, split[3])) ###write out the MATalpha strain, barcode, and clustered barcode count number
            else: ##if a barcode is not an exact match with one of the parental MATalpha strain barcode
                for check in checks:
                    bc_check = check.split('\t')
                    if(match(bc[13:], bc_check[1][7:])): ###check if the barcodes have a levenshtein distance of 1 or less with one of the MATalpha barcode
                        outf.write('{0},{1},{2}'.format(bc_check[0], bc, split[3])) ##if yes, write
                        
print("--- %s seconds ---" % (time.time() - start_time))                        
outf.close()