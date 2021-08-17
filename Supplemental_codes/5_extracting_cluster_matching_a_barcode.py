####Takes the MATalpha assigend clustered barcodes from 4_extracting_cluster_matching_alpha_barcode and assign them to a parental MATa strain. The output file has 4 columns separated by commas: MATa strain, MATalpha strain, double barcode, and the clustered barcode count
import glob
import os
import sys
from Levenshtein import distance 
import time

path = '/Volumes/TM/All_segregants/fitness_assay_barcodes/Glucose/Barcodes/' ###path to where the parental MATa strain barcodes are kept
path2 = '/Volumes/TM/All_segregants/fitness_assay_barcodes/Rapamycin/' ###path to where the output files from 4_extracting_cluster_matching_alpha_barcode are kept
barcodes = path + 'only_MATa_barcodes.txt' ###name of th parental MATa barcode file
start_time = time.time()

def match(s1, s2):
    return distance(s1, s2) <= 1
    
for filename in glob.iglob(path2 + 'alpha_matching/*.csv'):
    filename = os.path.basename(filename)
    inf = path2 + 'alpha_matching/' + filename
    outf = open(path2 + 'both_matching/' + filename.split('_')[0] + '_' + filename.split('_')[1] + '_both_matching_cluster.csv', 'w') ###output file name
    
    with open(inf) as f:
        with open(barcodes) as g:
            checks = g.readlines()
            bc_list = [i.split('\t')[1].split('\n')[0][7:] for i in checks] ###extract out the last 13bps of each parental MATa strain barcode
        
        for line in enumerate(f):
            split = line[1].split(',')
            bc = split[1]
             
            if(bc[13:] in bc_list): ##if a barcode is an exact match with one of the parental MATa strain barcode
                for check in checks:
                    bc_check = check.split('\t')
                    if(bc[0:13] in bc_check[1]):
                        outf.write('{0},{1},{2},{3}'.format(bc_check[0], split[0], bc, split[2])) ###write out the MATalpha strain, MATa strain, double barcode, and clustered barcode count number
            else: ##if a barcode is not an exact match with one of the parental MATa strain barcode
                for check in checks:
                    bc_check = check.split('\t')
                    if(match(bc[0:13], bc_check[1][7:])):  ###check if the barcodes have a levenshtein distance of 1 or less with one of the MATa barcode
                        outf.write('{0},{1},{2},{3}'.format(bc_check[0], split[0], bc, split[2])) ##if yes, write
                        
print("--- %s seconds ---" % (time.time() - start_time))                        
outf.close()