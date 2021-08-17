###written for python2
###This code reads in the raw read1 fastq file from a 150bp PE HiSeq/NextSeq sequence run and extracts out the double barcodes as well as the UMI/multiplex tag. See methods section for details on the library structure.
import sys
import os
import time
import multiprocessing as mp
    
def process_wrapper(chunkStart, chunkSize):
    with open(inf) as f:
        f.seek(chunkStart) ##reads in chunk of lines of defined size
        chunks = f.read(chunkSize).split('@') ##split chunk of lines by @ into individual fastq reads
        
        for chunk in chunks:
            if not chunk: continue
            
            split_read = chunk.split('\n') ##split each fastq reads by new line
        
            if len(split_read) == 5:
                fixed = split_read[1]  ###line which contains the sequence of the insert region
                quality = split_read[3] ###line which contains the ASCII qulaity score of the fastq read

                if substring not in fixed: continue

                if fixed.find(substring):  ###only consider fastq reads which contains the defined fixed sequence, used to separate out PhiX reads
                    temp = []
                    for char in quality:
                        scores = temp.append(ord(char)-33) ##convert ASCII quality code to a numeric Q-score
                        average = sum(temp)/len(temp) ##finds average quality score
    
                    if average < 25: continue
                    if average >= 25: ###only consider fastq reads with a average score above 25
                        umi1 = split_read[0].split(':')[-1].split('+')[0] ###for extracting i5 3bp Unique molecule identifier (UMI) and 5bp multiplex index
                        umi2 = split_read[0].split(':')[-1].split('+')[1] ###for extracting i7 3bp Unique molecule identifier (UMI) and 5bp multiplex index
                        bc1 = split_read[1][0:20] ##for extracting 20bp MATa barcode
                        bc2 = split_read[1][108:128] ##for extracting 20bp MATalpha barcode
        
                        data = '{0}{1},{2}{3}{4}{5}\n'.format(bc1[7:], bc2[7:], umi1, umi2, bc1[:7], bc2[:7]) ###creates output file with two columsn separated by a comma. First column contains the last 13bps of the MATa and MATalpha barcodes, which will be used to assign the reads to a diploid strain. The second columns contains the UMI, multiplex tag, and the first 7bps of the MATa and MATalpha barcodes. The first 7bps of the barcodes are used to pair the UMIs to each double barcode pair.
                        outf.write(data)
                        
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

####
if __name__ == '__main__':
    path = '/Volumes/TM/All_segregants/fitness_assay_barcodes/' ###path to input file
    inf = path + 'CRC_2_CKDL200166600-1a_HLLTTDSXY_L4_1.fq' ###name of the Read1 input file
    outf = open(path + 'crc2_barcode_extraction.csv', 'w') ###name of the output file
    substring = 'ATAACTTCGTATAATGTATGCTAT' ###fixed sequence present in the double barcode region = used to parse out the double barcode sequence
    start_time = time.time()  

    pool = mp.Pool(mp.cpu_count() - 1) ##run using parallel prcoessing

    for chunkStart,chunkSize in chunkify(inf):
        results = [pool.apply_async(process_wrapper, args = (chunkStart, chunkSize,))]
        output = [result.get() for result in results]
    
    print("--- %s seconds ---" % (time.time() - start_time))
    pool.close()
    outf.close()