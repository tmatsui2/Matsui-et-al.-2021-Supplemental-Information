###all codes here done in terminal using bash, this is not a python code
###to run this code, first download and install bartender; see (https://github.com/LaoZZZZZ/bartender-1.1)
cd /Volumes/TM/All_segregants/fitness_assay_barcodes/Zeocin/bartender_files ###path to output files from 2_fastq_parsing
for files in *.csv ###for running bartender of many files
do
    bartender_single_com -f $files -o $files -d 3 ###cluster all barcode reads that are within a hamming distance of 3 (-d parameter) or less. Barcodes with the same UMIs are omitted as they are most likely PCR duplicates
done


###for running bartender on one file
bartender_single_com -f zeo_T3.csv -o zeo_T3.csv -d 3