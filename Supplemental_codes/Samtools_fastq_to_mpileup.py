###Written for python 2.7. Takes whole-genome sequence reads (fastq format) and align them to the S288c reference genome using BWA. Aligned data is then used to create SAM and BAM files using SAMTOOLS. Finally, the sorted BAM files are converted into SAMTOOLS mpileup files for downstream analysis.
######################################################################################################################################################
import sys
import os
from subprocess import Popen, PIPE

def getBarcodes(path):
    barcodes = []
    for file in os.listdir(path):
	print file
        if ".fq" in file:
            if file[2:7] not in barcodes:
                barcodes.append(file[2:7])		
    return barcodes



if __name__ == '__main__':
	genome = "genome"
	path = "/Users/tmatsui2/Desktop/re-wgs/split/"
	os.chdir(path)
	
	names = getBarcodes(path)
	
	for name in names:
		if name + ".sam" not in os.listdir(path):
			step1 = "./bwa mem -t 16 {0}.fasta 1_{1}.fq 2_{1}.fq > {1}.sam".format(genome, name)
			step2 = "./samtools view -bS -o " + name + ".bam " + name + ".sam"
			step3 = "./samtools sort " + name + ".bam " + name + "_sorted"
			step4 = "./samtools rmdup -S " + name + "_sorted.bam " + name + "_rsorted.bam "
			step5 = "./samtools mpileup -f " + genome + ".fasta " + name + "_rsorted.bam > " + name + "_mpileup.txt"
			
			commands = [step1, step2, step3, step4, step5]
		
			for command in commands:
				Popen(command, stdout=PIPE, shell=1).wait()

######################################################################################################################################################
###Takes the mpileups files and extracts out genome positions in which a SNP between BY and 3S was previously identified. For each SNP position, the chromosome number, position
#number, number of reads that aligned to that SNP (coverage), and the counts of each allele was recorded (number of BY alleles vs number of 3S alleles).   

if __name__ == '__main__':
	
	path = "/Users/lab/Desktop/Chromatin_Modifier_Project/WT/WT_sequence/mpileup/"
	bases = ["A", "C", "G", "T"]
	
	snpdic = {}
	infile1 = open("/Users/lab/Desktop/Parental_flo11/Parental_mpileup/BY_snps.txt", 'r')
	infile2 = open("/Users/lab/Desktop/Parental_flo11/Parental_mpileup/3S_snps.txt", 'r')
	
	for line in infile1:
		lin = line.split()
		snpdic[lin[0] + "_" + lin[1]] = [lin[3], lin[2]]
	infile1.close()
    
	print "parent 1 accounted for"
   
	for line in infile2:
		lin = line.split()
		if lin[0] + '_' + lin[1] not in snpdic.keys():
			snpdic[lin[0] + "_" + lin[1]] = [lin[2], lin[3]]
		else:
			snpdic[lin[0] + "_" + lin[1]][1] = lin[3]
	infile2.close()
	snp_set = set(snpdic.keys())
    
	print "parent 2 accounted for"
    
	
	seglist = []
	for file in os.listdir(path):
		if ".mp" in file:
			if file.split('.')[0] not in seglist:
				seglist.append(file.split('.')[0]) 
	
	for seg in seglist:
		if str(seg + "_segsnps.txt") not in os.listdir(path):
			outfile = open(path +  "segsnps/" + seg + "_segsnps.txt", 'w')
			infile = open(path + seg + ".mp", 'r')
			outdata = []
			for line in infile:
				lin = line.split("\t")
				if lin[0] == "mitochondria": lin[0] = "17"
				base = lin[0] + "_" + lin[1]
				if base in snp_set:
					read = lin[4].upper()
					read = read.replace(",", lin[2]).replace(".", lin[2])
					p1count = read.count(snpdic[base][0])
					p2count = read.count(snpdic[base][1])
					outdata.append(lin[0] + "\t" + lin[1] + "\t" + lin[3] + "\t" + snpdic[base][0] + "\t" + str(p1count) + "\t" + snpdic[base][1] + "\t" + str(p2count))
			outfile.write('\n'.join(outdata))
			infile.close()
			outfile.close()
 
 

