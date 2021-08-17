###This code was used to a genome-wide scan for one-locus effects using FaST-LMM. 
###Please note that different parts of the code use different programming softwares.
####################################################################################################################################
###Run using R
#Reformatting genotype and phenotype table for plink format
load('/oak/stanford/groups/msalit/tmatsui2/Glucose/Glucose_geno.RData')
load('/oak/stanford/groups/msalit/tmatsui2/Glucose/Glucose_pheno.RData')
geno = as.matrix(geno)
geno = geno[,4:ncol(geno)]
geno = t(geno)
fid = 1:nrow(geno)
iid = fid
pid = rep(0, nrow(geno))
mid = pid
sex = pid

###
pheno_table = cbind(fid, iid, pheno$qnorm)
write.table(pheno_table, file = '/oak/stanford/groups/msalit/tmatsui2/Glucose/pheno.txt', quote = F, row.names = F, col.names = F, sep = '\t')

geno = cbind(fid, iid, pid, mid, sex, pheno$qnorm, geno)
write.table(geno, file = '/oak/stanford/groups/msalit/tmatsui2/Glucose/geno_plink.txt', quote = F, row.names = F, col.names = F, sep = '\t')

####################################################################################################################################
###Run using python
#Convert genotype table (each SNP position is denoted by the number of 3S alleles, eg, 0, 1, or 2) into biallelic genotype table (each SNP position
#is denoted as either AA, AC, or CC, corresponding the 0, 1, and 2 genotypes, respectively)
path = '/oak/stanford/groups/msalit/tmatsui2/Glucose/'
inf = path + 'geno_plink.txt'
outf = open(path + 'glu.ped', 'w')

with open(inf) as f:
    chunks = f.readlines()
    
    for chunk in chunks:
        if not chunk: continue
        data = chunk.split('\n')[0].split('\t')
        data1 = data[0:6]
        data2 = data[6:len(data)]
        
        for n, i in enumerate(data2):
            if i == '0':
                data2[n] = 'AA'
            elif i == '1':
                data2[n] = 'AC'
            elif i == '2':
                data2[n] = 'CC'
        data = '\t'.join(data1) + '\t' + '\t'.join(data2) + '\n'
        outf.write(data)

outf.close()

####################################################################################################################################
###Run in terminal with plink v1.07
#convert biallelic genotype table into plink binary biallelic genotype table (.bed format)
plink --file glu --compound-genotypes --make-bed --out glu --noweb

####################################################################################################################################
###Run in terminal or HPC cluster
#run FaST-LMM to scan for one-locus effects, first scan
import sys
import numpy as np
sys.path.insert(0,'/home/users/tmatsui2/anaconda3/lib/python3.8/site-packages/')
import logging
from fastlmm.association import single_snp
from pysnptools.snpreader import Bed

path = ('/oak/stanford/groups/msalit/tmatsui2/Glucose/')
logging.basicConfig(level = logging.INFO)
test_snps = Bed(path + 'glu', count_A1 = True)
pheno_fn = path + 'glu_pheno.txt'
out_fn = path + 'fastlmm/output_fastlmm/fast_lmm_glu1.txt'

results_dataframe = single_snp(test_snps = test_snps, pheno = pheno_fn, count_A1=True, output_file_name = out_fn)

####################################################################################################################################
###Run using R
#take results from the first FaST-LMM scan and identify the most significant SNP above the significance threshold from each chromosome
#make a separate file which contains the genotypes at the detected loci to be included as covariates in the second FaST-LMM scan
load('/oak/stanford/groups/msalit/tmatsui2/Glucose/Glucose_geno.RData')
geno = as.matrix(geno)
pheno = cbind(1:(ncol(geno) - 3), 1:(ncol(geno) - 3))
pval = read.table('/oak/stanford/groups/msalit/tmatsui2/Glucose/fastlmm/output_fastlmm/fast_lmm_glu1.txt', header = T, as.is = T)
pval$LOD = -log10(pval$PValue)

listy = list()
for(x in 1:16){
    temp = pval[which(pval$Chr == x),]
    if(max(temp$LOD) > 4.3){
        df = temp[which(abs(temp$SnpWeight) == max(abs(temp$SnpWeight))),]
        if(nrow(df) == 1){
            listy[[x]] = df
        } else if(nrow(df) > 1){
            listy[[x]] = df[1,]
        }
    }
}
tab = do.call(rbind, listy)
tab$row = paste(tab$Chr, tab$ChrPos, sep = '_')
tab$snp = as.numeric(unlist(strsplit(tab$SNP, split = 'rs'))[seq(2, nrow(tab) * 2, 2)])
covar = lapply(1:nrow(tab), function(x){
    geno[tab$snp[x], 4:ncol(geno)]
})
covar2 = do.call(cbind, covar)
covar3 = cbind(pheno, covar2)
rownames(covar3) = NULL
write.table(covar3, file = '/oak/stanford/groups/msalit/tmatsui2/Glucose/fastlmm/code/glu_cov_fr1.txt', quote = F, col.names = F, row.names = F, sep = '\t')

####################################################################################################################################
###Run in terminal or HPC cluster
#second scan of FaST-LMM with loci detected from the first scan included as covariates
import sys
import numpy as np
sys.path.insert(0,'/home/users/tmatsui2/anaconda3/lib/python3.8/site-packages/')
import logging
from fastlmm.association import single_snp
from pysnptools.snpreader import Bed

path = ('/oak/stanford/groups/msalit/tmatsui2/Glucose/')
logging.basicConfig(level = logging.INFO)
test_snps = Bed(path + 'glu', count_A1 = True)
pheno_fn = path + 'glu_pheno.txt'
out_fn = path + 'fastlmm/output_fastlmm/fast_lmm_glu2.txt'
cov_fn = path + 'fastlmm/code/glu_cov_fr1.txt'

results_dataframe = single_snp(test_snps = test_snps, pheno = pheno_fn, count_A1=True, output_file_name = out_fn, covar = cov_fn)

####################################################################################################################################
###Run using R
#take results from the second FaST-LMM scan and identify the most significant SNP above the significance threshold from each chromosome
#make a separate file which contains the genotypes at the loci detected from the first and second FaST-LMM to be included as covaraites for the third FaST-LMM scan
load('/oak/stanford/groups/msalit/tmatsui2/Glucose/Glucose_geno.RData')
geno = as.matrix(geno)
pheno = cbind(1:(ncol(geno) - 3), 1:(ncol(geno) - 3))
pval = read.table('/oak/stanford/groups/msalit/tmatsui2/Glucose/fastlmm/output_fastlmm/fast_lmm_glu2.txt', header = T, as.is = T)
pval$LOD = -log10(pval$PValue)

listy = list()
for(x in 1:16){
    temp = pval[which(pval$Chr == x),]
    if(max(temp$LOD) > 4.3){
        df = temp[which(abs(temp$SnpWeight) == max(abs(temp$SnpWeight))),]
        if(nrow(df) == 1){
            listy[[x]] = df
        } else if(nrow(df) > 1){
            listy[[x]] = df[1,]
        }
    }
}
tab = do.call(rbind, listy)
tab$row = paste(tab$Chr, tab$ChrPos, sep = '_')
tab$snp = as.numeric(unlist(strsplit(tab$SNP, split = 'rs'))[seq(2, nrow(tab) * 2, 2)])
covar = lapply(1:nrow(tab), function(x){
    geno[tab$snp[x], 4:ncol(geno)]
})
covar2 = do.call(cbind, covar)
rownames(covar2) = NULL
covar3 = read.table('/oak/stanford/groups/msalit/tmatsui2/Glucose/fastlmm/code/glu_cov_fr1.txt', header = F, as.is = T) #change for each iteration
covar4 = cbind(covar3, covar2)
write.table(covar4, file = '/oak/stanford/groups/msalit/tmatsui2/Glucose/fastlmm/code/glu_cov_fr2.txt', quote = F, col.names = F, 
    row.names = F, sep = '\t')

###############################################################################################################################
###the above 2 codes were repeated until no significant SNPs were identified