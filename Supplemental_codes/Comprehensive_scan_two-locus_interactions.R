###This code was used to do a comprehensive scan for two-locus interactions
###Written for R. Ran on a HPC cluster using job array (1-779) on SLURM
####################################################################################################################################
library('parallel')

###Get job array number as numeric values (1-779)
task_id = Sys.getenv('SLURM_ARRAY_TASK_ID')
x = as.numeric(task_id)
outpath = '/oak/stanford/groups/msalit/tmatsui2/CoCl2/output_2D_scan/' ##change

load('/oak/stanford/groups/msalit/tmatsui2/CoCl2/CoCl2_geno.RData')
geno = as.matrix(geno)
load('/oak/stanford/groups/msalit/tmatsui2/CoCl2/CoCl2_pheno.RData')

##Find the non-additive portion of each diploids' phenotype
p = lm(pheno$qnorm ~ pheno$midparent)$residuals

##Subset 7742 SNP markers into smaller groups of 200
start = seq(1, 7742, 200)
end = c(seq(200, 7742, 200), 7742)
row = 1:7742
names = rownames(geno)

##Take 2 subsets of SNP markers and test all pairwise combinations for interaction
if(x <= 741){
    combo = combn(1:length(start), 2)[, x]
    c1 = start[combo[1]]:end[combo[1]]
    c2 = start[combo[2]]:end[combo[2]]
    test = expand.grid(c1, c2)
    
    geno1 = geno[c1, 4:ncol(geno)]
    row1 = row[which(row %in% c1)]
    for(z in 1:nrow(geno1)){
        geno1[z,] = as.character(geno1[z,])
    }
     
    geno2 = geno[c2, 4:ncol(geno)]
    row2 = row[which(row %in% c2)]
    for(z in 1:nrow(geno2)){
        geno2[z,] = as.character(geno2[z,])
    }
    
    rm(geno) 
    
    output = lapply(1:nrow(test), function(y){
        l1 = geno1[which(row1 == test[y, 1]),]
        l2 = geno2[which(row2 == test[y, 2]),]
        c(names[test[y, 1]], names[test[y, 2]], -log10(anova(lm(p ~ l1 * l2))[[5]][3]))
    })
##Within a subset of SNP markers, test all unique pairs for interaction    
} else {
    c1 = start[x - 741]:end[x - 741]
    test = t(combn(c1, 2))
    
    geno1 = geno[c1, 4:ncol(geno)]
    row1 = row[which(row %in% c1)]
    for(z in 1:nrow(geno1)){
        geno1[z,] = as.character(geno1[z,])
    }
    
    rm(geno)

    output = lapply(1:nrow(test), function(y){
        l1 = geno1[which(row1 == test[y, 1]),]
        l2 = geno1[which(row1 == test[y, 2]),]
        c(names[test[y, 1]], names[test[y, 2]], -log10(anova(lm(p ~ l1 * l2))[[5]][3]))
    })
}

output = do.call(rbind, output)
output = data.frame(output, stringsAsFactors = F)
output[,3] = as.numeric(output[,3])
colnames(output) = c('l1', 'l2', 'pval')
save(output, file = paste(outpath, 'output_', x, '.RData', sep = ''))

