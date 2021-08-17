###This code was used to do a comprehensive scan for three-locus interactions
###Written for R. Ran on a HPC cluster using job array (1-500) on SLURM
####################################################################################################################################
library('parallel')

###Get job array number as numeric values (1-500)
task_id = Sys.getenv('SLURM_ARRAY_TASK_ID')
x = as.numeric(task_id)
path = '/oak/stanford/groups/msalit/tmatsui2/CoCl2/'
subDir = '3D'
dir.create(file.path(path, subDir), showWarnings = FALSE)
outpath = '/oak/stanford/groups/msalit/tmatsui2/CoCl2/3D/'

load(paste(path, 'CoCl2_geno.RData', sep = ''))
load(paste(path, 'CoCl2_pheno.RData', sep = ''))
load('/oak/stanford/groups/msalit/tmatsui2/5cM_subset_row.RData')

subset = as.matrix(geno)[subset, 4:ncol(geno)]
rm(geno)
##Find the non-additive portion of each diploids' phenotype
residual = lm(qnorm ~ midparent, data = pheno)$residuals

##Determine all unique 3-locus combinations and divide the combinations into 500 parts
combo = combn(nrow(subset), 3)
col = round(ncol(combo)/500, 0)
start = seq(from = 1, to = ncol(combo), col)[x]
end = (start + col - 1)
if(x == 500){
	end = ncol(combo)
}

##For each job array, test all 3-locus combinations within one of the 500 parts for interaction
s = Sys.time()
out = lapply(start:end, function(z) {
	l1 = as.character(subset[combo[,z][1],])
    l2 = as.character(subset[combo[,z][2],])
    l3 = as.character(subset[combo[,z][3],])
	c(rownames(subset)[combo[,z][1]], rownames(subset)[combo[,z][2]], rownames(subset)[combo[,z][3]], 
		anova(lm(residual ~ l1 * l2 * l3))[[5]][7])
})
e = Sys.time()
print(e - s) 

output = do.call('rbind',out)
output = data.frame(output, stringsAsFactors = F)
output[,4] = as.numeric(output[,4])
colnames(output) = c('l1', 'l2', 'l3', 'pval')
save(output, file = paste(outpath, 'CoCl2_3D_output_', x, '.RData', sep = ''))

