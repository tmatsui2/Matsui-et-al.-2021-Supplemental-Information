###Run using R. This code was used to do a genome-wide scan for one-locus effects with non-additive effects (e.g. dominance)
######################################################################################################################################################
library(parallel)

##Load genotype, phenotype, and snp info
load("/oak/stanford/groups/msalit/tmatsui2/CoCl2/CoCl2_geno.RData")
geno = as.matrix(geno)
load("/oak/stanford/groups/msalit/tmatsui2/CoCl2/CoCl2_pheno.RData")
snp_info = geno[,1:3]

###Find the non-additive portion of each diploids's phenotype, which also effectively accounts for family structure
residual = lm(qnorm ~ midparent, data = pheno)$residuals

##Genome-wide scan for one-locus effects with non-additive effects using the non-additive portion of each diploids' phenotype, first scan
dominance = mclapply(1:nrow(snp_info), function(y){
        l1 = as.character(geno[y, 4:ncol(geno)])
        model = anova(lm(residual ~ l1))
        c(model[[5]][1], model[[4]][1])
})

table = data.frame(matrix(unlist(dominance), ncol = 2, byrow = T))
table = cbind(snp_info, table)
names(table) = c(names(snp_info), 'pval', 'fstat')
table$LOD = -log10(table$pval)
save(table,  file = paste('/oak/stanford/groups/msalit/tmatsui2/CoCl2/1D_fr/FR0_CoCl2.RData', sep = ''))

###Find the most significant SNP above the significance threshold from each chromosome
sig = lapply(1:16, function(x){
    temp = table[which(table$c == x),]
    temp = temp[which(temp$fstat == max(temp$fstat, na.rm = T)),]
})
sig = do.call(rbind, sig)
sig = sig[!duplicated(sig$fstat),]
threshold = 5.3
sig = sig[which(sig$LOD > threshold),]
save(sig,  file = paste('/oak/stanford/groups/msalit/tmatsui2/CoCl2/1D_fr/Peaks_FR0_CoCl2.RData', sep = ''))

######################################################################################################################################################
###Forward linear regression to identify one-locus effects with non-additive effects. For each round of forward regression, identified loci were included as covariates
for(i in 1:10){
    ###make new formula where the identified loci were included as covariates
    formula = 'residual ~ '
    for(x in 1:nrow(sig)){
        covar = paste('c', x, sep = '')
        assign(covar, as.character(geno[which(rownames(geno) == rownames(sig)[x]), 4:ncol(geno)]))
        if(x < nrow(sig)){
            formula = paste(formula, covar, '+', sep = ' ')
        } else {
            formula = paste(formula, covar, sep = ' ')
        }
    }
    residual2 = lm(formula)$residuals
    
    ##Forward scan for one-locus effects with non-additive effects using phenotype residuals not explained by identified loci
    dominance1 = mclapply(1:nrow(snp_info), function(y){
            l1 = as.character(geno[y, 4:ncol(geno)])
            model = anova(lm(residual2 ~ l1))
            c(model[[5]][1], model[[4]][1])
    })
    table1 = data.frame(matrix(unlist(dominance1), ncol = 2, byrow = T))
    table1 = cbind(snp_info, table1)
    names(table1) = c(names(snp_info), 'pval', 'fstat')
    table1$LOD = -log10(table1$pval)
    save(table1,  file = paste('/oak/stanford/groups/msalit/tmatsui2/CoCl2/1D_fr/FR', i, '_CoCl2.RData', sep = ''))
    
    ##Find the most significant SNP above the significance threshold from each chromosome, and 
    ##add newly identified loci to the list of loci to be included as covariates in the next round of forward scan
    sig1 = lapply(1:16, function(x){
        temp = table1[which(table1$c == x),]
        temp = temp[which(temp$fstat == max(temp$fstat, na.rm = T)),]
    })
    sig1 = do.call(rbind, sig1)
    sig1 = sig1[!duplicated(sig1$fstat),]
    sig1 = sig1[which(sig1$LOD > threshold),]
    sig = rbind(sig, sig1)
    save(sig,  file = paste('/oak/stanford/groups/msalit/tmatsui2/CoCl2/1D_fr/Peaks_FR', i, '_CoCl2.RData', sep = ''))
}
