###Run using R. This code was used to do a forward scan to identify loci that interact with a hub locus.
######################################################################################################################################################
library(parallel)

##genotype, phenotype, snp info
load("/oak/stanford/groups/msalit/tmatsui2/CoCl2/CoCl2_geno.RData")
geno = as.matrix(geno)
load("/oak/stanford/groups/msalit/tmatsui2/CoCl2/CoCl2_pheno.RData")
snp_info = geno[,1:3]
hub = c('10_657858', '6_211652', '4_590826', '12_677900', '14_376313')

###Find the non-additive portion of each diploids's phenotype, which also effectively accounts for family structure
residual = lm(qnorm ~ midparent, data = pheno)$residuals

##FR with X lod threshold
for(z in 1:length(hub)){
    ##Genome-wide scan for loci that interact with a hub locus, first scan
    l1 = as.character(geno[which(rownames(geno) == hub[z]), 4:ncol(geno)])
    hotspotFR = mclapply(1:nrow(snp_info), function(y){
            l2 = as.character(geno[y, 4:ncol(geno)])
            model = anova(lm(residual ~ l1*l2))
            c(model[[5]][3], model[[4]][3])
    })

    table = data.frame(matrix(unlist(hotspotFR), ncol = 2, byrow = T))
    table = cbind(snp_info, table)
    names(table) = c(names(snp_info), 'pval', 'fstat')
    table$LOD = -log10(table$pval)
    save(table,  file = paste('/oak/stanford/groups/msalit/tmatsui2/CoCl2/2D_fr/FR0_', hub[z], '.RData', sep = ''))
    
    ###Find the most significant SNP above the significance threshold that interact with a hub locus from each chromosome
    sig = lapply(1:16, function(x){
        temp = table[which(table$c == x),]
        temp = temp[which(temp$fstat == max(temp$fstat, na.rm = T)),]
    })
    sig = do.call(rbind, sig)
    sig = sig[!duplicated(sig$fstat),]
    threshold = 5.3
    sig = sig[which(sig$LOD > threshold),]
    save(sig,  file = paste('/oak/stanford/groups/msalit/tmatsui2/CoCl2/2D_fr/Peaks_FR0_', hub[z], '.RData', sep = ''))
    
    ###Forward linear regression to identify hub modifiers. For each round of forward regression, identified loci were included as covariates
    for(i in 1:10){
        ###make new formula where the identified loci and interaction with the hub were included as covariates
        formula = 'residual ~ l1'
        for(x in 1:nrow(sig)){
            covar = paste('c', x, sep = '')
            assign(covar, as.character(geno[which(rownames(geno) == rownames(sig)[x]), 4:ncol(geno)]))
            formula = paste(formula, '+', covar, '+', paste(covar,':l1', sep = ''))
        }
        residual2 = lm(formula)$residuals
        
        ##Forward scan for hub modifiers using phenotype residuals not explained by identified loci and interaction with the hub
        hotspotFR1 = mclapply(1:nrow(snp_info), function(y){
                l2 = as.character(geno[y, 4:ncol(geno)])
                model = anova(lm(residual2 ~ l1*l2))
                c(model[[5]][3], model[[4]][3])
        })
        table1 = data.frame(matrix(unlist(hotspotFR1), ncol = 2, byrow = T))
        table1 = cbind(snp_info, table1)
        names(table1) = c(names(snp_info), 'pval', 'fstat')
        table1$LOD = -log10(table1$pval)
        save(table1,  file = paste('/oak/stanford/groups/msalit/tmatsui2/CoCl2/2D_fr/FR', i, '_', hub[z], '.RData', sep = ''))
        
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
        save(sig,  file = paste('/oak/stanford/groups/msalit/tmatsui2/CoCl2/2D_fr/Peaks_FR', i, '_', hub[z], '.RData', sep = ''))
    }
}
