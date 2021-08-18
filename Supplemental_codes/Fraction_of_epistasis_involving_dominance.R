###This code was used to calculate the fraction of epistasis involving dominance for each two-locus effects.
#######################################################################################################################################
###Each pairwise interaction was first partitioned into four epistasis types: additive-additive, dominance-additive, additive-dominance, and dominance-dominance
###Then the percent variance explained by each epistasis type was determined
###Using these PVE values, the fraction of epistasis involving dominance was calculated
listy = list()
for(p in 1:8){
    setwd('/Users/GIAB/Dropbox/Diploid_project/geno')
    load(dir()[p])
    
    setwd('/Users/GIAB/Dropbox/Diploid_project/pheno')
    load(dir()[p])
    
    setwd('/Users/GIAB/Dropbox/Diploid_project/2D')    
    load(dir()[p])
    all6 = all6[,1:15]
    
    residual = lm(qnorm ~ midparent, data = pheno)$residuals

    non = lapply(1:nrow(all6), function(x){
        print(x)
        temp = all6[x,]
        a1 = geno[which(rownames(geno) == temp$geno1), 4:ncol(geno)] ###locus1 additive
        a2 = geno[which(rownames(geno) == temp$geno2), 4:ncol(geno)] ###locus2 additive
        d1 = as.factor(geno[which(rownames(geno) == temp$geno1), 4:ncol(geno)]) ###locus1 dominance
        d2 = as.factor(geno[which(rownames(geno) == temp$geno2), 4:ncol(geno)]) ###locus2 dominance
        model = anova(lm(residual ~ d1 + d2 + a1:a2 + a1:d2 + d1:a2 + d1:d2))
        
        iaa = model[[5]][3]  ###PVE of additive-additive
        iad = model[[5]][4]  ###PVE of additive-dominance
        ida = model[[5]][5]  ###PVE of dominance-additive
        idd = model[[5]][6]  ###PVE of dominance-dominance
        
        aa = model[[2]][3] / sum(model[[2]])
        ad = model[[2]][4] / sum(model[[2]])
        da = model[[2]][5] / sum(model[[2]])
        dd = model[[2]][6] / sum(model[[2]])
        
        l1a = aa + ad
        l1d = da + dd
        l2a = aa + da
        l2d = ad + dd

        data.frame(iaa = iaa, iad = iad, ida = ida, idd = idd, l1a = l1a, l1d = l1d, l2a = l2a, l2d = l2d)
    })
    
    
    ####Determines if the involved locus is a hub (involved in >20 interactions)
    hubs = names(which(table(c(all6$congeno1, all6$congeno2)) > 20))
    all6$hub1 = sapply(1:nrow(all6), function(x){
        temp = all6[x,]
        if(temp$congeno1 %in% hubs){
            1
        } else {
            0
        }
    })
    all6$hub2 = sapply(1:nrow(all6), function(x){
        temp = all6[x,]
        if(temp$congeno2 %in% hubs){
            1
        } else {
            0
        }
    })
    all6$env = cond2[p]
    
    save(all6, file = dir()[p])
    
    listy[[p]] = all6 
}


