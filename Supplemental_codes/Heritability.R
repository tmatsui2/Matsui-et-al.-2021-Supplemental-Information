###For finding broad-sense and narrow-sense heritability
############################################################################################################
library(lme4)
library(foreach)
library(doMC)
library(rrBLUP)
library(qtl)
library(reshape2)
library(doParallel)
library(parallel)
library(sommer)

load('/oak/stanford/groups/msalit/tmatsui2/CoCl2/Corrected_estimated_fitness_of_replicates.RData')
p = replicates
p$geno = 1:nrow(p)

############################################################################################################
###calculate broad-sense heriability. Because of the large memory overhead required, 2500 strains were randomly chosen to calculate broad-sense heritability
###This process was repeated 1000 times
calcH2 = function(i) {
    mclapply(1:1000, function(x){
        sample = sample(1:nrow(i), 2500)
        i = i[sample, ]
    	test = na.omit(melt(i, id.vars = 'geno'))
    	geno = as.factor(test$geno)
        value = as.numeric(test$value)
    	var = lm(value ~ geno)
    	vcomp = anova(var)[[2]]
    	vcomp[1] / (vcomp[2] + vcomp[1])
    }, mc.cores = 4)
}

broad1 = calcH2(p)
save(broad1, file = '/oak/stanford/groups/msalit/tmatsui2/CoCl2/broad_heritaility2.RData')

############################################################################################################
###Calculate haritability due to additvity, dominance, and epistasis using the sommer R package
###Similar to broad-sense heriability, 2500 strains were randomly chosen each iteration
###This process was repeated 1000 times
load('/oak/stanford/groups/msalit/tmatsui2/CoCl2/Corrected_estimated_fitness_of_replicates.RData')
load('/oak/stanford/groups/msalit/tmatsui2/CoCl2/CoCl2_geno.RData')

p = replicates
p$mean_fitness = sapply(1:nrow(p), function(x){
    mean(unlist(p[x, 1:4]), na.rm = TRUE)
})
g = as.matrix(geno)
info = g[, 1:3]
g = g[,which(colnames(g) %in% p$geno)]
g = g-1

start_time = Sys.time()
sommer = mclapply(1:1000, function(x) {
    random = sample(4:ncol(g), 2500)
    g_random = g[, random]
    g_random = cbind(info, g_random[,order(colnames(g_random))])
    p_random = p[which(p$geno %in% colnames(g_random)),]
    p_random = p_random[order(p_random$geno),]

    #table(colnames(g_random[,4:ncol(g_random)]) == p_random$geno)

    pheno = data.frame(pheno = p_random$mean_fitness, id = p_random$geno, idd = p_random$geno, ide = p_random$geno)
    gdata = t(g_random[, 4:ncol(g_random)])

    A = A.mat(gdata) ###additivity
    D = D.mat(gdata) ###dominance
    E = E.mat(gdata) ###epistasis

    ans.ADE2 = mmer(pheno~1, random = ~vs(id, Gu = A) + vs(idd, Gu = D) + vs(ide, Gu = E), rcov = ~units, data = pheno)

    out = c(pin(ans.ADE2, h2 ~ (V1) / (V1 +  V2 + V3 + V4)), 
    pin(ans.ADE2, h2 ~ (V2) / (V1 +  V2 + V3 + V4)),
    pin(ans.ADE2, h2 ~ (V3) / (V1 +  V2 + V3 + V4)),
    pin(ans.ADE2, h2 ~ (V1 + V2 + V3) / (V1 + V2 + V3 + V4)))
    
    out = unlist(out)
}, mc.cores = 8)

save(sommer, file = '/oak/stanford/groups/msalit/tmatsui2/CoCl2/heritability.RData')
end_time = Sys.time()
end_time - start_time