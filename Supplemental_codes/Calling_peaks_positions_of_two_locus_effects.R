###Using data from the comprehensive scan for two-locus interactions, this code is used to identify two-locus interactions with significant effect on fitness.
###For each two-locus interaction, the peak position of each involved locus as well as the 95% confidence intervals of each involved locus were found.
###Written for R
####################################################################################################################################
library(parallel)
library(DescTools)

##Load genotype, phenotype, snp info
load("/Volumes/TM/All_segregants/fitness_assay_barcodes/CoCl2/CoCl2_geno.RData")
geno = as.matrix(geno)
load("/Volumes/TM/All_segregants/fitness_assay_barcodes/CoCl2/CoCl2_pheno.RData")
load("/Volumes/TM/All_segregants/fitness_assay_barcodes/CoCl2/snp_info.RData")

###############################################################################################################################
###Combine scan results from comprehenesive scan for two-locus effects
setwd('/Volumes/TM/All_segregants/fitness_assay_barcodes/CoCl2/2D_scan/output/')
dir = dir()

pvals = data.frame()
threshold = 5
for(x in dir){
    load(x)
    output = output[which(output$pval > threshold),]
    pvals = rbind(pvals, output)
}
data = t(sapply(1:nrow(pvals), function(x){
    one = strsplit(pvals$l1[x], split = '_')[[1]]
    c1 = one[1]
    p1 = one[2]
    
    two = strsplit(pvals$l2[x], split = '_')[[1]]
    c2 = two[1]
    p2 = two[2]
    as.numeric(c(c1, p1, c2, p2))
}))
pvals = cbind(data, pvals)
names(pvals) = c('c1', 'p1', 'c2', 'p2', 'locus1', 'locus2', 'pval')

residual = lm(qnorm ~ midparent, data = pheno)$residuals
inf = pvals[which(pvals$pval == 'Inf'),]
if(nrow(inf) > 0){
    pvals = pvals[-which(pvals$pval == 'Inf'),]
}

###############################################################################################################################
####Used to find the peak positions of each involved loci for each two-locus effects. This was achieved by first determining the pair of loci with the most significant p-value. 
###Then, all pairs of loci in which both loci are within 50kb of the identified two-locus interactions were removed to avoid counting the same two-locus 
###interactions multiple times. From the remaining sets of pairs of loci, the next pair of loci with the most significant p-value was found. Again, all pairs of loci 
###in which both loci are within 50kb of the identified two-locus interactions were removed. This process was repeated until no pairs of loci had p-values 
###above the significance threshold.
findpeaks = function(boundsize = 50000, d) { 
    threshold = 5.3
    peaks = data.frame()
    check = TRUE
    while (check) {
        if(max(d$pval, na.rm = T) > threshold){
            peaks = rbind(peaks, d[which(d$pval == max(d$pval, na.rm = T)),])
            d = d[((d$c1 == peaks[nrow(peaks),]$c1) & (d$p1 %in% (peaks[nrow(peaks),]$p1 - boundsize):(peaks[nrow(peaks),]$p1 + boundsize)) &
                  (d$c2 == peaks[nrow(peaks),]$c2) & (d$p2 %in% (peaks[nrow(peaks),]$p2 - boundsize):(peaks[nrow(peaks),]$p2 + boundsize))) == FALSE,]
        } else {
            check <- FALSE
        }
    }
	return(peaks)
}

hub = c('10_657858', '6_211652', '4_590826', '12_677900', '14_376313') ##change for each condition
#hub = c('10_657858', '6_211652')

peaks = findpeaks(d = pvals)
peaks2 = peaks[!duplicated(peaks$pval),]
peaks2 = peaks2[-which(peaks2$c1 == peaks2$c2 & abs(peaks2$p1 - peaks2$p2) < 100000),] ##remove all pairs of loci where the involved loci are 100kb of each other
sort(table(unlist(peaks2[,c('locus1','locus2')])))

###remove all pairs of loci involving a hub locus
remove = sapply(1:length(hub), function(x){
    q = hub[x]
    qc = as.numeric(strsplit(q, split = '_')[[1]][1])
    qp = as.numeric(strsplit(q, split = '_')[[1]][2])
    remove1 = -which(peaks2$c1 == qc & abs(peaks2$p1 - qp) < 100000)
    remove2 = -which(peaks2$c2 == qc & abs(peaks2$p2 - qp) < 100000)
    remove = c(remove1, remove2)
})
all_remove = unique(unlist(remove))
peaks3 = peaks2[all_remove,]

###############################################################################################################################
#####For each two-locus interactin, find the 95% confidence interval of each involved locus
pvals = data.frame()
threshold = 2
for(x in dir){
    load(x)
    output = output[which(output$pval > threshold),]
    pvals = rbind(pvals, output)
}
data = t(sapply(1:nrow(pvals), function(x){
    one = strsplit(pvals$l1[x], split = '_')[[1]]
    c1 = one[1]
    p1 = one[2]
    
    two = strsplit(pvals$l2[x], split = '_')[[1]]
    c2 = two[1]
    p2 = two[2]
    as.numeric(c(c1, p1, c2, p2))
}))
pvals = cbind(data, pvals)
names(pvals) = c('c1', 'p1', 'c2', 'p2', 'geno1', 'geno2', 'pval')

p1lod = sapply(1:nrow(peaks3), function(x){
    temp = peaks3[x,]
    query = pvals[which(pvals$c2 == temp$c2 & pvals$p2 == temp$p2 & pvals$c1 == temp$c1 & abs(pvals$p1 - temp$p1) < 50000),]
    peak = min(which(query$pval == temp$pval))
    
    peak2 = peak + 1
    check = TRUE
    while(check) {
        if(peak2 > nrow(query)) {
            check = FALSE
            p1max = query[peak2 - 1, 'p1']
        } else if(abs(query[peak, 'pval'] - query[peak2, 'pval']) < 3) {
            peak2 = peak2 + 1
        } else {
            check = FALSE
            p1max = query[peak2, 'p1']
        }
    }
    
    peak3 = peak - 1
    check = TRUE
    while(check) {
        if(peak3 == 0) {
            check = FALSE
            p1min = query[peak3 + 1, 'p1']
        } else if(abs(query[peak, 'pval'] - query[peak3, 'pval']) < 3) {
            peak3 = peak3 - 1
        } else {
            check = FALSE
            p1min = query[peak3, 'p1']
        }
    }
    c(p1min, p1max)
})
p2lod = sapply(1:nrow(peaks3), function(x){
    temp = peaks3[x,]
    query = pvals[which(pvals$c1 == temp$c1 & pvals$p1 == temp$p1 & pvals$c2 == temp$c2 & abs(pvals$p2 - temp$p2) < 50000),]
    peak = min(which(query$pval == temp$pval))
    
    peak2 = peak + 1
    check = TRUE
    while(check) {
        if(peak2 > nrow(query)) {
            check = FALSE
            p2max = query[peak2 - 1, 'p2']
        } else if(abs(query[peak, 'pval'] - query[peak2, 'pval']) < 3) {
            peak2 = peak2 + 1
        } else {
            check = FALSE
            p2max = query[peak2, 'p2']
        }
    }
    
    peak3 = peak - 1
    check = TRUE
    while(check) {
        if(peak3  == 0) {
            check = FALSE
            p2min = query[peak3 + 1, 'p2']
        } else if(abs(query[peak, 'pval'] - query[peak3, 'pval']) < 3) {
            peak3 = peak3 - 1
        } else {
            check = FALSE
            p2min = query[peak3, 'p2']
        }
    }
    c(p2min, p2max)
})
peaks4 = cbind(peaks3, t(p1lod), t(p2lod))
names(peaks4) = c(names(pvals), 'p1min', 'p1max', 'p2min', 'p2max')
save(peaks4, file = '../3LOD_2D_peaks.RData')
