###Using data from "FaST-LMM_for_one_locus_effects.R", this code is used to identify one-locus effects with significant effect on fitness.
###For each one-locus effect, the peak position as well as the 95% confidence intervals of the locus were found.
###Written for R
#######################################################################################################################################
setwd('/Volumes/TM/All_segregants/fitness_assay_barcodes/CoCl2/fastlmm/')
dir = dir()

##From each forward regression, find the most significant position on each chromosome that is above the significance threshold
qtl1 = data.frame()
for(z in 1:length(dir)){
    pval = read.table(dir[z], header = T, as.is = T)
    pval$LOD = -log10(pval$PValue)

    listy = list()
    for(x in 1:16){
        temp = pval[which(pval$Chr == x),]
        if(max(temp$LOD) > 4.3){
            df = temp[which(abs(temp$LOD) == max(abs(temp$LOD))),]
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
    tab$fr = z
    qtl1 = rbind(qtl1, tab)
}

###Find the 95% confidence interval for each one-locus effect
df = data.frame()
for(i in 1:length(dir)){
    print(i)
    infile1 = paste('./fast_lmm_cocl2', i, '.txt', sep = '')
    pval = read.table(infile1, header = T, as.is = T)
    pval$row = as.numeric(unlist(strsplit(pval$SNP, split = 'rs'))[seq(2, (nrow(pval)*2), 2)])
    pval = pval[order(pval$row),]
    row.names(pval) = paste(pval$Chr, pval$ChrPos, sep = '_')
    pval$LOD = -log10(pval$PValue)
    #pval = na.omit(pval)

    data = for(z in 1:nrow(qtl1)){
        print(z)
        if(round(pval[qtl1[z, 'row'], 'LOD'], 5) == round(qtl1[z, 'LOD'], 5)) {
            temp = qtl1[z,]
            query = pval[which(pval$Chr == temp$Chr & abs(pval$ChrPos - temp$ChrPos) < 50000),]
            peak = min(which(round(query$LOD, 5) == round(temp$LOD, 5)))
        
            peak2 = peak + 1
            check = TRUE
            while(check) {
                if(peak2 > nrow(query)) {
                    check = FALSE
                    p1max = query[peak2 - 1, 'ChrPos']
                } else if(abs(query[peak, 'LOD'] - query[peak2, 'LOD']) < 3) {
                    peak2 = peak2 + 1
                } else {
                    check = FALSE
                    p1max = query[peak2, 'ChrPos']
                }
            }

            peak3 = peak - 1
            check = TRUE
            while(check) {
                if(peak3 == 0) {
                    check = FALSE
                    p1min = query[peak3 + 1, 'ChrPos']
                } else if(abs(query[peak, 'LOD'] - query[peak3, 'LOD']) < 3) {
                    peak3 = peak3 - 1
                } else {
                    check = FALSE
                    p1min = query[peak3, 'ChrPos']
                }
            }
            tab = data.frame(c1 = qtl1[z, 'Chr'], p1 = qtl1[z, 'ChrPos'], pval = temp$LOD,                    
                geno1 = qtl1[z, 'row'], p1min = p1min, p1max = p1max)
            df = rbind(df, tab)
        }
    }
}


###If two one-locus effects have overlapping confidence intervals, remove locus with the less significance
df2 = data.frame()
for(x in unique(df$c1)){
    temp = df[which(df$c1 == x),]
    if(nrow(temp) == 1){
        df2 = rbind(df2, temp)
    } else if(nrow(temp) > 1){
        temp = temp[order(temp$pval, decreasing = T),]
        temp = temp[which(!duplicated(temp$p1)),]
        temp = temp[order(temp$pval, decreasing = F),]
        
        check = TRUE
        while(check){
            if(nrow(temp) == 1){
                df2 = rbind(df2, temp)
                check = FALSE
            } else {
                range = temp[2, 'p1min']:temp[2, 'p1max']
                for(z in 2:nrow(temp)){
                    range = c(range, temp[z, 'p1min']:temp[z, 'p1max'])
                }
                
                if(temp[1, 'p1min']:temp[1, 'p1max'] %overlaps% range){
                    temp = temp[-1,]
                } else {
                    df2 = rbind(df2, temp[1,])
                    temp = temp[-1,]
                }
            }
        }
    }
}

###For each one-locus effects, find the mean fitness values of diploids that are BY/BY, heterozygous, or 3S/3S at the given locus
avg = t(sapply(1:nrow(df2), function(x){
    locus = as.character(geno[as.character(df2$geno1[x]), 4:ncol(geno)])
    by = mean(pheno$qnorm[which(locus == '0')])
    het = mean(pheno$qnorm[which(locus == '1')])
    s = mean(pheno$qnorm[which(locus == '2')])
    c(by, het, s)
}))
df2 = cbind(df2, avg)
names(df2) = c(names(df), 'by', 'het', 's')

save(df2, file = paste('../1D_3LOD_cocl2.Rdata', sep = ''))