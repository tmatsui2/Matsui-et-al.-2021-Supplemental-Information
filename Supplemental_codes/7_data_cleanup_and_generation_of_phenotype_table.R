###written for R. The output files from FitSeq is filtered based on log likelihood scores and number of biological replicates (there mst be at least 2+). Biological replicates with extremely different fitness estimates in comparison to others (outliers) are also omitted, as this is most likely due to spontaneous mutations. The filtered data is then used to generate the phenotype and genotype data tables for each environment, which will used for all downstream analyses. The code must be run in this exact order. Each fucntion creates an output file, which is then used as the input file for the next function.
library(plotrix)
library(foreach)
library(doParallel)
library(stringr)
library(gtools)
library(ggplot2)
library(tidyr)
library(zoo)
library(outliers)
library(RColorBrewer)
library(parallel)
library(outliers)
library(preprocessCore)
library(bestNormalize)

###############################################################################################################
setwd('/Volumes/TM/All_segregants/fitness_assay_barcodes/Glucose2/pre-process/')

#Create histogram of likelihood scores and filter any double barcodes with outlier values.
likelihood = function(file, env){
    mainDir = '/Volumes/TM/All_segregants/fitness_assay_barcodes/CoCl2/'
    subDir = "cluster"
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    
    fit = read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = ',')
    name = read.table('Combined_PCR_corrected_fitseq_name.txt', header = F, as.is = T)
    names(name) = 'ID'

    fit = cbind(fit, name)
    write.table(fit, file = 'FitSeq_with_ID.csv', quote = F, sep = ',', col.names = T, row.names = F)
    hist(fit$Likelihood_Log, main = 'Histogram of log likelihood scores', xlab = 'log likelihood score', breaks = 100)
    quartz.save('./figures/histogram likelihood score.pdf', type = 'pdf')
    dev.off()
    
    cutoff1 = as.numeric(quantile(fit$Likelihood_Log)[2]-IQR(fit$Likelihood_Log)*1.5) ##determine which double barcodes have outlier likelihood score
    fit = fit[which(fit$Likelihood_Log > cutoff1),]
    seg = fit[which(str_count(fit$ID, 'many') == 0 & str_count(fit$ID, 'unique') == 0),]
    write.table(seg, file = './cluster/FitSeq_with_ID_good_seg_only.csv', quote = F, sep = ',', col.names = T, row.names = F)
    
	###create figure showing the lineage trajectory for each double barcode, colored by the estimated fitness
    adjusted = lapply(4:(3 + env), function(x){
        average = mean(fit[, x])
        (100/average) * fit[,x]
    })
    adjusted = do.call(cbind, adjusted)
    fit = cbind(fit, adjusted)
    fit2 = fit[, c(1, 8:(8 + env))]
    fit2 = gather(fit2, names(fit2)[3:(2 + env)], key = 'time_point', value = 'read_count')
    fit2$log_read = log10(fit2$read_count)

    ###reduced number for plotting
    id = unique(fit2$ID)
    set.seed(10)
    random = sample(1:length(id), 30000)
    fit3 = fit2[which(fit2$ID %in% id[random]),]

    mid = median(fit3$Estimated_Fitness)
    ggplot(fit3, aes(x = time_point, y = log_read, group = ID)) + 
        geom_line(aes(colour = Estimated_Fitness), size = 0.3) + 
        scale_colour_gradient2(midpoint = mid, low = 'blue', mid = 'white', high = 'red')  
    ggsave('./figures/lineage trajectory.pdf', height = 8, width = 15)
}
likelihood('output_FitSeq.csv', 4) ###output file from fitseq, and the number of time points in this environment

#Create histogram of number of replicates per diploid. Here, we only use genotypes with at least 2 replicates, to make sure the fitness estimates are reproducible across biologival replicates
histogram = function(file){
    seg = read.table(file, header = T, as.is = T, sep = ',')
    geno = lapply(1:nrow(seg), function(x){
        MATa = strsplit((strsplit(seg$ID[x], split = '\\+')[[1]][2]), split = '\\.')[[1]][2]
        MATalpha = strsplit((strsplit(seg$ID[x], split = '\\+')[[1]][2]), split = '\\.')[[1]][1]
        tab = cbind(MATa, MATalpha)
    })
    geno = do.call(rbind, geno)
    seg = cbind(seg, geno, stringsAsFactors = FALSE)
    seg$diploid = paste(seg$MATalpha, seg$MATa, sep = '.')
    write.table(seg, file = './cluster/FitSeq_with_ID_good_seg_only_detail.csv', quote = F, row.names = FALSE, col.names = T, sep = ',')
    rep = sapply(1:4, function(x){
        length(which(table(seg$diploid) == x))
    })
    
    legend = data.frame(replicates = 1:4, counts = rep, percent = round(rep/length(unique(seg$diploid)), digits = 3))

    hist(table(seg$diploid), xlab = '# of replicates', main = 'Number of replicates per genotype', xlim = c(1,4))
    addtable2plot(2.2, 55000, legend, xpad = 0.5, ypad = 1, bty = 'n', display.rownames = F, hlines = F, vlines = F, cex = 1)
    quartz.save('./figures/number of replicates per genotype.pdf', type = 'pdf')
    dev.off()
}
histogram('./cluster/FitSeq_with_ID_good_seg_only.csv')

#For each genotype, calculate standard error across replicates. Plot standard error against the mean of the replicates to show how fitness estimate accuracy improves with increasing fitness.
stderror = function(file){
    seg = read.table(file, header = T, as.is = T, sep = ',')
    one = names(which(table(seg$diploid) == 1))
    seg = seg[-which(seg$diploid %in% one),] ###remove all genotypes that only have 1 replicate
    write.table(seg, file = './cluster/Fitseq_2+_replicates.csv', quote = F, row.names = FALSE, col.names = T, sep = ',')
    reps_geno = unique(seg$diploid)
    error = t(sapply(1:length(reps_geno), function(x){
        if(x %in% seq(1, length(reps_geno), 1000)){
            print(x/length(reps_geno))
        }
        
        temp = seg[which(seg$diploid == reps_geno[x]),]
        c(mean(temp$Estimated_Fitness, na.rm = T), std.error(temp$Estimated_Fitness))
    }))

    error1 = data.frame(error, geno = reps_geno, stringsAsFactors = FALSE)
    names(error1) = c('mean', 'std', 'geno')
    save(error1, file = './cluster/error.RData')

    sliding = data.frame(se = rollmean(error1[, 2], 10), fit = rollmean(error1[, 1], 10))

    #standard error as a function of fitness
    #ggplot(error1, aes(x = error[, 1], y = error[, 2])) + geom_bin2d(bins = 300) + scale_fill_continuous(type = "viridis") +
      #theme_bw() + labs(x = 'Mean Estimated Fitness', y = 'Standard error', title = 'Standard error of fitness across barcode replicates')

    ggplot(sliding, aes(x = fit, y = se)) + geom_bin2d(bins = 300) + scale_fill_continuous(type = "viridis") +
        theme_bw() + labs(x = 'Mean Estimated Fitness', y = 'Standard error', title = 'Sliding window SE of fitness across replicates') + 
        geom_smooth()
    ggsave('../figures/sliding window SE of fitness across replicates.pdf', height = 8, width = 8)
}
stderror('./cluster/FitSeq_with_ID_good_seg_only.csv')

#For each genotype, examine all pairwise combinations of replicates and look at the absolute difference in fitness. Creates a table with the fitness of all replicates for all genotypes.
correlation = function(file){
    seg = read.table(file, header = T, as.is = T, sep = ',')
    reps_geno = unique(seg$diploid)
    replicates = data.frame(do.call(rbind, mclapply(1:length(reps_geno), function(x){
        temp = c(seg[which(seg$diploid == reps_geno[x]), 'Estimated_Fitness'])
        if(length(temp) == 1){
            c(temp, 0, 0, 0, reps_geno[x])
        } else if(length(temp) == 2){
            c(temp, 0, 0, reps_geno[x])
        } else if(length(temp) == 3){
            c(temp, 0, reps_geno[x])
        } else if(length(temp) == 4){
            c(temp, reps_geno[x])
        } else if(length(temp) > 4){
            c(temp[1:4], reps_geno[x])
        }
    })))
    names(replicates) = c('r1', 'r2', 'r3', 'r4', 'geno')
    for(x in 1:5){
        if (x != 5){
            replicates[,x] = as.numeric(as.character(replicates[,x]))
        } else {
            replicates[,x] = as.character(replicates[,x])
        }
    }
    save(replicates, file = './cluster/Estimated_fitness_of_replicates.RData')

    #correlation of all pairs of replicates per genotype
    rep12 = data.frame(r1 = replicates$r1, r2 = replicates$r2, stringsAsFactors = FALSE, geno = replicates$geno)
    rep13 = data.frame(r1 = replicates$r1, r2 = replicates$r3, stringsAsFactors = FALSE, geno = replicates$geno)
    rep14 = data.frame(r1 = replicates$r1, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep23 = data.frame(r1 = replicates$r2, r2 = replicates$r3, stringsAsFactors = FALSE, geno = replicates$geno)
    rep24 = data.frame(r1 = replicates$r2, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep34 = data.frame(r1 = replicates$r3, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep = rbind.data.frame(rep12, rep13, rep14, rep23, rep24, rep34, stringsAsFactors = FALSE)
    rep$r1 = as.numeric(rep$r1)
    rep$r2 = as.numeric(rep$r2)
    rep$p = rep$r1 * rep$r2
    rep = rep[-which(rep$p == 0),]

    ###looking at the range of estimated fitness across replicates
    rep1 = rep 
    rep1$range = abs(rep1$r1 - rep1$r2)

    hist(rep1$range, breaks = 100, xlab = 'Difference in estimated fitness', 
        main = 'Difference in estimated fitness between replicates')
    abline(v = 0.5, col = 'red', lty = 2)
    quartz.save('./figures/difference in estimated fitness between replicates.pdf', type = 'pdf')
    dev.off()
}
correlation('./cluster/Fitseq_2+_replicates.csv')

#For each genotype, remove outlier replicates with fitness values that differ from other replicates by more than 0.5 (cutoff). Re-create a new table with the fitness of the replicates without outliers for all genotypes.
outlier = function(cutoff){
    seg = read.table('./cluster/Fitseq_2+_replicates.csv', header = T, as.is = T, sep = ',')
    reps_geno = unique(seg$diploid)
    load('./cluster/Estimated_fitness_of_replicates.RData')
    
    rep12 = data.frame(r1 = replicates$r1, r2 = replicates$r2, stringsAsFactors = FALSE, geno = replicates$geno)
    rep13 = data.frame(r1 = replicates$r1, r2 = replicates$r3, stringsAsFactors = FALSE, geno = replicates$geno)
    rep14 = data.frame(r1 = replicates$r1, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep23 = data.frame(r1 = replicates$r2, r2 = replicates$r3, stringsAsFactors = FALSE, geno = replicates$geno)
    rep24 = data.frame(r1 = replicates$r2, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep34 = data.frame(r1 = replicates$r3, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep = rbind.data.frame(rep12, rep13, rep14, rep23, rep24, rep34, stringsAsFactors = FALSE)
    rep$r1 = as.numeric(rep$r1)
    rep$r2 = as.numeric(rep$r2)
    rep$p = rep$r1 * rep$r2
    rep = rep[-which(rep$p == 0),]

    ###looking at the range of estimated fitness across replicates
    rep1 = rep 
    rep1$range = abs(rep1$r1 - rep1$r2)
    rep1_sig = rep1[which(rep1$range > cutoff),]  ##change for each environment

    ##remove outlier replicates
    outlier = unique(rep1_sig$geno)
    outlier_id = mclapply(1:length(outlier), function(x){
        temp = seg[which(seg$diploid %in% rep1_sig$geno[x]),]
        if(nrow(temp) == 2){
            NA
        } else if(nrow(temp) > 2){
            check = TRUE
            while(check){
                md = median(temp$Estimated_Fitness)
                if(abs(max(temp$Estimated_Fitness) - md) > abs(min(temp$Estimated_Fitness) - md)){
                    temp = temp[-which(temp$Estimated_Fitness == max(temp$Estimated_Fitness)),]
                } else {
                    temp = temp[-which(temp$Estimated_Fitness == min(temp$Estimated_Fitness)),]
                }
                if(abs(max(temp$Estimated_Fitness) - min(temp$Estimated_Fitness)) < cutoff){ ##change for each environment
                    check = FALSE
                }
            }
            temp$ID
        } 
    })
    outlier_id = unlist(outlier_id)
    outlier_id = outlier_id[!is.na(outlier_id)]
    outlier = seg[which(seg$diploid %in% rep1_sig$geno),]
    remove = outlier$ID[-which(outlier$ID %in% outlier_id)]
    seg1 = seg[-which(seg$ID %in% remove),]
    write.table(seg1, file = './cluster/Corrected_Fitseq_2+_replicates.csv', row.names = F, col.names = T, sep = ',', quote = F)

    ###re-make correlation table after removing outlier replicates
    reps_geno = unique(seg1$diploid)
    replicates = data.frame(do.call(rbind, mclapply(1:length(reps_geno), function(x){
        temp = c(seg1[which(seg1$diploid == reps_geno[x]), 'Estimated_Fitness'])
        if(length(temp) == 1){
            c(temp, 0, 0, 0, reps_geno[x])
        } else if(length(temp) == 2){
            c(temp, 0, 0, reps_geno[x])
        } else if(length(temp) == 3){
            c(temp, 0, reps_geno[x])
        } else if(length(temp) == 4){
            c(temp, reps_geno[x])
        } else if(length(temp) > 4){
            c(temp[1:4], reps_geno[x])
        }
    })))
    names(replicates) = c('r1', 'r2', 'r3', 'r4', 'geno')
    for(x in 1:5){
        if (x != 5){
            replicates[,x] = as.numeric(as.character(replicates[,x]))
        } else {
            replicates[,x] = as.character(replicates[,x])
        }
    }
    save(replicates, file = './cluster/Corrected_estimated_fitness_of_replicates.RData')

    #Calculate correlation of all pairs of replicates per genotype
    rep12 = data.frame(r1 = replicates$r1, r2 = replicates$r2, stringsAsFactors = FALSE, geno = replicates$geno)
    rep13 = data.frame(r1 = replicates$r1, r2 = replicates$r3, stringsAsFactors = FALSE, geno = replicates$geno)
    rep14 = data.frame(r1 = replicates$r1, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep23 = data.frame(r1 = replicates$r2, r2 = replicates$r3, stringsAsFactors = FALSE, geno = replicates$geno)
    rep24 = data.frame(r1 = replicates$r2, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep34 = data.frame(r1 = replicates$r3, r2 = replicates$r4, stringsAsFactors = FALSE, geno = replicates$geno)
    rep = rbind.data.frame(rep12, rep13, rep14, rep23, rep24, rep34, stringsAsFactors = FALSE)
    rep$r1 = as.numeric(rep$r1)
    rep$r2 = as.numeric(rep$r2)
    rep$p = rep$r1 * rep$r2
    rep[rep == 0] = NA
    rep = rep[which(is.na(rep$p) == FALSE),]

    ###XY plot showing the fitnesses of all pairs of replicates for every genotype and calculating the correlation between them
    spear = cor(rep$r1, rep$r2, method = 'spearman')
    ggplot(rep, aes(x = r1, y = r2)) + geom_bin2d(bins = 300) + scale_fill_continuous(type = "viridis") +
      theme_bw() + labs(x = 'Replicate 1', y = 'Replicate 2', title = 'Correlation between replicates') +
      annotate("text", x = -0.55, y = 0.4, label = paste('Spearman Cor = ', round(spear, digits = 3), sep = '')) +
      geom_abline(slope = 1, intercept = 0, linetype = 'dashed', size = 0.2)
    ggsave('./figures/correlation between replicates.pdf', height = 8, width = 8)
}
outlier(0.5)

####Filter data based on likelihood scores, then plot the distribution of fitnesses of the diploids along with the fitnesses of the parental controls
distribution = function(file, file2, env){
	###determine double barcodes with significantly poor likelihood scores and remove from dataset
    fit = read.table(file, header = T, stringsAsFactors = FALSE, sep = ',')
    cutoff1 = as.numeric(quantile(fit$Likelihood_Log)[2]-IQR(fit$Likelihood_Log)*1.5) 
    fit = fit[which(fit$Likelihood_Log > cutoff1),]
    control = fit[c(which(str_count(fit$ID, 'unique') == 2), which(str_count(fit$ID, 'many') == 2)),]
    high = control[which(rowSums(control[,4:(3 + env)]) > 50),]
	
	###read in table with the barcodes of the parental controls
    a_parent = read.table('/Volumes/TM/Parental_barcode_control/parental_only_MATa_barcodes_detail.txt', header = F, as.is = F)
    alpha_parent = read.table('/Volumes/TM/Parental_barcode_control/parental_only_MATalpha_barcodes_detail.txt', header = F, as.is = F)
    a_parent$ID = paste(a_parent[,3], '_a', sep = '')
    alpha_parent$ID = paste(alpha_parent[,3], '_alpha', sep = '')
	
	
    listy = list()
    for(x in 1:nrow(high)){
        listy[[x]] = strsplit(strsplit(high$ID[x], split = '\\+')[[1]][2], split = '\\.')[[1]][2]
    }
    high$MATa = unlist(listy)

    listy = list()
    for(x in 1:nrow(high)){
        listy[[x]] = strsplit(strsplit(high$ID[x], split = '\\+')[[1]][2], split = '\\.')[[1]][1]
    }
    high$MATalpha = unlist(listy)

    listy = list()
    for(x in 1:nrow(high)){
        listy[[x]] = as.character(a_parent[which(a_parent$ID == high[x, 'MATa']), 2])
    }
    high$MATa_geno = unlist(listy)

    listy = list()
    for(x in 1:nrow(high)){
        listy[[x]] = as.character(alpha_parent[which(alpha_parent$ID == high[x, 'MATalpha']), 2])
    }
    high$MATalpha_geno = unlist(listy)

    high$diploid = paste(high$MATa_geno, high$MATalpha_geno, sep = '')
    diploid = c('BYBY', 'BY3S', '3SBY', '3S3S')

	###determine median fitness of all biological replicates of the parental controls
    listy = list()
    for(x in 1:length(diploid)){
       listy[[x]] = median(high[which(high$diploid == diploid[x]), 'Estimated_Fitness']) 
    }
    control_fitness = unlist(listy)
    
    seg = read.table(file2, header = T, as.is = T, sep = ',')
    mypal = brewer.pal(8,'Dark2')
    hist(seg$Estimated_Fitness, xlab = 'Estimated Fitness', main = 'Distribution of estimated fitness', breaks = 100)
    abline(v = control_fitness[1], col = mypal[1], lwd = 2)
    abline(v = control_fitness[2], col = mypal[2], lwd = 2)
    abline(v = control_fitness[3] + 0.005, col = mypal[3], lwd = 2)
    abline(v = control_fitness[4], col = mypal[4], lwd = 2)
    legend('topright', c('BY/BY', 'BY/3S', '3S/BY', '3S/3S'), col = mypal[1:4], pch = 20, box.lty = 0) 
    quartz.save('./figures/distribution of estimated fitness.pdf', type = 'pdf')
    dev.off()
}
distribution('FitSeq_with_ID.csv', './cluster/Corrected_Fitseq_2+_replicates.csv', 3)

#Make phenotype table with the mean fitness of the replicates for each genotype
phenotype = function(file){
    load(file)
    table = data.frame(t(sapply(1:nrow(replicates), function(x){
        temp = replicates[x,]
        temp[temp == 0] = NA
        avg = mean(unlist(temp[,1:4]), na.rm = T)
        parents = strsplit(temp$geno, split = '\\.')[[1]]
        c(pheno = avg, MATa = parents[2], MATalpha = parents[1], geno = paste(parents[1], parents[2], sep = '.'))
    })), stringsAsFactors = F)
    table[,1] = as.numeric(table[,1])
    save(table, file = './pheno.RData')
}
phenotype('./cluster/Corrected_estimated_fitness_of_replicates.RData')

#Make genotype table with only diploids with sufficient coverage (>5 per time point) and 2+ replicates for each condition
genotype = function(geno, snp_info, pheno){
    load(geno)
    load(pheno)
    load(snp_info)
    diploid_genotypes = diploid_genotypes[which(rownames(diploid_genotypes) %in% rownames(snp_info)),]
    diploid_genotypes = diploid_genotypes[,c(which(names(diploid_genotypes) %in% table$geno))]
    pheno = table[order(table$geno),]
    diploid_genotypes = diploid_genotypes[,order(names(diploid_genotypes))]
    geno = diploid_genotypes
    geno = cbind(snp_info, geno)
    geno = as.matrix(geno)
    save(pheno, file = 'CoCl2_pheno.RData')
    save(geno, file = 'CoCl2_geno.RData')
}
genotype('/Volumes/TM/All_segregants/fitness_assay_barcodes/diploid_genotypes.RData', 
    '/Volumes/TM/All_segregants/fitness_assay_barcodes/Glucose/snp_info.RData',
    'pheno.RData')

#Use the bestNormalize package to transform mean fitness values such that it is normalized. Normalized data is then used to estimate the fitnesses of the haploid parental strains.
midparent = function(file){
    load(file)
    un_a = unique(pheno$MATa)
    un_alpha = unique(pheno$MATalpha)
	
	###normalize fitness values using bestNormalize
	#norm = bestNormalize(pheno$pheno)
    pheno$qnorm = orderNorm(pheno$pheno)$x.t
    
	##calculate the fitness of each MATa haploid straina by taking the mean fitness of all diploids that originated from a haploid strain.
    mean_a = data.frame(geno = un_a, MATa_mean = sapply(1:length(un_a), function(x){
        mean_a = mean(pheno[which(pheno$MATa == un_a[x]), 'qnorm'])
    }), stringsAsFactors = F)
    mean_a = mean_a[order(mean_a$MATa_mean),]
    mean_a$rank = 1:nrow(mean_a)
   
    ##calculate the fitness of each MATalpha haploid straina by taking the mean fitness of all diploids that originated from a haploid strain. 
    mean_alpha = data.frame(geno = un_alpha, MATalpha_mean = sapply(1:length(un_alpha), function(x){
        mean_alpha = mean(pheno[which(pheno$MATalpha == un_alpha[x]), 'qnorm'])
    }), stringsAsFactors = F)
    mean_alpha = mean_alpha[order(mean_alpha$MATalpha_mean),]
    mean_alpha$rank = 1:nrow(mean_alpha)
    
    ##find the fitness rank of each haploid MATa and MATalpha parental strains
	data = data.frame(t(sapply(1:nrow(pheno), function(x){
        temp = pheno[x,]
        
        data = c(MATa_mid = mean_a$MATa_mean[which(mean_a$geno == temp$MATa[1])], MATa_rank = mean_a$rank[which(mean_a$geno == temp$MATa[1])],
        MATalpha_mid = mean_alpha$MATalpha_mean[which(mean_alpha$geno == temp$MATalpha[1])], 
        MATalpha_rank = mean_alpha$rank[which(mean_alpha$geno == temp$MATalpha[1])])
        
        return(data)
    })), stringsAsFactors = F)
    pheno = cbind(pheno, data)
	
	##calculate the expected fitness ('midparent') of each diploid assuming by taking the mean of the two haploid parental strain
    pheno$midparent = sapply(1:nrow(pheno), function(x){mean(c(pheno$MATa_mid[x], pheno$MATalpha_mid[x]))})
    save(pheno, file = 'CoCl2_pheno.RData')
    
    #plotting midparent
    data = pheno
    data$MATalpha_rank = data$MATalpha_rank + length(un_a)
    combined = length(un_a) + length(un_alpha)
    
    temp = data[c(which(data$MATa_rank == 1), which(data$MATalpha_rank == 1)),]
    
    plot(rep(1, nrow(temp)), temp$qnorm, pch = 20, cex = 0.5, col = 'black', xlim = c(1, combined), ylim = c(min(data$qnorm), max(pheno$qnorm)),
        main = 'fitness by parent', xlab = 'parent rank', ylab = 'fitness')
    points(temp$MATa_mid[1], col = 'red', pch = 20)
    
    for(x in 2:combined){
        temp = data[c(which(data$MATa_rank == x), which(data$MATalpha_rank == x)),]
        
        if(x %in% 2:length(un_a)){
            points(rep(x, nrow(temp)), temp$qnorm, pch = 20, cex = 0.5, col = 'black')
            points(x, temp$MATa_mid[1], col = 'red', pch = 20)
            
        } else if(x > length(un_a)){
            points(rep(x, nrow(temp)), temp$qnorm, pch = 20, cex = 0.5, col = 'black')
            points(x, temp$MATalpha_mid[1], col = 'red', pch = 20)
        }
    }
    
    quartz.save('./figures/midparent.pdf', type = 'pdf')
    dev.off()
}
midparent('CoCl2_pheno.RData')
