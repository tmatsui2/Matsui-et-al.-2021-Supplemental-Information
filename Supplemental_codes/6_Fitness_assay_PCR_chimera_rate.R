###written for R. The output files from 5_extracting_cluster_matching_a_barcode is used to calculate the PCR chimera rate by combining all time points for each environment. The calculated PCR chimera rate is then used to estimate the true barcode read count after correction for PCR chimeras. The output file is then fed into FitSeq to calculate fitness for each double barcode using the PCR corrected read counts from multiple time points. For more information about FitSeq see(https://github.com/FangfeiLi05/PyFitSeq).
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
env = c('CoCl2', 'CuSO4', 'Glucose', 'Glucose2', 'H2O2', 'NaCl', 'Rapamycin', 'Zeocin')

###############################################################################################################
####reading in tables and correcting for pcr chimera using all time points
pcr_correction = function(files){
    #mainDir = '/Volumes/TM/All_segregants/fitness_assay_barcodes/CoCl2/' 
    #subDir = "figures"
    #dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    
	###read in the first time point
    data = read.table(files[1], header = F, as.is = T, sep = ',')
    names(data) = c('MATa', 'MATalpha', 'double_barcode', strsplit(files[1], split = '_')[[1]][2])
    data$ID = paste(data$MATalpha, data$MATa, sep = '.')
    data = data[,c(1,2,3,5,4)]
    
	###read in subsequent time points, and merge barcode count number with the first time point for each diploid. This creates a data frame with the following columns: MATa strain, MATalpha strain, double barcode, diploid name, several columns with the read counts per time point
    for(file in files[-1]){
        print(file)
        merge = read.table(file, header = F, as.is = T, sep = ',')
        names(merge) = c('MATa', 'MATalpha', 'double_barcode', strsplit(file, split = '_')[[1]][2])
        merge$ID = paste(merge$MATalpha, merge$MATa, sep = '.')
        assign(strsplit(file, split = '_')[[1]][2], merge[, c('double_barcode', strsplit(file, split = '_')[[1]][2])])
        data = merge(data, get(strsplit(file, split = '_')[[1]][2]), by = 'double_barcode', all = T)
    }
    
    merged = data
    merged[is.na(merged)] = 0
    merged = merged[-which(merged$MATa == 0),]
    merged$sum = rowSums(merged[,5:ncol(merged)])

    ###Calculate the number of times each MATa barcode is present
    a = unique(merged$MATa)
    a_BC = data.frame(t(sapply(1:length(a), function(y){
        c(a[y], sum(merged[which(merged$MATa == a[y]), 'T3']), sum(merged[which(merged$MATa == a[y]), 'T6']), 
        sum(merged[which(merged$MATa == a[y]), 'T9']), sum(merged[which(merged$MATa == a[y]), 'T12']), 
        sum(merged[which(merged$MATa == a[y]), 'sum']))
    })), stringsAsFactors = FALSE) 
    
    for(x in 2:ncol(a_BC)){
        a_BC[,x] = as.numeric(as.character(a_BC[,x]))
    }
	
	###Extract double barcodes where the MATa barcode is from a control strain
    parent = c(which(grepl('unique', a_BC[, 1])), which(grepl('many', a_BC[, 1])))
    parental_a = a_BC[parent,]
    seg_a = a_BC[-parent,]

    ###Calculate the number of times each MATalpha barcode is present
    alpha = unique(merged$MATalpha)
    alpha_BC = data.frame(t(sapply(1:length(alpha), function(y){
        c(alpha[y], sum(merged[which(merged$MATalpha == alpha[y]), 'T3']), sum(merged[which(merged$MATalpha == alpha[y]), 'T6']), 
        sum(merged[which(merged$MATalpha == alpha[y]), 'T9']), sum(merged[which(merged$MATalpha == alpha[y]), 'T12']), 
        sum(merged[which(merged$MATalpha == alpha[y]), 'sum']))
    })), stringsAsFactors = FALSE) 
    
    for(x in 2:ncol(alpha_BC)){
        alpha_BC[,x] = as.numeric(as.character(alpha_BC[,x]))
    }
	###Extract double barcodes where the MATalpha barcode is from a control strain
    parent = c(which(grepl('unique', alpha_BC[, 1])), which(grepl('many', alpha_BC[, 1])))
    parental_alpha = alpha_BC[parent,]
    seg_alpha = alpha_BC[-parent,]
    
    ###For each double barcode, add columns showing the number of time the associated MATa barcode was present
    listy4 = lapply(1:nrow(merged), function(i){
        a_BC[which(a_BC[, 1] == merged[i, 'MATa']), 2:ncol(a_BC)] 
    })
    listy4 = do.call(rbind, listy4)
    merged = cbind(merged, listy4)
    colnames(merged)[ncol(merged)] = 'MATa_read_count'
	
	###For each double barcode, add columns showing the number of time the associated MATalpha barcode was present
    listy5 = lapply(1:nrow(merged), function(w){
        alpha_BC[which(alpha_BC[, 1] == merged[w, 'MATalpha']), 2:ncol(alpha_BC)] 
    })
    listy5 = do.call(rbind, listy5)
    merged = cbind(merged, listy5)
    colnames(merged)[ncol(merged)] = 'MATalpha_read_count'
	
	###For each double barcode, multiply the barcode counts of the MATa and MATalpha barcodes
    merged$readxread = merged$MATa_read_count * merged$MATalpha_read_count
	
	###Extract out double barcode pairs between a control strain and a segregant strain. These pairing should not exist as they were never mated to each other. Thus, these must be PCR chimeras.
    chimera = which(str_count(merged$ID, 'unique') + str_count(merged$ID, 'many') == 1)
    parental_chimera = merged[chimera,]
    
	###For double barcodes between a MATa control and MATalpha segregant, find the number of PCR chimera counts and the product of the number of MATa and MAT barcodes present. We expect more PCR chimeras involving barcodes that are present at a high frequency.
    no_read_chimera1 = expand.grid(seg_a[,1], parental_alpha[,1])
    no_read_chimera1$ID = paste(no_read_chimera1[,2], no_read_chimera1[,1], sep = '.')
    reads1 = t(sapply(1:nrow(no_read_chimera1), function(i){
       MATa = a_BC[which(a_BC[, 1] == no_read_chimera1[i, 1]),2]
       MATalpha = alpha_BC[which(alpha_BC[, 1] == no_read_chimera1[i, 2]),2]   
       c(MATa, MATalpha)
    }))
    no_read_chimera1 = cbind(no_read_chimera1, reads1)
    no_read_chimera1 = no_read_chimera1[-which(no_read_chimera1$ID %in% parental_chimera$ID),]
    no_read_chimera1$sum = 0
    names(no_read_chimera1) = c('MATa', 'MATalpha', 'ID', 'MATa_read_count', 'MATalpha_read_count', 'sum')
    no_read_chimera1$readxread = no_read_chimera1$MATa_read_count * no_read_chimera1$MATalpha_read_count
    
	###For double barcodes between a MATalpha control and MATa segregant, find the number of PCR chimera counts and the product of the number of MATa and MAT barcodes present
    no_read_chimera2 = expand.grid(parental_a[,1], seg_alpha[,1])
    no_read_chimera2$ID = paste(no_read_chimera2[,2], no_read_chimera2[,1], sep = '.')
    reads2 = t(sapply(1:nrow(no_read_chimera2), function(i){
       MATa = a_BC[which(a_BC[, 1] == no_read_chimera2[i, 1]),2]
       MATalpha = alpha_BC[which(alpha_BC[, 1] == no_read_chimera2[i, 2]),2]   
       c(MATa, MATalpha)
    }))
    no_read_chimera2 = cbind(no_read_chimera2, reads2)
    no_read_chimera2 = no_read_chimera2[-which(no_read_chimera2$ID %in% parental_chimera$ID),]
    no_read_chimera2$sum = 0
    names(no_read_chimera2) = c('MATa', 'MATalpha', 'ID', 'MATa_read_count', 'MATalpha_read_count', 'sum')
    no_read_chimera2$readxread = no_read_chimera2$MATa_read_count * no_read_chimera2$MATalpha_read_count
    
    total_chimera = data.frame(readxread = c(parental_chimera$readxread, no_read_chimera1$readxread, no_read_chimera2$readxread),
        sum = c(parental_chimera$sum, no_read_chimera1$sum, no_read_chimera2$sum))
    
    ###Calculate the relationship between the number of PCR chimeras relatice to the product of the number of MATa and MATalpha barcodes present using linear regression.
    model = lm(total_chimera$sum ~ total_chimera$readxread) 
    intercept = summary(model)[[4]][1,1]
    slope = summary(model)[[4]][2,1]

    plot(total_chimera$readxread, total_chimera$sum, xlab = 'Product of Parental BC and Segregant BC', 
            ylab = '# of PCR chimera', main = paste('CoCl2 - All time points', ': R-squared = ', 
            round(summary(model)$r.squared, digits = 3), ' Mean coverage = ', round(mean(merged$sum)),
            sep = ''), cex = 0.1)
    abline(model, col = 'red', lwd = 2)  
    #quartz.save('../figures/pcr chimera all time point.pdf', type = 'pdf')
    #dev.off()
    
    ###Using the linear regression model above, subtract out the expected number of PCR chimeras based on the frequency of the involved barcodes for each double barcode. Any double barcodes that become negative are set as 0.
    corr = merged
    cnames = c('t3_a', 't6_a', 't9_a', 't12_a', 't3_alpha', 't6_alpha', 't9_alpha', 't12_alpha')
    colnames(corr)[c(10:13, 15:18)] = cnames 
    product = sapply(1:length(files), function(x){
        corr[, cnames[x]] * corr[, cnames[x + 4]]
    })
    corr = cbind(corr, product)
    
    no_chimera = sapply(1:length(files), function(x){
        chi = intercept + slope * corr[, 8 + length(files) * 3 + x]
    })
    corr = cbind(corr, no_chimera)
    
    remove = sapply(1:length(files), function(x){
        chi = corr[, 4 + x] - corr[, 8 + length(files) * 4 + x]
        chi[chi < 0] = 0
        round(chi)
    })
    temp = cbind(corr, remove)
    save(temp, file = '../All_info_table_PCR_chimera.RData')
    table = cbind(corr[,c('double_barcode', 'MATa', 'MATalpha', 'ID')], remove)
    save(table, file = '../Combined_PCR_corrected_table_V2.RData')
    
    ###Take the PCR chimera corrected read count and format it for FitSeq. removed any double barcodes with an average read count of less than 4 per time point, as FitSeq cannot accurately estimate fitness for double barcodes with very low read counts. Also removed all PCR chimeras involving a control straina and a segregant.
    table$average = sapply(1:nrow(table), function(x){
        temp = table[x,]
        mean(unlist(temp[5:(4 + length(files))]))
    })
    good = table[which(table$average > 4),]
    good = good[-c(unique(which(str_count(good$ID, 'unique') == 1), which(str_count(good$ID, 'many') == 1))),]
    good$name = paste(good$double_barcode, '+', good$ID, sep = '')
    identifier = good[, 'name']
    abundance = good[, 5:(4 + length(files))]
    names(abundance) = NULL

    #write.table(abundance, file = '../Combined_PCR_corrected_fitseq_table.txt', quote = F, row.names = F, col.names = F, sep = ',')
    #write.table(identifier, file = '../Combined_PCR_corrected_fitseq_name.txt', quote = F, row.names = F, col.names = F, sep = ',')
}

###getting all files in directory
for(p in 1:8){
    d = paste('/Volumes/TM/All_segregants/fitness_assay_barcodes/', env[p], '/pre-process/both_matching/', sep = '')
    setwd(d) ##change directory to where the output files from 5_extracting_cluster_matching_a_barcode are kept for each environment.
    files = dir()
    files = mixedsort(files[grep('csv', files)])
    pcr_correction(files)
}

###############################################################################################################
###Running FitSeq on the output file Combined_PCR_corrected_fitseq_table: run in terminal using python3
cd /Volumes/TM/All_segregants/fitness_assay_barcodes/CoCl2/
python3 /Users/tmatsui2/PyFitSeq/pyfitseq.py -i Combined_PCR_corrected_fitseq_table.txt -t 0 3 6 9 -o output








