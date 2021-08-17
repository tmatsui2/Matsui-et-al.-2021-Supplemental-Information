###Takes the output files from Samtools_fastq_to_mpileup.py to create a HMM corrected genotype table for each strain. The files should have 7 columns, where the columns correspond to
#chromosome, position, coverage, BY allele, number of BY allele. 3S allele, and number of 3S alleles.
########################################################################################################################################################################
###Getting files  
library("HMM")
setwd("~/allsnps/test/")
geno = read.table("ctk1_3_H9allpos.txt", header = F, as.is = T)
geno[,1] = as.numeric(geno[,1])
geno[,3] = (geno[,5] + geno[,7])

c_l = sapply(1:16, function(i) {max(geno[,2][geno[,1] == i])})
g_l = c(0, cumsum(c_l))

files = dir()
letters = c("txt")
headers = c(paste(letters, sep = ""))
files = files[unlist(lapply(headers, function(i) {which(grepl(i, files))}))]
files2 = gsub('.{4}$', '', files)

###Finding files with low coverage (<1 calls/base)
smallsnps = lapply(files, function(g) {
	data = read.table(g, header = F, as.is = T)
	data[,1] = as.numeric(data[,1])
	data$gp = data[,2] + g_l[data[,1]]
	data = data[order(data$gp),]
	cov =  sum(data[,3]/nrow(data))
	cov1 = (nrow(data[data[,3] == 0,])/nrow(data))*100
	cov3 = g[cov < 1] #average coverage less than 1
	cov6 = as.data.frame(matrix(data = c(cov3)))
})
	
cov7 = as.data.frame(do.call("rbind", lapply(smallsnps, "[",)))
	
###For finding average coverage of segsnps files

avgcov = lapply(files, function(g) {
	data = read.table(g, header = F, as.is = T)
	data[,1] = as.numeric(data[,1])
	data$gp = data[,2] + g_l[data[,1]]
	data = data[order(data$gp),]
	covo = unique(data[,1] > 0)
	cov = sum(data[,3])/length(data[,3])
	cov1 = as.data.frame(matrix(data= cov)
)})

cov2 = as.data.frame(do.call("rbind", lapply(avgcov, "[",)))	
cov3 = mean(cov2[,1])

###For extracting rows where the number of counts for both alleles are the same
z = lapply(files, function(g) {
	print(g)
	data = read.table(g, header = F, as.is = T)
	data[,1] = as.numeric(data[,1])
	data$gp = data[,2] + g_l[data[,1]]
	data = data[order(data$gp),]
 	data = data[data[,5] == data[,7],]	
	data = data[(data[,4] == 0) == FALSE,]
	data = data[order(data$gp),]
})
 			
gp = "gp"
l = as.matrix(unlist(lapply(z, '[', gp )))
n = as.vector(l[duplicated(l[1:nrow(l),]) == FALSE,])

###Function for running HMM on genotype calls
flipstate = function(t) {
	t = split(t, t[,1])
	t = lapply(t, function(d) {
		nas = which(is.na(d$h))
		newvals = sapply(nas, function(x) {
			up = na.omit(rev(d$h[1:(x-1)]))[1]
			down = na.omit(d$h[(x+1):length(d$h)])[1]
			if (x != 1 & x != length(d$h)) {
				if (is.na(up) == FALSE & is.na(down) == FALSE) {
					if (up == down) {return(up)}
					else {return(NA)}					
				} else if (is.na(up) == FALSE & is.na(down) == TRUE) {
					return(up)
				} else if (is.na(up) == TRUE & is.na(down) == FALSE) {
					return(down)
				}
			} else if (x == 1) {
				return(na.omit(d$h[(x+1):length(d$h)])[1])
			} else if (x == length(d$h)) {
				return(na.omit(rev(d$h[1:(x-1)]))[1])
			}		
		})
		if (length(nas) > 0) {d$h[nas] = newvals} 
		return(d)
	})
	t = do.call("rbind", t)
	return(t)
}


callgeno2 = function(t) {
	na.data = t[is.na(t$sign),]   	
	t = na.omit(t)
	model = initHMM(as.character(0:1), as.character(0:1), c(.25,.75), matrix(c(.999,.001,.001,.999),2),matrix(c(.25,.75,.75,.25),2))
	t = split(t, t[,1])
	t = lapply(t, function(d) {
		d$h = as.numeric(viterbi(model, as.character(d$sign)))
		d
	})
	t= do.call("rbind", t)
	if(nrow(na.data) > 0) {
		na.data$h = NA
		t = rbind(t, na.data)
	}
	t = t[order(t$gp),]
	t = flipstate(t)
	return(t)
}

###Running HMM on files

w = lapply(files, function(f) {
	data = read.table(f, header = F, as.is = T)
	data[,1] = as.numeric(data[,1]) 
	data$gp = data[,2] + g_l[data[,1]]
	data = data[order(data$gp),]
	data[,3] = (data[,5] + data[,7])   
	data = data[(data$gp %in% n) == FALSE,]
	data$f = data[,5]/data[,3]
	data = data[order(data$gp),]
	data$sign = as.numeric(data$f >= 0.5)  
	t= callgeno2(data)
}) 

###Organizing output into a single data frame
results3 = cbind(cbind(w[[1]][,1:2], w[[1]]$gp), data.frame(matrix(unlist(lapply(w, function(w) {w$h})), ncol= length(w))))
results3 = na.omit(results3)
results3$f = rowSums(results3[,4:ncol(results3)])/(ncol(results3)-3)
names(results3)[1:3] = c("c", "p", "gp")

### Parsing results into unique genotypes
results4 = results3[duplicated(results3[,4:(ncol(results3)-1)]) == FALSE,]


