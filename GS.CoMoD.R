
# enrich.test: performs over-representation test with Fisher's exact test
# x: disease gene vector, bg: reference gene vector, gs: a functional gene set, 
# pval.only: if true, the function returns p value; otherwise, returns p value and odds ratio. 

enrich.test <- function(x, bg, gs, pval.only = T){	
	a1 <- length(intersect(x, gs))
	a2 <- length(x) - a1
	a3 <- length(intersect(setdiff(bg, x), gs))
	a4 <- length(setdiff(bg, x)) - a3
	tmp <- fisher.test(matrix(c(a1, a2, a3, a4), nrow=2), alternative='greater')
	if(pval.only) {
		return(c(tmp[[1]]))
	} else {
		return(c(tmp[[3]], tmp[[1]]))
	}
}


# enrich.vector: returns a p value vector
# x: disease gene vector, bg: reference gene vector, gs.list: a list containing functional gene sets

enrich.vector <- function(x, bg, gs.list){
	tmp.re <- (unlist(lapply(gs.list, enrich.test, x=x, bg=bg)))
	names(tmp.re) <- names(gs.list)
	return(tmp.re)
}



# gs.comod.perm: performs simulation analysis with sampled disease genes
# x: disease gene vector hg.ref: human reference gene vector, gs.list: a list containing functional gene sets
# proportion: a parameter for splitting ratio. For example, 0.3 divide total disease genes into two groups, 
# and the number of the groups have a ratio of 3:7.

gs.comod.perm <- function(x, hg.ref, gs.list, proportion = 0.5){
	temp.dg1 <- sample(x, proportion*length(x))
	temp.dg2 <- setdiff(x, temp.dg1)
	pvec1 <- enrich.vector(temp.dg1, hg.ref, gs.list)
	pvec2 <- enrich.vector(temp.dg2, hg.ref, gs.list) 
	re <- cor(exp(pvec1), exp(pvec2), use='complete.obs')
	return(re)
}


#gs.comod: returns GS.CoM (gene-set-based comorbidity) score
#x, y: disease gene vectors, bg: reference gene vector, gs.list: a list having functional gene sets. 

gs.comod <- function(x, y, bg, gs.list){
	px <- enrich.vector(x, bg, gs.list)
	py <- enrich.vector(y, bg, gs.list)
	re <- cor(exp(px), exp(py))
	return(re)
}




1. Estimating GS.CoM with gs.comod function
# disA: a vector having genes related to disease A, disB: a vector having genes related to disease B
# hg.ref: reference gene vector containing human genes. 

	re <- gs.comod(disA, disB, hg.ref, gobp)


2. Estimating GS.CoM (batch mode)

# When multiple disease gene sets are applied, gs.comod function is not relevant. 
# After obtaining p value vector for each disease gene set, GS.CoM scores can be easily determined. 

disease.gene.pvec <- list()
for(i in 1:length(disease.gene.list)){
	disease.gene.pvec[[i]] <- enrich.vector(disease.gene.list[[i]], bg, gs.list)
}

disease.gene.pmat <- as.data.frame(disease.gene.pvec)
gs.com.mat <- cor(exp(disease.gene.pmat))

# gs.com.mat: matrix having all disease pair-wise GS.CoM scores
# Using mat2vec function, GS.CoM score vector for all disease pairs is obtained. 
# Row and column names should be filled in the gs.com.mat matrix. 

mat2vec <- function(x){
	n.col <- ncol(x)
	ind <- 1
	col.names <- colnames(x)
	pnames <- list()
	re <- list()

	for(i in 1:(n.col-1)){
		for(j in (i+1):n.col){
			re[[ind]] <- x[i, j] 
			pnames[[ind]] <- paste(col.names[i], col.names[j], sep='_')
			ind <- ind+1
		}
	}

	names(re) <- unlist(pnames)	
	return(unlist(re))
}


re <- mat2vec(gs.com.mat) 



3.  Simulation analysis 

# Simulation analysis with splitting known disease genes
# dm.genes: type 2 diabetes genes from DisGeNET database.
# In this case, proportion is 0.3, which indicates that the dm.genes are divided 
# into two groups, and the ration of gene numbers in two groups is 3:7. 

re <- gs.com.perm(dm.genes, hg.ref, gobp, proportion = 0.3)


4. Logistic regression analysis. 

# Using GS.CoM score, logistic models are determined according to area under curvie (AUC) value. 

library(pROC)
# ind can be replaced by rr.ind5 (RR >5), and rr.ind10 (RR>10) in each dataset (menche data and rubio data). 

ind <- rr.ind
m <- glm(ind~ gobp.gs.com+gomf.gs.com, family='binomial', data=menche.data)
roc(ind, predict(m, type='response'))
m <- glm(ind~ gobp.gs.com+reactome.gs.com, family='binomial', data=menche.data)
roc(ind, predict(m, type='response'))
m <- glm(ind~ gomf.gs.com+reactome.gs.com, family='binomial', data=menche.data)
roc(ind, predict(m, type='response'))
m <- glm(ind~ gobp.gs.com+gomf.gs.com+reactome.gs.com, family='binomial', data=menche.data)
roc(ind, predict(m, type='response'))


5. Identification of core subsets



library(flexmix)

# core.subset returns an indicator for functional gene sets having higher GS.CoM score. 
# disA: disease A genes, disB: disease B genes, bg: reference gene vector, gs.list: a list having functional gene sets. 
core.subset <- function(disA, disB, bg, gs.list){
	exp.px <- exp(enrich.vector(disA, bg, gs.list))
	exp.py <- exp(enrich.vector(disB, bg, gs.list))
	mix.re <- flexmix(exp.py ~ exp.px-1, k = 2, model = FLXMRglm(family = "gaussian"))
	ind <- clusters(mix.re)
	re1 <- cor(exp.px[ind==1], exp.py[ind==1])
	re2 <- cor(exp.px[ind==1], exp.py[ind==1])
	if(re1 > re2) return(which(ind==1)) else return(which(ind==2))
}

re <- core.subset(dm.genes, htn.genes, hg.ref, gobp)
cor(exp.pa[re], exp.pb[re])






	
	










