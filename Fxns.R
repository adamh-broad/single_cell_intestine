
library(RColorBrewer)
brewer16 = c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
brewer16[6] = "khaki2"
brewer16[8] = "lightskyblue2"

### Compute the group-wise mean of a dataset.
group.means <- function(counts, groups, fn=mean, use.data.table=F)
{
	counts <- aggregate(t(counts), by=list(groups), FUN=fn)
	rownames(counts) = counts$Group.1
	counts$Group.1 = NULL
	r = t(counts)
	return(r)
}

# Logging utility function
info <- function(text, ...)
{
	cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

# Logging utility function
warn <- function(text, ...)
{
	cat(sprintf(paste(Sys.time(),"WARN:", text,"\n")))
}

### Compute TPM expression values from raw UMI counts
tpm <- function(counts, mult=10000)
{
	info("Running TPM normalisation")
	total.counts = colSums(counts)
	scaled.counts = t(t(counts) / total.counts) 
	scaled.counts * mult
}


### Run ComBat batch correction from the SVA package
batch.normalise.comBat <- function(counts, batch.groups, max.val=6)
{
	library(sva)
    batch.groups = factor(batch.groups) ## drop zero levels
    batch.id = 1:length(unique(batch.groups))
    names(batch.id) = unique(batch.groups)
    batch.ids = batch.id[batch.groups]
    correct.data = ComBat(counts,batch.ids, prior.plots=FALSE, par.prior=TRUE)
    correct.data[correct.data > max.val] = max.val
    as.data.frame(correct.data)
}


### Get variable genes. Code adapted from:
### | Brennecke et al, Accounting for technical noise in single-cell RNA-seq experiments
### | Nature Methods 10, 1093–1095 (2013), doi:10.1038/nmeth.2645
### 	See: https://images.nature.com/original/nature-assets/nmeth/journal/v10/n11/extref/nmeth.2645-S2.pdf
### 	and: http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
get.variable.genes <- function(ed, min.cv2=2, pdf=NULL, width=9, height=8, do.plot=T, p.thresh=0.05)
{
	library(statmod)
	means <- rowMeans(ed)
	vars <- apply(ed,1,var)
	cv2 <- vars/means^2
	minMeanForFit <- unname( quantile( means[ which( cv2 > min.cv2 ) ], .95 ) )
	useForFit <- means >= minMeanForFit # & spikeins
	info(sprintf("Fitting only the %s genes with mean expression > %s", sum(useForFit), minMeanForFit))
	fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ), cv2[useForFit] )
	a0 <- unname( fit$coefficients["a0"] )
	a1 <- unname( fit$coefficients["a1tilde"])
	if(do.plot){par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2))}
	xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
	vfit <- a1/xg + a0
	if(do.plot){lines( log(xg), log(vfit), col="black", lwd=3 )}
	
	df <- ncol(ed) - 1
	# add confidence interval
	if(do.plot){
		lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
		lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")
	}
	afit <- a1/means+a0
	varFitRatio <- vars/(afit*means^2)
	varorder <- order(varFitRatio, decreasing=T)
	oed <- ed[varorder,]
	pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
	adj.pval <- p.adjust(pval,"fdr")
	r = data.frame(rownames(ed), varFitRatio, pval, adj.pval)
	colnames(r) = c("Gene", "VarianceFitRatio", "p", "p.adj")
	v = r[!is.na(r$p.adj),]
	n.sig = sum(v$p.adj<p.thresh)
  	info(sprintf("Found %s variable genes (p<0.05)", n.sig))

	# add top 100 genes
	if(do.plot){
		points(log(means[varorder[1:n.sig]]),log(cv2[varorder[1:n.sig]]),col=2)
	}
	r = r[order(r$VarianceFitRatio, decreasing=T), ]
	r$Rank = 1:nrow(r)
	return(r)
}


# Test for significant PCs adapted from: 
#	' Permutation Parallel Analysis
#	'
#	' Estimate a number of significant principal components from a permutation test
#   B is the number of permutations
#   threshold is p-value for significance
#'
sig.pcs.perm <- function (dat, B = 100, threshold = 0.05,
                                        randomized=F,  
                                        verbose=TRUE, seed = NULL,
                                        max.pc=100, n.cores=1, 
                                        center=T, scale=T) {
    ptm <- proc.time()
    if(B %% n.cores != 0){stop("Permutations must be an integer multiple of n.cores")}
    cat(sprintf("Scaling input matrix [center=%s, scale=%s]\n", center, scale))
    dat = t(dat)
    dat = as.matrix(t(scale(t(dat), center=center, scale=scale)))
    if (!is.null(seed)) set.seed(seed)
    n <- min(max.pc, ncol(dat))
    m <- nrow(dat)
    print(paste0("Considering only the top ", n, " PCs. Supply max.pc if you wish to change"))
    cat(sprintf("Running initial PCA\n"))
    if(randomized){
        library(rsvd)
        uu <- rsvd(as.matrix(dat), k=max.pc)
    }else{
        uu <- corpcor::fast.svd(dat, tol = 0)
    }
    ndf <- n - 1
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    dstat0 <- matrix(0, nrow = B, ncol = ndf)
    if(verbose==TRUE) message("Estimating number of significant principal components. Permutation: ")
    #permutations
    if(n.cores==1){
        for (i in 1:B) {
            if(verbose==TRUE) cat(paste(i," "))
            dat0 <- t(apply(dat, 1, sample, replace = FALSE))
            if(randomized){
                library(rsvd)
                uu0 <- rsvd(as.matrix(dat0), k=max.pc)
            }else{
                uu0 <- corpcor::fast.svd(dat0, tol = 0)
            }
            dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
        }
    }else{
        library(parallel)
        library(foreach)
        library(doParallel)
        cl<-makePSOCKcluster(n.cores, outfile="")
        registerDoParallel(cl, n.cores)
        chunksize = B/n.cores
        vals = split(1:B, ceiling(seq_along(1:B)/chunksize))
        dstat0 = foreach(run.id=1:n.cores, .packages="corpcor", .combine=cbind) %dopar% {
            v = vals[[run.id]]
            #cat(sprintf("Core %s will run perms: %s \n", run.id, paste(v, collapse=",")))
            do.call(rbind, lapply(v, function(i) {
                if(verbose==TRUE) cat(paste(i," "))
                dat0 <- t(apply(dat, 1, sample, replace = FALSE))
                
                if(randomized){
                    library(rsvd)
                    uu0 <- rsvd(as.matrix(dat0), k=max.pc)
                }else{
                    uu0 <- corpcor::fast.svd(dat0, tol = 0)
                }
                uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
            }))
        }
        cat("\nUnregistering parallel backend..")
        stopCluster(cl)
        registerDoSEQ()
        cat(" done\n");
    }
    p <- rep(1, n)
    for (i in 1:ndf) {
      p[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:ndf) {
      p[i] <- max(p[(i - 1)], p[i])
    }
    r <- sum(p <= threshold)
    y = proc.time() - ptm
    cat(sprintf("\n\n PC permutation test completed. \n %s PCS significant (p<%s, %s bootstraps)\n Runtime: %s s\n ", r,  threshold, B,signif(y[["elapsed"]], 3)))
    return(list(r = r, p = p))
}


plot_graph <- function(g, 
	node.names=NULL, 
	node.clustering=NULL, 
	node.label.size=2, 
	use.cols=NULL,
	layout=NULL,
	layout.type=NULL,
	grid=F, 
	simplify=F, 
	community=F,
	node.size=500,
	max.edges.for.plot=100000,
	edge.width=1)
{
	if(is.null(node.names))
	{
		#node.names = c(rep("", vcount(g)))
		node.names = V(g)$name
	}
	V(g)$label.cex = node.label.size
	
	if(!is.null(node.clustering))
	{
		if(length(node.clustering) != vcount(g))
		{
			error("Length of cluster vector must match number of nodes in graph!")
			return (FALSE)
		}
		n.colors = length(unique(node.clustering))
		if(is.null(use.cols))
		{
			use.cols = intense.cols(n.colors)	
		}
		if(community) # simplified graph showing community structure only
		{
			g <- contract.vertices(g, node.clustering)
			g <- igraph::simplify(g, remove.loops=TRUE)
			V(g)$color = use.cols
		}else
		{
			V(g)$color = mapvalues(node.clustering, unique(node.clustering), use.cols)
		}
		
		#E(g)$weight <- 1
		
	}else
	{
		if(community)
		{
			error("Cannot draw a community graph without being passed a graph clustering!")
			return (FALSE)
		}
	}

	if(simplify)
	{
		info("Simplifying graph..")
		g = igraph::simplify(g)
	}

	if(ecount(g) > max.edges.for.plot)
	{
		#TODO: should trim graph to the maximum edge num.
		warn(sprintf("Graph contains more than %s edges. Trimming some..",max.edges.for.plot))
		before = ecount(g)
		g = delete.edges(g, which(E(g)$weight <=.35))
		after = ecount(g)
		info(sprintf("Removed %s edges for plotting [%s remain]..", before-after, after))
	}

	if(is.null(layout))
	{
		info(sprintf("Laying out graph [%s nodes]..", vcount(g)))
		if(is.null(layout.type))
		{
			info("Using autolayout. Try: fruchterman.reingold, reingoild.tilford, kamada, or spring for better results")
			l = layout.auto(g)
		}else{
			if(layout.type=="fruchterman.reingold"){l = layout.fruchterman.reingold(g)}
			else{if(layout.type=="reingold.tilford"){
					g = simplify(g)
					l <- layout.reingold.tilford(g, circular=T)}else{
					if(layout.type=="kamada"){l = layout.kamada.kawai(g)}else{
						if(layout.type == "spring"){l = layout.spring(g)}else{
							error("unknown layout type")
						}
					}
				}
			}
		}
	}else{
		l = layout
	}
	
	info("Plotting..")
	plot(g, layout=l, vertex.size=node.size, edge.width=edge.width, edge.arrow.size=0,
			vertex.label=node.names, rescale=FALSE, xlim=range(l[,1]), ylim=range(l[,2]),
		 	vertex.label.dist=1)
	return(l)
}



build_knn_graph <- function(dm, k=200, verbose=F)
{
	library(cccd) 
	if(k==0)
	{
		k = floor(sqrt(nrow(dm))/2)
	}
	if(verbose)
	{
		info(sprintf("Building %s-nearest [%s] neighbor graph..", k, dist.type))
	}
	g <- nng(dx=dm,k=k)
	V(g)$name = rownames(dm)
	if(verbose)
	{
		info(sprintf("%s %s-NN computed. Average degree: %s", dist.type, k, mean(degree(g))))
	}
	return(g)
}



# graph.type can be jaccard, invlogweighted or dice, community detect
# can be louvain, infomap or markov. 
cluster_graph <- function(	g, 
							graph.type="knn", # can be threshold (binarise the distance matrix), jaccard or knn.
							dm=NULL,
							community.detect="infomap", 
							distance.method="euclidean",
							k=0)
{
	if(identical(toupper(community.detect), toupper("markov")))
	{
		r = igraph::cluster.markov(g)
		clusters = r$Cluster
		
	}else{
		if(identical(toupper(community.detect), toupper("louvain")))
		{
			r = igraph::multilevel.community(as.undirected(g))
			clusters = r$membership
		}else{
			if(identical(toupper(community.detect), toupper("infomap")))
			{
				r = igraph::infomap.community(g, modularity=TRUE)
				clusters = r$membership
			}else{
				error(sprintf("Unknown community detection method: %s", community.detect))
				return (FALSE)
			}
		}
	}
	n.clusters =length(unique(clusters))
	
	f = function(i){as.vector(clusters==i)}
	clist= lapply(1:n.clusters, f)
	m = igraph::modularity(g, clusters)
	return (list("result"=r,
		"clustermethod"=paste(graph.type, "-graph clustering [", community.detect,"]", sep=""), 
		"nc"=n.clusters, 
		"modularity"=m, 
		"clusterlist"=clist,		
		"partition"=clusters))
}


merge_clusters <- function(clustering, clusters.to.merge, new.name=NULL)
{
	if(length(clustering) < 2)
	{
		cat("ERROR: Must provide 2 or more cluster ID's to merge!")
		return (clusterings)
	}

	i = 1
	if(!is.null(new.name)){
		use.id = new.name
		levels(clustering) = c(levels(clustering), use.id)
		clustering[which(clustering == clusters.to.merge[1])] = use.id
	}else
	{
		use.id = clusters.to.merge[1]
	}
	for(id in clusters.to.merge)
	{
		if(i > 1)
		{
			cat(sprintf("Merging cluster %s into %s ..\n", id, use.id))
			clustering[which(clustering == id)] = use.id
		}
		i = i + 1 
		
	} 
	return (factor(clustering))
}






