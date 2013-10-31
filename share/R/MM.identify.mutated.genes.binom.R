###############################################################################
# Identification of significantly mutated genes
# With binomial test
#
# Magali Michaut, m.michaut@nki.nl
# NKI, 23/01/2013
#
###############################################################################

prepare.data.for.binomial.test <- function(m,bp.covered){
  # get mutation counts
  counts <- data.frame(
    gene = colnames(m),
    count = apply(m,2,sum),
    stringsAsFactors = F
    );
  
  # add bp covered
  df <- merge(counts,bp.covered);
  
  # add mutation rates
  nb.samples <- 1;#nrow(m)# version used so far does not use the nb of samples
  df$mutation.rate <- df$count / (df$bp.covered  * nb.samples);
  global.mutation.rate <- sum(df$count) / ( sum(df$bp.covered) * nb.samples);
  
  data <- data.frame(
    #eg = df$gene,
    gene = df$gene,
    x = df$count,
    n = df$bp.covered,
    p = global.mutation.rate,
    stringsAsFactors = F
  );
  
  data
}

# data isa data.frame with a test on each line (x.n.p)
compute.binom.signif <- function(data){
  results <- data.frame();
  for (i in 1:nrow(data)){
    test <- binom.test(
      x = data$x[[i]],
      n = data$n[[i]],
      p = data$p[[i]],
      alternative = 'greater'
    );
    res <- data.frame(
      #eg = data$gene[[i]],
      gene = data$gene[[i]],
      test = test$method,
      x = data$x[[i]],
      n = data$n[[i]],
      p = data$p[[i]],
      alternative = test$alternative,
      p.value = test$p.value,
      stringsAsFactors = F
    );
    results <- rbind(results,res);
  }
  results$q <- p.adjust(results$p.value,method='fdr');
  results <- results[order(results$p.value),];
  results
}

run.identification.signif.genes <- function(m,bp.covered){
 data <- prepare.data.for.binomial.test(m,bp.covered);
 res <- compute.binom.signif(data);
}

get.sizes.default <- function(){
  library("org.Hs.eg.db");
  x <- org.Hs.egCHRLOC;
  eg2start <- as.list(x[mappedkeys(x)]);
  x <- org.Hs.egCHRLOCEND;
  eg2end <- as.list(x[mappedkeys(x)]);
  
  genes <- sort(intersect(names(eg2start),names(eg2end)));
  
  data <- data.frame(
    gene = genes,
    start = sapply(genes,function(g) min(eg2start[[g]])),
    end = sapply(genes,function(g) min(eg2end[[g]])),
    stringsAsFactors = F
    );
  
  data$bp.covered <- abs(data$end - data$start) + 1;
  data
}

# bed is a tab delimited file with required column names
# - chr
# - start
# - end
# - gene
# extension gives the nb of base pairs to add around the probe definition
# if you want to add some buffer around
get.sizes.from.BED <- function(bedfile,extension=0){
  bed <- read.delim(bedfile,as.is=T);
  
  library('intervals');
  
  # get all baits associated to a given gene
  genes <- sort(unique(bed$gene));
  
  baits <- lapply(genes,function(g){
    df <- subset(bed,gene==g);
    ch <- unique(df$chr);
    stopifnot( length(ch)==1 );
    df
  });
  
  # compute the number of bp covered
  bp.cov <- lapply(baits,function(df){
    df$start <- df$start - extension;
    df$start[df$start<0] <- 0;
    df$end <- df$end + extension;
    m <- as.matrix(df[,c('start','end')]);
    intervals <- Intervals(m);
    region <- reduce(intervals);
    size <- size(region);
    bp.covered <- sum(size);
    bp.covered
  });
  
  res <- data.frame(
    gene = genes,
    bp.covered = unlist(bp.cov),
    stringsAsFactors = F
  );
  res
}

example <- function(){
  # random mutation matrix
  m <- matrix(0,ncol=60,nrow=10);
  rownames(m) <- paste('Sample',1:nrow(m),sep='');
  colnames(m) <- paste('Gene',1:ncol(m),sep='');
  m[sample.int(nrow(m),size=4),sample.int(ncol(m),size=30)] <- 1;
  
  # random gene sizes
  # you can use the function get.sizes.default for gene size
  bp.covered <- data.frame(
    gene = paste('Gene',1:ncol(m),sep=''),
    bp.covered = ceiling(rnorm(ncol(m),mean=1000,sd=10)),
    stringsAsFactors = F
    );
  
  # run function
  res <- run.identification.signif.genes(m,bp.covered);
}