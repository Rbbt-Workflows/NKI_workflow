###############################################################################
# Study mutation patterns in genes 
# e.g. Identification of hotspot variants
#
# Magali Michaut, m.michaut@nki.nl
# NKI, 23/01/2013
#
###############################################################################

# variants is a data.frame with all variants
# required columns are:
# - gene: in which gene is the variant
# - pos: position of the variant
# - sample: in which sample is this variant present
get.variant.pattern <- function(variants){
  genes <- sort(unique(variants$gene));
  
  mutated.samples <- sapply(genes,function(g){
    s <- subset(variants,gene==g);
    samples <- sort(unique(s$sample));
    samples.freq <- sapply(samples,function(sa){
      nrow(subset(s,sample==sa));
    });
    names(samples.freq) <- samples;
    samples.freq
  });
  
  nb.samples <- sapply(mutated.samples,length);
  sample.freq.mean <- sapply(mutated.samples,mean);
  
  mutated.sites <- sapply(genes, function(g){
    s <- subset(variants,gene==g);
    sites <- sort(unique(s$pos));
    sites.freq <- sapply(sites,function(si){
      nrow(subset(s,pos==si));
    });
    names(sites.freq) <- sites;
    sites.freq
  });
  
  nb.variants <- sapply(mutated.sites,sum);
  nb.sites <- sapply(mutated.sites,length);
  site.freq.mean <- sapply(mutated.sites,mean);
  
  pattern <- data.frame(
    gene = genes,
    nb.samples = nb.samples,
    sample.freq.mean = sample.freq.mean,
    nb.variants = nb.variants,
    nb.sites = nb.sites,
    site.freq.mean = site.freq.mean,
    stringsAsFactors = F
  );
  
  pattern <- pattern[order(pattern$nb.variants,decreasing=T),];
  
  pattern$sites.variants.ratio <- pattern$nb.sites / pattern$nb.variants;
  pattern
}

example <- function(){
  nb.genes <- 100;
  # we assume they all have the same size
  size <- 15;
  # how many variants has each gene
  max.nb.variants <- 20;
  nb.variants <- sample.int(n=max.nb.variants,size=nb.genes,replace=T);
  # what are the positions of the variants
  positions <- sapply(nb.variants,function(x){
    sample.int(size,size=x,replace=T);
  });
  gene.index <- rep(1:nb.genes,sapply(positions,length));
  # samples
  nb.samples <- 10;
  samples <- sample.int(nb.samples,size=length(gene.index),replace=T);
  samples <- paste('Sample',samples);
  variants <- data.frame(
    gene = paste('Gene',gene.index),
    pos = unlist(positions),
    sample = samples,
    stringsAsFactors = F
    );
  
  pattern <- get.variant.pattern(variants);
  pattern
}
