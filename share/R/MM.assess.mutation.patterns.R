###############################################################################
# Script to identify ME
#
# Magali Michaut, m.michaut@nki.nl
# NKI, 05/06/2012
#
###############################################################################

is.binary <- function(m){
  length(table(m))==2 & all(c(0,1) %in% m);
}

# a and b are binary vectors of the same size
# and ordered the same way
contingency.table <- function(a,b){
  conting <- matrix(
    c(
      sum(a & b), sum(a & !b), sum(!a & b), sum(!a & !b)
    ),
    ncol = 2
  );
  rownames(conting) <- paste('b',c('Mut','WT'),sep='.');
  colnames(conting) <- paste('a',c('Mut','WT'),sep='.');
  conting
}


# m should be a binary matrix
# in rows the samples
# in cols the genes
# 1/0 if event or not
assess.me.event.matrix <- function(m,max.genes=100,alternative='less'){
  stopifnot(is.binary(m));
  
  if (ncol(m)>max.genes){
    cat('The matrix has',ncol(m),'genes, more than the max',max.genes,'\n');
    counts <- apply(m,2,sum);
    counts <- sort(counts,decreasing=T);
    keep <- names(counts)[1:max.genes];
    cat('We keep the',max.genes,'most frequently altered genes\n');
    m <- m[,keep];
  }
  
  res <- data.frame();
  n <- ncol(m);
  for (i in 1:(n-1)){
    va <- m[,i];
    for (j in (i+1):n){
      vb <- m[,j];
      conting <- contingency.table(va,vb);
      test.fisher <- fisher.test(conting,alternative=alternative);
      test.chisq <- chisq.test(conting);
      type <- ifelse(test.fisher$estimate>1,'CO','ME');
      
      df <- data.frame(
        a = colnames(m)[[i]],
        b = colnames(m)[[j]],
        samples = nrow(m),
        a.mut = sum(va),
        b.mut = sum(vb),
        ab.mut = sum(va&vb),
        odds = test.fisher$estimate,
        fisher.p = test.fisher$p.value,
        type = type,
        fisher.alternative = test.fisher$alternative,
        chisq.estimate = test.chisq$statistic,
        chisq.p = test.chisq$p.value,
        stringsAsFactors = F
      );
      
      res <- rbind(res,df);
    }
  }
  res$fisher.q <- p.adjust(res$fisher.p,method='fdr');
  res <- res[order(res$fisher.p),];
  res
}