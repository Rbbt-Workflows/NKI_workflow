###############################################################################
# Study gene mutations in a dataset
# Generate visual summary
#
# Magali Michaut, m.michaut@nki.nl
# NKI, 12/02/2013
#
###############################################################################

FONT <- 'Times'; #'Helvetica' 'Courier'
BLACK <- rgb(0,0,0,0.99);

legend.col <- function(col, lev){
  
  opar <- par
  
  n <- length(col)
  
  bx <- par("usr")
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}

replace.na <- function(v,by=0){
  v[is.na(v)] <- by;
  v
}

get.color.variant.ratio <- function(variant){
  cols <- colorpanel(100,'red','#DDDDDD');
  #values <- variant[variant>0];
  #start <- min(values);
  #indices <- 1+ceiling((values-start) * length(cols));
  indices <- ceiling(variant * length(cols));
  assigned <- rep(NA,length(indices));
  assigned[variant>0] <- cols[indices];
  assigned[variant==0] <- 'grey';
  cuts <- cut(variant,length(cols),label=F);
  #lev <- info$sites.variants.ratio[cuts] + start;
  #lev <- variant[cuts] + start;
  lev <- variant[cuts];
  res <- list(cols=cols,assigned=assigned,lev=lev);
  res
}

# m: mutation matrix to compute mutation frequencies
# bin: results of the binomial test
# patterns: results of the mutation pattern in gene
get.data.info <- function(m,bin,patterns){
  count <- apply(m,2,sum);
  freq <- count / nrow(m);
  df <- data.frame(
    gene = colnames(m),
    mut.count = count,
    mut.freq = freq,
    stringsAsFactors = F
  );
  
  # add results binomial
  df <- merge(df,bin,all=T);
  # add info mutation patterns in gene
  df <- merge(df,patterns,all=T);
  
  df$mut.count <- replace.na(df$mut.count);
  df$mut.freq <- replace.na(df$mut.freq);
  df$nb.samples <- replace.na(df$nb.samples);
  df$sample.freq.mean <- replace.na(df$sample.freq.mean);
  df$nb.variants <- replace.na(df$nb.variants);
  df$nb.sites <- replace.na(df$nb.sites);
  df$site.freq.mean <- replace.na(df$site.freq.mean);
  df$sites.variants.ratio <- replace.na(df$sites.variants.ratio);
  
  stopifnot(all(df$mut.count==df$x));
  stopifnot(all(df$mut.count==df$nb.samples));
  
  df
}

# info is a data.frame with required columns
# - mut.freq
# - score (e.g. from binomial test)
# - n
# - gene
# - sites.variants.ratio
# - p.value
create.mutation.bubble.plot <- function(info,use.names,title = "Study",text.size=3){
  freq <- info$mut.freq;
  score <- - log10(info$p.value);
  coverage <- info$n;
  label <- info$gene;
  
  pal <- get.color.variant.ratio(info$sites.variants.ratio);
  # remove Inf to take a numeric max
  finite <- score[is.finite(score)];
  score <- score + max(finite)/10;
  symbols(
    x = log10(coverage), xlab='Gene coverage',
    y = freq, ylab = 'Mutation frequency',
    ylim = c(0,max(freq)+0.15),
    circles = score, inch=1,
    cex.axis = 1.4,cex.lab = 1.4,
    bty = 'n', lwd = 2,las=1,
    fg = pal$assigned,
    main = paste('Mutated genes in',title)
  );
    #bg = pal$assigned, fg = BLACK,
  
  if (use.names){
    keep <- info$p.value<0.05 | info$mut.freq>0.1;
    label[!keep] <- '';
    text(
      x = log10(coverage),
      y = jitter(freq),
      labels = label,
      col = 'black',
      cex = text.size
    );
  }
  return()
  #legend.col(pal$cols, pal$lev); 
}

# genes can be used to restrict to a set of genes
generate.bubble.plot <- function(
  out='./',file='bubble.plot.png',genes=NULL,
  m,bin,patterns
  ){
  
  cat('\t- Mutation bubble plots\n');
  png(file=paste(out,file,sep='/'));
  info <- get.data.info(m,bin,patterns);
  
  if (!is.null(genes)){
    info <- subset(info,gene %in% genes);
  }
  op <- par(oma=c(1,1,1,1));
  create.mutation.bubble.plot(info,use.names=T,text.size=1);
  par(op);
  dev.off(); 
}


