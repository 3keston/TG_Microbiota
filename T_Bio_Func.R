# functions

fix.ttest.const <- function(x, y) {
  obj = try(t.test(x, y), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj)      
}

fit.g <- function(mat, gene.idvec, probe.id, test.name) {
  # fits t-test for groups assuming groups are equal
  # TODO: 

  # assign rownames to mat
  rownames(mat) = probe.id

  # ranges
  range.split.lower = dim(mat)[2] / 2
  range.split.upper = range.split.lower + 1
  range.lower = 1
  range.upper = dim(mat)[2]

  output.df = data.frame()
  loop.range = 1:dim(mat)[1]
  #loop.range = 10300:10400 # debug
  for (i in loop.range) {
    cat("Running >", test.name, "<>", i, "\n")
    # gene
    slice = log2(as.numeric(mat[i, ])+1)

    # fold change for upper / lower
    fc.upper = slice[range.upper:range.split.upper]
    fc.lower = slice[range.lower:range.split.lower]
    
    up.mean = mean(fc.upper)
    lo.mean = mean(fc.lower)

    fold.change = up.mean / lo.mean

    fit = fix.ttest.const(scale(fc.lower), scale(fc.upper))
    
    # addressing the issue of constant vectors
    if (is.na(fit)[1]) {
      temp.out.df = data.frame(gene=gene.idvec[i],
                               tac=lo.mean,
                               com=up.mean,
                               fold=fold.change,
                               tttstat=NA,
                               ttpval=NA,
                               lmbeta=NA,
                               lmse=NA,
                               lmpval=NA)
      output.df = rbind(output.df, temp.out.df)
    } else {
      fit.tstat = fit$statistic
      fit.pvalu = fit$p.value
      
      sm.lm =
      summary(lm(c(rep(0,range.split.lower),rep(1,range.split.lower)) 
                 ~ scale(as.numeric(slice))))
      lm.beta = sm.lm$coefficients[2,1]
      lm.stdr = sm.lm$coefficients[2,2]
      lm.pval = sm.lm$coefficients[2,4]

      temp.out.df = data.frame(gene=gene.idvec[i],
                               tac=lo.mean,
                               com=up.mean,
                               fold=fold.change,
                               tttstat=fit.tstat,
                               ttpval=fit.pvalu,
                               lmbeta=lm.beta,
                               lmse=lm.stdr,
                               lmpval=lm.pval)
      output.df = rbind(output.df, temp.out.df) 
    }
  }
  # calculate qvalues
  rownames(output.df) = rownames(mat)[loop.range]
    
  return(output.df)
}

anno.qv <- function(output.df) {
  # annotates fit.g with qvalues
  require(qvalue)
  if (length(is.na(fitted$tttstat)) == length(fitted$fold)) {
    return(fitted)
  } else {
    output.qv = output.df[!is.na(output.df$ttpval),] # remove NA for qval
    qv.tt = qvalue(output.qv$ttpval)
    qv.lm = qvalue(output.qv$lmpval)

    qv.temp = data.frame(qv.lm$qvalues, qv.tt$qvalues)
    rownames(qv.temp) = rownames(output.qv)
    output.df = merge(output.df, qv.temp, by = "row.names", all.x = TRUE)
    colnames(output.df)[1] = "probe"
    colnames(output.df)[c(length(colnames(output.df))  )] = "ttqv"
    colnames(output.df)[c(length(colnames(output.df))-1)] = "lmqv"
    
    return(output.df)
 }
}

outliers.func <- function(x, y, nclust, spec.x, spec.y) {
  # this was taking x and y vectors and clustering. might not be the best
  # choice and will likely be rewritten
  # TODO rewrite or trash
  temp.resid = lm(y ~ x)$residuals
  temp.km = kmeans(temp.resid, nclust)$cluster
  temp.res = data.frame(c1=x , c2=y, c3=temp.km)
  p = ggplot(temp.res, aes(c1, c2))
  p = p + geom_point(aes(colour = factor(c3)))
  p = p + xlab(spec.x) + ylab(spec.y) 
  p = p + theme_bw() + coord_fixed()
  return(p)
}

sd.func <- function(x1, x2, x, y) {
  # take replicates x and y, regress y on x and calculate the variance which we
  # will use as a cutoff for the next function
  # TODO this function is so massively wrong
  main.resid = lm(scale(x1) ~ scale(x2))$residuals
  main.cut = sd(main.resid)
  xy.resid = lm(scale(y) ~ scale(x))$residuals
  xy.resid[xy.resid <= main.cut] = NA
  temp.df = data.frame(x, y, xy.resid)
  return(temp.df)
}

calc.fold.change <- function(over, undr, names.vec, cutfc=NULL) {
  # calculates fold change
  # TODO
  names(over) = names.vec
  names(undr) = names.vec
  fc = over / undr
  names(fc) = names.vec

  if(!is.null(cutfc)) {
    fc[fc < cutfc] = NA
  }

  return(fc)
}

clean.func <- function(df.clean, rownames.df, vfloor=TRUE, vcut=NULL, 
                       rsum=FALSE, vscale=FALSE, logt=FALSE) {
  # this was taking a data frame of observations, averaging of columns, then
  # flooring the to 0, then log2 transforming when necessary. 
  # TODO
  rownames(df.clean) = rownames.df

  if (rsum == TRUE) {
    df.clean = scale(df.clean, center=TRUE, scale=TRUE)
    df.clean = rowSums(df.clean)
  }
  
  if (!is.null(vcut)) {
    df.clean[df.clean < vcut] = NA
  }
 
  if (vfloor == TRUE) {
    df.clean[df.clean < 0] = 0
  }

  if (logt == TRUE) {
    df.clean = log2(df.clean + 1)
  }

  if (vscale == TRUE){
    df.clean = scale(df.clean, scale=TRUE, center=TRUE)
  }
  
  return(df.clean)
}

#
