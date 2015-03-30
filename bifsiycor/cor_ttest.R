# cor

'%&%' <- function(a, b) paste(a, b, sep="")

day.vec = c(3, 7, 10, 14)

for (x in 1:length(day.vec)) {
  day = day.vec[x]

  data = read.csv('../TG_noGit/Day_'%&% day %&%'.csv',
                  header=T,sep=";")

  dat2 = data[,2:dim(data)[2]]
  dat2 = t(dat2)
  colnames(dat2) = data[,1]

  cut.perc = dim(dat2)[1] * 0.35

  id.df = data.frame()
  for (i in 1:dim(dat2)[2]) {
    taxon.temp = dat2[,i]
    tax.id = colnames(dat2)[i]
    taxon.temp.zero = taxon.temp[taxon.temp <= 0.0]
    cat(i, length(taxon.temp.zero),"\n")

    if (length(taxon.temp.zero) < cut.perc) {
      tempdf = data.frame(c1=tax.id)
      id.df = rbind(id.df, tempdf)
    } 
  }

  rownames(id.df) = id.df[,1]

  dat.cut = dat2[,intersect(rownames(id.df), colnames(dat2))]

  out.fit = data.frame()
  for (i in 1:(dim(dat.cut)[2]-1)) {
    cat(i, "\n")
    tax.test = dat.cut[,i]
    g1 = tax.test[1:6]
    g2 = tax.test[7:13]
    fit.t = t.test(g1, g2)
    fit.w = wilcox.test(g1,g2, exact=F)
    name = colnames(dat.cut)[i]
    t.stat = fit.t$statistic
    t.pv = fit.t$p.value
    w.stat = fit.w$statistic
    w.pv = fit.w$p.value
    temp = data.frame(id=name,
                      wstat=w.stat,
                      wpv=w.pv,
                      tstat=t.stat,
                      tpv=t.pv)
    out.fit = rbind(out.fit, temp)

  }

  library(qvalue)

  qv.w = qvalue(out.fit$wpv)
  qv.t = qvalue(out.fit$tpv)

  out.fit$wqv = qv.w$qvalues
  out.fit$tqv = qv.t$qvalues

  out.fit = out.fit[order(out.fit$wqv),]

  write.table(out.fit, file = "../TG_noGit/L6_D"%&% day %&%"_ttest.txt", 
              quote=F, col.names=T, row.names=F)
}
