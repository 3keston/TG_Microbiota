# cor

peer.calc <- function(expr, n.fac) {
  # calculate PEER factors for m x n matrix
  # n > m
  cat("Calculating PEER factors... \n")
  require(peer)
  model = PEER()
  PEER_setPhenoMean(model,as.matrix(expr))
  dim(PEER_getPhenoMean(model))
  PEER_setNk(model,n.fac)
  PEER_getNk(model)
  PEER_update(model)

  factors = PEER_getX(model)
  rownames(factors) = rownames(expr)
  return(factors)
}

data = read.csv('../TG_noGit/Final_OTU__rarefied_forSIY.csv',
                header=T,sep=";")
rownames(data) = data$X.OTU.ID

dat2 = data[,2:dim(data)[2]]
dat2 = t(dat2)
response = dat2[,dim(dat2)[2]]
dat2 = dat2[,-c(dim(dat2)[2])]

cut.perc = 20 * 0.35

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
x.peer = peer.calc(dat.cut, 2)

out.fit = data.frame()
for (i in 1:(dim(dat.cut)[2]-1)) {
  cat(i, "\n")
  tax.test = dat.cut[,i]
  fit = coef(summary(lm(scale(as.vector(response)) ~
                        scale(as.vector(tax.test)) + x.peer)))
  name = colnames(dat.cut)[i]
  est = fit[2,1]
  sd = fit[2,2]
  pv = fit[2,4]
  temp = data.frame(id=name,
                    est=est,
                    sd=sd,
                    pv=pv)
  out.fit = rbind(out.fit, temp)

}
out.fit = out.fit[order(out.fit[,4]),]

library(qvalue)

qv = qvalue(out.fit$pv)

out.fit$qv = qv$qvalues

write.table(out.fit, file = "tax_cor_RA_rarefied_PEER2.txt", 
            quote=F, col.names=T, row.names=F)
