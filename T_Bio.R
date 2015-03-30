# file analyzes TG microbiota data

source('T_Bio_Func.R')
library(gplots)
library(RColorBrewer)

dat1 = read.table('../TG_noGit/quan_norm_bgsub.txt', 
                  sep='\t', 
                  header=T, 
                  as.is=T)
rownames(dat1) = dat1[,2] # probeid

# making data frames for the groups
df.g1.tac = data.frame(dat1$TG1.AVG_Signal, dat1$TG7.AVG_Signal)
df.g1.jax = data.frame(dat1$TG3.AVG_Signal, dat1$TG9.AVG_Signal) 
df.g1.bif = data.frame(dat1$TG5.AVG_Signal, dat1$TG11.AVG_Signal) 

df.g2.tac = data.frame(dat1$TG2.AVG_Signal, dat1$TG8.AVG_Signal)
df.g2.jax = data.frame(dat1$TG4.AVG_Signal, dat1$TG10.AVG_Signal)
df.g2.bif = data.frame(dat1$TG6.AVG_Signal, dat1$TG12.AVG_Signal)

df.g3.tac = data.frame(dat1$TG14.AVG_Signal)
df.g3.jax = data.frame(dat1$TG13.AVG_Signal)
df.g3.bif = data.frame(dat1$TG15a.AVG_Signal) 

# cleaning tac to measure outliers
rn = rownames(dat1)
df.g1.tac.cl = clean.func(df.g1.tac, rn)
df.g2.tac.cl = clean.func(df.g2.tac, rn)
df.g3.tac.cl = clean.func(df.g3.tac, rn)

vcut=10
df.g1.jax.cl = clean.func(df.g1.jax, rn, vcut=vcut)
df.g2.jax.cl = clean.func(df.g2.jax, rn, vcut=vcut)
df.g3.jax.cl = clean.func(df.g3.jax, rn, vcut=vcut)

df.g1.bif.cl = clean.func(df.g1.bif, rn, vcut=vcut)
df.g2.bif.cl = clean.func(df.g2.bif, rn, vcut=vcut)
df.g3.bif.cl = clean.func(df.g3.bif, rn, vcut=vcut)

# calculating fold change 04/07/15
cutfc = NULL
# g1
fc.g1.exp1.jax.tac = calc.fold.change(df.g1.jax.cl[,1], 
                                      df.g1.tac.cl[,1],
                                      rownames(dat1),
                                      cutfc=cutfc)

fc.g1.exp2.jax.tac = calc.fold.change(df.g1.jax.cl[,2],
                                      df.g1.tac.cl[,2],
                                      rownames(dat1),
                                      cutfc=cutfc)

fc.g1.exp1.bif.tac = calc.fold.change(df.g1.bif.cl[,1],
                                      df.g1.tac.cl[,1],
                                      rownames(dat1),
                                      cutfc=cutfc)

fc.g1.exp2.bif.tac = calc.fold.change(df.g1.bif.cl[,2],
                                      df.g1.tac.cl[,2],
                                      rownames(dat1),
                                      cutfc=cutfc)

# g2 
fc.g2.exp1.jax.tac = calc.fold.change(df.g2.jax.cl[,1], 
                                      df.g2.tac.cl[,1],
                                      rownames(dat1),
                                      cutfc=cutfc)


fc.g2.exp1.bif.tac = calc.fold.change(df.g2.bif.cl[,1],
                                      df.g2.tac.cl[,1],
                                      rownames(dat1),
                                      cutfc=cutfc)

# g3
fc.g3.exp1.jax.tac = calc.fold.change(df.g3.jax.cl[,1],
                                      df.g3.tac.cl[,1],
                                      rownames(dat1),
                                      cutfc=cutfc)

fc.g3.exp1.bif.tac = calc.fold.change(df.g3.bif.cl[,1],
                                      df.g3.tac.cl[,1],
                                      rownames(dat1),
                                      cutfc=cutfc)

# g1
results.fc.g1 = data.frame(probeid=rownames(dat1),
                           gene=dat1$TargetID,
                           jax_tac1=fc.g1.exp1.jax.tac,
                           jax_tac2=fc.g1.exp2.jax.tac,
                           bif_tac1=fc.g1.exp1.bif.tac,
                           bif_tac2=fc.g1.exp2.bif.tac)
mean.fc.g1.1 = rowMeans(results.fc.g1[,3:4], na.rm=T)
mean.fc.g1.2 = rowMeans(results.fc.g1[,5:6], na.rm=T)
results.fc.g1$four_mean_1 = mean.fc.g1.1
results.fc.g1$four_mean_2 = mean.fc.g1.2

results.fc.g1.cut.50 = results.fc.g1[which(results.fc.g1[,7] >= 1.5 &
                                           results.fc.g1[,8] >= 1.5),]
results.fc.g1.cut.25 = results.fc.g1[which(results.fc.g1[,7] >= 1.25 &
                                           results.fc.g1[,8] >= 1.25),]

results.fc.g1.cut.jax = results.fc.g1[which(results.fc.g1[,7] >= 1.50),]
results.fc.g1.cut.bif = results.fc.g1[which(results.fc.g1[,8] >= 1.50),]

write.table(results.fc.g1.cut.50, file = "results.fc.g1.cut.50.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g1.cut.25, file = "results.fc.g1.cut.25.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g1.cut.jax, file = "results.fc.g1.cut.jax.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g1.cut.bif, file = "results.fc.g1.cut.bif.txt", quote=F,
            col.names=T, row.names=F)

# g2
results.fc.g2 = data.frame(probeid=rownames(dat1),
                           gene=dat1$TargetID,
                           jax_tac1=fc.g2.exp1.jax.tac,
                           jax_tac2=fc.g2.exp2.jax.tac,
                           bif_tac1=fc.g2.exp1.bif.tac,
                           bif_tac2=fc.g2.exp2.bif.tac)
mean.fc.g2.1 = rowMeans(results.fc.g2[,3:4], na.rm=T)
mean.fc.g2.2 = rowMeans(results.fc.g2[,5:6], na.rm=T)
results.fc.g2$four_mean_1 = mean.fc.g2.1
results.fc.g2$four_mean_2 = mean.fc.g2.2

results.fc.g2.cut.50 = results.fc.g2[which(results.fc.g2[,7] >= 1.5 &
                                           results.fc.g2[,8] >= 1.5),]
results.fc.g2.cut.25 = results.fc.g2[which(results.fc.g2[,7] >= 1.25 &
                                           results.fc.g2[,8] >= 1.25),]

results.fc.g2.cut.jax = results.fc.g2[which(results.fc.g2[,7] >= 1.50),]
results.fc.g2.cut.bif = results.fc.g2[which(results.fc.g2[,8] >= 1.50),]

write.table(results.fc.g2.cut.50, file = "results.fc.g2.cut.50.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g2.cut.25, file = "results.fc.g2.cut.25.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g2.cut.jax, file = "results.fc.g2.cut.jax.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g2.cut.bif, file = "results.fc.g2.cut.bif.txt", quote=F,
            col.names=T, row.names=F)

# g3
results.fc.g3 = data.frame(probeid=rownames(dat1),
                           gene=dat1$TargetID,
                           jax_tac1=fc.g3.exp1.jax.tac,
                           bif_tac1=fc.g3.exp1.bif.tac)

results.fc.g3.cut.50 = results.fc.g3[which(results.fc.g3[,3] >= 1.5 &
                                           results.fc.g3[,4] >= 1.5),]
results.fc.g3.cut.25 = results.fc.g3[which(results.fc.g3[,3] >= 1.25 &
                                           results.fc.g3[,4] >= 1.25),]

results.fc.g3.cut.jax = results.fc.g3[which(results.fc.g3[,3] >= 1.50),]
results.fc.g3.cut.bif = results.fc.g3[which(results.fc.g3[,4] >= 1.50),]

write.table(results.fc.g3.cut.50, file = "results.fc.g3.cut.50.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g3.cut.25, file = "results.fc.g3.cut.25.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g3.cut.jax, file = "results.fc.g3.cut.jax.txt", quote=F,
            col.names=T, row.names=F)
write.table(results.fc.g3.cut.bif, file = "results.fc.g3.cut.bif.txt", quote=F,
            col.names=T, row.names=F)

# read choice genes
gene.list = read.table('../TG_noGit/targets3.txt')
rownames(gene.list) = gene.list[,1]

dat1.list = dat1[intersect(rownames(gene.list), rownames(dat1)),]
gnames = dat1.list$TargetID

# tac tac jax jax bif bif for g1
heat.df = data.frame(df.g1.tac.cl, df.g1.jax.cl, df.g1.bif.cl)
heat.dm = data.matrix(heat.df)
heat.dm = heat.dm[intersect(rownames(gene.list), rownames(heat.dm)),]
rownames(heat.dm) = gnames
colnames(heat.dm) = c("Tac_1", "Tac_2", "Jax_1", "Jax_2", "Bif_1", "Bif_2")
heat.dm = log2(heat.dm+1)
heat.dm = t(scale(t(heat.dm), center=TRUE, scale=TRUE))

pdf('~/TG_noGit/heatrg_t1.pdf')
heatmap.2(as.matrix(heat.dm),
          lhei=c(0.2,1.0),
          lwid=c(0.2,0.1),
          Rowv=NULL,
          Colv=NULL,
          dendrogram="none", 
          notecol="black",
          col=bluered(256),
          scale="none",
          key=TRUE,
          keysize=1.5,
          density.info="none",
          trace="none",
          sepwidth=c(0.05, 0.05),
          sepcolor='black',
          rowsep=1:nrow(heat.dm),
          colsep=1:ncol(heat.dm),
          cexRow=0.7,
          cexCol=1.2
          )
dev.off()

sig.taxa =
read.table('../TG_noGit/Final_table_for_heatmap_miseq3_072915.csv',
           header=T,sep=";")
rownames(sig.taxa) = sig.taxa[,1]
sig.taxa = sig.taxa[,2:dim(sig.taxa)[2]]
sig.taxa = data.matrix(sig.taxa)
sig.taxa = log2(sig.taxa+1)
sig.taxa = t(scale(t(sig.taxa), center=TRUE, scale=TRUE))

pdf('../TG_noGit/SigTaxa_bluered.pdf')
heatmap.2(as.matrix(sig.taxa),
          lhei=c(0.2,0.5),
          lwid=c(0.01,0.1),
          margin=c(5, 10),
          Rowv=NULL,
          Colv=NULL,
          dendrogram="none",
          notecol="black",
          #col=rev(brewer.pal(8, "RdBu")),
          col=bluered(256),
          scale="none",
          key=FALSE,
          keysize=0.001,
          density.info="none",
          trace="none",
          sepwidth=c(0.005, 0.005),
          sepcolor='black',
          rowsep=1:nrow(sig.taxa),
          colsep=1:ncol(sig.taxa),
          cexRow=0.7,
          cexCol=1.2
          )
dev.off()


sig.taxa =
read.table('../TG_noGit/Final_table_for_heatmap_miseq2_FDR.csv',
           header=T,sep=",")
rownames(sig.taxa) = sig.taxa[,1]
sig.taxa = sig.taxa[,2:dim(sig.taxa)[2]]
sig.taxa = data.matrix(sig.taxa)
sig.taxa = log2(sig.taxa+1)
sig.taxa = t(scale(t(sig.taxa), center=TRUE, scale=TRUE))

pdf('../TG_noGit/miSeq2_bluered_FDR.pdf')
heatmap.2(as.matrix(sig.taxa),
          lhei=c(0.2,0.1),
          lwid=c(0.1,0.1),
          margin=c(5, 10),
          Rowv=NULL,
          Colv=NULL,
          dendrogram="none",
          notecol="black",
          #col=rev(brewer.pal(8, "RdBu")),
          col=bluered(256),
          scale="none",
          key=FALSE,
          keysize=0.001,
          density.info="none",
          trace="none",
          sepwidth=c(0.005, 0.005),
          sepcolor='black',
          rowsep=1:nrow(sig.taxa),
          colsep=1:ncol(sig.taxa),
          cexRow=0.7,
          cexCol=1.2
          )
dev.off()


#
