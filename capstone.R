# Capstone project
library(ggplot2)
library(randomForest)
library(ggfortify)
library(cowplot)
library(matrixStats)

# Methylation
meth=readRDS("/data/compgen2016/day10_projectDay/methylation.rds")
meth$dat[1:5,1:5]
meth$cpgs[1:5,]

# Expression
exp=readRDS("/data/compgen2016/day10_projectDay/rnaseq.rds")
exp$dat[1:5,1:5]

# CNV
cna=readRDS("/data/compgen2016/day10_projectDay/cna.rds")
cna$dat[1:5,1:5]

# Annotation
pat=readRDS("/data/compgen2016/day10_projectDay/patient2subtypes.rds")
head(pat)
pat[1:5,,drop=FALSE]

# PCA
expr <- log2(exp$dat)
expr <- expr[rowSums(expr[, -1])>0, ]

meth.pca <- prcomp(t(meth$dat))
expr.pca <- prcomp(t(expr))
cnv.pca <- prcomp(t(cna$dat))
  
title <- ggdraw() + draw_label("Methylation data PCA", fontface='bold')
meth.pl <- autoplot(meth.pca, data = pat, colour="subtype")
pl1 <- plot_grid(title, meth.pl, ncol=1, rel_heights=c(0.1, 1))
title <- ggdraw() + draw_label("Expression data PCA", fontface='bold')
expr.pl <- autoplot(expr.pca, data = pat, colour="subtype")
pl2 <- plot_grid(title, expr.pl, ncol=1, rel_heights=c(0.1, 1))
title <- ggdraw() + draw_label("CNV data PCA", fontface='bold')
cnv.pl <- autoplot(cnv.pca, data = pat, colour="subtype")
pl3 <- plot_grid(title, cnv.pl, ncol=1, rel_heights=c(0.1, 1))

pl <- plot_grid(pl1, pl2, pl3, nrow=1, ncol=3, align="hv")
save_plot("plots/PCA.pdf", pl, base_width=15, nrow=1)


# Most variable data
library(matrixStats)

meth.vars <- rowVars(meth$dat)
expr.vars <- rowVars(expr)
cnv.vars <- rowVars(as.matrix(cna$dat))

meth.varmean <- mean(meth.vars)
expr.varmean <- mean(expr.vars)
cnv.varmean <- mean(cnv.vars)

cnv.N <- 50
meth.N <- 13
expr.N <- 50
meth.ndx <- order(meth.vars, decreasing = T)[1:meth.N]
expr.ndx <- order(expr.vars, decreasing = T)[1:expr.N]
cnv.ndx <- order(cnv.vars, decreasing = T)[1:cnv.N]

meth.f <- meth$dat[meth.ndx, ]
expr.f <- expr[expr.ndx, ]
cnv.f <- cna$dat[cnv.ndx, ]


# Random forest

# CNV
cnv.df = data.frame(pat, t(cnv.f))
cnv.rf <- randomForest(subtype ~ ., data=cnv.df, importance=TRUE,
                        proximity=TRUE)
print(cnv.rf)
cnv.genes <- round(importance(cnv.rf), 2)

# Methylation
meth.df = data.frame(pat, t(meth.f))
meth.rf <- randomForest(subtype ~ ., data=meth.df, importance=TRUE,
                       proximity=TRUE)
print(meth.rf)
meth.genes <- round(importance(meth.rf), 2)

# Expression
expr.df = data.frame(pat, t(expr.f))
expr.rf <- randomForest(subtype ~ ., data=expr.df, importance=TRUE,
                        proximity=TRUE)
print(expr.rf)
expr.genes <- round(importance(expr.rf), 2)


## Look at variable importance:
#round(importance(iris.rf), 2)
## Do MDS on 1 - proximity:
#iris.mds <- cmdscale(1 - iris.rf$proximity, eig=TRUE)
#op <- par(pty="s")
#pairs(cbind(iris[,1:4], iris.mds$points), cex=0.6, gap=0,
#      col=c("red", "green", "blue")[as.numeric(iris$Species)],
#      main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
#par(op)
#print(iris.mds$GOF)
