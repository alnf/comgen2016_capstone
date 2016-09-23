# Capstone project

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


