#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
data <- read.table(args[1], col.names=c("x","y"), header=F)
attach(data)
maximum1 = which.max(y)
print(x[maximum1])
