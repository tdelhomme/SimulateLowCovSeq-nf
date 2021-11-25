#!/usr/bin/env Rscript

meancov=as.numeric(commandArgs(TRUE)[1])
targetcov=as.numeric(commandArgs(TRUE)[2])

res = round(100 / (meancov / targetcov), 1)

if(res * 10 < 100){ # this is numbers lower than 10 percent
  res = gsub("\\.", "", paste(0, res, sep=""))
}

cat(res)
