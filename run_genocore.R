#!/bin/env Rscript
## 150507
## Seongmun Jeong
## lovemun@hanyang.ac.kr

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("file", nargs=1, help="Input txt file")
parser$add_argument("-p", "--predefined", default = "NN")
parser$add_argument("-cv", "--coverage", default=99,
    help = "User defined coverage")
parser$add_argument("-d", "--delta", default=0.001,
    help = "user defined delta")
parser$add_argument("-o", "--output", default="Run",
    help = "output name")

args <- parser$parse_args()
file <- args$file
cv <- args$coverage
delta <- args$delta
pfile <- args$predefined
output <- args$output
cat("Reading input argument", "\n")
source(paste0(getwd(), "/", "genocore.R")

if (grepl("csv", file)){
    tdata <- read.csv(file, header=T)
} else {
    tdata <- read.table(file, header = T, sep='\t')
}
data.set <- tdata[,-1]
rownames(data.set) <- tdata[,1]
rm(tdata)
for (i in 1:ncol(data.set)){
    idx <- which(data.set[,i] == -1)
    data.set[idx,i] <- NA
}
cat(getwd(), "\n")
fixedfile <- paste0(output, "_fixed.csv")
write.csv(data.set, file = fixedfile)
if (pfile != "NN"){
    preset <- scan(pfile, what = "character")
    preset <- make.names(preset)
    cat(head(preset), "\n")
} else {
    preset <- NULL
}
Temp_file = paste0(output, "_Temp.csv")
Cover_file = paste0(output, "_Coverage.csv")
Coreset_file = paste0(output, "_Coreset.csv")
core.set(data.set, preset = preset, coverage = cv, delta = delta, Temp_file = Temp_file, coverage_filename = Cover_file, Coreset = Coreset_file)
