if("argparse" %in% rownames(installed.packages()) == FALSE){
    install.package("argparse")
}
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
parser$add_argument("-m", "--maf", default = 0, help = "Remove Minor allele frequency")

args <- parser$parse_args()
file <- args$file
cv <- args$coverage
delta <- args$delta
pfile <- args$predefined
output <- args$output
maf <- args$maf
source("genocore.R")


if (grepl("csv", file)){
    tdata <- read.csv(file, header=T, check.names=F)
} else {
    tdata <- read.table(file, header = T, sep='\t', check.names = F)
}
data.set <- tdata[,-1]
rownames(data.set) <- tdata[,1]
rm(tdata)
data.set[data.set == -1] <- NA
if (pfile != "NN"){
    preset <- scan(pfile, what = "character")
} else {
    preset <- NULL
}
if (maf != 0){
    source("/data/lovemun/src_packages/GenoCore/calc_maf.R")
    cm <- calc_maf(data.set)
    rm.ix <- which(cm$maf < maf)
    write(rownames(data.set)[rm.ix], file = paste0(output, "Remove_markers.txt"))
    data.set <- data.set[-rm.ix,]
    cat("Remove ", length(rm.ix), " markers from dataset", "\n")
}

Temp_file = paste0(output, "_Temp.csv")
Cover_file = paste0(output, "_Coverage.csv")
Coreset_file = paste0(output, "_Coreset.csv")
core.set(data.set, preset = preset, coverage = cv, delta = delta, Temp_file = Temp_file, coverage_filename = Cover_file, Coreset = Coreset_file)
