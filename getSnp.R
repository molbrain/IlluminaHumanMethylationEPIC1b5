processUCSCsnp <- function(snpfile) {
    require(GenomicRanges)
    cat("Reading file\n")
    df <- read.delim(gzfile(snpfile), header = FALSE,
                     stringsAsFactors = FALSE)
    names(df) <- c("chr", "start", "end", "name", "strand",
                   "refNCBI", "class", "alleleFreqs")
    print(table(df$chr))
    cat("Only keeping chrs 1-22, X, Y\n")
    df <- df[df$chr %in% paste0("chr", c(1:22, "X", "Y")),]
    print(table(df$class))
    cat("Only keeping class 'single'\n")
    df <- df[df$class == "single",]
    cat("Computing MAF\n")
    df$alleleFreqs <- sub(",$", "", df$alleleFreqs)
    sp <- strsplit(df$alleleFreqs, ",")
    minFreq <- sapply(sp, function(xx) min(as.numeric(xx)))
    cat("Instantiating object\n")
    grSNP <- GRanges(seqnames = df$chr, strand = df$strand,
                     ranges = IRanges(start = df$start + 1, end = df$end),
                     MAF = minFreq, ref = df$refNCBI)
    names(grSNP) <- df$name
    grSNP
}


grSnp138CommonSingle <- processUCSCsnp("./snp138Common_small.txt.gz")
save(grSnp138CommonSingle, file = "./grSnp138CommonSingle.rda")

grSnp142CommonSingle <- processUCSCsnp("./snp142Common_small.txt.gz")
save(grSnp142CommonSingle, file = "./grSnp142CommonSingle.rda")

grSnp144CommonSingle <- processUCSCsnp("./snp144Common_small.txt.gz")
save(grSnp144CommonSingle, file = "./grSnp144CommonSingle.rda")

grSnp146CommonSingle <- processUCSCsnp("./snp146Common_small.txt.gz")
save(grSnp146CommonSingle, file = "./grSnp146CommonSingle.rda")

grSnp147CommonSingle <- processUCSCsnp("./snp147Common_small.txt.gz")
save(grSnp147CommonSingle, file = "./grSnp147CommonSingle.rda")

#grSnp150CommonSingle <- processUCSCsnp("./snp150Common_small.txt.gz")
#save(grSnp150CommonSingle, file = "./grSnp150CommonSingle.rda")

#grSnp151CommonSingle <- processUCSCsnp("./snp151Common_small.txt.gz")
#save(grSnp151CommonSingle, file = "./grSnp151CommonSingle.rda")
