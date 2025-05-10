print(c("START",date()))

library(minfi)

manifestFile <- "./sources/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
maniTmp <- minfi:::read.manifest.EPIC(manifestFile)
anno <- maniTmp$manifest
manifestList <- maniTmp$manifestList

IlluminaHumanMethylationEPIC1b5manifest <- do.call(
    IlluminaMethylationManifest,
    list(TypeI = manifestList$TypeI,
    TypeII = manifestList$TypeII,
    TypeControl = manifestList$TypeControl,
    TypeSnpI = manifestList$TypeSnpI,
    TypeSnpII = manifestList$TypeSnpII,
    annotation = "IlluminaHumanMethylationEPIC1b5"))

anno$IlmnID <- NULL
nam <- names(anno)
names(nam) <- nam
nam[c("AddressA_ID", "AddressB_ID",
      "AlleleA_ProbeSeq", "AlleleB_ProbeSeq",
      "Infinium_Design_Type",
      "Next_Base", "Color_Channel")] <-  c("AddressA", "AddressB",
                                           "ProbeSeqA", "ProbeSeqB",
                                           "Type", "NextBase", "Color")

names(nam) <- NULL
names(anno) <- nam
rownames(anno) <- anno$Name
anno <- anno[getManifestInfo(IlluminaHumanMethylationEPIC1b5manifest, type = "locusNames"),]


## Locations
print("## Locations")

Locations <- anno[, c("CHR", "MAPINFO")]
names(Locations) <- c("chr", "pos")
Locations$pos <- as.integer(Locations$pos)
Locations$chr <- paste("chr", Locations$chr, sep = "")
Locations$strand <- ifelse(anno$Strand == "F", "+", "-")
table(Locations$chr, exclude = NULL)
rownames(Locations) <- anno$Name
Locations <- as(Locations, "DataFrame")


## Manifest
print("## Manifest")

Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")


## Islands.UCSC
print("## Islands.UCSC")

Islands.UCSC <- anno[, c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")]
names(Islands.UCSC) <- c("Islands_Name", "Relation_to_Island")
Islands.UCSC <- as(Islands.UCSC, "DataFrame")
Islands.UCSC$Relation_to_Island[Islands.UCSC$Relation_to_Island == ""] <- "OpenSea"
table(Islands.UCSC$Relation_to_Island, exclude = NULL)


## SNPs.Illumina
print("## SNPs.Illumina")

SNPs.Illumina <- anno[, c("SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]
SNPs.Illumina <- as(SNPs.Illumina, "DataFrame")


## Other
print("## Other")

usedColumns <- c(names(Manifest), names(SNPs.Illumina), 
                 c("CHR", "MAPINFO", "Strand",
                   "Chromosome_36", "Coordinate_36", "Genome_Build"),
                 c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"))
Other <- anno[, setdiff(names(anno), usedColumns)]
nam <- names(Other)
nam <- sub("_NAME", "_Name", nam)
nam[nam == "X450k_Enhancer"] <- "Methyl450_Enhancer"
nam
Other <- as(Other, "DataFrame")


## grSNP
print("## grSNP")
## We now use an exisitng grSnp object containing a GRanges of relevant SNPs.
processUCSCsnp <- function(snpfile) {
    require(GenomicRanges)
    cat("Reading file\n")
    df <- read.delim(gzfile(snpfile), header = FALSE,
                     stringsAsFactors = FALSE)
    df <- df[,c(2,3,4,5,7,8,12,25)]
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

print("## grSnp138CommonSingle")
grSnp138CommonSingle <- processUCSCsnp("./sources/snp138Common.txt.gz")
print("## grSnp141CommonSingle")
grSnp141CommonSingle <- processUCSCsnp("./sources/snp141Common.txt.gz")
print("## grSnp142CommonSingle")
grSnp142CommonSingle <- processUCSCsnp("./sources/snp142Common.txt.gz")
print("## grSnp144CommonSingle")
grSnp144CommonSingle <- processUCSCsnp("./sources/snp144Common.txt.gz")
print("## grSnp146CommonSingle")
grSnp146CommonSingle <- processUCSCsnp("./sources/snp146Common.txt.gz")
print("## grSnp147CommonSingle")
grSnp147CommonSingle <- processUCSCsnp("./sources/snp147Common.txt.gz")
print("## grSnp150CommonSingle")
grSnp150CommonSingle <- processUCSCsnp("./sources/snp150Common.txt.gz")
print("## grSnp151CommonSingle")
grSnp151CommonSingle <- processUCSCsnp("./sources/snp151Common.txt.gz")


## SNP overlap
print("## SNP overlap")
map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos, width = 1),
               Strand = map$strand, Type = map$Type)
map <- minfi:::.getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

## dbSNP
print("## SNPs.151")
SNPs.151CommonSingle <- minfi:::.doSnpOverlap(map, grSnp151CommonSingle)
print("## SNPs.150")
SNPs.150CommonSingle <- minfi:::.doSnpOverlap(map, grSnp150CommonSingle)
print("## SNPs.147")
SNPs.147CommonSingle <- minfi:::.doSnpOverlap(map, grSnp147CommonSingle)
print("## SNPs.146")
SNPs.146CommonSingle <- minfi:::.doSnpOverlap(map, grSnp146CommonSingle)
print("## SNPs.144")
SNPs.144CommonSingle <- minfi:::.doSnpOverlap(map, grSnp144CommonSingle)
print("## SNPs.142")
SNPs.142CommonSingle <- minfi:::.doSnpOverlap(map, grSnp142CommonSingle)
print("## SNPs.141")
SNPs.141CommonSingle <- minfi:::.doSnpOverlap(map, grSnp141CommonSingle)
print("## SNPs.138")
SNPs.138CommonSingle <- minfi:::.doSnpOverlap(map, grSnp138CommonSingle)


## create annoObj
print("## create annoObj")

annoNames <- c(
    "Locations", "Manifest", "SNPs.Illumina",
    "SNPs.151CommonSingle", "SNPs.150CommonSingle", "SNPs.147CommonSingle", "SNPs.146CommonSingle",
    "SNPs.144CommonSingle", "SNPs.142CommonSingle", "SNPs.141CommonSingle",
    "SNPs.138CommonSingle",
    "Islands.UCSC", "Other")
for(nam in annoNames) {
    cat(nam, "\n")
    save(list = nam, file = file.path("./data", paste(nam, "rda", sep = ".")), compress = "xz")
}
annoStr <- c(array = "IlluminaHumanMethylationEPIC1b5",
             annotation = "ilm10b5",
             genomeBuild = "hg19")
defaults <- c("Locations", "Manifest", "SNPs.138CommonSingle", "Islands.UCSC", "Other")
pkgName <- sprintf("%sanno.%s.%s", annoStr["array"], annoStr["annotation"],
                    annoStr["genomeBuild"])

annoObj <- IlluminaMethylationAnnotation(objectNames = annoNames, annotation = annoStr,
                              defaults = defaults, packageName = pkgName)


## create package
print("## create package")

assign(pkgName, annoObj)
save(list = pkgName,
     file = file.path("./data", paste(pkgName, "rda", sep = ".")), compress = "xz")

print(c("DONE",date()))

