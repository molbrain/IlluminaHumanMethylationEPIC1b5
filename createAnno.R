library(minfi)
manifestFile <- "./infinium-methylationepic-v-1-0-b5-manifest-file.csv"
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

Locations <- anno[, c("CHR", "MAPINFO")]
names(Locations) <- c("chr", "pos")
Locations$pos <- as.integer(Locations$pos)
Locations$chr <- paste("chr", Locations$chr, sep = "")
Locations$strand <- ifelse(anno$Strand == "F", "+", "-")
table(Locations$chr, exclude = NULL)
rownames(Locations) <- anno$Name
Locations <- as(Locations, "DataFrame")

Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")

Islands.UCSC <- anno[, c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")]
names(Islands.UCSC) <- c("Islands_Name", "Relation_to_Island")
Islands.UCSC <- as(Islands.UCSC, "DataFrame")
Islands.UCSC$Relation_to_Island[Islands.UCSC$Relation_to_Island == ""] <- "OpenSea"
table(Islands.UCSC$Relation_to_Island, exclude = NULL)

SNPs.Illumina <- anno[, c("SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]
SNPs.Illumina <- as(SNPs.Illumina, "DataFrame")

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

## We now use an exisitng grSnp object containing a GRanges of relevant SNPs.
## This is created in a separate script

##
## SNP overlap
##

map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos, width = 1),
               Strand = map$strand, Type = map$Type)
map <- minfi:::.getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

## dbSNP
load("./grSnp151CommonSingle.rda")
SNPs.151CommonSingle <- minfi:::.doSnpOverlap(map, grSnp151CommonSingle)
load("./grSnp150CommonSingle.rda")
SNPs.150CommonSingle <- minfi:::.doSnpOverlap(map, grSnp150CommonSingle)
load("./grSnp147CommonSingle.rda")
SNPs.147CommonSingle <- minfi:::.doSnpOverlap(map, grSnp147CommonSingle)
load("./grSnp146CommonSingle.rda")
SNPs.146CommonSingle <- minfi:::.doSnpOverlap(map, grSnp146CommonSingle)
load("./grSnp144CommonSingle.rda")
SNPs.144CommonSingle <- minfi:::.doSnpOverlap(map, grSnp144CommonSingle)
load("./grSnp142CommonSingle.rda")
SNPs.142CommonSingle <- minfi:::.doSnpOverlap(map, grSnp142CommonSingle)
load("./grSnp141CommonSingle.rda")
SNPs.141CommonSingle <- minfi:::.doSnpOverlap(map, grSnp141CommonSingle)
load("./grSnp138CommonSingle.rda")
SNPs.138CommonSingle <- minfi:::.doSnpOverlap(map, grSnp138CommonSingle)

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

assign(pkgName, annoObj)
save(list = pkgName,
     file = file.path("./data", paste(pkgName, "rda", sep = ".")), compress = "xz")
