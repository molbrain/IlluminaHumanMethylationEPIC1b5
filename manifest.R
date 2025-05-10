library(minfi)

manifestFile <- "./EPIC/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
maniTmp <- minfi:::read.manifest.EPIC(manifestFile)
maniList <- maniTmp$manifestList

IlluminaHumanMethylationEPIC1b5manifest <- IlluminaMethylationManifest(
    TypeI = maniList$TypeI,
    TypeII = maniList$TypeII,
    TypeControl = maniList$TypeControl,
    TypeSnpI = maniList$TypeSnpI,
    TypeSnpII = maniList$TypeSnpII,
    annotation = "IlluminaHumanMethylationEPIC1b5"
)

save(IlluminaHumanMethylationEPIC1b5manifest,
     compress = "xz",
     file = "./IlluminaHumanMethylationEPIC1b5manifest.rda"
)
