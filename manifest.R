print(c("START",date()))

library(minfi)

manifestFile <- "./sources/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
maniTmp <- minfi:::read.manifest.EPIC(manifestFile)
manifestList <- maniTmp$manifestList

IlluminaHumanMethylationEPIC1b5manifest <- IlluminaMethylationManifest(
    TypeI = manifestList$TypeI,
    TypeII = manifestList$TypeII,
    TypeControl = manifestList$TypeControl,
    TypeSnpI = manifestList$TypeSnpI,
    TypeSnpII = manifestList$TypeSnpII,
    annotation = "IlluminaHumanMethylationEPIC1b5"
)

save(IlluminaHumanMethylationEPIC1b5manifest,
     compress = "xz",
     file = "./data/IlluminaHumanMethylationEPIC1b5manifest.rda"
)


print(c("DONE",date()))
