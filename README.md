# IlluminaHumanMethylationEPIC1b5
The manifest/anno Library for the IlluminaHumanMethylationEPIC 1.0b5 microarray with `minfi`.

## Install
```bash
unzip IlluminaHumanMethylationEPIC1b5manifest.v0.1.0.zip
unzip IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19.v0.1.0.zip
```
```R
library(minfi)
install.packages("IlluminaHumanMethylationEPIC1b5manifest",repos=NULL,type="source")
install.packages("IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19",repos=NULL,type="source")
```

## Usage example
```R
library(minfi)
library(IlluminaHumanMethylationEPIC1b5manifest)
library(IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19)

Targets <- read.metharray.sheet(base="./YOURDATADIR",pattern = "SampleSheet")
RGset <- read.metharray.exp(base="./YOURDATADIR",targets=Targets)
RGset@annotation<-c(array="IlluminaHumanMethylationEPIC1b5",annotation="ilm10b5.hg19")
```
Other usages are same as `IlluminaHumanMethylationEPICmanifest` and `IlluminaHumanMethylationEPICanno.ilm10b4.hg19` (used for EPICv1 as default)


## Downloads
The source file for installation is a little large, so it is stored in an external repository. Please download it from the link below.

### version 1.0.1 (250510)
[IlluminaHumanMethylationEPIC1b5manifest.v0.1.0.zip](https://www.dropbox.com/scl/fi/jt3w30e2c4k2p60hglsao/IlluminaHumanMethylationEPIC1b5manifest.v0.1.0.zip?rlkey=wxcl540jibp3oo0jn4y8g83r6&dl=0 "External link for IlluminaHumanMethylationEPIC1b5manifest.v0.1.0.zip")

[IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19.v0.1.0.zip](https://www.dropbox.com/scl/fi/j7gagvkkm2fgdp8gi3fuj/IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19.v0.1.0.zip?rlkey=j4gxtlbbww9pkmxgt6jh0deka&dl=0 "External link for IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19.v0.1.0.zip")

### version 0.0.1 (221215)
[IlluminaHumanMethylationEPIC1b5manifest.v0.0.1.zip](https://www.dropbox.com/scl/fi/18fss70jgmgedd7abu6a0/IlluminaHumanMethylationEPIC1b5manifest.v0.0.1.zip?rlkey=wxq0cmjmnhaibb0wydmuyuwki&dl=0 "External link for IlluminaHumanMethylationEPIC1b5manifest.v0.0.1.zip")

[IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19.v0.0.1.zip](https://www.dropbox.com/scl/fi/7ytno3inr4esamobdaw8s/IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19.v0.0.1.zip?rlkey=q25vswk42lqooeq9utl78kyiz&dl=0 "External link for IlluminaHumanMethylationEPIC1b5anno.ilm10b5.hg19.v0.0.1.zip")
