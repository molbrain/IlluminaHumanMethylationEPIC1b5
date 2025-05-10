#!/bin/zsh

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp150Common.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp147Common.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp146Common.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp144Common.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp142Common.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp138Common.txt.gz

gunzip -c snp151Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp151Common_small.txt.gz
gunzip -c snp150Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp150Common_small.txt.gz
gunzip -c snp147Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp147Common_small.txt.gz
gunzip -c snp146Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp146Common_small.txt.gz
gunzip -c snp144Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp144Common_small.txt.gz
gunzip -c snp142Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp142Common_small.txt.gz
gunzip -c snp138Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp138Common_small.txt.gz

