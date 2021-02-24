# GEO platform information is downloaded from https://www.ncbi.nlm.nih.gov/geo/browse/?view=platforms
# in five `.csv` files and combined/saved as `platforms.csv`. As of Feb 24, 2021,
# there are 21,888 platforms exist in GEO and 5,272 of them are Homo sapiens.

download.file(url = "https://www.dropbox.com/s/o3q0h05jeh2fiyw/GEOplatforms_2021-02-23.zip?dl=0",
              method = "wget",
              destfile = "gpl.zip")
unzip("gpl.zip")
platformlist <-
  lapply(1:5, function(n)
    read.csv(file = paste0("platform_", n, ".csv"), as.is = TRUE))
all.platforms <- Reduce(rbind, platformlist)
dim(all.platforms) #All platforms

rm(platformlist)

hs_ind <- which(all.platforms$Taxonomy == "Homo sapiens")
all.platforms <- all.platforms[hs_ind, ]
dim(all.platforms) #Homo sapiens platforms only

# human Entrez Gene identifier symbols
library(org.Hs.eg.db)
hgnc.reference <- as.character(org.Hs.egSYMBOL)
names(hgnc.reference) <- NULL

library(HGNChelper)

raw.platform.dir <- "platforms"
dir.create(raw.platform.dir)

hgnc.vec.dir <- "hgnc.vecs"
dir.create(hgnc.vec.dir)


##### Create `gpls_already_tested.csv` file ####################################
if (file.exists("gpls_already_tested.csv")) {
  gpls.already.tested <- read.csv("gpls_already_tested.csv", header=TRUE)
} else {
  write.table(t(c("platform", "colname", "frac.hgnc", "nrow",
                  "valid.frac", "valid.after.hgnchelper.frac",
                  "distribution", "submission_date")),
              file="gpls_already_tested.csv",
              row.names=FALSE, col.names=FALSE, sep=",")
}

currentHumanMap = getCurrentHumanMap()

library(GEOquery)
for (i in seq(nrow(all.platforms))){
  gpl <- all.platforms[i, 1]
  hgnc.vec.file <- file.path(hgnc.vec.dir, paste0(gpl, "_hgnc.vec.RData"))

  if(exists("gpls.already.tested") && gpl %in% gpls.already.tested$platform) next
  gpldat <- try(getGEO(gpl, destdir=raw.platform.dir))
  if(class(gpldat) == "try-error") next
  gpltable <- try(Table(gpldat))
  if(class(gpltable) == "try-error") next
  hgnc.frac <- apply(gpltable, 2, function(x) sum(unique(x) %in% hgnc.reference) / length(unique(x)))

  # any column containing gene symbol will have >0 value
  if(any(hgnc.frac > 0.5)) {   # updated from `hgnc.frac > 0`

    # assumed that there is only one 'symbol' column
    hgnc.vec <- unique(as.character(gpltable[, which.max(hgnc.frac)]))
    # get rid of anything after a space
    hgnc.vec <- gsub("[ ].+", "", hgnc.vec)
    # convert to ascii
    HGNChelper.output <- checkGeneSymbols(iconv(hgnc.vec, "latin1", "ASCII", ""),
                                          map = currentHumanMap)

    # fraction of valid gene symbols BEFORE HGNChelper
    valid.frac <- sum(HGNChelper.output$Approved) / length(hgnc.vec)
    # fraction of valid gene symbols AFTER HGNChelper
    after.HGNChelper.valid.frac <- sum(!is.na(HGNChelper.output$Suggested.Symbol)) / length(hgnc.vec)

    hgnc.vec <- c(gpldat@header$distribution, gpldat@header$submission_date, hgnc.vec)
    save(hgnc.vec, file=hgnc.vec.file, compress="bzip2")

    info.for.file <- t(c(gpl, colnames(gpltable)[which.max(hgnc.frac)],
                         max(hgnc.frac), nrow(gpltable),
                         valid.frac, after.HGNChelper.valid.frac,
                         gpldat@header$distribution,
                         gpldat@header$submission_date))

    print(paste(i, gpl, length(hgnc.vec), "HGNC symbols found and saved."))

  } else {
    info.for.file <- t(c(gpl, colnames(gpltable)[which.max(hgnc.frac)],
                         max(hgnc.frac), nrow(gpltable), NA, NA,
                         gpldat@header$distribution,
                         gpldat@header$submission_date))
    print(paste(i, gpl, ": No HGNC symbols."))
  }

  write.table(info.for.file, file="gpls_already_tested.csv",
              append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
}

