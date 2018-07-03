##Script to fetch HGNC symbols from GEO for all GPL in column 1 of
##platforms.csv.

library(GEOquery)
all.platforms <- read.csv("platform.csv", as.is=TRUE)

library(org.Hs.eg.db)
hgnc.reference <- as.character(org.Hs.egSYMBOL)
names(hgnc.reference) <- NULL

library(HGNChelper)



raw.platform.dir <- "platforms"
dir.create(raw.platform.dir)

hgnc.vec.dir <- "hgnc.vecs"
dir.create(hgnc.vec.dir)

if(file.exists("gpls_already_tested.csv")){
    gpls.already.tested <- read.csv("gpls_already_tested.csv", header=TRUE)
}else{
    write.table(t(c("platform", "colname", "frac.hgnc", "nrow", "valid.frac", "valid.after.hgnchelper.frac", "distribution", "submission_date")), file="gpls_already_tested.csv", row.names=FALSE, col.names=FALSE, sep=",")
}

for (i in 1:nrow(all.platforms)){
    gpl <- all.platforms[i, 1]
    hgnc.vec.file <- paste(hgnc.vec.dir, "/", gpl, "_hgnc.vec.RData", sep="")
    if(exists("gpls.already.tested") && gpl %in% gpls.already.tested$platform) next
    gpldat <- try(getGEO(gpl, destdir=raw.platform.dir))
    if(class(gpldat) == "try-error") next
    gpltable <- try(Table(gpldat))
    if(class(gpltable) == "try-error") next
    hgnc.frac <- apply(gpltable, 2, function(x) sum(unique(x) %in% hgnc.reference) / length(unique(x)))
    if(any(hgnc.frac > 0)){
        hgnc.vec <- unique(as.character(gpltable[, which.max(hgnc.frac)]))
        hgnc.vec <- gsub("[ ].+", "", hgnc.vec)  ##get rid of anything after a space
        HGNChelper.output <- checkGeneSymbols(iconv(hgnc.vec, "latin1", "ASCII", "")) #convert to ascii
        valid.frac <- sum(HGNChelper.output$Approved) / length(hgnc.vec)
        after.HGNChelper.valid.frac <- sum(!is.na(HGNChelper.output$Suggested.Symbol)) / length(hgnc.vec)
        hgnc.vec <- c(gpldat@header$distribution, gpldat@header$submission_date, hgnc.vec)
        save(hgnc.vec, file=hgnc.vec.file, compress="bzip2")
        info.for.file <- t(c(gpl, colnames(gpltable)[which.max(hgnc.frac)], max(hgnc.frac), nrow(gpltable), valid.frac, after.HGNChelper.valid.frac, gpldat@header$distribution, gpldat@header$submission_date))
        print(paste(i, gpl, length(hgnc.vec), "HGNC symbols found and saved."))
    }else{
       info.for.file <- t(c(gpl, colnames(gpltable)[which.max(hgnc.frac)], max(hgnc.frac), nrow(gpltable), NA, NA, gpldat@header$distribution, gpldat@header$submission_date))
       print(paste(i, gpl, ": No HGNC symbols."))
    }
    write.table(info.for.file, file="gpls_already_tested.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
}

