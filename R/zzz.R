.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(
      "Please cite our software :) \n \n",
      "Sehyun Oh et al.",
      "HGNChelper: identification and correction of invalid gene symbols for human and mouse.",
      "F1000Research 2020, 9:1493.",
      "DOI: https://doi.org/10.12688/f1000research.28033.1 \n \n",
      "Type `citation('HGNChelper')` for a BibTeX entry."
    )
  )
  invisible()
}