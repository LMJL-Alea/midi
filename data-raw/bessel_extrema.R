# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8581507/
# Mecholsky, Nicholas (2021), “Extrema for Bessel Functions of the First Kind,
# J”, Mendeley Data, V2, doi: 10.17632/cdhw4wn5sy.2

utils::download.file(
  url = "https://prod-dcd-datasets-cache-zipfiles.s3.eu-west-1.amazonaws.com/cdhw4wn5sy-2.zip",
  destfile = "data-raw/BesselMaxAndMin.zip"
)

utils::unzip(
  zipfile = "data-raw/BesselMaxAndMin.zip",
  exdir = "data-raw"
)

bessel_extrema <- readr::read_tsv(
  file = "data-raw/BesselMaxAndMin.tsv",
  col_names = FALSE
)

bessel_extrema <- as.matrix(bessel_extrema)
colnames(bessel_extrema) <- NULL
rownames(bessel_extrema) <- NULL
bessel_extrema <- bessel_extrema[, 1:1000]

fs::file_delete("data-raw/BesselMaxAndMin.zip")
fs::file_delete("data-raw/BesselMaxAndMin.tsv")

usethis::use_data(
  bessel_extrema,
  internal = TRUE,
  compress = "xz",
  version = 3,
  overwrite = TRUE
)
