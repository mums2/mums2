# download the PSU-MSMLS database from GNPS using massdatabase
massdatabase::download_gnps_spectral_library(gnps_library = "PSU-MSMLS", path = here::here("tests/testthat/exttestdata/database_data"))
# read the database msp file into R as S3 class "gnps_ref"
psu_msmls <- massdatabase::read_msp_data_gnps(file = here::here("tests/testthat/exttestdata/database_data/PSU-MSMLS.msp"))
db2 <- massdatabase::read_msp_data(file = here::here("tests/testthat/exttestdata/database_data/PSU-MSMLS.msp"), source = "gnps")

# length(psu_msmls)
# psu_msmls[[1]]
# names(psu_msmls[[1]])

saveRDS(psu_msmls[[1]], file = here::here("tests/testthat/exttestdata/database_data/psu_msmls_oneref"))

massdatabase:: download_massbank_compound(url = "https://github.com/MassBank/MassBank-data/releases/download/2021.12",
  source = c("nist", "riken"), path = here::here("tests/testthat/exttestdata/database_data")
)
massbank <- massdatabase::read_msp_data(file = here::here("tests/testthat/exttestdata/database_data/massbank_compound/MassBank_NIST.msp"))
# length(massbank)
# massbank[[1]]
# names(massbank[[1]])

# psu_msmls[[1]]$info$key
# massbank[[1]]$info$key

# library(data.table)
# psu_msmls[[1]]$info[psu_msmls[[1]]$info %like% "precursorMZ", ]

# class(psu_msmls[[1]]$info)
# like(psu_msmls[[1]]$info$key, "precursorMZ")

# tolower(psu_msmls[[1]]$info$key)
# tolower(massbank[[1]]$info$key)

# pmz_col <- grep("precursorMZ", psu_msmls[[1]]$info$key, ignore.case = TRUE)
# pmz_col <- grep("precursorMZ", massbank[[1]]$info$key, ignore.case = TRUE)

# usethis::use_data(psu_msmls, overwrite = TRUE)
