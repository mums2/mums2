test_that("annotate_ms_featrues returns the correct annotations in the
           first n ms features", {
            dir <- "exttestdata"
            q_file <- "matched_data.RDS"
            r_file <- "database_data/PSU-MSMLS.msp"

            dat <- readRDS(test_path(dir, q_file))

            psu_msmls <- read_msp(test_path(dir, r_file))

            annotations <- annotate_ms2(dat, psu_msmls,
                                        gnps_params(0.5), 2, .2, min_peaks = 0)

            colnames <- c("query_ms1_id", "query_ms2_id", "query_mz",
                          "query_rt", "ref_idx", "query_formula", "chemical_similarity", "score", "NAME", "PRECURSORMZ",
                          "PRECURSORTYPE", "FORMULA", "Ontology", "INCHIKEY",
                          "INCHI", "SMILES", "RETENTIONTIME", "IONMODE",
                          "INSTRUMENTTYPE", "INSTRUMENT",
                          "COLLISIONENERGY", "Comment", "Num Peaks")
            expect_equal(colnames(annotations), colnames)
            expect_true(nrow(annotations) > 0)
            expect_s3_class(annotations, "data.frame")
})

test_that("annotate_ms_featrues returns the omu where the query is present", {
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  cluster <- cluster_data(distances, dat,  0.3, "opticlust")
  psu_msmls <- read_msp(test_path(dir, r_file))
  annotations <- annotate_ms2(dat, psu_msmls,
    gnps_params(0.5), 2, .2, min_peaks = 0, cluster_data = cluster)
  
  expect_true("OMU" %in% colnames(annotations))
})

test_that("annotate_ms_featrues returns the correct amount of rows and columns", {
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  psu_msmls <- read_msp(test_path(dir, r_file))
  annotations <- annotate_ms2(dat, psu_msmls,
    gnps_params(0.5), 2, .2, min_peaks = 0)
  
  expect_true(nrow(annotations) == 50)
  expect_true(ncol(annotations) == 23)
})
