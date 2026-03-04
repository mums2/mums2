test_that("annotate_ms_features returns the correct annotations in the
           first n ms features", {
            dir <- "exttestdata"
            q_file <- "matched_data.RDS"
            r_file <- "database_data/PSU-MSMLS.msp"

            dat <- readRDS(test_path(dir, q_file))

            psu_msmls <- read_msp(test_path(dir, r_file))

            annotations <- annotate_ms2(dat, psu_msmls,
                                        modified_cosine_params(0.5), 20,
                                        .2, 0, min_peaks = 0)

            colnames <- c("query_ms1_id", "query_ms2_id", "query_mz",
                          "query_rt", "ref_idx", "query_formula",
                          "chemical_similarity", "score",
                          "num.peaks", "name", "precursormz",
                          "precursortype", "formula", "ontology", "inchikey",
                          "inchi", "smiles", "retentiontime", "ionmode",
                          "instrumenttype", "instrument",
                          "collisionenergy", "comment")
            expect_true(all(colnames %in% colnames(annotations)))
            expect_true(nrow(annotations) > 0)
            expect_s3_class(annotations, "data.frame")
          })

test_that("annotate_ms_features returns the omu where the query is present", {
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  cluster <- cluster_data(distances, dat,  0.3, "opticlust")
  psu_msmls <- read_msp(test_path(dir, r_file))
  annotations <- annotate_ms2(dat, psu_msmls,
                              modified_cosine_params(0.5), 20, .2, 0,
                              min_peaks = 0, cluster_data = cluster)

  expect_true("omu" %in% colnames(annotations))
})

test_that("annotate_ms_features returns the correct 
          of rows and columns", {
            dir <- "exttestdata"
            r_file <- "database_data/PSU-MSMLS.msp"
            dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
            distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
            psu_msmls <- read_msp(test_path(dir, r_file))
            annotations <- annotate_ms2(dat, psu_msmls,
                                        modified_cosine_params(0.5), 20,
                                        .2, 0, min_peaks = 0)

            expect_true(nrow(annotations) == 3)
            expect_true(ncol(annotations) == 23)

            annotations <- annotate_ms2(dat, psu_msmls,
                                        modified_cosine_params(0.5), -1,
                                        .7, 0, min_peaks = 0)

            expect_true(nrow(annotations) == 209)
            expect_true(ncol(annotations) == 23)
          })

test_that("annotate_ms_features works with predicted molecular formulas", {
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat <- compute_molecular_formulas(dat)
  psu_msmls <- read_msp(test_path(dir, r_file))
  annotations <- annotate_ms2(dat, psu_msmls,
                              modified_cosine_params(0.5),
                              1000, .1, 0, min_peaks = 0)
  expect_true("query_formula" %in% colnames(annotations))
})

test_that("annotate_ms2 will return a warning if the annotations are empty", {
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat <- compute_molecular_formulas(dat)
  psu_msmls <- read_msp(test_path(dir, r_file))
  expect_warning(annotate_ms2(dat, psu_msmls,
                              modified_cosine_params(0.5),
                              100, .1, 0, min_peaks = 0))
})

test_that("annotate_ms2 will fail if not supplied a mass_data object", {
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  psu_msmls <- read_msp(test_path(dir, r_file))

  expect_error(annotate_ms2(c(), psu_msmls,
                              modified_cosine_params(0.5),
                              100, .1, 0, min_peaks = 0),
              "The mass_data object must be created using")
})

test_that("annotate_ms2 will fail if not supplied a correct scoring method", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  psu_msmls <- read_msp(test_path(dir, r_file))
  err <-  "The mass_data object must be created using the `ms2_ms1_compare()`"

  expect_error(annotate_ms2(dat, psu_msmls,
                              c(),
                              100, .1, 0, min_peaks = 0),
              "score_params must be created using the")
})

test_that("annotate_ms2 will fail if not supplied correct clustering data", {
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  psu_msmls <- read_msp(test_path(dir, r_file))

  expect_error(annotate_ms2(dat, psu_msmls,
              modified_cosine_params(0.5), 20,
              .2, 0, min_peaks = 0, cluster_data = data.frame()),
              "cluster_data must be created using")
})

test_that("annotate_ms2 will fail if not supplied the correct numeric data", {
  dir <- "exttestdata"
  r_file <- "database_data/PSU-MSMLS.msp"
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  psu_msmls <- read_msp(test_path(dir, r_file))

  expect_error(annotate_ms2(dat, psu_msmls,
                            modified_cosine_params(0.5), "20",
                            .2, 0, min_peaks = 0),
              "ppm")
  
  expect_error(annotate_ms2(dat, psu_msmls,
                            modified_cosine_params(0.5), 20,
                            ".2", 0, min_peaks = 0),
              "min_score")
  
  expect_error(annotate_ms2(dat, psu_msmls,
                            modified_cosine_params(0.5), 20,
                            .2, "0", min_peaks = 0),
              "chemical_min_score")
  
  expect_error(annotate_ms2(dat, psu_msmls,
                            modified_cosine_params(0.5), 20,
                            .2, 0, min_peaks = "0"),
              "min_peaks")
  
  expect_error(annotate_ms2(dat, psu_msmls,
                            modified_cosine_params(0.5), 20,
                            .2, 0, min_peaks = 0,
                            number_of_threads = "1"),
              "number_of_threads")
})

