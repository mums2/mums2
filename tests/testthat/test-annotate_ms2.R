test_that("annotate_ms_featrues returns the correct annotations in the
           first n ms features", {
            dir <- "exttestdata"
            q_file <- "matched_data.RDS"
            r_file <- "database_data/PSU-MSMLS.msp"

            dat <- readRDS(test_path(dir, q_file))

            psu_msmls <- read_msp(test_path(dir, r_file))

            annotations <- annotate_ms2(dat, psu_msmls,
                                        gnps_params(0.5), 2, .2)

            colnames <- c("query_ms1_id", "query_ms2_id", "query_mz",
                          "query_rt", "ref_idx", "score", "NAME", "PRECURSORMZ",
                          "PRECURSORTYPE", "FORMULA", "Ontology", "INCHIKEY",
                          "INCHI", "SMILES", "RETENTIONTIME", "IONMODE",
                          "INSTRUMENTTYPE", "INSTRUMENT",
                          "COLLISIONENERGY", "Comment", "Num Peaks")
            expect_equal(colnames(annotations), colnames)
            expect_true(nrow(annotations) > 0)
            expect_s3_class(annotations, "data.frame")
})

############################################
###             Helper funs              ###
############################################

test_that("get_ref_precursor works correctly", {
  dir <- "exttestdata/database_data"
  file <- "psu_msmls_oneref"
  single_ref <- readRDS(test_path(dir, file))

  pmz <- get_ref_precursor(single_ref)

  expect_equal(pmz, 148.061)
})
