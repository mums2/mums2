limit_cores <- function() {
  Sys.setenv("OMP_THREAD_LIMIT" = 1L)
  options("Ncpu" = 1L)
  data.table::setDTthreads(threads = 1L)
}
