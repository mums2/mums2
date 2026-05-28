limit_cores <- function() {
  Sys.setenv("OMP_THREAD_LIMIT" = 1)
  data.table::setDTthreads(threads = 1)
}
