library(testthat)
library(modstatespace)

if (Sys.getenv("NOT_CRAN") == "true") {
  test_check("modstatespace")
}
