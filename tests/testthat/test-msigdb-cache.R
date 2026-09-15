test_that("bundled MSigDB cache resolves by name and by explicit version", {
  sets <- DeepCC:::resolveGeneSets("MSigDB")
  expect_type(sets, "list")
  expect_gt(length(sets), 20000)

  version <- attr(sets, "db_version")
  expect_type(version, "character")
  expect_match(version, "^[0-9]{4}\\.[0-9]+\\.Hs$")

  expect_identical(DeepCC:::resolveGeneSets(paste0("MSigDB_", version)), sets)
})

test_that("bundled cache records provenance and collection mapping", {
  sets <- DeepCC:::resolveGeneSets("MSigDB")
  expect_identical(attr(sets, "source"), "msigdbr")
  collections <- attr(sets, "collections")
  expect_length(collections, length(sets))
  expect_true(all(c("H", "C2", "C5") %in% collections))
})

test_that("unknown geneSets names fail with an actionable error", {
  expect_error(DeepCC:::resolveGeneSets("MSigDBv7"), "Unknown geneSets value")
  expect_error(DeepCC:::resolveGeneSets("MSigDB_1999.1.Hs"), "Unknown geneSets value")
})

test_that("functional spectra work offline with the bundled cache", {
  set.seed(7)
  eps <- matrix(rnorm(4 * 500), nrow = 4, ncol = 500)
  colnames(eps) <- as.character(sample(1:30000, 500))
  geneSets <- list(setA = colnames(eps)[1:20], setB = colnames(eps)[50:80])

  fs <- getFunctionalSpectra(as.data.frame(eps), geneSets = geneSets, scale = FALSE)
  expect_equal(dim(fs), c(4L, 2L))
  expect_identical(colnames(fs), c("setA", "setB"))
})

test_that("get_msigdbr reports how to proceed when msigdbr is missing", {
  testthat::local_mocked_bindings(
    msigdbrAvailable = function() FALSE,
    .package = "DeepCC"
  )
  expect_error(get_msigdbr(), "bundled snapshot")
})
