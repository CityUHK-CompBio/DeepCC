test_that("batch kernel matches legacy per-sample calculation", {
  set.seed(42)
  S <- 12; G <- 40
  eps <- matrix(rnorm(S * G), nrow = S, ncol = G)
  colnames(eps) <- paste0("G", seq_len(G))
  rownames(eps) <- paste0("S", seq_len(S))
  eps_df <- as.data.frame(eps)

  # Small synthetic gene sets with varying sizes
  geneSets <- list(
    setA = c("G1", "G5", "G20"),
    setB = c("G2", "G2", "G7", "G9"),
    setC = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10"),
    setD = c("G99", "G100"),
    setE = character(0)
  )

  # Reference: old per-sample path
  ref <- do.call(rbind, lapply(seq_len(S), function(i) {
    geneList <- DeepCC:::preprocessGeneList(eps_df[i, ])
    vapply(geneSets, function(gs) {
      DeepCC::calcEnrichmentScore(geneList, gs)
    }, numeric(1))
  }))
  rownames(ref) <- rownames(eps)

  got <- DeepCC::getFunctionalSpectra(eps_df, geneSets = geneSets, scale = FALSE, cores = NULL)
  expect_equal(as.matrix(got), ref, tolerance = 1e-10)
})

test_that("batch kernel matches legacy per-sample calculation with scaling", {
  set.seed(43)
  S <- 8; G <- 30
  eps <- matrix(rnorm(S * G, mean = 5, sd = 2), nrow = S, ncol = G)
  colnames(eps) <- paste0("G", seq_len(G))
  eps_df <- as.data.frame(eps)

  geneSets <- list(setA = c("G1", "G10", "G20"), setB = c("G5", "G15", "G25"))
  ref <- do.call(rbind, lapply(seq_len(S), function(i) {
    centered <- scale(eps_df, scale = FALSE)
    geneList <- DeepCC:::preprocessGeneList(as.data.frame(centered)[i, ])
    vapply(geneSets, function(gs) DeepCC::calcEnrichmentScore(geneList, gs), numeric(1))
  }))

  got <- DeepCC::getFunctionalSpectra(eps_df, geneSets = geneSets, scale = TRUE, cores = NULL)
  expect_equal(as.matrix(got), ref, tolerance = 1e-10)
})

test_that("geneSets string shortcuts are removed with clear error", {
  expect_error(DeepCC:::resolveGeneSets("MSigDBv7"), "have been removed")
  expect_error(DeepCC:::resolveGeneSets("MSigDBv5"), "have been removed")
  expect_error(DeepCC:::resolveGeneSets("MSigDBv6"), "have been removed")
})

test_that("single sample path uses same scoring kernel", {
  set.seed(44)
  ep <- rnorm(50)
  names(ep) <- paste0("G", seq_len(50))
  geneSets <- list(setA = c("G1", "G10", "G20"), setB = c("G5", "G15"))
  gs <- DeepCC:::resolveGeneSets(geneSets)
  geneList <- DeepCC:::preprocessGeneList(ep)
  ref <- vapply(seq_along(gs), function(i) DeepCC::calcEnrichmentScore(geneList, gs[[i]]), numeric(1))
  names(ref) <- names(gs)
  got <- DeepCC::getFunctionalSpectrum(ep, geneSets = geneSets, logChange = TRUE)
  expect_equal(got, ref, tolerance = 1e-12)
})
