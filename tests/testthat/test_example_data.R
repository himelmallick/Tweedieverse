test_that("bundled real-data containers retain source groups and full-library factors", {
  examples <- new.env()
  load(
    system.file("extdata", "tweedieverse_data.RData", package = "Tweedieverse"),
    envir = examples
  )
  for (object_name in ls(examples)) {
    rds_file <- system.file(
      "extdata", paste0(object_name, ".rds"), package = "Tweedieverse"
    )
    expect_identical(readRDS(rds_file), examples[[object_name]])
  }
  microbiome <- examples$microbiome_tse
  umi <- examples$umi_sce
  nonumi <- examples$nonumi_sce
  expect_s4_class(microbiome, "TreeSummarizedExperiment")
  expect_s4_class(umi, "SingleCellExperiment")
  expect_s4_class(nonumi, "SingleCellExperiment")
  expect_equal(dim(SummarizedExperiment::assay(microbiome)), c(20L, 48L))
  expect_equal(dim(umi), c(20L, 60L))
  expect_equal(dim(nonumi), c(20L, 60L))
  expect_equal(as.integer(table(microbiome$diagnosis)), rep(16L, 3))
  expect_equal(levels(microbiome$diagnosis), c("nonIBD", "CD", "UC"))
  expect_equal(as.integer(table(umi$CellType)), c(30L, 30L))
  expect_equal(levels(umi$CellType), c("CD8", "CD4"))
  expect_equal(as.integer(table(nonumi$celltype)), c(30L, 30L))
  expect_setequal(levels(nonumi$celltype), c("astrocytes", "oligodendrocytes"))
  expect_s4_class(SummarizedExperiment::assay(umi), "sparseMatrix")
  expect_true(all(
    umi$library_size >= Matrix::colSums(SummarizedExperiment::assay(umi))
  ))
  expect_true(all(is.finite(nonumi$source_size_factor) & nonumi$source_size_factor > 0))
})
