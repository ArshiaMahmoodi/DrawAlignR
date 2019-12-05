context("Visualization functions")
library(DrawAlignR)

test_that("Plot data type is null", {
  p <- NULL
  #p <- plot_chrom_reference(chrom12_GR.mzML_ChromSelected, 0, 0, 0, 0, 0, 0, "test")
  expect_true(is.null(p))
})

test_that("Test that input values are valid", {
  #chrom <- chrom13_GR.mzML_ChromSelected
  precursor <- 960
  Run_ID <- 13
  RT <- 2420.42
  Left_width <- 2404.61
  Right_width <- 2440.96
  mz <- 630.2936
  sequence <- "VGEGTY(Phospho)GVVYK(Label:13C(6)15N(2))"

  #p <- plot_chrom_reference(chrom, precursor, Run_ID, RT, Left_width, Right_width, mz, sequence)

  expect_true(is.numeric(precursor))
  expect_true(precursor > 0)

  expect_true(is.numeric(Run_ID))
  expect_true(Run_ID > 0)

  expect_true(is.numeric(RT))

  expect_true(is.numeric(Left_width))
  expect_true(is.numeric(Right_width))

  expect_true(Right_width > RT)
  expect_true(Left_width < RT)
  expect_true(Left_width < Right_width)


  expect_true(is.numeric(mz))
  expect_true(mz > 0)

  expect_true(is.character(sequence))
})
