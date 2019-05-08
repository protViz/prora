context("Run scripts")

test_that("run script for sigora_longformat works", {
  source("../../inst/RunScripts/sigora_longformat.R")
})

test_that("run script for webGestaltR_ORA_longformat works", {
  source("../../inst/RunScripts/webGestatR_ORA_longformat.R")
})

test_that("run script for GSEA works", {
  source("../../inst/RunScripts/GSEA_longformat.R")
})
