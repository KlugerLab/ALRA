A <- b_nk_example
labels <- labels_example

test_that("alra works", {
  # Library and log normalize the data
  A_norm <- normalize_data(A)

  # Choose k. 
  k_choice <- choose_k(A_norm)

  # run ALPRA
  A_norm_completed <- alra(A_norm,k=k_choice$k)[[3]]

  #NCAM1 should be expressed in all NK cells, but is only expressed in 4% of of NK cells in the original data.  
  expect_lte(mean(A_norm[labels == "cd56_nk", "NCAM1", drop = FALSE] > 0), .1)
  expect_gte(mean(A_norm_completed[labels == "cd56_nk", "NCAM1", drop = FALSE] > 0), .5)

  #CR2 should be expressed in all B cells, but is only expressed in 1% of of B cells in the original data.
  expect_lte(mean(A_norm[labels == "b_cells", "CR2", drop = FALSE] > 0), .1)
  expect_gte(mean(A_norm_completed[labels == "b_cells", "CR2", drop = FALSE] > 0), .5)
})
