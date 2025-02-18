#' Generate simulated data

cells <- simulate_cells()
write.csv(cells, "data-raw/simulated.csv", row.names = FALSE, quote = FALSE)
spe <- make_spe_from_expr_data(expression_file = "data-raw/simulated.csv",
                               hierarchy_file = "inst/extdata/hierarchy.yaml")
usethis::use_data(spe, overwrite = TRUE)
