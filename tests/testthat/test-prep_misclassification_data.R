# Pins current correct-path behavior of prep_misclassification_data():
# balanced-table invariants, zero-padding, weights, controls, and binning.

test_that("balanced table invariants hold without controls", {
  md = make_categorical_microdata()

  out = prep_misclassification_data(
    data = md$micro,
    outcome_1 = "y1var",
    outcome_2 = "y2var",
    regressor = "xvar",
    X_names = md$levels,
    Y1_names = md$levels,
    Y2_names = md$levels,
    record_vals = FALSE
  )

  # Returned object structure
  expect_named(out, c("tab", "J", "K", "X_names", "Y1_names", "Y2_names",
                      "X_vals", "Y_vals", "W_names"))
  expect_equal(out$J, 3)
  expect_equal(out$K, 3)
  expect_equal(out$X_names, md$levels)
  expect_equal(out$Y1_names, md$levels)
  expect_equal(out$Y2_names, md$levels)
  expect_identical(out$X_vals, NA)
  expect_identical(out$Y_vals, NA)
  expect_identical(out$W_names, NA)

  # Balanced J^2 * K table with no duplicated cells
  expect_equal(nrow(out$tab), 3^2 * 3)
  expect_named(out$tab, c("X", "Y1", "Y2", "n"))
  expect_false(any(duplicated(out$tab[, c("X", "Y1", "Y2")])))

  # Counts match the microdata exactly (unit weights)
  expect_equal(sum(out$tab$n), nrow(md$micro))
  expect_equal(out$tab$n, expected_combo_counts(out$tab, md$combos))

  # Rows are ordered by (Y2, Y1, X) with X varying fastest
  expect_equal(as.character(out$tab$X), rep(md$levels, times = 9))
  expect_equal(as.character(out$tab$Y1), rep(rep(md$levels, each = 3), times = 3))
  expect_equal(as.character(out$tab$Y2), rep(md$levels, each = 9))
})

test_that("row order follows the user-supplied (non-alphabetical) name order", {
  lev = c("c", "a", "b")
  md = make_categorical_microdata(levels = lev)

  out = prep_misclassification_data(
    data = md$micro,
    outcome_1 = "y1var",
    outcome_2 = "y2var",
    regressor = "xvar",
    X_names = lev,
    Y1_names = lev,
    Y2_names = lev,
    record_vals = FALSE
  )

  expect_equal(as.character(out$tab$X), rep(lev, times = 9))
  expect_equal(as.character(out$tab$Y1), rep(rep(lev, each = 3), times = 3))
  expect_equal(as.character(out$tab$Y2), rep(lev, each = 9))
  expect_equal(out$tab$n, expected_combo_counts(out$tab, md$combos))
})

test_that("empty cells are zero-padded so the table stays balanced", {
  md = make_categorical_microdata()
  drop = md$micro$xvar == "a" & md$micro$y1var == "b" & md$micro$y2var == "c"
  micro = md$micro[!drop, ]

  out = prep_misclassification_data(
    data = micro,
    outcome_1 = "y1var",
    outcome_2 = "y2var",
    regressor = "xvar",
    X_names = md$levels,
    Y1_names = md$levels,
    Y2_names = md$levels,
    record_vals = FALSE
  )

  expect_equal(nrow(out$tab), 27)
  hole = out$tab$n[as.character(out$tab$X) == "a" &
                     as.character(out$tab$Y1) == "b" &
                     as.character(out$tab$Y2) == "c"]
  expect_equal(hole, 0)
  expect_equal(sum(out$tab$n), nrow(md$micro) - sum(drop))
})

test_that("weights are summed within cells", {
  md = make_categorical_microdata()
  micro = md$micro
  # Alternating weights that sum exactly to nrow(micro): no normalization warning
  micro$w = rep(c(0.5, 1.5), length.out = nrow(micro))

  out = prep_misclassification_data(
    data = micro,
    outcome_1 = "y1var",
    outcome_2 = "y2var",
    regressor = "xvar",
    weights = "w",
    X_names = md$levels,
    Y1_names = md$levels,
    Y2_names = md$levels,
    record_vals = FALSE
  )

  expect_equal(sum(out$tab$n), sum(micro$w))
  expected = aggregate(w ~ xvar + y1var + y2var, data = micro, FUN = sum)
  key_tab = paste(as.character(out$tab$X), as.character(out$tab$Y1),
                  as.character(out$tab$Y2))
  key_exp = paste(expected$xvar, expected$y1var, expected$y2var)
  expect_equal(out$tab$n, expected$w[match(key_tab, key_exp)])
})

test_that("controls split the data into per-cell balanced tables", {
  md = make_categorical_microdata()
  micro_u = md$micro
  micro_u$ctrl = "u"
  micro_v = md$micro
  micro_v$ctrl = "v"
  micro = rbind(micro_u, micro_v)

  # Cells have 378 observations each, so the small-cell warning must fire
  expect_warning(
    out <- prep_misclassification_data(
      data = micro,
      outcome_1 = "y1var",
      outcome_2 = "y2var",
      regressor = "xvar",
      controls = "ctrl",
      X_names = md$levels,
      Y1_names = md$levels,
      Y2_names = md$levels,
      record_vals = FALSE
    ),
    "less than 1000"
  )

  expect_equal(out$W_names, c("u", "v"))
  expect_true(is.list(out$tab))
  expect_named(out$tab, c("u", "v"))
  expect_equal(unlist(out$J), c(u = 3, v = 3))
  expect_equal(unlist(out$K), c(u = 3, v = 3))
  for (cell in out$tab) {
    expect_equal(nrow(cell), 27)
    expect_false(any(duplicated(cell[, c("X", "Y1", "Y2")])))
    expect_equal(sum(cell$n), nrow(md$micro))
    expect_equal(cell$n, expected_combo_counts(cell, md$combos))
  }
})

test_that("multiple control columns are collapsed with underscores", {
  md = make_categorical_microdata()
  micro_u = md$micro
  micro_u$c1 = "u"
  micro_v = md$micro
  micro_v$c1 = "v"
  micro = rbind(micro_u, micro_v)
  micro$c2 = "z"

  expect_warning(
    out <- prep_misclassification_data(
      data = micro,
      outcome_1 = "y1var",
      outcome_2 = "y2var",
      regressor = "xvar",
      controls = c("c1", "c2"),
      X_names = md$levels,
      Y1_names = md$levels,
      Y2_names = md$levels,
      record_vals = FALSE
    ),
    "less than 1000"
  )

  expect_equal(out$W_names, c("u_z", "v_z"))
})

test_that("continuous variables are weighted-averaged within bins", {
  combos = expand.grid(xb = 1:3, y1b = 1:3, y2b = 1:3, KEEP.OUT.ATTRS = FALSE)
  micro = combos[rep(seq_len(nrow(combos)), each = 3), ]
  rownames(micro) = NULL
  micro$offset = rep(c(-5, 0, 5), times = nrow(combos))
  # Weights sum to nrow(micro) = 81; weighted mean offset within each X bin is
  # (1.8*(-5) + 0.6*0 + 0.6*5) / 3 = -2
  micro$w = rep(c(1.8, 0.6, 0.6), times = nrow(combos))
  micro$xc = 100 * micro$xb + micro$offset
  micro$y1c = 10 * micro$y1b
  micro$y2c = 10 * micro$y2b

  out = prep_misclassification_data(
    data = micro,
    outcome_1 = "y1c",
    outcome_2 = "y2c",
    regressor = "xc",
    outcome_1_bin = "y1b",
    outcome_2_bin = "y2b",
    regressor_bin = "xb",
    weights = "w",
    record_vals = TRUE
  )

  expect_equal(out$J, 3)
  expect_equal(out$K, 3)
  expect_equal(nrow(out$tab), 27)

  # Binned values: X = 100*b - 2 (weighted), Y1 = Y2 = 10*b (exact)
  expect_equal(out$X_vals, c(98, 198, 298))
  expect_equal(out$Y_vals, c(10, 20, 30))
  expect_equal(out$X_names, c("98", "198", "298"))
  expect_equal(out$Y1_names, c("10", "20", "30"))
  expect_equal(out$Y2_names, c("10", "20", "30"))

  # Each combination carries its summed weight
  expect_equal(out$tab$n, rep(3, 27))
  expect_equal(sum(out$tab$n), sum(micro$w))

  # Order: X fastest, then Y1, then Y2 (numeric values here)
  expect_equal(out$tab$X, rep(c(98, 198, 298), times = 9))
  expect_equal(out$tab$Y1, rep(rep(c(10, 20, 30), each = 3), times = 3))
  expect_equal(out$tab$Y2, rep(c(10, 20, 30), each = 9))
})

test_that("input validation catches malformed calls", {
  md = make_categorical_microdata()

  # NAs in a relevant column
  bad = md$micro
  bad$y1var[1] = NA
  expect_error(
    prep_misclassification_data(
      data = bad, outcome_1 = "y1var", outcome_2 = "y2var", regressor = "xvar",
      X_names = md$levels, Y1_names = md$levels, Y2_names = md$levels,
      record_vals = FALSE),
    "NA values"
  )

  # record_vals = FALSE requires names
  expect_error(
    prep_misclassification_data(
      data = md$micro, outcome_1 = "y1var", outcome_2 = "y2var",
      regressor = "xvar", record_vals = FALSE),
    "should be provided"
  )

  # record_vals = TRUE forbids names
  expect_error(
    prep_misclassification_data(
      data = md$micro, outcome_1 = "y1var", outcome_2 = "y2var",
      regressor = "xvar", X_names = md$levels, Y1_names = md$levels,
      Y2_names = md$levels, record_vals = TRUE),
    "should not be provided"
  )

  # record_vals = TRUE requires numeric variables
  expect_error(
    prep_misclassification_data(
      data = md$micro, outcome_1 = "y1var", outcome_2 = "y2var",
      regressor = "xvar", record_vals = TRUE),
    "numeric"
  )

  # Partial binning arguments
  expect_error(
    prep_misclassification_data(
      data = md$micro, outcome_1 = "y1var", outcome_2 = "y2var",
      regressor = "xvar", outcome_1_bin = "y1var",
      X_names = md$levels, Y1_names = md$levels, Y2_names = md$levels,
      record_vals = FALSE),
    "all should be provided"
  )

  # Values absent from the supplied names
  bad = md$micro
  bad$xvar[1] = "d"
  expect_error(
    prep_misclassification_data(
      data = bad, outcome_1 = "y1var", outcome_2 = "y2var", regressor = "xvar",
      X_names = md$levels, Y1_names = md$levels, Y2_names = md$levels,
      record_vals = FALSE),
    "NAs introduced"
  )
})

test_that("reserved column names trigger a warning but still tabulate", {
  md = make_categorical_microdata()
  micro = md$micro
  colnames(micro)[colnames(micro) == "xvar"] = "X"

  expect_warning(
    out <- prep_misclassification_data(
      data = micro, outcome_1 = "y1var", outcome_2 = "y2var", regressor = "X",
      X_names = md$levels, Y1_names = md$levels, Y2_names = md$levels,
      record_vals = FALSE),
    "will be overwritten"
  )
  expect_equal(nrow(out$tab), 27)
  expect_equal(sum(out$tab$n), nrow(micro))
})
