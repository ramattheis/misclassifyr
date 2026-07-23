# Deterministic microdata builders shared across tests.
# No RNG: every (X, Y1, Y2) combination i (in expand.grid order) appears
# exactly `i` times, so expected tabulations are known in closed form.

make_categorical_microdata = function(levels = c("a", "b", "c")) {
  combos = expand.grid(xvar = levels, y1var = levels, y2var = levels,
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  combos$count = seq_len(nrow(combos))
  micro = combos[rep(seq_len(nrow(combos)), combos$count),
                 c("xvar", "y1var", "y2var")]
  rownames(micro) = NULL
  list(micro = micro, combos = combos, levels = levels)
}

# Expected counts for the rows of a prep_misclassification_data() tabulation,
# looked up from the combo table above.
expected_combo_counts = function(tab, combos) {
  key = paste(as.character(tab$X), as.character(tab$Y1), as.character(tab$Y2))
  combos$count[match(key, paste(combos$xvar, combos$y1var, combos$y2var))]
}
