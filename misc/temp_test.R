# Temporary script to debug corner case

test = synthetic_data(dgp_delta = "No error", J=5, K=5, sample_size = 1e6)
test_tab = test$tab

test_out = misclassifyr(
  test_tab, J = 5, K = 5,
  X_names = c("a","b","c","d","e"),
  Y1_names = c("A","B","C","D","E"),
  Y2_names = c("A","B","C","D","E"),
  W_names = c("w1","w2"),
  bayesian = T
)

