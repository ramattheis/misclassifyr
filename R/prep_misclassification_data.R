#' prep_misclassification_data
#'
#' This package provides a menu of options for estimation and inference of misclassification models in which the analyst has access to two noisy measures, `Y1` and `Y2` of a latent outcome `Y*`, a correctly measured covariate `X`, and discrete controls `W`.
#'
#' @param data A data.frame containing the outcome variable,
#' @param outcome_1 A character string denoting the variable in the dataframe to be used as the first measure of an outcome, Y_1.
#' @param outcome_2 A character string denoting the variable in the dataframe to be used as the second measure of an outcome, Y_2.
#' @param regressor A character string denoting the variable in the dataframe to be used as a regressor, X.
#' @param controls A character string or vector of character strings denoting the variable/variables to be used as non-parametric controls, W.
#' @param weights A character string denoting a variable containing individual level weights
#' @param record_vals A logical value indicating whether to record the unique values of the outcomes and the regressor.
#' @return A list of objects including tabulated data to be used in misclassifyr()
#' @export
prep_misclassification_data <- function(
    data,
    outcome_1,
    outcome_2,
    regressor,
    controls = NA,
    weights = NA,
    record_vals = F) {

  #------------------------------------------------------------
  # Catching input errors
  #------------------------------------------------------------

  # Data should be a data.frame
  if(!is.data.frame(data)){stop("Eror: data should be a data.frame.")}

  # Throwing an error if any of the entered column names are redundant
  keep = c(regressor, outcome_1, outcome_2)
  if(!identical(controls, NA)){ keep = c(keep, controls) }
  if(!is.na(weights)){ keep = c(keep, weights)}
  if(any(duplicated(keep))){stop("Error: column name arguments should not be repeated.")}

  # Variable names should be present in 'data'
  if(!(outcome_1 %in% colnames(data))){stop("Error: outcome_1 should be a column in data.")}
  if(!(outcome_2 %in% colnames(data))){stop("Error: outcome_2 should be a column in data.")}
  if(!(regressor %in% colnames(data))){stop("Error: regressor should be a column in data.")}
  if(!identical(controls,NA)){if(!all(controls %in% colnames(data))){stop("Error: controls should be a column / columns in data.")}}
  if(!is.na(weights)){if(!(weights %in% colnames(data))){stop("Error: weights should be a column in data.")}}

  # Dropping data to the required columns
  data = data[,keep]

  # Passing a warning if any of the reserved column names (generated in this function) are in use
  if(("weight" %in% colnames(data))){warning("A column named 'weight' is present in data and will be overwritten")}
  if(("W" %in% colnames(data))){warning("A column named 'W' is present in data and will be overwritten")}
  if(("X" %in% colnames(data))){warning("A column named 'X' is present in data and will be overwritten")}
  if(("Y1" %in% colnames(data))){warning("A column named 'Y1' is present in data and will be overwritten")}
  if(("Y2" %in% colnames(data))){warning("A column named 'Y2' is present in data and will be overwritten")}

  # record_vals should be logical
  if(!is.logical(record_vals)){stop("Error: record_vals should be logical.")}

  # If record_vals is true, the outcome and regressor variables should be numeric in data
  if(record_vals){
    if(!is.numeric(unlist(data[,outcome_1]))){stop("Error: if record_vals == T, ourcome_1 in data should be numeric.")}
    if(!is.numeric(unlist(data[,outcome_2]))){stop("Error: if record_vals == T, ourcome_2 in data should be numeric.")}
    if(!is.numeric(unlist(data[,regressor]))){stop("Error: if record_vals == T, regressor in data should be numeric.")}
  }

  # If provided, weights should always be numeric
  if(!is.na(weights)){if(any(is.na(as.numeric(data[,weights])))){stop("Error: if provided, weights should be a numeric vector in data.")} }

  # Checking whether the sum of weights is equal to the number of rows in the data, sending a warning if not
  if(!is.na(weights)){if(abs(sum(data[,weights]) - nrow(data)) > 1){warning("The sum of weights is not equal to the number of rows in data. Are you sure the weights are normalized correctly?") } }

  #------------------------------------------------------------
  # Defining weight, remaming columns, and collapsing controls to a single column
  #------------------------------------------------------------

  # Generating a new weights column
  if(!is.na(weights)){
    data$weight = data[,weights]
  } else {
    data$weight = 1
  }

  # Defining a control vector W if controls isn't NA
  if(!identical(controls, NA)){
    if(length(controls) > 1){
      data$W = apply(data[,controls],1,paste, collapse = "_")
    } else {
      data$W = as.character(data[,controls])
    }
  }

  # Renaming the outcomes and the regressor
  data$X = data[,regressor]
  data$Y1 = data[,outcome_1]
  data$Y2 = data[,outcome_2]

  # Dropping to necessary columns
  keep = c("X","Y1","Y2","weight")
  if(!identical(controls,NA)){keep = c(keep, "W")}
  data = data[,c(keep)]

  #------------------------------------------------------------
  # Tabulating
  #------------------------------------------------------------

  # Defining a function to tabulate one cell of data
  tabulate = function(cell){

    # Throwing an error if the number of unique values in Y1 and Y2 is different
    if(length(unique(cell$Y1)) != length(unique(cell$Y2))){stop("Error: Y1 and Y2 must be the same dimension.")}

    # Recording the names of the outcome, the regressor, and the controls
    X_names = unique(cell$X)
    Y_names = unique(cell$Y1)

    # Recording the values of the outcome and the regressor
    if(record_vals){
      Y_vals = cell$Y1 |> unique() |> sort(decreasing = F)
      X_vals = cell$X |> unique() |> sort(decreasing = F)
    } else {
      Y_vals = NA
      X_vals = NA
    }

    # Ensuring that the tabulation is balanced across X, Y1, Y2
    empty_cell = expand.grid(unique(cell$X), unique(cell$Y1), unique(cell$Y2))
    colnames(empty_cell) = c("X","Y1","Y2")
    empty_cell$weight = 0
    empty_cell$W = NA
    empty_cell = empty_cell[,colnames(cell)]
    cell = rbind(cell,empty_cell)

    # Tabulating the cell by regressor and outcome
    tab = cell |>
      dplyr::arrange(X,Y1,Y2) |>
      dplyr::group_by(X,Y1,Y2) |>
      dplyr::summarise(n = sum(weight),.groups = "drop") |>
      as.data.frame()

    # Returning tabulations, names, and values
    return(list(
     tab = tab,
     J = length(Y_names),
     K = length(X_names),
     X_names = X_names,
     Y_names = Y_names,
     X_vals = X_vals,
     Y_vals = Y_vals
    ))
  }

  if(identical(controls,NA)){

    # No controls provided, tabulating the full population
    out = tabulate(data)

    # Adding empty entry for W_names
    out = append(out, list(W_names = NA))

  } else {

    # Splitting the data by control cell
    cells = split(data, data$W)

    # Recording W_names
    W_names = names(cells)

    # Tabulating within each cell
    out = lapply(cells, tabulate)

    # Splitting the output into separate lists
    tab_list = lapply(out, "[[", 1)
    J_list = lapply(out,"[[",2)
    K_list = lapply(out,"[[",3)
    X_names_list = lapply(out,"[[",4)
    Y_names_list = lapply(out,"[[",5)
    X_vals_list = lapply(out,"[[",6)
    Y_vals_list = lapply(out,"[[",7)

    # Returning warnings if control cell sizes are small
    cell_sizes = sapply(tab_list, function(tab) sum(tab$n))
    if(min(cell_sizes) < 25){
      warning("Smallest control cell size is less than 25; within-cell estimates will be unreliable.")
    } else if(min(cell_sizes) < 100){
      warning("Smallest control cell size is less than 100; *strongly* consider coarsening controls.")
    } else if(min(cell_sizes) < 1000){
      warning("Smallest control cell size is less than 1000; consider coarsening controls.")
    }

    # Returning warnings if J and K are inconsistent across cells
    if(length(unique(unlist(J_list))) != 1 | length(unique(unlist(J_list))) != 1){
      warning("The dimension of X or Y is inconsistent across cells. Either coarsen cells or allow the dimension of phi / psi to vary across cells.")
      }

    # Recombining into a list of lists
    out = list(
      tab = tab_list,
      J = J_list,
      K = K_list,
      X_names = X_names_list,
      Y_names = Y_names_list,
      X_vals = X_vals_list,
      Y_vals = Y_vals_list,
      W_names = W_names
    )
  }

  # Returning a list including the tabulated data and the names/values of the regressor, outcome, and controls
  return(out)

}
