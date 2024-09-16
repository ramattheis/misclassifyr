#' prep_misclassification_data
#'
#' This function tabulates data and generates metadata in a format to be used with the misclassifyr() function.
#'
#' @import dplyr
#'
#' @param data A data.frame containing the outcome variable,
#' @param outcome_1 A character string denoting the variable in the dataframe to be used as the first measure of an outcome, Y_1.
#' @param outcome_2 A character string denoting the variable in the dataframe to be used as the second measure of an outcome, Y_2.
#' @param regressor A character string denoting the variable in the dataframe to be used as the regressor, X.
#' @param outcome_1_bin A character string denoting the variable in the dataframe to be used to group / bin and average the first measure of an outcome, Y_1.
#' @param outcome_2_bin A character string denoting the variable in the dataframe to be used to group / bin and average the second measure of an outcome, Y_2.
#' @param regressor_bin A character string denoting the variable in the dataframe to be used to group / bin and average the regressor, X.
#' @param controls A character string or vector of character strings denoting the variable/variables to be used as non-parametric controls, W.
#' @param weights A character string denoting a variable containing individual level weights
#' @param X_names A vector of character strings denoting the values of the regressor in the desired order. If NA, as is default, names will be inferred from the data.
#' @param Y1_names A vector of character strings denoting the values of the outcome in the desired order. If NA, as is default, names will be inferred from the data.
#' @param Y2_names A vector of character strings denoting the values of the instrument in the desired order. If NA, as is default, names will be inferred from the data.
#' @param record_vals A logical value indicating whether to record the unique values of the outcomes and the regressor. If record_vals = F, you likely want to order the data by the regressor and outcomes before applying prep_misclassification_data.
#' @return A list of objects including tabulated data to be used in misclassifyr()
#' @export
prep_misclassification_data <- function(
    data,
    outcome_1,
    outcome_2,
    regressor,
    outcome_1_bin = NA,
    outcome_2_bin = NA,
    regressor_bin = NA,
    controls = NA,
    weights = NA,
    X_names = NA,
    Y1_names = NA,
    Y2_names = NA,
    record_vals = F) {

  #------------------------------------------------------------
  # Catching input errors
  #------------------------------------------------------------

  # Data should be a data.frame
  if(!is.data.frame(data)){stop("Data should be a data.frame.")}

  # If any _bin variables are present, they all should be
  if(any(c(!identical(outcome_1_bin,NA),
           !identical(outcome_2_bin,NA),
           !identical(regressor_bin,NA))) &
     !all(c(!identical(outcome_1_bin,NA),
           !identical(outcome_2_bin,NA),
           !identical(regressor_bin,NA)))){
    stop("If any `..._bin` variables provided, all should be provided.")
  }

  # Throwing an error if any of the entered column names are redundant
  keep = c(regressor, outcome_1, outcome_2)
  if(!identical(controls, NA)){ keep = c(keep, controls) }
  if(!identical(weights, NA)){ keep = c(keep, weights)}
  if(!identical(regressor_bin, NA)){ keep = c(keep, regressor_bin, outcome_1_bin, outcome_2_bin)}
  if(any(duplicated(keep))){stop("Column name arguments should not be repeated.")}

  # Variable names should be present in 'data'
  if(!(outcome_1 %in% colnames(data))){stop("`outcome_1` should be a column in data.")}
  if(!(outcome_2 %in% colnames(data))){stop("`outcome_2` should be a column in data.")}
  if(!(regressor %in% colnames(data))){stop("`regressor` should be a column in data.")}
  if(!identical(controls,NA)){if(!all(controls %in% colnames(data))){stop("`controls` should be a column / columns in data.")}}
  if(!identical(weights,NA)){if(!(weights %in% colnames(data))){stop("`weights` should be a column in data.")}}
  if(!identical(regressor_bin,NA)){if(!(regressor_bin %in% colnames(data))){stop("`regressor_bin` should be a column in data.")}}
  if(!identical(outcome_1_bin,NA)){if(!(outcome_1_bin %in% colnames(data))){stop("`outcome_1_bin` should be a column in data.")}}
  if(!identical(outcome_2_bin,NA)){if(!(outcome_2_bin %in% colnames(data))){stop("`outcome_2_bin` should be a column in data.")}}

  # Dropping data to the required columns
  data = data[,keep]

  # Passing a warning if any of the reserved column names (generated in this function) are in use
  if(("weight" %in% colnames(data))){warning("A column named 'weight' is present in `data` and will be overwritten")}
  if(("W" %in% colnames(data))){warning("A column named 'W' is present in `data` and will be overwritten")}
  if(("X" %in% colnames(data))){warning("A column named 'X' is present in `data` and will be overwritten")}
  if(("Y1" %in% colnames(data))){warning("A column named 'Y1' is present in `data` and will be overwritten")}
  if(("Y2" %in% colnames(data))){warning("A column named 'Y2' is present in `data` and will be overwritten")}
  if(("X_bin" %in% colnames(data))){warning("A column named 'X_bin' is present in `data` and will be overwritten")}
  if(("Y1_bin" %in% colnames(data))){warning("A column named 'Y1_bin' is present in `data` and will be overwritten")}
  if(("Y2_bin" %in% colnames(data))){warning("A column named 'Y2_bin' is present in `data` and will be overwritten")}

  # record_vals should be logical
  if(!is.logical(record_vals)){stop("`record_vals` should be logical.")}

  # Either record_vals = T OR X_names and Y_names should be provided
  if((record_vals)){
    if(!identical(X_names, NA) & !identical(Y1_names, NA) & !identical(Y2_names, NA)  ){
      stop("If `record_vals`==T, `X_names`, `Y1_names`, and `Y_2names` should not be provided.")
      }
  } else {
    if(identical(X_names, NA) | identical(Y1_names, NA) | identical(Y2_names,NA)){
      stop("If `record_vals`==F, `X_names`, `Y1_names`, and `Y2_names` should be provided.")
    }
  }

  # If record_vals is true, the outcome and regressor variables should be numeric in data
  if(record_vals){
    if(!is.numeric(unlist(data[,outcome_1]))){stop("If `record_vals` == T, `ourcome_1` in `data` should be numeric.")}
    if(!is.numeric(unlist(data[,outcome_2]))){stop("If `record_vals` == T, `ourcome_2` in `data` should be numeric.")}
    if(!is.numeric(unlist(data[,regressor]))){stop("If `record_vals` == T, `regressor` in `data` should be numeric.")}
  }

  # If provided, weights should always be numeric
  if(!identical(weights,NA)){if(any(is.na(as.numeric(data[,weights])))){stop("If provided, `weights` should be a numeric vector in `data`.")} }

  # Checking whether the sum of weights is equal to the number of rows in the data, sending a warning if not
  if(!identical(weights,NA)){if(abs(sum(data[,weights]) - nrow(data)) > 1){warning("The sum of `weights` is not equal to the number of rows in `data`. Are you sure weights are normalized correctly?") } }

  #------------------------------------------------------------
  # Defining weight, renaming columns, and collapsing controls to a single column
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

  # Defining factor variables to encode the desired ordering
  if(!identical(X_names, NA) & !identical(Y1_names, NA) & !identical(Y2_names, NA)){
    data[,regressor] = factor(data[,regressor], levels = X_names)
    data[,outcome_1] = factor(data[,outcome_1], levels = Y1_names)
    data[,outcome_2] = factor(data[,outcome_2], levels = Y2_names)
  }

  # Renaming the outcomes and the regressor
  data$X = data[,regressor]
  data$Y1 = data[,outcome_1]
  data$Y2 = data[,outcome_2]

  # recording the outcome and regressor groups if provided
  if(!identical(outcome_1_bin,NA) &
     !identical(outcome_2_bin,NA) &
     !identical(regressor_bin,NA)){
    data$X_bin = data[,regressor_bin]
    data$Y1_bin = data[,outcome_1_bin]
    data$Y2_bin = data[,outcome_2_bin]
  }

  # Dropping to necessary columns
  keep = c("X","Y1","Y2","weight")
  if(!identical(controls,NA)){keep = c(keep, "W")}
  if(!identical(regressor_bin,NA)){keep = c(keep, "X_bin", "Y1_bin", "Y2_bin")}
  data = data[,c(keep)]

  #------------------------------------------------------------
  # Tabulating
  #------------------------------------------------------------

  # First, ordering data by X, Y1, Y2
  data = data[order(data$X,data$Y1,data$Y2),]

  # Defining a function to tabulate one cell of data
  tabulate = function(cell){

    # Averaging X, Y1, Y2 within cell if _bins are provided
    if(!identical(regressor_bin, NA) & (record_vals)){
      # averaging within bins for
      cell = cell |>
        dplyr::group_by(X_bin) |>
        dplyr::mutate(X = mean(X)) |>
        dplyr::group_by(Y1_bin) |>
        dplyr::mutate(Y1 = mean(Y1)) |>
        dplyr::group_by(Y2_bin) |>
        dplyr::mutate(Y2 = mean(Y2)) |>
        as.data.frame()

      # Dropping binning variables
      cell$X_bin = NULL
      cell$Y1_bin = NULL
      cell$Y2_bin = NULL
    }

    # Throwing an error if the number of unique values in Y1 and Y2 is different
    if(length(unique(cell$Y1)) != length(unique(cell$Y2))){stop("`Y1` and `Y2` or `Y1_bin` and `Y2_bin` must be the same dimension.")}

    # Recording the names of the outcome and the regressor
    if(identical(X_names,NA)){
      X_names = unique(cell$X) |> as.character()
      Y1_names = unique(cell$Y1) |> as.character()
      Y2_names = unique(cell$Y2) |> as.character()
    }

    # Recording the values of the outcome and the regressor
    if(record_vals){
      Y_vals = cell$Y1 |> unique() |> sort(decreasing = F)
      X_vals = cell$X |> unique() |> sort(decreasing = F)
    } else {
      Y_vals = NA
      X_vals = NA
    }

    # Ensuring that the tabulation is balanced across X, Y1, Y2
    empty_cell = expand.grid(X_names, Y1_names, Y2_names)
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
     J = length(Y1_names),
     K = length(X_names),
     X_names = X_names,
     Y1_names = Y1_names,
     Y2_names = Y2_names,
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
    Y1_names_list = lapply(out,"[[",5)
    Y2_names_list = lapply(out,"[[",6)
    X_vals_list = lapply(out,"[[",7)
    Y_vals_list = lapply(out,"[[",8)

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
      Y1_names = Y1_names_list,
      Y2_names = Y2_names_list,
      X_vals = X_vals_list,
      Y_vals = Y_vals_list,
      W_names = W_names
    )
  }

  # Returning a list including the tabulated data and the names/values of the regressor, outcome, and controls
  return(out)

}
