#' Synthetic linked father-son census records from Ancien Regime France
#'
#' A synthetic dataset created to demonstrate the use of the `misclassifyr`
#' package. It contains father-son pairs in the Third Estate in Ancien Regime
#' France, linked across three censuses (that never actually happened) in
#' 1750, 1770, and 1780. Because the son is observed twice --- in 1770 and
#' again in 1780 --- the data supply the two noisy measures of a single latent
#' outcome that identification of a misclassification model requires.
#'
#' Occupations are ordered from lowest to highest income score:
#' `"Vagabond"`, `"Metayer"`, `"Journalier"`, `"Petit Metiers"`,
#' `"Petite Bourgeoisie"`, `"Haute Bourgeoisie"`.
#'
#' @format A data frame with 100,000 rows and 9 variables:
#' \describe{
#'   \item{birthplace}{The region of birth for the son.}
#'   \item{birthyear}{The year of birth for the son.}
#'   \item{linked_weight}{Weights intended to correct for selection into the linked sample based on birthplace and birthyear.}
#'   \item{father_occupation_1750}{The occupation of the father observed in the (fictional) 1750 census.}
#'   \item{father_income_1750}{The income 'score' of the father's occupation in 1750, where the score varies across birthyear and birthplace.}
#'   \item{son_occupation_1770}{The occupation of the son observed in the (fictional) 1770 census.}
#'   \item{son_income_1770}{The income 'score' of the son's occupation in 1770, where the score varies across birthyear and birthplace.}
#'   \item{son_occupation_1780}{The occupation of the son observed in the (fictional) 1780 census.}
#'   \item{son_income_1780}{The income 'score' of the son's occupation in 1780, where the score varies across birthyear and birthplace.}
#' }
#' @source Synthetic data generated for the package.
#' @examples
#' data(ancienregime)
#' str(ancienregime)
#'
#' # The two measures of the son's occupation disagree for a sizeable share of
#' # pairs: that disagreement is what identifies the misclassification model.
#' mean(ancienregime$son_occupation_1770 != ancienregime$son_occupation_1780)
#'
#' # Naive intergenerational elasticity, uncorrected for measurement error
#' coef(lm(son_income_1780 ~ father_income_1750, data = ancienregime,
#'         weights = ancienregime$linked_weight))
"ancienregime"
