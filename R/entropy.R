#' Entropy - Shannon Information
#'
#' @param test A Data frame, must contain numerical columns
#' @description Calculates the Entropy value for each row across numerical columns.
#'
#'
#' @return
#' Returns an extra column to the dataframe with a value related to the Entropy value. :)
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr
#' @references
#' https://en.wikipedia.org/wiki/Entropy_(information_theory)
#' @examples
#' data(rnapol_score)
#' test =  entropy(rnapol_score)
#'
#'
#'
entropy = function(test){
  test = test |>  mutate(suma = rowSums(across(where(is.numeric))))
  test = test |> mutate(suma = ifelse(suma > 0 , suma, 1))
  test <- test %>%
    mutate(across(where(is.numeric), ~ . / suma))

  test = test |> select(-suma)

  test <- test %>%
    rowwise() %>%
    mutate(Entropy = -sum(c_across(where(is.numeric)) * ifelse(c_across(where(is.numeric)) == 0,0, log2(c_across(where(is.numeric))))))

  return(test)
}
