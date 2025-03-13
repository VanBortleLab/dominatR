
#' Categorical Entropy
#'
#' @param test A Data frame, must contain numerical columns and Entropy should be calculated already
#'
#' @description Calculates the Categorical Entropy value related to each column. Entropy calculation should be performed before using this function
#'
#'
#'
#' @return
#' Returns the data frame with modified values for numerical columns.
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr
#' @examples
#' data(rnapol_score)
#' test =  entropy(rnapol_score)
#' test = Qentropy(test)
#'
#' @seealso [entropy()]
Qentropy = function(test){
  test <- test %>%
    rowwise() %>%
    mutate(across(where(is.numeric), ~ ifelse(. == 0, Inf, Entropy - log2(.))))

  test = test |> select(-Entropy)

  return(test)
}
