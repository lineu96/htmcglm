#' @name mc_multcomp
#'
#' @author Lineu Alberto Cavazani de Freitas,
#' \email{lineuacf@@gmail.com}
#'
#' @export
#'
#' @title Multiple comparisons test for each response.
#'
#' @description Performs a multiple comparisons test to compare
#' diferences between treatment levels for each response for model
#' objects produced by mcglm.
#'
#' @param object An object of \code{mcglm} class.
#'
#' @param effect A list of vector of variables. For each configuration
#' of these the estimate will be calculated.
#'
#' @param data Data frame with the dataset used in the model.
#' 
#' @param verbose a logical if TRUE print some information about the 
#' tests performed. Default verbose = TRUE.
#'
#' @return Table of multiple comparisons.
#'
#' @seealso \code{mc_mult_multcomp}.
#'
#' @examples
#'
#' library(mcglm)
#' library(Matrix)
#' library(htmcglm)
#'
#' data("soya", package = "mcglm")
#'
#' form.grain <- grain ~ water * pot
#' form.seed <- seeds ~ water * pot
#'
#' soya$viablepeasP <- soya$viablepeas / soya$totalpeas
#' form.peas <- viablepeasP ~ water * pot
#'
#' Z0 <- mc_id(soya)
#' Z1 <- mc_mixed(~0 + factor(block), data = soya)
#'
#' fit_joint <- mcglm(linear_pred = c(form.grain,
#'                                    form.seed,
#'                                    form.peas),
#'                    matrix_pred = list(c(Z0, Z1),
#'                                       c(Z0, Z1),
#'                                       c(Z0, Z1)),
#'                    link = c("identity",
#'                             "log",
#'                             "logit"),
#'                    variance = c("constant",
#'                                 "tweedie",
#'                                 "binomialP"),
#'                    Ntrial = list(NULL,
#'                                  NULL,
#'                                  soya$totalpeas),
#'                    power_fixed = c(TRUE,TRUE,TRUE),
#'                    data = soya)
#'
#' mc_multcomp(object = fit_joint,
#'             effect = list(c('water'),
#'                           c('water'),
#'                           c('water')),
#'             data = soya)
#'
#' mc_multcomp(object = fit_joint,
#'             effect = list(c('pot'),
#'                           c('pot'),
#'                           c('pot')),
#'             data = soya)
#'
#'
#' mc_multcomp(object = fit_joint,
#'             effect = list(c('water', 'pot'),
#'                           c('water', 'pot'),
#'                           c('water', 'pot')),
#'             data = soya)
#'

mc_multcomp <- function(object, effect, data, verbose = TRUE){

  # Obter a matriz de combinações lineares dos parâmetros dos
  # modelos que resultam nas médias ajustadas
  # (geralmente denotada por L)
  #-------------------------------------------------------------------

  m_glm <- list()

  for (i in 1:length(object$linear_pred)) {
    m_glm[[i]] <- glm(formula = object$linear_pred[[i]], data = data)
  }

  mm <- list()

  for (i in 1:length(m_glm)) {
    mm[[i]] <- doBy::LE_matrix(m_glm[[i]], effect = effect[[i]])
  }

  #-------------------------------------------------------------------

  # Para testar os contrastes de uma média ajustada contra a outra
  # deve-se subtrair as linhas da primeira matriz duas a duas
  # (geralmente denotada por K)

  K1 <- list()

  for (i in 1:length(mm)) {
    K1[[i]] <- apc(mm[[i]])
  }

  K2 <- list()

  for (i in 1:length(mm)) {
    K2[[i]] <- by(K1[[i]], INDICES = row.names(K1[[i]]),
                  FUN = as.matrix)
  }

  #-------------------------------------------------------------------

  # Vetor beta chapeu e indice de resposta
  beta <- coef(object, type = "beta")[,c(1, 4)]

  #----------------------------------------------------------------

  # Número de betas por resposta
  n_beta <- as.vector(table(beta$Response))

  #----------------------------------------------------------------

  # Número de respostas
  n_resp <- length(n_beta)

  #----------------------------------------------------------------

  # Lista vcov por resposta desconsiderando parametros de dispersao e
  # potencia

  vcov_betas <- list()

  if (n_resp == 1) {
    vcov_betas[[1]] <- vcov(object)[1:n_beta[1], 1:n_beta[1]]
  } else {
    vcov_betas[[1]] <- vcov(object)[1:n_beta[1], 1:n_beta[1]]
    for (i in 2:n_resp) {
      vcov_betas[[i]] <-
        vcov(object)[(cumsum(n_beta)[i-1]+1):(cumsum(n_beta)[i]),
                     (cumsum(n_beta)[i-1]+1):(cumsum(n_beta)[i])]

    }
  }

  #----------------------------------------------------------------

  ## Tabela

  tabela <- list()

  for (j in 1:n_resp) {

    W <- vector() # Vetor para a estatística de teste
    gl <- vector() # Vetor para graus de liberdade
    p_val <- vector() # Vetor para p-valor


    for (i in 1:dim(K2[[j]])) {
      W[i] <- as.numeric((t(K2[[j]][[i]] %*%
                              subset(beta,
                                     beta$Response == j)$Estimates)) %*%
                           (solve(K2[[j]][[i]]%*%
                                    vcov_betas[[j]]%*%
                                    t(K2[[j]][[i]])))%*%
                           (K2[[j]][[i]] %*%
                              subset(beta,
                                     beta$Response == j)$Estimates))
      gl[i] <- nrow(K2[[j]][[i]])
      p_val[i] <- pchisq(W[i], df = gl[i], lower.tail = FALSE)

    }
    tabela[[j]] <-
      data.frame(Contrast = names(K2[[j]]),
                 Df = gl,
                 Chi = round(W, 4),
                 'Pr(>Chi)' = round(p.adjust(p_val,
                                             method = 'bonferroni'), 4),
                 check.names = F)
  }

  #----------------------------------------------------------------

  if (verbose == TRUE) {
    cat("Multiple comparisons test for each outcome using Wald statistic\n\n")
    for (i in 1:n_resp) {
      cat("Call: ")
      print(object$linear_pred[[i]])
      cat("\n")
      print(tabela[[i]])
      cat("\n")
    }
    
    return(invisible(tabela))
  } else {
    return(tabela)
  }
  
}
