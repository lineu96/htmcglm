#' @name mc_mult_multcomp
#'
#' @author Lineu Alberto Cavazani de Freitas,
#' \email{lineuacf@@gmail.com}
#'
#' @export
#'
#' @title Multivariate multiple comparisons test.
#'
#' @description Performs a multiple comparisons test to compare
#' diferences between treatment levels for all responses of model
#' objects produced by mcglm.
#'
#' @param object An object of \code{mcglm} class.
#'
#' @param effect A vector of variables. For each configuration of these
#' the estimate will be calculated.
#'
#' @param data Data frame with the dataset used in the model.
#' 
#' @param verbose a logical if TRUE print some information about the 
#' tests performed. Default verbose = TRUE.
#'
#' @return Table of multiple comparisons.
#'
#' @seealso \code{mc_multcomp}.
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
#' mc_mult_multcomp(object = fit_joint,
#'                  effect = c('water'),
#'                  data = soya)
#'
#' mc_mult_multcomp(object = fit_joint,
#'                  effect = c('pot'),
#'                  data = soya)
#'
#' mc_mult_multcomp(object = fit_joint,
#'                  effect = c('water', 'pot'),
#'                  data = soya)
#'

mc_mult_multcomp <- function(object, effect, data, verbose = TRUE){

  # Obter a matriz de combinações lineares dos parâmetros dos
  # modelos que resultam nas médias ajustadas (geralmente denotada
  # por L)

  m_glm <- glm(formula = object$linear_pred[[1]], data = data)
  mm <- doBy::LE_matrix(m_glm, effect = effect)

  #-------------------------------------------------------------------

  # Para testar os contrastes de uma média ajustada contra a outra
  # deve-se subtrair as linhas da primeira matriz duas a duas
  # (geralmente denotada por K)

  K1 <- apc(mm)
  K2 <- by(K1, INDICES = row.names(K1), FUN = as.matrix)

  #-------------------------------------------------------------------

  # Aplicando K no teste Wald

  #----------------------------------------------------------------

  # Vetor beta chapeu
  beta <- coef(object, type = "beta")[,c(1, 4)]

  #----------------------------------------------------------------

  # Número de betas
  n_beta <- sum(as.vector(table(beta$Response)))

  #----------------------------------------------------------------

  # Número de respostas
  n_resp <- length(as.vector(table(beta$Response)))

  #----------------------------------------------------------------

  # ERROS E AVISOS

  preds <- c()

  for (i in 1:n_resp) {
    preds[i] <- sub(".*~", "",gsub(" ", "",
                                   as.character(object$linear_pred)[i]))
  }

  if(length(unique(preds)) != 1) stop("The predictors must be the same for all outcomes.")

  if(n_resp == 1) warning("You are applying a multivariate function to a univariate problem.")

  #----------------------------------------------------------------

  # vcov desconsiderando parametros de dispersao e potencia
  vcov_betas <- vcov(object)[1:n_beta, 1:n_beta]

  #----------------------------------------------------------------

  # Matriz G
  G <- diag(n_resp)

  #----------------------------------------------------------------

  # Matriz L

  L_par <- list()

  for (i in 1:length(K2)) {
    L_par[[i]] <- kronecker(G, K2[[i]])
  }


  #----------------------------------------------------------------

  ## Tabela

  W <- vector() # Vetor para a estatística de teste
  gl <- vector() # Vetor para graus de liberdade
  p_val <- vector() # Vetor para p-valor

  ### Estatística de teste:
  #### t(L*beta) x (L*vcov*t(L))^-1 x (L*beta) ~ Qui-quadrado(numero de parametros testados)

  for (i in 1:length(L_par)) {

    W[i] <- as.numeric((t(L_par[[i]]%*%beta$Estimates)) %*%
                         (solve(L_par[[i]]%*%vcov_betas%*%
                                  t(L_par[[i]]))) %*%
                         (L_par[[i]]%*%beta$Estimates))
    gl[i] <- nrow(L_par[[i]])
    p_val[i] <- pchisq(W[i], df = gl[i], lower.tail = FALSE)
  }

  tabela <- data.frame(Contrast = names(K2),
                       Df = gl,
                       Chi = round(W, 4),
                       'Pr(>Chi)' = round(
                         p.adjust(p_val, method = 'bonferroni'), 4),
                       check.names = F)

  #----------------------------------------------------------------
  
  if (verbose == TRUE) {
    cat("Multivariate multiple comparisons test using Wald statistic\n\n")
    cat("Call: ")
    cat(paste0('~ ', preds[[1]]))
    cat("\n")
    print(tabela)
    
    return(invisible(tabela))
  } else {
    return(tabela)
  }

}
