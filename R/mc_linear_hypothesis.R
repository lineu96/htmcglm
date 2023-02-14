#' @name mc_linear_hypothesis
#'
#' @author Lineu Alberto Cavazani de Freitas,
#' \email{lineuacf@@gmail.com}
#'
#' @export
#'
#' @import mcglm
#' @import Matrix
#'
#' @importFrom stats as.formula binomial coef dist fitted glm make.link
#' model.frame model.matrix na.exclude pchisq qchisq qnorm pnorm quasi
#' residuals vcov model.response terms p.adjust
#'
#' @title Test Linear Hypothesis
#'
#' @description Performs Wald tests for testing a linear hypothesis for
#' model objects produced by mcglm.
#'
#' @param object An object of \code{mcglm} class.
#'
#' @param hypothesis A vector of strings with the hypotheses to be
#' tested.
#' 
#' @param verbose a logical if TRUE print some information about the 
#' tests performed. Default verbose = TRUE.
#'
#' @return Table result of the hypothesis test specified.
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
#' mc_linear_hypothesis(object =  fit_joint,
#'                      hypothesis = c('beta11 = 0'))
#'
#' mc_linear_hypothesis(object =  fit_joint,
#'                      hypothesis = c('beta11 = 0',
#'                                     'beta12 = 0'))
#'
#' mc_linear_hypothesis(object =  fit_joint,
#'                      hypothesis = c('beta11 = 0',
#'                                     'beta12 = 0',
#'                                     'beta21 = 0',
#'                                     'beta22 = 0',
#'                                     'beta31 = 0',
#'                                     'beta32 = 0'))
#'
#' mc_linear_hypothesis(object =  fit_joint,
#'                      hypothesis = c('beta11 = beta21'))
#'
#' mc_linear_hypothesis(object =  fit_joint,
#'                      hypothesis = c('tau11 = 0'))
#'
#' mc_linear_hypothesis(object =  fit_joint,
#'                      hypothesis = c('tau11 = 0',
#'                                     'tau21 = 0',
#'                                     'tau31 = 0'))
#'
#' mc_linear_hypothesis(object =  fit_joint,
#'                      hypothesis = c('tau12 = tau22'))

mc_linear_hypothesis <- function(object, hypothesis, verbose = TRUE){

  # Vetor beta chapeu
  coefs <- coef(object, type = c("beta", "tau", "power"))

  #----------------------------------------------------------------

  # Número de parametros
  n_coefs <- sum(as.vector(table(coefs$Response)))

  #----------------------------------------------------------------

  # Número de respostas
  n_resp <- length(as.vector(table(coefs$Response)))

  #----------------------------------------------------------------

  # vcov
  vcov_coefs <- vcov(object)[as.vector(coefs$Parameters),
                             as.vector(coefs$Parameters)]

  #----------------------------------------------------------------

  # Quebrando as strings recebidas como argumento (hipóteses)
  hypothesis2 <- stringr::str_split(hypothesis,
                                    pattern = c('='),
                                    simplify = T)

  #----------------------------------------------------------------

  # Gerando um data frame
  hypothesis3 <- as.data.frame(hypothesis2)
  names(hypothesis3) <- c('parameters', 'null_hyp')

  #----------------------------------------------------------------

  # Retirando espaços
  hypothesis3$parameters <- stringr::str_replace(
    hypothesis3$parameters, " ",  "")
  hypothesis3$null_hyp <- stringr::str_replace(
    hypothesis3$null_hyp, " ",  "")

  #----------------------------------------------------------------

  if(sum(hypothesis3$parameters %in% coefs$Parameters) !=
     length(hypothesis3$parameters)) stop("You specified hypothesis about parameters that does not exist in the model.")

  #----------------------------------------------------------------

  # Gerando a matriz L
  L_user <- matrix(nrow = nrow(hypothesis3), ncol = n_coefs)

  colnames(L_user) <- as.vector(coefs$Parameters)

  for (i in 1:nrow(L_user)) {
    L_user[i,] <- ifelse(colnames(L_user) == hypothesis3$parameters[i],
                         1,0)
  }

  for (i in 1:nrow(L_user)) {
    L_user[i,] <- ifelse(colnames(L_user) == hypothesis3$null_hyp[i],
                         -1,L_user[i,])
  }

  #----------------------------------------------------------------

  # Substitui strings por 0 (hipóteses de igualdade)
  hypothesis3$null_hyp <- ifelse(stringr::str_detect(
    hypothesis3$null_hyp, c('beta')) == T,0,hypothesis3$null_hyp)
  hypothesis3$null_hyp <- ifelse(stringr::str_detect(
    hypothesis3$null_hyp, c('tau')) == T,0,hypothesis3$null_hyp)
  hypothesis3$null_hyp <- ifelse(stringr::str_detect(
    hypothesis3$null_hyp, c('power')) == T,0,hypothesis3$null_hyp)

  #----------------------------------------------------------------

  # Converte valores da hipótese nula para numericas
  hypothesis3$null_hyp <- as.numeric(hypothesis3$null_hyp)

  #----------------------------------------------------------------

  # Efetua o teste

  W <- vector() # Vetor para a estatística de teste
  gl <- vector() # Vetor para graus de liberdade
  p_val <- vector() # Vetor para p-valor

  W <- as.numeric((t((L_user%*%coefs$Estimates) - hypothesis3$null_hyp))
                  %*% (solve(L_user%*%vcov_coefs%*%t(L_user))) %*%
                    ((L_user%*%coefs$Estimates) - hypothesis3$null_hyp))
  gl <- nrow(L_user)
  p_val <- pchisq(W, df = gl, lower.tail = FALSE)

  tabela <- data.frame(Df = gl,
                       Chi = round(W, 4),
                       'Pr(>Chi)' = round(p_val, 4),
                       check.names = F)

  #----------------------------------------------------------------

  if (verbose == TRUE) {
    cat("Linear hypothesis test")
    cat("\n\n")
    cat("Hypothesis:")
    
    
    print_hyp <- as.data.frame(hypothesis)
    names(print_hyp) <- NULL
    print(print_hyp)
    
    cat("\n")
    
    cat("Results:\n")
    
    print(tabela)
    cat("\n")
    
    return(invisible(tabela))
  } else {
    return(tabela)
  }

}
