#' @name mc_manova_dispersion
#'
#' @author Lineu Alberto Cavazani de Freitas,
#' \email{lineuacf@@gmail.com}
#'
#' @export
#'
#' @title MANOVA tables for dispersion components.
#'
#' @description Performs Wald tests to generate multivariate analysis-of
#' -variance tables of the significance for the dispersion components
#' for model objects produced by mcglm.
#'
#' @param object An object of \code{mcglm} class.
#'
#' @param p_var A vector of indices that indicate how the dispersion
#' parameters are related. Parameters with the same index are tested
#' together.
#'
#' @param names Names to be shown in the table.
#' 
#' @param verbose a logical if TRUE print some information about the 
#' tests performed. Default verbose = TRUE.
#'
#' @return Type III MANOVA table for dispersion components of mcglm
#' objects.
#'
#' @seealso \code{mc_manova_I}, \code{mc_manova_II} and
#' \code{mc_manova_III}.
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
#' mc_manova_dispersion(fit_joint,
#'                p_var = c(0,1),
#'                names = c('tau11', 'tau21'))
#'

mc_manova_dispersion <- function(object, p_var, names, verbose = TRUE){

  # Vetor tau
  tau <- coef(object, type = "tau")[,c(1,2,4)]

  #----------------------------------------------------------------

  # Número de taus
  n_tau <- sum(as.vector(table(tau$Response)))

  #----------------------------------------------------------------

  # Número de respostas
  n_resp <- length(as.vector(table(tau$Response)))

  #----------------------------------------------------------------

  # ERROS E AVISOS

  preds <- c()

  for (i in 1:n_resp) {
    preds[i] <- sub(".*~", "",gsub(" ", "",
                                   as.character(object$linear_pred)[i]))
  }

  if(length(unique(preds)) != 1) stop("For MANOVA functions, the predictors must be the same for all outcomes.")

  if(n_resp == 1) warning("You are applying a MANOVA function to a univariate problem.")

  #----------------------------------------------------------------

  # Nomes dos parâmetros
  tau_names <- as.vector(unique(tau$Parameters))

  #----------------------------------------------------------------

  # vcov desconsiderando parametros de regressao e potencia
  vcov_taus <- vcov(object)[tau_names, tau_names]

  #----------------------------------------------------------------

  # Índice que associa tau a matriz Z
  p_var <- p_var

  #----------------------------------------------------------------

  # Matriz F para todos os parâmetros (Hypothesis matrix)
  F_all <- diag(length(p_var))

  #----------------------------------------------------------------

  # Matriz F por variável (Hypothesis matrix)
  F_par <- by(data = F_all,
              INDICES = p_var,
              FUN = as.matrix)

  #----------------------------------------------------------------

  # Matriz G
  G <- diag(n_resp)

  #----------------------------------------------------------------

  # Matriz L

  L_par <- list()

  for (i in 1:length(F_par)) {
    L_par[[i]] <- kronecker(G, F_par[[i]])
  }


  #----------------------------------------------------------------

  ## Tabela

  W <- vector() # Vetor para a estatística de teste
  gl <- vector() # Vetor para graus de liberdade
  p_val <- vector() # Vetor para p-valor

  ### Estatística de teste:
  #### t(L*beta) x (L*vcov*t(L))^-1 x (L*beta) ~ Qui-quadrado(numero de
  #### parametros testados)

  for (i in 1:length(L_par)) {

    W[i] <- as.numeric((t(L_par[[i]]%*%tau$Estimates)) %*%
                         (solve(L_par[[i]]%*%
                                  vcov_taus%*%
                                  t(L_par[[i]]))) %*%
                         (L_par[[i]]%*%tau$Estimates))
    gl[i] <- nrow(L_par[[i]])
    p_val[i] <- pchisq(W[i], df = gl[i], lower.tail = FALSE)
  }

  tabela <- data.frame(Covariate = names,
                       Df = gl,
                       Chi = round(W, 4),
                       'Pr(>Chi)' = round(p_val, 4),
                       check.names = F)

  #----------------------------------------------------------------
  
  if (verbose == TRUE) {
    cat("MANOVA type III using Wald statistic for dispersion parameters\n\n")
    cat("Call: ")
    cat(paste0('~ ', preds[[1]]))
    cat("\n")
    print(tabela)
    
    return(invisible(tabela))
  } else {
    return(tabela)
  }
  
}
