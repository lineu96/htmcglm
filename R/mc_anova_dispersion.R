#' @name mc_anova_dispersion
#'
#' @author Lineu Alberto Cavazani de Freitas,
#' \email{lineuacf@@gmail.com}
#'
#' @export
#'
#' @title ANOVA tables for dispersion components.
#'
#' @description Performs Wald tests to generate analysis-of-variance
#' tables of the significance for the dispersion components by response
#' variables for model objects produced by mcglm.
#'
#' @param object An object of \code{mcglm} class.
#'
#' @param p_var A list of indices that indicate how the dispersion
#' parameters are related. Parameters with the same index are tested
#' together.
#'
#' @param names Names to be shown in the table.
#' 
#' @param verbose a logical if TRUE print some information about the 
#' tests performed. Default verbose = TRUE.
#'
#' @return Type III ANOVA table for dispersion components of mcglm
#' objects.
#'
#' @seealso \code{mc_anova_I}, \code{mc_anova_II} and
#' \code{mc_anova_III}.
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
#' mc_anova_dispersion(fit_joint,
#' p_var = list(c(0,1), c(0,1), c(0,1)),
#' names = list(c('tau10', 'tau11'),
#'              c('tau20', 'tau21'),
#'              c('tau30', 'tau31')))
#'

mc_anova_dispersion <- function(object, p_var, names, verbose = TRUE){

  # Vetor tau e indice de resposta
  tau <- coef(object, type = "tau")[,c(1,2, 4)]

  #----------------------------------------------------------------

  # Número de taus por resposta
  n_tau <- as.vector(table(tau$Response))

  #----------------------------------------------------------------

  # Número de respostas
  n_resp <- length(n_tau)

  #----------------------------------------------------------------

  # Lista vcov por resposta desconsiderando parametros de
  # regressao e potencia

  vcov_taus <- list()

  padrao <- vector()

  for (j in 1:n_resp) {
    for (i in 1:length(row.names(vcov(object)))) {
      padrao[i] <- sjmisc::str_contains(rownames(vcov(object))[i],
                                        pattern = paste0('tau',j))
    }

    id <- NULL

    names2 <- data.frame(row_names = row.names(vcov(object)),
                         id = padrao)

    names2 <- as.vector(subset(names2, id == TRUE)$row_names)

    vcov_taus[[j]] <- vcov(object)[names2, names2]
  }

  #----------------------------------------------------------------

  # Matriz L para todos os parâmetros (Hypothesis matrix), por resposta
  L_all <- list()

  for (i in 1:n_resp) {
    L_all[[i]] <- diag(length(p_var[[i]]))
  }

  #----------------------------------------------------------------

  # Matriz L por variável (Hypothesis matrix), por resposta

  L_par <- list()

  for (i in 1:n_resp) {
    L_par[[i]] <- by(data = L_all[[i]],
                     INDICES = p_var[[i]],
                     FUN = as.matrix)
  }

  #----------------------------------------------------------------

  ## Tabela

  tabela <- list()

  for (j in 1:n_resp) {

    W <- vector() # Vetor para a estatística de teste
    gl <- vector() # Vetor para graus de liberdade
    p_val <- vector() # Vetor para p-valor


    for (i in 1:dim(L_par[[j]])) {
      W[i] <- as.numeric((t(L_par[[j]][[i]] %*%
                              subset(tau,
                                     tau$Response == j)$Estimates)) %*%
                           (solve(L_par[[j]][[i]]%*%
                                    vcov_taus[[j]]%*%
                                    t(L_par[[j]][[i]]))) %*%
                           (L_par[[j]][[i]] %*%
                              subset(tau, tau$Response == j)$Estimates))
      gl[i] <- ifelse(is.null(nrow(L_par[[j]][[i]])) == TRUE,
                      1,nrow(L_par[[j]][[i]]))
      p_val[i] <- pchisq(W[i], df = gl[i], lower.tail = FALSE)

    }
    tabela[[j]] <-
      data.frame(Dispersion = names[[j]],
                 Df = gl,
                 Chi = round(W, 4),
                 'Pr(>Chi)' = round(p_val, 4),
                 check.names = F)
  }

  #----------------------------------------------------------------
  
  if (verbose == TRUE) {
    cat(
      "ANOVA type III using Wald statistic for dispersion parameters\n\n")
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
