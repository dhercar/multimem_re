#' Fit a modular glmmTMB model with two-column multi-membership random effects
#'
#' @param col_i Name of re (col 1) (character)
#' @param col_j Name of re (col 2) (character)
#' @param data Data frame containing columns col_i, col_j, response, and predictors
#' @param control_nlminb A list of control parameters passed directly to optimised (nlminb) (e.g., list(eval.max = 300, iter.max = 300))
#' @param formula A two-sided model formula (without the fake term), passed to glmmTMB
#' @param parallel An integer specifying the number of cores to use in parallel model estimation (passed to TMB::openmp)
#' @param ... Additional arguments passed to glmmTMB() (e.g., family = gaussian())
#' @return A fitted glmmTMB model with multimembership RE
#' 
fit_mm_glmmTMB <- function(formula, 
                           col_i, 
                           col_j, 
                           data,
                           control_nlminb = list(), 
                           parallel = 1, # 
                           ...) {
  
  TMB::openmp(parallel)
  
  cat(format(Sys.time(), "%X"), ': getting data ready\n')
  
  stopifnot(is.character(col_i), is.character(col_j),
            all(c(col_i, col_j, 'date_i', 'date_j') %in% names(data)))
  
  # membership identifiers
  s_i <- data[[col_i]]
  s_j <- data[[col_j]]
  unique_s <- unique(c(s_i, s_j))
  
  # parse W (each row has two 1's)
  W <- Matrix(0, nrow = nrow(data), ncol = length(unique_s), sparse = TRUE)
  W[cbind(seq_len(nrow(data)), match(s_i, unique_s))] <- 1
  W[cbind(seq_len(nrow(data)), match(s_j, unique_s))] <- 1
  stopifnot(all(rowSums(W) == 2))
  colnames(W) <- unique_s
  
  # dummy factor 'fake' as last random term
  data$fake <- rep(unique_s, length.out = nrow(data))
  
  # modify formula
  formula_mod <- update(formula,  ~ . + (1|fake))
  
  # build model spec without fitting
  cat(format(Sys.time(), "%X"), ': building model components\n')
  m0 <- glmmTMB(formula_mod, data = data, doFit = FALSE, ..., control =  glmmTMBControl(parallel = parallel))
  mt0 <- fitTMB(m0, doOptim = FALSE)

  # Replace col in Z corresponding to mm re with W
  Z_cols <- ncol(mt0$env$data$Z)
  mt0$env$data$Z <- as(Matrix(cbind(mt0$env$data$Z[,-c((Z_cols - ncol(W) + 1):Z_cols)], W)), 'TsparseMatrix')
  
  cat(format(Sys.time(), "%X"), ': starting optimization\n')
  
  # finish fit
  fit <- with(mt0, nlminb(par, fn, gr, control = control_nlminb))
  
  cat(format(Sys.time(), "%X"), ': finalizing model\n')
  m1  <- finalizeTMB(m0, mt0, fit)
  
  cat(format(Sys.time(), "%X"), ': done\n')
  return(m1)
}
