library(spectralGraphTopology)

#' @title Laplacian matrix of a balanced k-component graph with heavy-tailed data
#'
#' Computes the Laplacian matrix of a balanced graph on the basis of an observed data matrix,
#' where we assume the data to be Student-t distributed.
#'
#' @param X an n x p data matrix, where n is the number of observations and p is
#'        the number of nodes in the graph.
#' @param k the number of components of the graph.
#' @param heavy_type a string which selects the statistical distribution of the data    .
#'        Valid values are "gaussian" or "student".
#' @param nu the degrees of freedom of the Student-t distribution.
#'        Must be a real number greater than 2.
#' @param w0 initial vector of graph weights. Either a vector of length p(p-1)/2 or
#'        a string indicating the method to compute an initial value.
#' @param beta hyperparameter that controls the regularization to obtain a
#'        k-component graph
#' @param update_beta whether to update beta during the optimization.
#' @param alpha spare regularization term parameter
#' @param t upperbound scale parameter
#' @param d the nodes' degrees. Either a vector or a single value.
#' @param rho ADMM hyperparameter.
#' @param update_rho whether or not to update rho during the optimization.
#' @param maxiter maximum number of iterations.
#' @param reltol relative tolerance as a convergence criteria.
#' @param verbose whether or not to show a progress bar during the iterations.
#' @export
#' @import spectralGraphTopology
learn_kcomp_heavytail_graph_balanced <- function(X,
                                                 alpha = 0.1,
                                                 k = 1,
                                                 heavy_type = "gaussian",
                                                 nu = NULL,
                                                 w0 = "naive",
                                                 d = 1,
                                                 t = 1.2,
                                                 beta = 1e-8,
                                                 update_beta = TRUE,
                                                 early_stopping = FALSE,
                                                 rho = 1,
                                                 update_rho = FALSE,
                                                 maxiter = 10000,
                                                 reltol = 1e-5,
                                                 verbose = TRUE,
                                                 record_objective = FALSE) {
  
  X <- scale(as.matrix(X))
  # number of nodes
  p <- ncol(X)
  
  t <- t*d/sqrt(p-1)
  # number of observations
  n <- nrow(X)
  LstarSq <- vector(mode = "list", length = n)
  for (i in 1:n)
    LstarSq[[i]] <- Lstar(X[i, ] %*% t(X[i, ])) / n
  # w-initialization
  if (assertthat::is.string(w0)) {
    w <- spectralGraphTopology:::w_init(w0, MASS::ginv(cor(X)))
    A0 <- A(w)
    A0 <- A0 / rowSums(A0)
    w <- spectralGraphTopology:::Ainv(A0)
  }
  else {
    w <-w0
  }
  # Theta-initilization
  Lw <- L(w)
  Aw <- A(w)
  Theta <- Lw
  U <- eigen(Lw, symmetric = TRUE)$vectors[, (p - k + 1):p]
  Y <- matrix(0, p, p)
  y <- rep(0, p)
  # Z-initialization
  Lambda <- matrix(0, p, p)
  Z <- A0
  
  # ADMM constants
  mu <- 2
  tau <- 2
  # residual vectors
  primal_lap_residual <- c()
  primal_deg_residual <- c()
  dual_residual <- c()
  # augmented lagrangian vector
  lagrangian <- c()
  beta_seq <- c()
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  elapsed_time <- c()
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    
    for (j in 1:1){
      # update w
      LstarLw <- Lstar(Lw)
      DstarDw <- Dstar(diag(Lw))
      LstarSweighted <- rep(0, .5*p*(p-1))
      if (heavy_type == "student") {
        for (q in 1:n)
          LstarSweighted <- LstarSweighted + LstarSq[[q]] * compute_student_weights(w, LstarSq[[q]], p, nu)
      } else if (heavy_type == "gaussian") {
        for (q in 1:n)
          LstarSweighted <- LstarSweighted + LstarSq[[q]]
      }
      grad <- LstarSweighted + Lstar(beta * crossprod(t(U)) - Y - rho * Theta) + Dstar(y - rho * d) + rho * (LstarLw + DstarDw)
      grad <- grad + Astar( rho * (Aw-Z) + Lambda )
      eta <- 1 / (2*rho * (2*p - 1)+ 2*rho)
      wi <- w - eta *  grad
      thr <- sqrt(2*alpha *eta )
      wi[wi< thr] <- 0
      Lwi <- L(wi)
      Awi <- A(wi)
    }
    
    # Update Z
    Z <- Awi + Lambda/rho
    for (i in 1:p){
      norm_i <- norm(Z[i,], type="2")
      if (norm_i>t) {
        Z[i,] <- Z[i,]/norm_i *t
      }
    }
    
    
    # update U
    U <- eigen(Lwi, symmetric = TRUE)$vectors[, (p - k + 1):p]
    
    # update Theta
    eig <- eigen(rho * Lwi - Y, symmetric = TRUE)
    V <- eig$vectors[,1:(p-k)]
    gamma <- eig$values[1:(p-k)]
    Thetai <- V %*% diag((gamma + sqrt(gamma^2 + 4 * rho)) / (2 * rho)) %*% t(V)
    
    # update Y
    R1 <- Thetai - Lwi
    Y <- Y + rho * R1
    
    # update y
    R2 <- diag(Lwi) - d
    y <- y + rho * R2
    
    R3 = Awi - Z
    Lambda <- Lambda + rho * R3
    # compute primal, dual residuals, & lagrangian
    primal_lap_residual <- c(primal_lap_residual, norm(R1, "F"))
    primal_deg_residual <- c(primal_deg_residual, norm(R2, "2"))
    dual_residual <- c(dual_residual, rho*norm(Lstar(Theta - Thetai), "2"))
    lagrangian <- c(lagrangian, compute_augmented_lagrangian_kcomp_ht_balanced(wi, LstarSq, Thetai, U, Y, y, d, heavy_type, n, p, k, rho, beta, nu, alpha, Z, Lambda))
    
    # update rho
    if (update_rho) {
      eig_vals <- spectralGraphTopology:::eigval_sym(Theta)
      n_zero_eigenvalues <- sum(eig_vals < 1e-9)
      if (k < n_zero_eigenvalues)
        rho <- .5 * rho
      else if (k > n_zero_eigenvalues)
        rho <- 2 * rho
      else {
        if (early_stopping) {
          has_converged <- TRUE
          break
        }
      }
    }
    if (update_beta) {
      eig_vals <- spectralGraphTopology:::eigval_sym(L(wi))
      n_zero_eigenvalues <- sum(eig_vals < 1e-9)
      if (k < n_zero_eigenvalues)
        beta <- .5 * beta
      else if (k > n_zero_eigenvalues)
        beta <- 2 * beta
      else {
        if (early_stopping) {
          has_converged <- TRUE
          break
        }
      }
      beta_seq <- c(beta_seq, beta)
    }
    if (verbose)
      pb$tick()
    
    elapsed_time <- c(elapsed_time, proc.time()[3] - start_time)
    has_converged <- (norm(Lwi - Lw, 'F') / norm(Lw, 'F') < reltol) && (i > 1)
    if (has_converged)
      break
    w <- wi
    Lw <- Lwi
    Aw <- Awi
    
    Theta <- Thetai
  }
  results <- list(laplacian = L(wi), adjacency = A(wi), theta = Thetai, maxiter = i,
                  convergence = has_converged, beta_seq = beta_seq,
                  primal_lap_residual = primal_lap_residual,
                  primal_deg_residual = primal_deg_residual,
                  dual_residual = dual_residual,
                  lagrangian = lagrangian,
                  elapsed_time = elapsed_time)
  return(results)
}

compute_augmented_lagrangian_kcomp_ht_balanced <- function(w, LstarSq, Theta, U, Y, y, d, heavy_type, n, p, k, rho, beta, nu, alpha, Z, Lambda) {
  eig <- eigen(Theta, symmetric = TRUE, only.values = TRUE)$values[1:(p-k)]
  Lw <- L(w)
  Aw <- A(w)
  Dw <- diag(Lw)
  u_func <- 0
  if (heavy_type == "student") {
    for (q in 1:n)
      u_func <- u_func + (p + nu) * log(1 + n * sum(w * LstarSq[[q]]) / nu)
  } else if (heavy_type == "gaussian"){
    for (q in 1:n)
      u_func <- u_func + sum(n * w * LstarSq[[q]])
  }
  u_func <- u_func / n
  return(u_func - sum(log(eig)) + sum(y * (Dw - d)) + sum(diag(Y %*% (Theta - Lw)))
         + .5 * rho * (norm(Dw - d, "2")^2 + norm(Lw - Theta, "F")^2) + beta * sum(w * Lstar(crossprod(t(U))))
         + .5 * rho * norm(Aw - Z, "2")^2 + sum(diag(Lambda %*% (Aw - Z)))
        + alpha * sum(w>0)  )
}

hardThresh <- function(v, thr){
  
  return( v * (abs(v) > thr) )
  
}
