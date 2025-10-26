Mode <- function(v) {
  uniqv <- unique(v)
  return ( uniqv[which.max(tabulate(match(v, uniqv)))] )
}


Balancedness <- function(netobj, p, q) {
  clust_size <- igraph::components(netobj)$csize
  k <- length(clust_size)
  if (k<q)
  {
    clust_size <- c(rep(0,q-k), clust_size)
  }
  return( sum(abs(clust_size - p/q))/q )
}



Balancedness_norm <- function(netobj, p, q) {
  clust_size <- igraph::components(netobj)$csize
  k <- length(clust_size)
  if (k<q)
  {
    clust_size <- c(rep(0,q-k), clust_size)
  }
  num <- sum(abs(clust_size - p/q))
  clust_size <- c(rep(1,(q-1)), p-(q-1))
  dennum <- sum(abs(clust_size - p/q))
  return( 1 - num/dennum )
}


GINI <- function(netobj) {
  clust_size <- igraph::components(netobj)$csize
  k <- length(clust_size)
  clust_size_avg <- mean(clust_size)
  num <- 0 
  for (i in 1:k)
  {
    for (j in 1:k) {
      num <- num + abs(clust_size[i]-clust_size[j]) }
  }
  return(num/ (2*k^2*clust_size_avg))
}

simplex_project <- function(x0, lb) {
  n <- length(x0)
  x0_pos <- x0
  x0_pos[x0_pos<0] <- 0
  
  if (sum(x0_pos) >= lb) {
    alpha <- 0
  }
  
  else if (sum(x0>=0) == n) {
    alpha <- (lb - sum(x0)) / n  
  }
  
  else{
    
    x0_sorted <- sort(x0, decreasing = TRUE)
    id_neg <- which(x0_sorted<0)
    if (id_neg[1]>1) {
      id_neg <- c(id_neg[1]-1,id_neg)
    }
    
    x0_csum <- cumsum(x0_sorted)
    
    for (j in 1:length(id_neg)) {
      id <- id_neg[j]
      alpha <- (lb - x0_csum[id])/id
      
      if (id<n) {
        if  (alpha>=0 & -alpha < x0_sorted[id] & -alpha >= x0_sorted[id+1]) {
          break
        }    
      }
      
    }
    
  }
  
  x <- x0  + alpha
  x[x<0] = 0
  
  return(x)
}






simplex_project_lb_ub <- function(x0, lb = 0, ub = Inf) {
  # Ensure x0 is a column vector
  x0 <- as.vector(x0)
  
  # Define the length of x0
  n <- length(x0)
  
  # Set negative elements in x0 to 0
  x0_pos <- pmax(x0, 0)
  
  # Sort x0 in descending order
  x0_sorted <- sort(x0, decreasing = TRUE)
  
  # Cumulative sum of sorted x0
  x0_csum <- cumsum(x0_sorted)
  
  # Initialize beta and alpha
  beta <- 0
  alpha <- 0
  
  # Calculate beta for the upper bound constraint
  if (sum(x0_pos) <= ub) {
    beta <- 0
  } else if (sum(x0 >= 0) == n) {
    beta <- (sum(x0) - ub) / n
  } else {
    for (j in 1:n) {
      beta <- (x0_csum[j] - ub) / j
      if (j < n) {
        if (beta >= 0 && beta < x0_sorted[j] && beta >= x0_sorted[j + 1]) {
          break
        }
      } else {
        if (beta >= 0 && beta < x0_sorted[j]) {
          break
        }
      }
    }
  }
  
  # Calculate alpha for the lower bound constraint
  if (sum(x0_pos) >= lb) {
    alpha <- 0
  } else if (sum(x0 >= 0) == n) {
    alpha <- (lb - sum(x0)) / n
  } else {
    id_neg <- which(x0_sorted < 0)
    if (length(id_neg) > 0 && id_neg[1] > 1) {
      id_neg <- c(id_neg[1] - 1, id_neg)
    }
    
    for (j in seq_along(id_neg)) {
      id <- id_neg[j]
      alpha <- (lb - x0_csum[id]) / id
      if (id < n) {
        if (alpha >= 0 && -alpha < x0_sorted[id] && -alpha >= x0_sorted[id + 1]) {
          break
        }
      } else {
        if (alpha >= 0 && -alpha < x0_sorted[id]) {
          break
        }
      }
    }
  }
  
  # Project x0 using the calculated alpha and beta
  x <- x0 - beta + alpha
  x <- pmax(x, 0)
  # x <- pmax(x0 - beta, lb)
  # x <- pmin(x + alpha, ub)
  
  return(x)
}




project_circle <- function(z0, d, u) {
  beta <- min(z0)
  itermax <- 1000
  reltol <- 1e-6
  p <- length(z0)
  c <- (sum(z0)-d)/p
  pu <- p*u
  
  for (j in 1:itermax) {
    z0_sub <- z0 - rep(beta,p)
    z0_sub_pos <- z0_sub
    z0_sub_pos[z0_sub_pos<0] <- 0;
    z0_sub_neg <- -z0_sub;
    z0_sub_neg[z0_sub_neg<0] <- 0;
    
    beta_new <- ( u*sum(z0_sub_neg) - d*max(0,norm(z0_sub_pos,'2')-u)  )/pu + c;
    if ( norm(beta - beta_new,'2')/norm(beta,'2') < reltol ) {
      break
    }
    
    beta <- beta_new
  }
  
  z <- z0_sub_pos
  z_norm <- norm(z,'2');
  if (z_norm>u) {
    z <- z/z_norm * u;
  }
  
  return( z )
}



evaluate_clustering <- function(net, true_labels, p, q) {
  labels_pred <- rep(0, p)
  memberships <- igraph::components(net)$membership
  for (j in 1:q){
    idx <- memberships %in% c(j)
    labels_pred[idx] <- Mode(true_labels[idx])
  }
  
  mask <- labels_pred != true_labels
  purity <- 1- sum(mask)/length(mask)
  
  
  labels_pred <- memberships
  labels_pred_sorted <- labels_pred 
  perms <- permn(c(1:q))
  acc_max <- 0
  for (k in 1:length(perms)){
    perm <- perms[[k]]
    for (j in 1:q){
      idx <- memberships %in% j
      labels_pred_sorted[idx] <- perm[j]
      
    }
    mask <- labels_pred_sorted != true_labels
    acc <- 1- sum(mask)/length(mask)
    if (acc>= acc_max) {
      acc_max <- acc
      ind_max <- k
    }
  }
  perm <- perms[[ind_max]]
  for (j in 1:q){
    idx <- labels_pred %in% j
    labels_pred_sorted[idx] <- perm[j]
  }
  
  NMI <- randnet::NMI(labels_pred_sorted, true_labels)
  
  ARI <- mclust::adjustedRandIndex(labels_pred_sorted, true_labels)
  mask <- labels_pred_sorted != true_labels
  accuracy <- 1- sum(mask)/length(mask)
  
  
  

  mod_gt <- modularity(net, true_labels)
  
  balanced_norm <- Balancedness_norm(net, p, q)
  GINI_metric <- GINI(net)
  
  metrics <- list( memberships = memberships,  labels_pred = labels_pred_sorted,  
                   accuracy = accuracy, purity = purity, 
                   mod = mod_gt,
                   balanced = balanced_norm, GINI = GINI_metric,
                   NMI = NMI, ARI = ARI)
  return( metrics )
}


