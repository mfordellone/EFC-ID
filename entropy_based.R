#------------------------------------------------------------------#
#                                                                  #
#' Entropy-based Fuzzy C-means algorithm for interval-valued data  #
#'                                                                 #
#' ----------------------------------------------------------------#
#' 
#' Authors: Fordellone M. (2022)
#'
#' @param X A data frame or a matrix with n rows (the units) and J columns (from 1 to J/2 are the centers 
#' and from J/2+1 to J are the radii). Then, J/2 variables there are. 
#' @param K Number of clusters (a priori fixed).
#' @param p Parameter of entropy (a priori fixed).
#' @param RS Number of random starts of the algorithm (default: 100).
#' @param maxiter Maximum number of iterations per each random start (default: 100).
#' @param prep Type of data preprocessing: one among "none" (default) and "stand".
#' @param mode Method: "centroids" (default), "medoids" or "trimmed"

entropy_based <- function(X, K, p, RS = 100, maxiter = 100, prep = c("none", "stand"), 
                          mode = c("centroids", "medoids", "trimmed")){
  conv <- 1e-09
  n <- nrow(X)
  J <- ncol(X)
  if (prep == "none") {
    X <- na.omit(X)
  } 
  else {
    X <- scale(na.omit(X), center = FALSE, scale = TRUE)[, ]
  }
  Xc <- X[,c(1:(J/2))]
  Xr <- X[,((J/2)+1):J]
  values  <- vector(length(RS), mode = "numeric")
  n_iter  <- vector(length(RS), mode = "numeric")
  iter_ev <- vector(length = length(maxiter), mode = "numeric")
  func_opt <- 10 ^ 10 * sum(X ^ 2)
  for (rs in 1:RS) {
    set.seed(rs)
    U <- matrix(runif(n * K, 0, 1), nrow = n, ncol = K)
    U <- U / apply(U, 1, sum)
    U_old <- U + 1
    iter  <- 0
    d_ic1  <- as.matrix(dist(Xc) ^ 2)[, 1:K]
    d_ic2  <- as.matrix(dist(Xr) ^ 2)[, 1:K]
    while ((sum(abs(U_old - U)) > conv) && (iter < maxiter)) {
      iter  <- iter + 1
      U_old <- U
      if (mode == "centroids"){
        Cc    <- centroids_entr(X = Xc, U = U)$center
        Cr    <- centroids_entr(X = Xr, U = U)$center
      }
      if (mode == "medoids"){
          Cc <- Xc[medoids_update(X = Xc, U = U, k = K, n = n, distance = d_ic1),]
          Cr <- Xr[medoids_update(X = Xr, U = U, k = K, n = n, distance = d_ic2),]
      }
      if (mode == "trimmed"){
          pt     = 0.5
          Cc    <- trimmed(X = Xc, U = U, p = pt)$center
          Cr    <- trimmed(X = Xr, U = U, p = pt)$center
      }
      d_ic1 <- as.matrix(dist(rbind(Cc, Xc))^2)[- (1:K), 1:K]
      d_ic2 <- as.matrix(dist(rbind(Cr, Xr))^2)[- (1:K), 1:K]
      w     <- min(sum(U*d_ic1)/sum(U*(d_ic1+d_ic2)),0.5)
      U     <- memberships_entr(Xc = Xc, Xr = Xr, Cc = Cc, Cr = Cr, d_ic1 = d_ic1, d_ic2 = d_ic2, p = p, w = w)$memb
      U     <- U / apply(U, 1, sum)
      entropy <- p*U*log(x = U, base = exp(1))
      for (i in 1:n){
        for (k in 1:K){
          if (entropy[i,k]=="NaN"){
            entropy[i,k]=0
          }
        }
      }
      iter_ev[iter] <- sum((U)*((1-w)^2*d_ic1+(w^2)*d_ic2))-sum(entropy)
    }
    entropy <- p*U*log(x = U, base = exp(1))
    for (i in 1:n){
      for (k in 1:K){
        if (entropy[i,k]=="NaN"){
          entropy[i,k]=0
        }
      }
    }
    #U     <- U / apply(U, 1, sum)
    func <- sum((U)*((1-w)^2*d_ic1+(w^2)*d_ic2))-sum(entropy)
    values[rs] <- func
    n_iter[rs] <- iter
    if (func < func_opt){
      func_opt  <- func
      U_opt     <- U
      Cr_opt    <- Cr
      Cc_opt    <- Cc
      w_opt     <- w
      iter_opt  <- iter
      rownames(Cr_opt) <- paste("Cluster", 1:K)
      rownames(Cc_opt) <- paste("Cluster", 1:K)
      colnames(Cr_opt) <- colnames(Xr)
      colnames(Cc_opt) <- colnames(Xc)
    }
    print(paste("Entropy-based fuzzy c-means:", "Random start", rs, "->", "iter", n_iter[rs], "|", "Obj. Function", round(func_opt, 3)))
  }
  
  # Plot
  par(mfrow=c(1,1))
  plot(x = 1:iter_opt, y = iter_ev[1:iter_opt], lwd = 2, type = "b", lty = 2, xlab = "Iterations", ylab = "Objective function")
  cl <- cluster(U=U_opt)$cl
  
  par(mfrow=c(K,1))
  for (k in 1:K){
    barplot(U_opt[,k], col = round(k), xlab = "Observations", ylab = expression(u[k]), names = c(1:n))
    abline(h = 0.6, lty = 2, lwd = 2, col = "gray")
  }
  
  if (J/2 == 2){
    library(graphics)
    par(mfrow=c(1,1))
    plot(x = X[,1], y = X[,2], lwd = 2, pch = 16, col = cl, 
         main = "Interval-valued data", xlab = "Variable 1", ylab = "Variable 2", xlim = c(min(X[,1])-5,max(X[,1])+5), ylim = c(min(X[,2])-5,max(X[,2])+5))
    for (i in 1:n){
      rect(xleft = (X[i,1]-X[i,3]), xright = (X[i,1]+X[i,3]), ybottom = (X[i,2]-X[i,4]), ytop = (X[i,2]+X[i,4]), density = NULL, border = cl[i])
      legend("topleft", legend = c(paste("Group", 1:3), "No Grouped"), col = c(1:3, 8), pch = 16, cex = 1.3, bg='lightblue')
    }
    par(mfrow=c(1,1))
    plot(x = Xc[,1], y = Xc[,2], lwd = 2, pch = cl, col = cl, xlab = "Variable 1", ylab = "Variable 2", main = "Centers")
    legend("topleft", legend = c(paste("Group", 1:K), "No Grouped"), col = c(1:K, 8), pch = c(1:K, 8), cex = 1.3, bg='lightblue')
    par(mfrow=c(1,1))
    plot(x = Xr[,1], y = Xr[,2], lwd = 2, pch = cl, col = cl, xlab = "Variable 1", ylab = "Variable 2", main = "Radii")
    legend("topleft", legend = c(paste("Group", 1:K), "No Grouped"), col = c(1:K, 8), pch = c(1:K, 8), cex = 1.3, bg='lightblue')
  }
  
  names(values)    <- paste("Random start", 1:RS, sep = " ")
  names(n_iter)    <- names(values)
  results          <- list()
  results$U        <- round(U_opt, 4)
  results$Cc       <- round(Cc_opt, 4)
  results$Cr       <- round(Cr_opt, 4)
  results$w        <- round(w_opt, 4)
  results$values   <- values
  results$func_obj <- min(values)
  results$n_iter   <- iter_opt
  results$cluster  <- cl
  results$iter_ev  <- iter_ev[1:iter_opt]
  results$mode     <- mode
  return(results)
}  


# Trimmed function --------------------------------------------------------

trimmed <- function(X, U, p){
  J = ncol(X)
  K = ncol(U)
  center <- matrix(data = 0, nrow = K, ncol = J)
  for(j in 1:J){
    Q         <- quantile(X[,1], probs=c(p/2, 1-p/2), na.rm = FALSE)
    iqr       <- IQR(X[,j])
    up        <- Q[2]+1.5*iqr # Upper Range  
    low       <- Q[1]-1.5*iqr # Lower Range
    new_data  <- subset(X, X[,j] > (Q[1] - 1.5*iqr) & X[,j] < (Q[2]+1.5*iqr))
    new_U     <- subset(U, X[,j] > (Q[1] - 1.5*iqr) & X[,j] < (Q[2]+1.5*iqr))
  }
  for (j in 1:J){
    for (k in 1:K){
      center[k,j] <- weighted.mean(x = new_data[,j], w = new_U[,k])
    }
  }
  results <- list()
  results$center <- center
  return(results)
}

# Medoids function --------------------------------------------------------

medoids_update <- function(X, U, k, n, distance) {
  medoids <- vector(mode = 'numeric', length = k)
  for (c in 1:k) {
    min_dist_i <- 10^5 * sum(X ^ 2, na.rm = TRUE)
    for (i in 1:n) {
      dist_i <- sum(U[,c]*distance[i,])
      if ((dist_i < min_dist_i) &
          match(i, medoids, nomatch = 0) == 0) {
        min_dist_i <- dist_i
        medoids[c] <- i
      }
    }
  }
  return(medoids)
}

# Centroids function ------------------------------------------------------

centroids_entr <- function(X, U){
  K <- ncol(U)
  J <- ncol(X)
  center <- matrix(data = 0, nrow = K, ncol = J)
  for (j in 1:J){
    for (k in 1:K){
      center[k,j] <- weighted.mean(x = X[,j], w = U[,k])
    }
  }
  results <- list()
  results$center <- center
  return(results)
}

# Memberships function ----------------------------------------------------

memberships_entr <- function(Xc, Xr, Cc, Cr, d_ic1, d_ic2, p, w){
  n   <- nrow(Xc)
  K   <- nrow(Cc)
  Uc  <- matrix(data = 0, nrow = n, ncol = K)
  for (i in 1:n){
    Uc[i, ] <- ((1/(exp((p^(-1))*((1-w)^2*d_ic1[i, ]+(w^2)*d_ic2[i, ])))))/(sum(1/(exp((p^(-1))*(1-w)^2*d_ic1[i,]+(w^2)*d_ic2[i,]))))
    Uc[is.na(Uc)] <- 1
  }
  results <- list()
  results$memb  <- Uc
  return(results)
}  

# Clusters vector ---------------------------------------------------------

cluster <- function(U){
  n <- nrow(U)
  K <- ncol(U)
  Ucc <- matrix(data = 8, nrow = n, ncol = 1)
  for (i in 1:n){
    for (k in 1:K){
      if (U[i,k]>0.7){
        Ucc[i,1]=k
      }
    }
  }
  output    <- list()
  output$cl <- Ucc
  return(output)
}