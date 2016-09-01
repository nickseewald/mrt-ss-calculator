# Copyright 2016 Nicholas J. Seewald and Peng Liao
# This file is part of MRT-SS Calculator.

# MRT-SS Calculator is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MRT-SS Calculator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with MRT-SS Calculator.  If not, see <http://www.gnu.org/licenses/>.

# Contact the authors at nseewald@umich.edu or pengliao@umich.edu

SampleSize = function(beta_t, tau, delta, alpha0, beta0, setup, p, q, Nmax=1000){

  ### MRT Sample Size Calculation (Liao et. al 2015) ###

  ### INPUT:
  ### standardized treatment effect (beta_t), expected availability (tau),
  ### randomization probability (delta), significance level (alpha0)
  ### desired power (beta0), study setup (setup = list(days, occ.per.days))
  ### p, q = number of parameters in Z and B (assume Z are quadratic, linear or constant)
  ### Nmax is the maximum number of participants to acheive the power

  ### OUTPUT: Sample size

  para = function(order, Total){

    if(order == 0){

      return(matrix(1, Total, 1))

    }else{

      Z <- matrix(0, Total, order + 1)
      Z[,1] <- 1;
      for(i in 1:order){

        Z[, i+1] <- (rep(c(0:(setup$days-1)), each = setup$occ.per.day))^i

      }

    }

    return(Z)

  }


  Total <- setup$days * setup$occ.per.day;
  Z <- para(order = p-1, Total);

  stopifnot(length(beta_t) == Total)
  stopifnot(length(tau) == Total);

  if(length(delta) == 1){

    delta <- rep(delta, Total);

  }else if(length(delta) == setup$days){

    delta <- rep(delta, each = setup$occ.per.day);

  }else{

    stopifnot(length(delta) == Total)
  }


  sigma.beta <- solve(t(Z) %*% diag(tau * delta * (1-delta)) %*% Z)

  d <- solve(t(Z) %*% diag(tau) %*% Z) %*% t(Z)  %*% diag(tau) %*% beta_t


  N.all <- c(10:Nmax);


  PowerCal <- function(p, q, N, d, Sigma_beta, alpha){

    C <- N * t(d) %*% solve(Sigma_beta) %*% d ## noncentral para.
    df2 <- N-(p+q) ## degrees of freedom
    adj_c <- qf(1-alpha, df1 = p, df2 = df2) ### critical value

    # output the power
    power <- 1 - pf(q = adj_c, df1 = p, df2 = df2 , ncp = C)

    return(power)
  }

  power.all <- sapply(N.all, function(n) PowerCal(p,q,n,d,sigma.beta,alpha0))

  if(all(power.all < beta0)){

    stop(paste("Cannot attain",beta0,"power when sample size is below", N.try));

  }else{

    opt.index <- which(power.all>=beta0)[1]
    N <- N.all[opt.index]

  }


  return(N)
}

PowerCalculation = function(N, beta_t, tau, delta, alpha0, setup, p, q){

  ### MRT Power Calculation (Liao et. al 2015) ###

  ### INPUT:
  ### sample size (N)
  ### standardized treatment effect (beta_t), expected availability (tau),
  ### randomization probability (delta), significance level (alpha0)
  ### study setup (setup = list(days, occ.per.days))
  ### p, q = number of parameters in Z and B (assume Z are quadratic, linear or constant)

  ### OUTPUT: (Estimated) power

  para = function(order, Total){

    if(order == 0){

      return(matrix(1, Total, 1))

    }else{

      Z <- matrix(0, Total, order + 1)
      Z[,1] <- 1;
      for(i in 1:order){

        Z[, i+1] <- (rep(c(0:(setup$days-1)), each = setup$occ.per.day))^i

      }

    }

    return(Z)

  }


  Total <- setup$days * setup$occ.per.day;
  Z <- para(order = p-1, Total);

  stopifnot(length(beta_t) == Total)
  stopifnot(length(tau) == Total);


  if(length(delta) == 1){

    delta <- rep(delta, Total);

  }else if(length(delta) == setup$days){

    delta <- rep(delta, each = setup$occ.per.day);

  }else{

    stopifnot(length(delta) == Total)
  }



  sigma.beta <- solve(t(Z) %*% diag(tau * delta * (1-delta)) %*% Z)

  d <- solve(t(Z) %*% diag(tau) %*% Z) %*% t(Z)  %*% diag(tau) %*% beta_t


  PowerCal <- function(p, q, N, d, Sigma_beta, alpha){

    C <- N * t(d) %*% solve(Sigma_beta) %*% d ## noncentral para.
    df2 <- N-(p+q) ## degrees of freedom
    adj_c <- qf(1-alpha, df1 = p, df2 = df2) ### critical value

    # output the power
    power <- 1 - pf(q = adj_c, df1 = p, df2 = df2 , ncp = C)

    return(power)
  }

  power <- PowerCal(p,q,N,d,sigma.beta,alpha0)


  return(power)

}
