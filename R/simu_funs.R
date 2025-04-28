# Simulate functional data (fd)
# Generate eigen functions for fd
eigen_fun <- function(K0, tt)
{
  # tt: time points
  nT = length(tt)
  Phi = matrix(0, K0, nT)

  if(max(tt) > 1){
    # type 1
    for(k in 1:K0){
      Phi[k,] = cos(pi*k*tt/5) / sqrt(5)
    }

  }else{
    # type 2
    for(k in 1:K0){
      Phi[k,] = sqrt(2)*sin(pi*k*tt)
    }

  }
  return(Phi)     # K0*nT matrix
}


# Simulate sample data of predictors
Zt_gen <- function(n, K0, xi01, tt, v)
{
  nT = length(tt)

  if(max(tt) > 1){

    lamb = 15*(1:K0)^(-2)
    xi0 = matrix(rnorm(n*K0), n, K0)
    xi0 = cbind(xi01, xi0[,-1])
    xi0 = sweep(xi0, 2, sqrt(lamb), FUN = "*")
    Phi = eigen_fun(K0, tt)       # K0*nT, eigen functions

    e = matrix(sqrt(v)*rnorm(n*nT), n, nT)
    Zt = xi0%*%Phi + e
    xi = pnorm(xi0, sd = rep(sqrt(lamb), each = n))

  }else{

    lamb = (1:K0)^(-1.5)                   # 1*K0, variance sequence / eigen values
    xi0 = matrix(rnorm(n*K0), n, K0)       # n*K0 matrix, loadings / coeficients
    if(!is.null(xi01)){
      xi0 = cbind(xi01, xi0[,-1])
    }
    xi0 = sweep(xi0, 2, sqrt(lamb), FUN = '*')
    Phi = eigen_fun(K0, tt)     # K0*nT, eigen functions

    e = matrix(sqrt(v)*rnorm(n*nT), n, nT)
    Zt = xi0%*%Phi + e
    xi = pnorm(xi0, sd = rep(sqrt(lamb), each = n))    # n*K0 matrix, transformed pc scores
  }

  return(list(xi = xi, Zt = Zt, xi0 = xi0))
}


# additive or linear functions for fd effects
gg <- function(u, typ)
{
  # u: transformed fpc scores, n*K0 matrix
  if(typ == 1){

    gg = 1.5*((u[,1]-1/2)^2-1/12) + u[,2]-1/2 + 1.5*(sin(pi*u[,3])-2/pi) + c( (u[,4:(dim(u)[2])] - 1/2)%*%(1/(4:(dim(u)[2]))))

  }else if(typ == 2){
    # linear effect
    # u: Zt, n*nT matrix, t in [0,1]
    tt = seq(0, 1, by = 0.01)[-1]
    gg = c(u %*% (1+log(1 + tt))) * 0.01

  }

  return(c(gg)) # 1*n vector

}
