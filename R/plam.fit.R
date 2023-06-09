# Fitting partial linear functional additive model
plam.fit <- function(X.train, Z.train, y.train, folds, X.pred, Z.pred, y.pred)
{
  X.train = as.matrix(X.train)
  Z.train = as.matrix(Z.train)
  n.train = length(y.train)
  nx = ncol(X.train)
  nz = ncol(Z.train)

  traindata = data.frame(cbind(y.train, X.train, Z.train))
  colnames(traindata) <- xznam <- c("y", paste0("x", 1:nx), paste0("z", 1:nz))
  preddata = data.frame(cbind(y.pred, X.pred, Z.pred))
  colnames(preddata) <- xznam

  plam <- as.formula(paste("y ~ ", paste(xznam[2:(1+nx)], collapse = " + "), " + ", paste0( "s(", paste(xznam[(2+nx):(1+nx+nz)], collapse = ", k = 5) + s(" ), ", k = 5)")))
  fit0 <- gam(plam, data = traindata, family = gaussian)
  pred <- predict(fit0, type = "response", newdata = preddata)

  # cv residuals
  cvresi = rep(NA, n.train)
  for(i in 1:length(folds)){
    fit1 <- gam(plam, data = traindata[-folds[[i]],], family = gaussian)
    cvpred <- predict(fit1, type = "response", newdata = traindata[folds[[i]],])
    cvresi[folds[[i]]] <- y.train[folds[[i]]] - as.vector(cvpred)
    # cat("cv", i, "\n")
  }

  return(list(muhat.train = fit0$fitted.values, ehat.train = fit0$residuals,
              muhat.pred = as.vector(pred), edf = sum(fit0$edf),
              cvehat.train = cvresi))
}

