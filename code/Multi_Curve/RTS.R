### Kalman Filter and RTS for multi-curve
ks_multi<-function(y,Phi,mu,Lambda,V, W_t,N,i){
  ### run the kf first
    n<-nrow(y)
    #Phi<-t(NS_loading(tau, lambda))

    ### initialize beta_0|0 and sigma_0|0
    x0<-mu/(1-diag(Lambda))
    xt_filter<-x0
    Pt_filter<-diag(diag(V) / (1 - diag(Lambda)^2) )

    ###
    nll=0
    x_filter<-matrix(NA, nrow=n, ncol=3*i)
    x_prediction<-matrix(NA, nrow=n, ncol=3*i)
    y_prediction<-matrix(NA, nrow=n, ncol=N)
    Pt_filter_table<-array(NA, dim=c(3*i, 3*i, n))
    Pt_prediction_table<-array(NA, dim=c(3*i, 3*i, n))
    x_smoothing<-matrix(NA, nrow=n, ncol=3*i)
    Pt_smoothing_table<-array(NA, dim=c(3*i, 3*i, n))
    ### predicting and filtering
    for(t in 1:n){
      ### predicting
      xt_prediction  <- mu + Lambda %*% xt_filter # a_t+1|t
      Pt_prediction <- Lambda %*% Pt_filter %*% t(Lambda) + V # P_t+1|t
      yt_prediction <- Phi_tilde %*% xt_prediction # ytilde_t|t-1 = d_t + F_t a_t|t-1

      ### filtering
      et <- y[t, ] - yt_prediction # e_t = y_t - ytilde_t|t-1
      F_t <- Phi_tilde %*% Pt_prediction %*% t(Phi_tilde) + W_t[t,,] # Covariance matrix of et
      invF_t <- solve(F_t) # inverse of F
      K <- Pt_prediction %*% t(Phi_tilde) %*% invF_t # Kalman gain matrix: K_t
      xt_filter <- xt_prediction + K %*% et # a_t
      Pt_filter <- (diag(3*i) - K %*% Phi_tilde) %*% Pt_prediction # P_t

      ### nll, store the filter and prediction
      x_filter[t, ]<-xt_filter
      x_prediction[t, ]<-xt_prediction
      y_prediction[t, ]<-yt_prediction
      Pt_filter_table[,,t]<-Pt_filter
      Pt_prediction_table[, ,t]<-Pt_prediction
      nll=nll+(log(det(F_t))+t(et)%*%invF_t%*%et)
    }
    ### based on kf, run the smoothing
    ## we need x_prediction here
    ## also x_filter stands for beta_T|T
    x_smoothing[n, ]<-xt_filter
    xt_smoothing<-xt_filter
    Pt_smoothing<-Pt_filter
    Pt_smoothing_table[,,n]<-Pt_filter_table[,,n]
    for(t in (n-1):1){
      J<-Pt_filter_table[,,t]%*%t(Lambda)%*%Pt_prediction_table[,,(t+1)]
      xt_smoothing<-x_filter[t,]+J%*%(xt_smoothing-x_prediction[(t+1),])
      Pt_smoothing<-Pt_filter_table[,,t]+J%*%(Pt_smoothing-Pt_prediction_table[,,(t+1)])%*%t(J)
      x_smoothing[t,]<-xt_smoothing
      Pt_smoothing_table[,,t]<-Pt_smoothing
    }
  return(list(nll=nll,beta_filter=x_filter,Pt_filter=Pt_filter_table,
              beta_prediction=x_prediction,Pt_prediction=Pt_prediction_table,
              beta_smoother=x_smoothing,Pt_smoothing=Pt_smoothing_table))

}
