EM_KM<-function(y, Phi, X,mu,Lambda,V,n,N,p,max_iter,B_initial){
  ### initialize B and Psi
  B<-B_initial
  Psi<-diag(rep(1e-3,N))
  W_t<-array(NA, c(n, N, N))
  M<-array(NA, c(n, 1, N))
  v_t<-rep(NA, n)
  m_t<-rep(NA, n)
  nll_table<-rep(NA,max_iter)
  B_trace<-array(NA, c(N,p,max_iter))
  Psi_trace<-array(NA, c(N,N,max_iter))
  for(iter in 1:max_iter){
    #### E-step
    for(t in 1:n){
      W_t[t,,]<-B%*%X[t,]%*%t(B%*%X[t,])+Psi
    }
    ### run kalman smoother first
    ks_result<-ks_ns(y,Phi,mu,Lambda,V, W_t)
    beta_smoother<-ks_result$beta_smoother
    #beta_smoother<-beta
    Sigma_smoothing<-ks_result$Pt_smoothing
    nll_table[iter]<-ks_result$nll
    B_trace[,,iter]<-B
    Psi_trace[,,iter]<-Psi
    Pt_smoother<-array(NA, c(n, N, N))
    ### update beta_t|y_1:n and gamma_t|y_1:n
    for(t in 1:n){
      Pt_smoother[t,,]<-Phi%*%Sigma_smoothing[,,t]%*%t(Phi)
      #Pt_smoother[t,,]<-matrix(rep(0, N^2), N,N)
      v_t[t]<-(1+t(B%*%X[t,])%*%solve(Psi)%*%B%*%X[t,])^(-1)
      M[t,,]<-t(B%*%X[t,])%*%solve(Psi)*v_t[t]
      m_t[t]<-M[t,,]%*%(y[t,]-Phi%*%beta_smoother[t,])
    }
    #### M-step
    ## B
    B_term1<-0
    B_term2<-0
    for(t in 1:n){
      M_t<-matrix(M[t,,], nrow=1, ncol=N)
      B_term1<-B_term1+(m_t[t]*(y[t,]-Phi%*%beta_smoother[t,])+Pt_smoother[t,,]%*%t(M_t))%*%t(X[t,])
      #B_term1<-B_term1+(m_t[t]*(y[t,]-Phi%*%beta[t,])+Pt_smoother[t,,]%*%t(M_t))%*%t(X[t,])
      B_term2<-B_term2+(m_t[t]^2+M_t%*%Pt_smoother[t,,]%*%t(M_t)+
        v_t[t])[1,1]*X[t,]%*%t(X[t,])
    }
    B<-B_term1%*%solve(B_term2)
    B<-abs(B)
    ## Psi
    Psi_term<-0
    for(t in 1:n){
      M_t<-matrix(M[t,,], nrow=1, ncol=N)
      Psi_term<-Psi_term+(y[t,]-Phi%*%beta_smoother[t,]-m_t[t]*B%*%X[t,])%*%
                t(y[t,]-Phi%*%beta_smoother[t,]-m_t[t]*B%*%X[t,])+
                (diag(rep(1, N))-B%*%X[t,]%*%M_t)%*%Pt_smoother[t,,]%*%
                t(diag(rep(1, N))-B%*%X[t,]%*%M_t)+v_t[t]*B%*%X[t,]%*%t(B%*%X[t,])
    }
    Psi<-diag(diag(Psi_term/n))
    #Psi<-diag(rep(,N))
  }
  return(list(B=B,Psi=Psi, nll=nll_table))
}
