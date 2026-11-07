EM_KS_2<-function(y, Phi_tilde, X,mu_mul,Lambda_mul,V_mul,n,N,p,max_iter,B_initial, B_est,Psi_est,L,i, delta, beta){
  ### initialize B and Psi
  B_est[[i]]<-B_initial
  Psi_tilde_pre<-delta*diag(rep(1, N))
  for(j in 1:(i-1)){
    Psi_tilde_pre<-Psi_tilde_pre+L[[j]]%*%Psi_est[[j]]%*%t(L[[j]])
  }
  Psi_est[[i]]<-diag(1e-3, N)
  W_update<-array(NA, c(n, N, N))
  M<-list()
  for(t in 1:n){
    M[[t]]<-array(NA, c(i, 1, N))
  }
  V_t<-array(NA, c(n, i, i))
  m_t<-matrix(NA, nrow=n, ncol=i)
  nll_table<-rep(NA,max_iter)
  #B_trace<-array(NA, c(N,p,max_iter))
  #Psi_trace<-array(NA, c(N,N,max_iter))
  W_initial<-array(NA, c(n, N,N))
  for(t in 1:n){
    W_initial[t,,]<-0
    for(j in 1:(i-1)){
      W_initial[t,,]<-W_initial[t,,]+L[[j]]%*%B_est[[j]]%*%X[[j]][t,]%*%t(L[[j]]%*%B_est[[j]]%*%X[[j]][t,])
    }
  }
  for(iter in 1:max_iter){
    Psi_tilde<-Psi_tilde_pre+L[[i]]%*%Psi_est[[i]]%*%t(L[[i]])
    #### E-step
    for(t in 1:n){
      W_update[t,,]<-W_initial[t,,]+(L[[i]]%*%B_est[[i]])%*%X[[i]][t,]%*%t(L[[i]]%*%B_est[[i]]%*%X[[i]][t,])+Psi_tilde
    }
    ### run kalman smoother first
    ks_result<-ks_multi(y,Phi_tilde,mu_mul,Lambda_mul,V_mul, W_update,N,i)
    beta_smoother<-ks_result$beta_smoother
    #beta_smoother<-beta
    #beta_smoother<-simulates$xt
    Sigma_smoothing<-ks_result$Pt_smoothing
    #Sima_smoothing<-matrix(0, nrow=3*i, ncol=3*i)
    #Sigma_smoothing<-0
    nll_table[iter]<-ks_result$nll
    #B_trace[,,iter]<-B
    #Psi_trace[,,iter]<-Psi
    Pt_smoother<-array(NA, c(n, N, N))
    ### update beta_t|y_1:n and gamma_t|y_1:n
    for(t in 1:n){
      Pt_smoother[t,,]<-Phi_tilde%*%Sigma_smoothing[,,t]%*%t(Phi_tilde)
      #Pt_smoother[t,,]<-matrix(rep(0, N^2), N,N)
      A<-matrix(NA, nrow=i, ncol=i)
      for(j in 1:i){
        for(k in 1:i){
          A[j,k]<-t(L[[j]]%*%B_est[[j]]%*%X[[j]][t,])%*%solve(Psi_tilde)%*%
                  L[[k]]%*%B_est[[k]]%*%X[[k]][t,]
        }
      }
      V_t[t,,]<-solve(diag(rep(1, i))+A)
      for(j in 1:i){
        M[[t]][j,,]<-0
        for(k in 1:i){
          M[[t]][j,,]<-M[[t]][j,,]+V_t[t,,][j, k]*t(L[[k]]%*%B_est[[k]]%*%X[[k]][t,])%*%solve(Psi_tilde)
        }
        m_t[t,j]<-M[[t]][j,,]%*%(y[t,]-Phi_tilde%*%beta_smoother[t,])
      }
    }

    #### M-step
    ## B
    B_term1<-0
    B_term2<-0
    for(t in 1:n){
      #M_ti<-matrix(M[[t]][i,,], nrow=1, ncol=N)
      #B_term1<-B_term1+(m_t[t]*(y[t,]-Phi_tilde%*%beta_smoother[t,])+Pt_smoother[t,,]%*%t(M_t))%*%t(X[t,])
      #B_term1<-B_term1+(m_t[t]*(y[t,]-Phi%*%beta[t,])+Pt_smoother[t,,]%*%t(M_t))%*%t(X[t,])
      B_sub_term1<-0
      B_sub_term2<-0
      B_sub_term3<-0
      for(j in 1:(i-1)){
        B_sub_term1<-B_sub_term1+L[[j]]%*%B_est[[j]]%*%X[[j]][t,]*m_t[t,j]
        B_sub_term2<-B_sub_term2+L[[j]]%*%B_est[[j]]%*%X[[j]][t,]%*%M[[t]][j,,]
        B_sub_term3<-B_sub_term3+L[[j]]%*%B_est[[j]]%*%X[[j]][t,]*V_t[t,i,j]
      }
      B_term1<-B_term1+(y[t,]-Phi_tilde%*%beta_smoother[t,]-B_sub_term1)%*%t(m_t[t,i]*X[[i]][t,])+
        (diag(rep(1, N))-B_sub_term2)%*%Pt_smoother[t,,]%*%t(X[[i]][t,]%*%M[[t]][i,,])-
        B_sub_term3%*%t(X[[i]][t,])
      B_term2<-B_term2+(m_t[t,i]^2+t(M[[t]][i,,])%*%Pt_smoother[t,,]%*%M[[t]][i,,]+
        V_t[t,i,i])*X[[i]][t,]%*%t(X[[i]][t,])
    }
    B_tilde<-B_term1%*%solve(B_term2)
    B_tilde<-B_tilde
    B_est[[i]]<-solve(L[[i]])%*%B_tilde
    ## Psi
    Psi_term_1<-0
    Psi_term_2<-0
    Psi_term_3<-0
    Psi_term_4<-0
    for(t in 1:n){
      Psi_sub_term1<-0
      Psi_sub_term2<-0
      #Psi_sub_term3<-0
      for(j in 1:i){
        Psi_sub_term1<-Psi_sub_term1+m_t[t, j]*L[[j]]%*%B_est[[j]]%*%X[[j]][t,]
        Psi_sub_term2<-Psi_sub_term2+L[[j]]%*%B_est[[j]]%*%X[[j]][t,]%*%M[[t]][j,,]
        for(k in 1:i){
          Psi_term_3<-Psi_term_3+V_t[t,j,k]*L[[j]]%*%B_est[[j]]%*%X[[j]][t,]%*%
             t(L[[k]]%*%B_est[[k]]%*%X[[k]][t,])
        }
      }
      Psi_term_1<-Psi_term_1+(y[t,]-Phi_tilde%*%beta_smoother[t,]-Psi_sub_term1)%*%
          t(y[t,]-Phi_tilde%*%beta_smoother[t,]-Psi_sub_term1)
      Psi_term_2<-Psi_term_2+(diag(rep(1,N))-Psi_sub_term2)%*%Pt_smoother[t,,]%*%
          t(diag(rep(1,N))-Psi_sub_term2)
    }
    Psi_tilde<-(Psi_term_1+Psi_term_2+Psi_term_3)/n
    for(j in 1:(i-1)){
       Psi_term_4<-Psi_term_4+L[[j]]%*%Psi_est[[j]]%*%t(L[[j]])
    }
    Psi_est[[i]]<-diag(diag(solve(L[[i]])%*%(Psi_tilde-Psi_term_4-delta*diag(rep(1,N)))%*%t(solve(L[[i]]))))
  }
  B<-solve(L[[i]])%*%B_tilde
  Psi<-diag(solve(L[[i]])%*%(Psi_tilde-Psi_term_4-delta*diag(rep(1,N)))%*%t(solve(L[[i]])))
  return(list(B=B,Psi=Psi, nll=nll_table))
}
