remotes::install_github("HMJiangGatech/huge")
library(huge)
## B_est, Psi_est: Estimated value of Reference curves and initial value for the target curve
## beta_smoother: Estimated value of Reference curves and initial value for the target curve
## In this function, instead of checking convergence, we repeatedly run 100 times
GLS_EM_RTS<-function(B_est, Psi_reg, Psi_est, beta_smoother, max_iter,delta,d,N){
  R<-diag(rep(1, d*N))
  for(rep in 1:100){
  L_matrix<-t(chol(R))
  for(i in 2:d){
      y_tilde<-0
      Phi_tilde<-list()
    for(j in 1:i){
        L[[j]]<-solve(L_matrix)[((i-1)*N+1):(i*N), ((j-1)*N+1):(j*N)]
        Phi_tilde[[j]]<-L[[j]]%*%Phi[[j]]
        y_tilde<-y_tilde+L[[j]]%*%t(y[, ((j-1)*N+1):(j*N)])
      }
      y_tilde<-t(y_tilde)
     for(j in 1:(i-1)){
        y_tilde<-y_tilde-beta_smoother[[j]]%*%t(Phi_tilde[[j]])
     }
    em_result2<-EM_KS_2_abo_reg(y_tilde, Phi_tilde[[i]], X_com,mu[[i]],Lambda[[i]],V[[i]],n,N,p,max_iter,
                B_est[[i]], B_est,Psi_est,L,i, delta, Psi_reg)
    B_est[[i]]<-em_result2[[1]]
    Psi_est[[i]]<-em_result2[[2]]
    beta_smoother[[i]]<-em_result2[[3]]
    }

  y_stand<-y
  for(i in 1:d){
    W_t<-array(NA,c(n, N,N))
    for(t in 1:n){
      W_t[t,,]<-B_est[[i]]%*%X_com[[i]][t, ]%*%t(B_est[[i]]%*%X_com[[i]][t, ])+Psi_est[[i]]
      y_stand[t,((i-1)*N+1):(i*N)]<-solve(t(chol(W_t[t,,])))%*%(y[t,((i-1)*N+1):(i*N)]-
      Phi[[i]]%*%beta_smoother[[i]][t,])
   }
  }
  huge_result<-huge(y_stand, method="glasso")
  Precision<-huge_result$icov[[4]]
  R<-solve(Precision)
  R[1:20,1:20]<-diag(1, N)
  R[21:40, 21:40]<-diag(1,N)
}
return(list(R=R, B=em_result2[[1]], Psi=em_result2[[2]], beta=em_result2[[3]]))
}
