rm(list=ls())
set.seed(12345)

###############
# Library:
###############
library(INLA)
library(MASS)

#########################
# Length of series:
#########################
N <- 30
T <- 500
d <- 3

##################
# Parameters:
#################
tau <- 300
beta <- 0.4

# Covariance matrix \Sigma :
sigma11 <- 1/2
sigma22 <- 1/4
sigma33 <- 1/3
rho_12 <- 0.8
rho_13 <- 0.7
rho_23 <- 0.6

sigma12 <- rho_12*sqrt(sigma11)*sqrt(sigma22)
sigma13 <- rho_13*sqrt(sigma11)*sqrt(sigma33)
sigma23 <- rho_23*sqrt(sigma22)*sqrt(sigma33)

level.sigma.true <- matrix(c(sigma11,sigma12,sigma13
                             ,sigma12,sigma22,sigma23
                             ,sigma13,sigma23,sigma33),nrow=d)
# AR latent states:
phi1 <- 0.7
phi2 <- 0.7
phi3 <- 0.6

# Covariance matrix for the evolution of gamma_j,t
sigma_w1_sq <- 1/20
sigma_w2_sq <- 1/25
sigma_w3_sq <- 1/20

########################
# Simulation of AR(1):
########################
ar_sim <- function(rho_AR,sigma_w_sq,T){
  x <- as.vector(arima.sim(list(order = c(1,0,0), ar = rho_AR), n = T,sd=sqrt(sigma_w_sq)))
  return(x)
}

#########################
# Simulate x_{j,t}
#########################
rho_AR <- c(phi1,phi2,phi3)
sigma_w_sq <- c(sigma_w1_sq,sigma_w2_sq,sigma_w3_sq)
x_jt <- list()
for(j in 1:d){
  x_jt[[j]] <- ar_sim(rho_AR[j],sigma_w_sq[j],T=T)
}

x_jt_vec <- matrix(unlist(x_jt))
x_jt_df <- data.frame(j.index= rep(1:d,each=T),t.index=rep(1:T,d),x_jt_vec)


x_jit <- list()
for(i in 1:N){
  x_jit[[i]] <- data.frame(i.index=i,x_jt_df)
}

x_jt <- do.call(rbind,x_jit)

#alpha_it: (Matrix of dimension N*T ** d)
alpha_it_list <- list()
for(i in 1:N){
  alpha_it_list[[i]] <- data.frame(i.index=i,j.index=rep(1:d,each=T),t.index= rep(1:T,d)
                                   ,alpha=as.vector(mvrnorm(n = T,mu = rep(0,d),Sigma = level.sigma.true)))
}

alpha_it <- do.call(rbind,alpha_it_list)

# Static Regressor:
S_jit_mat <- data.frame(i.index=rep(1:N,each=T*d),j.index=rep(1:d,each=T*N)
                        ,t.index=rep(1:T,N*d),S_jit = rnorm(d*N*T,sd=1))

theta_jit <- data.frame(i.index=alpha_it$i.index,j.index= alpha_it$j.index
                        ,t.index=alpha_it$t.index
                        , theta=exp(alpha_it$alpha + x_jt$x_jt_vec + beta*S_jit_mat$S_jit))

# Formatting responses:
response_jit <- data.frame(i.index=theta_jit$i.index,j.index=theta_jit$j.index
                           ,t.index=theta_jit$t.index,theta=theta_jit$theta,response=NA)

for(k in 1:nrow(response_jit)){
  response_jit$response[k] <- rgamma(n = 1,shape = tau,scale = response_jit$theta[k]/tau)
}

# Response:
all_response <- response_jit
new_comp1 <- response_jit[response_jit$j.index==1,]$response
new_comp2 <- response_jit[response_jit$j.index==2,]$response
new_comp3 <- response_jit[response_jit$j.index==3,]$response

dat_res <- data.frame(new_comp1,new_comp2,new_comp3)
N_1 <- d*T

# Index for random effects (There are N replicated d-variate random effects)
alpha_index <- rep(c(1:N_1),N)

# Index for each component 
ind_1 <- rep(c(1:T, rep(NA,T),rep(NA,T)),N)
ind_2 <- rep(c(rep(NA,T), 1:T,rep(NA,T)),N)
ind_3 <- rep(c(rep(NA,T), rep(NA,T),1:T),N)

# Replicates:
rep_alpha_index <- rep(1:N,each=N_1)
rep_ind1 <-  rep(1:N, each=N_1)
rep_ind2 <-  rep(1:N, each=N_1)
rep_ind3 <-  rep(1:N, each=N_1)

response_df_new <- data.frame(all_response,alpha_index
                              ,rep_alpha_index,ind_1,ind_2,ind_3,
                              rep_ind1,rep_ind2,rep_ind3)
####################
## R-INLA setup:
####################

# Priors:
prior1 <- list(prec = list(prior="loggamma",param=c(1,0.1)))
prior2 <- list(prec = list(prior="loggamma",param=c(1,0.1)))
prior3 <- list(prec = list(prior="loggamma",param=c(1,0.1)))

formula.inla= response ~ -1+f(alpha_index, model="iid3d", n=N_1
                              ,replicate = rep_alpha_index)+ S_jit_mat$S_jit +
  f(ind_1,model="ar1", replicate = rep_ind1,hyper = prior1) +
  f(ind_2,model="ar1", replicate = rep_ind2,hyper = prior2) +
  f(ind_3,model="ar1", replicate = rep_ind3,hyper = prior3)

result= inla(formula.inla,family="gamma",
             data =response_df_new,
             control.compute=list(dic=FALSE),
             control.family = list(hyper = list(prec = list(prior="loggamma"
                                                            ,param=c(300,1)))),
             verbose = TRUE,control.inla=list(control.vb=list(enable=FALSE)))

## Results:

all_parameters_hyper <- result$summary.hyperpar
all_parameters_fixed <- result$summary.fixed
tab <- rbind.data.frame(all_parameters_fixed[,-7],all_parameters_hyper)
tab$true <- c(beta,tau,1/sigma11,1/sigma22,1/sigma33
              ,rho_12,rho_13,rho_23
              ,(1/sigma_w1_sq)*(1 - phi1^2)
              ,phi1,(1/sigma_w2_sq)*(1 - phi2^2)
              ,phi2,(1/sigma_w3_sq)*(1 - phi3^2)
              ,phi3)

row_id_marginal_prec <-  grep("Precision for ind",rownames(tab))
row_id_rho <-  grep("Rho for ind",rownames(tab))

tab[row_id_marginal_prec,] <- tab[row_id_marginal_prec,]/(1-tab[row_id_rho,]^2)
rownames(tab)[row_id_marginal_prec] <- c("1/W1",
                                         "1/W2",
                                         "1/W3")
# Renaming the rows:
rownames(tab)[1:8] <- c("beta","tau","1/sigma1^2","1/sigma2^2","1/sigma3^2"
                        ,"rho12","rho13","rho23")
rownames(tab)[c(10,12,14)] <- c("phi11","phi22","phi33")

# Rearranging the rows to make it similar to tale output:
tab[c(2,1,3:8,10,12,14,9,11,13),c(1,2,3,5,7)]
