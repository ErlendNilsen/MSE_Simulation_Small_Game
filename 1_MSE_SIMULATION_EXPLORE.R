

### SETTING WD AND LOADING LIBRARIES;

setwd("C:/Users/erlend.nilsen/OneDrive - NINA/Prosjekter/Arbeidsbord/Ryper/Lagopus MSE/MSE_RFF")
library(MSEtools) # devtools::install_git("https://github.com/ErlendNilsen/MSEtools")
library(tidyverse)
library(MDPtoolbox)


#### Using the TL-model and 
#### vizualising the relationship between growth rate
#### and density

r_low <- 0.7
r_high <- 0.7

x <- seq(1,50, length.out = 100)
lam1 <- lam2 <- numeric()

for(i in 1:length(x)){
  
  lam1[i] <- PopMod_TLpre(X_t0=x[i], sigma2_e=0, N_Harv=0, K=20, theta=0.3, r=r_low)$r
  lam2[i] <- PopMod_TLpre(X_t0=x[i], sigma2_e=0, N_Harv=0, K=20, theta=0.9, r=r_high)$r
}

par(mfrow=c(1,2))

plot(x, lam2, ylim=c(-1, 1), xlab="Population size (X)", ylab="r", type="l")
lines(x, lam1, col="red")
mtext("A)", adj=0.05)

plot(x, exp(lam2), xlab="log(X)", ylab=expression(lambda), type="l")
lines(x, exp(lam1), col="red")
mtext("B)", adj=0.05)


####

cpue <- 
  
  Impl3b(Effort=E_sim, alpha=alpha_sim, beta=(log10(qb)-(log10(alpha_sim)))*-1, Pop_state=x_sim[t,i], Est_Pop_state=Est_PopState[t,i], c=c_sim)  
  
#####################################################################################################################
#####################################################################################################################
### Assuming perfect controllability; 
# Control parameters - population model

K_km2 <- 20
A <- 1000

sim_K <- K_km2*A
sim_theta <- c(0.3, 1)
r_sim <- c(0.5, 0.75, 1.0)
s2e_sim <- c(0.015, 0.15, 0.3)
cv_monitor <- c(0, 0.125, 0.25, 0.5)
p_sim <- c(0.5, 0.75, 1)

# Control parameter - harvest decision model
q <- seq(0, 0.8, by=0.05)
n_sim <- 1000
t_step <- 25


sim_TAC1 <- X_star1 <- lam <- obs_xt <- array(dim=c(t_step+1, n_sim, length(q), length(s2e_sim), length(r_sim), length(sim_theta), length(cv_monitor), length(p_sim)))
x_sim1 <- H1 <- array(0, dim=c(t_step+1, n_sim, length(q), length(s2e_sim), length(r_sim), length(sim_theta), length(cv_monitor), length(p_sim)))

for(l in length(p_sim)){
  for(m in 1:length(cv_monitor)){
    for(n in 1:length(sim_theta)){
      for (k in 1:length(r_sim)){  ## growth rate
        for(s in 1:length(s2e_sim)){  ## Environ. stochasticity
          for(o in 1:length(q)){
            
            for(i in 1:n_sim){
              k_scale <- runif(1, min=0.2, max=1.2)
              x_sim1[1,i,o,s,k,n,m,l] <- k_scale*sim_K  ## Initial population size
              
              
              for (t in 1:t_step){
                
                obs_xt[t,i,o,s,k,n,m,l] <- obs_mod1(scale="Abund", value=x_sim1[t,i,o,s,k,n,m,l], bias=1, cv=cv_monitor[m], LogNorm="LND")
                
                sim_TAC1[t,i,o,s,k,n,m,l] <- min(x_sim1[t,i,o,s,k,n,m,l], HarvDec1 (HD_type="A", qu=q[o], PopState_est=obs_xt[t,i,o,s,k,n,m,l]))
                H1[t,i,o,s,k,n,m,l] <-  min(x_sim1[t,i,o,s,k,n,m,l]*0.9, Impl_rd(round(sim_TAC1[t,i,o,s,k,n,m,l]), ModType="A", p=p_sim[l]))
                
                temp <- PopMod_TLpre(X_t0=x_sim1[t,i,o,s,k,n,m,l], sigma2_e=s2e_sim[s], theta=sim_theta[n], N_Harv=H1[t,i,o,s,k,n,m,l], K=sim_K, r=r_sim[k])
                x_sim1[t+1,i,o,s,k,n,m,l] <- temp$X_t1
                X_star1[t+1,i,o,s,k,n,m,l] <- temp$X_star
                lam[t,i,o,s,k,n,m,l] <- temp$r
                if(x_sim1[t+1,i,o,s,k,n,m,l]<2*A)
                  break
                
              }
            } 
          } 
        }
      }}}}

####################################
####################################

x_sim1["time", "rep", "harvest rate", "e_var", "growth rate", "theta", "monitor", "p_q"]


### DIFFERENT GROWTH RATES; 

Harv1 <- Harv2 <- numeric()

  for(i in 1:length(q)){
  
Harv1[i] <- median(H1[20:25,,i,1,2,1,1,3])
Harv2[i] <- median(H1[20:25,,i,1,3,1,1,3])


}

par(bty="l", lwd=2, cex=1.5)
plot(q, Harv2, type="l", col="dark orange", lwd=3, ylab="Harvest bag", xlab="Harvest rate")
lines(q, Harv1, type="l", col="dark green", lwd=3)
legend(x=0.5, y=1500, col=c("dark orange", "dark green"), lwd=5, legend=c("r=1.00", "r=0.75"), bty="n", cex=1)


### Figure for "NO MSY -APPROACH"
plot(q, Harv2, type="l", col="dark orange", lwd=3, ylab="Harvest bag", xlab="Harvest rate", xlim=c(0, 0.43), main="MSY")


### DIFFERENT theta; 

Harv3 <- Harv4 <- numeric()

for(i in 1:length(q)){
  
  Harv3[i] <- median(H1[20:25,,i,1,2,1,1,3])
  Harv4[i] <- median(H1[20:25,,i,1,2,2,1,3])
  
  
}

par(bty="l", lwd=2, cex=1.5)
plot(q, Harv4, type="l", col="dark orange", lwd=3, ylab="Harvest bag", xlab="Harvest rate")
lines(q, Harv3, type="l", col="dark green", lwd=3)
legend(x=0.55, y=1500, col=c("dark orange", "dark green"), lwd=5, legend=c("theta=1.00", "theta=0.30"), bty="n", cex=0.7)


### DIFFERENT e_var; 

Harv5 <- Harv6 <- numeric()

for(i in 1:length(q)){
  
  Harv5[i] <- median(H1[20:25,,i,1,2,1,1,3])
  Harv6[i] <- median(H1[20:25,,i,3,2,1,1,3])
  
  
}

par(bty="l", lwd=2, cex=1.5)
plot(q, Harv5, type="l", col="dark orange", lwd=3, ylab="Harvest bag", xlab="Harvest rate")
lines(q, Harv6, type="l", col="dark green", lwd=3)
legend(x=0.55, y=1500, col=c("dark orange", "dark green"), lwd=5, legend=c("e_varr=0.02", "e_var=0.15"), bty="n", cex=0.7)


### DIFFERENT CV_Monitor; 

Harv7 <- Harv8 <- numeric()

for(i in 1:length(q)){
  
  Harv7[i] <- median(H1[20:25,,i,2,2,1,1,3])
  Harv8[i] <- median(H1[20:25,,i,2,2,1,4,3])
  
  
}

par(bty="l", lwd=2, cex=1.5)
plot(q, Harv7, type="l", col="dark orange", lwd=3, ylab="Harvest bag", xlab="Harvest rate")
lines(q, Harv8, type="l", col="dark green", lwd=3)
legend(x=0.5, y=1300, col=c("dark orange", "dark green"), lwd=5, legend=c("cv_mon=0.00", "cv_mon=0.50"), bty="n", cex=0.7)


####################################################
### ; 

temp1 <- temp2 <- Via1 <- Via2 <- Harv10 <- numeric()

for(i in 1:length(q)){
  for(j in 1:n_sim){
    
  temp1[j] <- ifelse(min(x_sim1[,j,i,2,3,1,2,3])==0, 1, 0)
  temp2[j] <- ifelse(min(x_sim1[,j,i,2,3,1,2,3])<6000, 1, 0)
  
  }
  Via1[i] <- mean(temp1)
  Via2[i] <- mean(temp2)
  Harv10[i] <- median(H1[20:25,,i,2,3,1,2,3])
  }
  

  par(bty="l", lwd=2, cex=1.5)
  plot(q, Via1, type="l", col="dark orange", lwd=3, ylab="P(quasi-extinction)", xlab="Harvest rate", ylim=c(0,1))
  lines(q, Via2, type="l", col="dark green", lwd=3)
  legend(x=0.3, y=0.3, col=c("dark orange", "dark green"), lwd=5, 
         legend=c("threshold=2 ptarmigan/sq.km", "threshold=6 ptarmigan/sq.km"), bty="n", cex=0.7)
  abline(v=0.25, lty=3)


  
  ###########################################################################################################################
  ###########################################################################################################################
  
  EffDec_L1 <- function(slo_solgte=1.1403888, se_solgte=0.1442594,
                        slo_dag=0.7923, se_dag=0.0684,
                        slo_periode=0.3564, se_periode=0.0824,
                        slo_kombi=0.6326, se_kombi=0.1341,
                        solgte=0.75,
                        dag=1,
                        periode=0,
                        kombi=0,
                        pop_dens=PopDensSample){
    beta1 <- rnorm(1, mean=slo_solgte, sd=se_solgte)
    beta2 <- rnorm(1, mean=slo_dag, sd=se_dag)
    beta3 <- rnorm(1, mean=slo_periode, sd=se_periode)
    beta4 <- rnorm(1, mean=slo_kombi, sd=se_kombi)
    Effort <- exp(solgte*beta1 + dag*beta2 + periode*beta3 + kombi*beta4)
    Effort
  }
  
  
  Impl_L1 <- function(Effort="Effort", inte=-0.507206, se_i=0.133139, 
                      slo_dens=0.0553671, se_dens=0.0084397, 
                      pop_dens=PopDensSample){
    
    beta <- rnorm(1, mean=slo_dens, sd=se_dens)
    alpha <- rnorm(1, mean=inte, sd=se_i)
    
    H_I <- exp(alpha + (pop_dens*beta))*Effort
    H_I
    
  }
  
  
  
  
  
  ###########################################################################################################################
  ###########################################################################################################################
  
  
  K_km2 <- 20
  A <- 1000
  
  sim_K <- K_km2*A
 
  cv_monitor <- c(0.125, 0.5)
  threshold <- c(0, 4, 8, 10, 12)
  
  # Control parameter - harvest decision model
  sold <- seq(0, 3, by=0.05)
  n_sim <- 1000
  t_step <- 25
  
  
 X_star2 <- lam2 <- obs_xt2 <- array(dim=c(t_step+1, n_sim, length(sold), length(cv_monitor), length(threshold)))
x_sim2 <- H2 <- sim_Effort <- array(0, dim=c(t_step+1, n_sim, length(sold), length(cv_monitor), length(threshold)))
  
 
  sim_theta <- 0.3
  r_sim <- 0.75
  s2e_sim <- 0.15
            
            
            for(l in 1:length(threshold)){
            for(m in 1:length(cv_monitor)){
            for(n in 1:length(sold)){
              
              for(i in 1:n_sim){
                k_scale <- runif(1, min=0.2, max=1.2)
                x_sim2[1,i,n,m,l] <- k_scale*sim_K  ## Initial population size
                
                
                for (t in 1:t_step){
                  
                  obs_xt2[t,i,n,m,l] <- obs_mod1(scale="Abund", value=x_sim2[t,i,n,m,l], bias=1, cv=cv_monitor[m], LogNorm="LND")
                  
                  sim_Effort[t,i,n,m,l] <- ifelse(threshold[l]*A>obs_xt2[t,i,n,m,l], 0, EffDec_L1(slo_solgte=1.1403888, se_solgte=0.1442594, slo_dag=0.7923, se_dag=0.0684,
                                              slo_periode=0.3564, se_periode=0.0824, slo_kombi=0.6326, se_kombi=0.1341,
                                              solgte=sold[n]-0.6965233,
                                              dag=0,
                                              periode=1,
                                              kombi=0,
                                              pop_dens=x_sim2[t,i,n,m,l]-(12.1*A)))
                  H2[t,i,n,m,l] <- min(x_sim2[t,i,n,m,l]*0.9, (Impl_L1(Effort=sim_Effort[t,i,n,m,l], inte=-0.507206, se_i=0.1, 
                                                                      slo_dens=0.0553671, se_dens=0.001, pop_dens=(x_sim2[t,i,n,m,l]/A)-12.1))*A)
                  H2[t,i,n,m,l] <- ifelse(sold[n]==0, 0, H2[t,i,n,m,l])
                  
                  temp <- PopMod_TLpre(X_t0=x_sim2[t,i,n,m,l], sigma2_e=s2e_sim, theta=sim_theta, N_Harv=H2[t,i,n,m,l], K=sim_K, r=r_sim)
                  x_sim2[t+1,i,n,m,l] <- temp$X_t1
                  X_star2[t+1,i,n,m,l] <- temp$X_star
                  lam2[t,i,n,m,l] <- temp$r
                  if(x_sim2[t+1,i,n,m,l]<2*A)
                    break
                  
                }
              } 
            } 
          }
        }
  

 ##########################################
### Relationship sold vs effort;   
  x <- seq(0,3, by=0.1)
  y <- exp(x-0.6965233*1.14 + 0.3564)
  
  plot(x,y, type="l", xlab="Sold licences", ylab="Hunting effort", xlim=c(0,2.5), col="dark orange", yaxt="n")
  axis(side = 2, labels = F)
  
  ### DIFFERENT sold; 
  
  temp1 <- temp2 <- Via1 <- Via2 <- Harv10 <- numeric()
  
  for(i in 1:length(sold)){
    for(j in 1:n_sim){
      
      temp1[j] <- ifelse(min(x_sim2[,j,i,1,1])==0, 1, 0)
      temp2[j] <- ifelse(min(x_sim2[,j,i,1,1])<6000, 1, 0)
      
    }
    Via1[i] <- mean(temp1)
    Via2[i] <- mean(temp2)
    #Harv10[i] <- median(H1[20:25,,i,2,3,1,2,3])
  }
  
  
  par(bty="l", lwd=2, cex=1.5)
  plot(sold, Via1, type="l", col="dark orange", lwd=3, ylab="P(quasi-extinction)", xlab="Sold licences/sq.km", ylim=c(0,1), xlim=c(0, 3.5))
  lines(sold, Via2, type="l", col="dark green", lwd=3)
  legend(x=1, y=0.3, col=c("dark orange", "dark green"), lwd=5, 
         legend=c("Quasi-extiction=2 ptarmigan/sq.km", "Quasi-extiction=6 ptarmigan/sq.km"), bty="n", cex=0.7)

  
  ###################################
  ### Different sold & threshold; 
  
  ### DIFFERENT sold; 
  
  temp1 <- temp2 <- temp3 <- temp4 <- temp5 <- Via1 <- Via2 <- Via3 <- Via4 <- Via5 <- numeric()
  
  for(i in 1:length(sold)){
    for(j in 1:n_sim){
      
      temp1[j] <- ifelse(min(x_sim2[,j,i,1,1])==0, 1, 0)
      temp2[j] <- ifelse(min(x_sim2[,j,i,1,2])==0, 1, 0)
      temp3[j] <- ifelse(min(x_sim2[,j,i,1,3])==0, 1, 0)
      temp4[j] <- ifelse(min(x_sim2[,j,i,2,3])==0, 1, 0)
      temp5[j] <- ifelse(min(x_sim2[,j,i,2,4])==0, 1, 0)
      
      
    }
    Via1[i] <- mean(temp1)
    Via2[i] <- mean(temp2)
    Via3[i] <- mean(temp3) 
    Via4[i] <- mean(temp4)
    Via5[i] <- mean(temp5)
    #Harv10[i] <- median(H1[20:25,,i,2,3,1,2,3])
  }
  
  
  par(bty="l", lwd=2, cex=1.5)
  plot(sold, Via1, type="l", col="dark orange", lwd=3, ylab="P(quasi-extinction)", xlab="Sold licences/sq.km", ylim=c(0,1), xlim=c(0, 3.5))
  #lines(sold, Via2, type="l", col="dark green", lwd=3)
  lines(sold, Via3, type="l", col="dark green", lwd=3)
  lines(sold, Via4, type="l", col="dark green", lwd=3)
  #lines(sold, Via5, type="l", col="dark green", lwd=3)
  
  legend(x=0.05, y=0.95, col=c("dark orange", "dark green"), lwd=5, 
         legend=c("No harvest threshold", "Threshold harvest"), bty="n", cex=0.7)
  
  
  
    
  
