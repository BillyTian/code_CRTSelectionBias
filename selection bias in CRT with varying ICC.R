library(nlme)

m_bar <- 500
n <- 20
cv <- 0
rho <- 0.1 #0.01

type <- c(rep(0,2), rep(1,4), rep(2,2), rep(4,2))
mu_a <- c(rep(c(1,2), 3), rep(1,4))
delta_a <- rep(c(0.2,0.8), 5)
lambda_a1 <- rep(c(0.1,0.2), 5)
lambda_a2 <- rep(c(0.1,0.3), 5)
mu_c <- c(1,2,2,1,1,2,2,2,2,2)
delta_c <- c(rep(c(0.2,0.8),2), rep(c(0.8,0.2),3))
lambda_c1 <- c(rep(c(0.1,0.2), 4), 0.2,0.1)
lambda_c2 <- c(rep(c(0.1,0.3), 4), 0.3,0.1)
table <- cbind(type, mu_a, delta_a, lambda_a1, lambda_a2,
               mu_c, delta_c, lambda_c1, lambda_c2)


empirical <- function(parameter, nsims=2000){
  #True parameters
  mu_a <- parameter[2]
  delta_a <- parameter[3]
  lambda_a <- parameter[4:5]
  mu_c <- parameter[6]
  delta_c <- parameter[7]
  lambda_c <- parameter[8:9]
  
  tau_R <- NULL
  truth <- NULL
  bias <- NULL
  beta <- NULL
  std.dev <- NULL
  LCI <- NULL
  UCI <- NULL
  
  count <- NULL
  for (j in 1:nsims){
    #j=1
    set.seed(2021+520*j)
    
    #Varying cluster sizes
    if (cv==0){
      m_vector <- rep(m_bar, n)
    } else{
      m_vector <- round(rgamma(n, shape=cv^(-2), rate=m_bar^(-1)*cv^(-2)))
      m_vector[m_vector<2] <- 2
    }
    
    cluster <- rep(1:n, m_vector)
    ind <- NULL
    for (k in 1:length(m_vector)){
      ind <- c(ind, 1:m_vector[k])
    }
    
    N <- sum(m_vector)
    #Generate two baseline covariates, infinite population regimen
    #Generate an individual covariate with cluster-specific population mean
    mu1 <- rnorm(n,0,1)
    X1 <- rnorm(N, rep(mu1, m_vector), 1)
    
    #Within each cluster, generate binary X2
    X2 <- NULL
    for (k in 1:length(m_vector)){
      X2 <- c(X2, rbinom(m_vector[k], 1, 0.4))
    }
    
    #Randomization
    Z <- rep(0,n)
    Z[sample(1:n, n/2)] <- 1
    Z <- rep(Z, m_vector)
    
    #Tune the mechanism parameters to reach approximately a:c:n=4:3:3
    vector_a <- c(0.3, 0.2, 0.1) 
    vector_c <- c(0.1, 0.2, -0.1)
    vector_n <- c(0, 0, 0)#for identifiability
    
    num_a <- exp(vector_a[1]+vector_a[2]*X1+vector_a[3]*X2)
    num_c <- exp(vector_c[1]+vector_c[2]*X1+vector_c[3]*X2)
    num_n <- exp(vector_n[1]+vector_n[2]*X1+vector_n[3]*X2)
    
    pa_vector <- num_a/(num_a + num_c + num_n)
    pc_vector <- num_c/(num_a + num_c + num_n)
    pn_vector <- num_n/(num_a + num_c + num_n)
    
    membership.matrix <- NULL
    for (i in 1:N){
      membership.matrix <- cbind(membership.matrix, as.vector(rmultinom(1,1,c(pa_vector[i], pc_vector[i], pn_vector[i]))))
    }
    a_membership <- membership.matrix[1,]
    c_membership <- membership.matrix[2,]
    n_membership <- membership.matrix[3,]
    
    pi_a <- mean(a_membership)
    pi_c <- mean(c_membership)
    pi_n <- mean(n_membership)
    
    G <- rep("a", N)
    G[c_membership==1] <- "c"
    G[n_membership==1] <- "n"
    
    pop_data <- data.frame(cbind(cluster, ind, X1, X2, Z))
    pop_data$G <- G
    trt_data <- pop_data[pop_data$Z==1,]
    ctr_data <- pop_data[pop_data$Z==0,]
    trt_data$R <- ifelse(trt_data$G=="a" | trt_data$G=="c", 1, 0)
    ctr_data$R <- ifelse(ctr_data$G=="a", 1, 0)
    trt_R <- trt_data[trt_data$R==1,]
    ctr_R <- ctr_data[ctr_data$R==1,]
    
    trt_cluster50_data <- NULL
    for (i in unique(trt_R$cluster)){
      trt_cluster50_data <- rbind(trt_cluster50_data, 
                                  trt_R[trt_R$cluster==i,][1:50,])
    }
    
    ctr_cluster50_data <- NULL
    for (i in unique(ctr_R$cluster)){
      ctr_cluster50_data <- rbind(ctr_cluster50_data, 
                                  ctr_R[ctr_R$cluster==i,][1:50,])
    }
    
    enroll_data <- rbind(trt_cluster50_data, ctr_cluster50_data)
    
    #Cluster-level random effect
    gamma <- rep(rnorm(n, 0, sqrt(rho)), each=50)
    epsilon <- rnorm(1000, 0, sqrt(1-rho))
    
    #Generate (potential) outcomes
    enroll_data$Y_a1 <- mu_a + delta_a + lambda_a[1]*enroll_data$X1 + lambda_a[2]*enroll_data$X2 + gamma + epsilon
    enroll_data$Y_a0 <- mu_a + lambda_a[1]*enroll_data$X1 + lambda_a[2]*enroll_data$X2 + gamma + epsilon
    enroll_data$Y_a <- mu_a + delta_a*enroll_data$Z + lambda_a[1]*enroll_data$X1 + lambda_a[2]*enroll_data$X2 + gamma + epsilon
    
    enroll_data$Y_c1 <- mu_c + delta_c + lambda_c[1]*enroll_data$X1 + lambda_c[2]*enroll_data$X2 + gamma + epsilon
    enroll_data$Y_c0 <- mu_c + lambda_c[1]*enroll_data$X1 + lambda_c[2]*enroll_data$X2 + gamma + epsilon
    enroll_data$Y_c <- mu_c + delta_c*enroll_data$Z + lambda_c[1]*enroll_data$X1 + lambda_c[2]*enroll_data$X2 + gamma + epsilon
    
    enroll_data$Y1[enroll_data$G=="a"] <- enroll_data$Y_a1[enroll_data$G=="a"]
    enroll_data$Y1[enroll_data$G=="c"] <- enroll_data$Y_c1[enroll_data$G=="c"]
    enroll_data$Y0[enroll_data$G=="a"] <- enroll_data$Y_a0[enroll_data$G=="a"]
    enroll_data$Y0[enroll_data$G=="c"] <- enroll_data$Y_c0[enroll_data$G=="c"]
    enroll_data$Y[enroll_data$G=="a"] <- enroll_data$Y_a[enroll_data$G=="a"]
    enroll_data$Y[enroll_data$G=="c"] <- enroll_data$Y_c[enroll_data$G=="c"]
    
    #True causal parameters
    tau_R[j] <- (2*pi_a+pi_c)/(2*pi_a+2*pi_c) * delta_a + pi_c/(2*pi_a+2*pi_c) * delta_c
    #tau_O <- pi_a/(pi_a+pi_c+pi_n) * delta_a + pi_c/(pi_a+pi_c+pi_n) * delta_c + pi_n/(pi_a+pi_c+pi_n) * delta_n
    #Causal effect by PO
    truth[j] <- mean(enroll_data$Y1)-mean(enroll_data$Y0) 
    
    
    fit <- try(lme(Y ~ Z + X1 + X2, data=enroll_data, random= ~ 1 | cluster), silent=T)
    if(class(fit)=="try-error"){next}
    count[j] <- 1
    
    beta[j] <- summary(fit)$tTable[2,1]
    std.dev[j] <- summary(fit)$tTable[2,2]
    bias[j] <- beta[j]-truth[j]
    LCI[j] <- beta[j]-qnorm(0.975)*std.dev[j]
    UCI[j] <- beta[j]+qnorm(0.975)*std.dev[j]
  }
  avg_tau_R <- mean(tau_R, na.rm=T)
  avg_truth <- mean(truth, na.rm=T)
  
  avg_bias <- mean(bias, na.rm=T)
  bias_prop <- abs(avg_bias/avg_truth)*100
  avg_sd <- mean(std.dev, na.rm=T)
  beta_sd <- sd(beta, na.rm=T)
  
  coverage_truth <- mean( truth<=UCI & truth>=LCI, na.rm=T )*100
  
  #nconv <- 1-sum(count, na.rm=T)/nsims
  return(c(avg_tau_R, avg_truth, 
           avg_bias, bias_prop, avg_sd, beta_sd,
           coverage_truth))
}

result <- NULL
for (i in 1:nrow(table)){
  result <- rbind(result, empirical(parameter=as.numeric(table[i,])))
}

result <- data.frame(cbind(table, result))
names(result)[10:16] <- c("avg.tau_R", "avg.truth",
                          "bias", "bias.prop (%)", "mc.sd", "est.sd",
                          "coverage.truth")

mainDir = '' #change the directory for outputs
write.table(result, paste0(mainDir, "/sb_cluster_ICC", rho,".csv"), sep =",", row.names=F)