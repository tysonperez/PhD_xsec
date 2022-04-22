#Objective 6: 
#To explore the normal reference range for PRO scores and neurophysiological 
#measures. An additional study will be carried out to recruit age- and 
#sex-matched non-ID participants in a similar manner as our clinical trial.

#Load libraries
library(rjags)  
library(readxl)   
library(tidyverse)
library(bayesplot)
library(bayestestR)
library(see)
library(patchwork)

##CLINICAL##

#Load & visualise clinical data
clinical_data <- read_excel("DESS.xlsx")#, sheet="Obj6")
clinical_data

#Gather relevant data (baseline 2=t1) & add variables (tx)
get_data_obj6.1 <- function(data){
  #update last number in data[c(...)] to reflect column of interest
  tmp = data[c(1,2,25)]
  tmp$value = tmp[,3][[1]]
  tmp$tx <- 1*(tmp$Group=="ID") 
  return(tmp)
}
data_obj6.1 <- get_data_obj6.1(clinical_data)
data_obj6.1

##HEALTHY##

#Load & visualise healthy data
healthy_data <- read_excel("Healthy.xlsx",sheet="Obj6")
healthy_data

#Gather relevant data (baseline 2=t1) & add variables (tx)
get_data_obj6.2 <- function(data){
  #update last number in data[c(...)] to reflect column of interest
  tmp = data[c(1,2,180)]
  tmp$tx = 2*(tmp$Group=="non-ID")
  tmp$value = tmp[,3][[1]]
  return(tmp)
}
data_obj6.2 <- get_data_obj6.2(healthy_data)
data_obj6.2

#bind clinical & healthy data
data_obj6 <- rbind(data_obj6.1,data_obj6.2)
data_obj6

#Define the model
m6 <- "model{
        for(i in 1:n.p){
          y[i] ~ dnorm(mu.y[i],tau[tx[i]])
          mu.y[i] <- mu[tx[i]]
          
          #predicted data 
          y.pred[i] ~ dnorm(mu.y.pred[i],tau[tx[i]])
          mu.y.pred[i] <- mu[tx[i]]
          
          #calculate Pearson residuals^2 for observed & predicted data
          pearson[i] <- (y[i] - mu.y[i])/sqrt(tau[tx[i]])
          pearson.pred[i] <- (y.pred[i] - mu.y.pred[i])/sqrt(tau[tx[i]])
          D[i] <- pow(pearson[i],2)
          D.pred[i] <- pow(pearson.pred[i],2)
        }
        
        #vague prior
        for(i in 1:2){
        mu[i] ~ dnorm(15,0.0001)
        tau[i] <- 1/(sd[i]*sd[i])
        sd[i] ~ dt(0,0.033,3)T(0,)   #half t distribution 
        }
        
        #nodes to monitor
        ID <- mu[1]
        non_ID <- mu[2]
        ID_vs_non_ID <- mu[1]-mu[2]
        # std_ID <- mu[1]/sd[1]
        # std_non_ID <- mu[2]/sd[2]
        std_ID_v_non_ID <- (mu[1]-mu[2])/sqrt((sd[1]^2+sd[2]^2)/2)
    
        
        #summary stat (T) = sum the Pearson res^2 for observed & predicted data
        T <- sum(D[])
        T.pred <- sum(D.pred[])
        
        #calculate probability T.pred >= T (Bayesian p-value)
        Bayes.p <- step(T.pred - T)
    }"

#Generate initial values
inits_func <- function(){
  mu <- runif(2,0,15)
  sd <- runif(2,0,5)
  return(list(mu=mu,sd=sd))
}

#Compile the model
m6_jags <- jags.model(textConnection(m6),
                      data=list(y=data_obj6$value,
                                tx=data_obj6$tx,
                                n.p=length(data_obj6$Participant)),
                      inits = inits_func,
                      n.chains=3,n.adapt=5000
)
#Burn in
update(m6_jags,10000)

#Simulate the posterior
m6_out <- coda.samples(model=m6_jags,
                       variable.names=c("ID","non_ID","ID_vs_non_ID",
                                        #"std_ID","std_non_ID",
                                        "std_ID_v_non_ID",
                                        "sd","T","T.pred","Bayes.p"),
                       n.iter=25000)

######Posteriors######

#plot posteriors
post_plot <- function(data,parameter){
  color_scheme_set("brightblue")
  mcmc_areas(data, 
             prob=0.95, 
             point_est="mean",
             pars=parameter) #+
  # labs(title="Posterior distributions", 
  #      subtitle="with means & 95% HDI")
}
#post_plot(m6_out,c("ID_vs_non_ID"))
# post_plot(m6_out,c("non_ID"))
# post_plot(m6_out,c("ID"))
post_plot(m6_out,c("std_ID_v_non_ID"))

#summary of posterior
summary(m6_out)                           

#describe the posterior
describe_posterior(m6_out[,c(9:11)],         # [,c(#,#)] parameters of interest
                   centrality = "mean",
                   ci_method="hdi",           
                   rope_ci = 1,
                   rope_range = c(-0.1,0.1))    

#convert model object to matrix
m6_out_matrix <- as.matrix(m6_out[,c(9:11)]) # [,c(#,#)] parameters of interest
head(m6_out_matrix)

#plot probability of direction
plot(p_direction(m6_out_matrix[,3], ci=1))      # [,#] parameter of interest = difference score

#plot % eff size in rope
plot(rope(m6_out_matrix[,3],                    # [,#] parameter of interest = difference score
          ci=1,
          ci_method="hdi",
          range=c(-0.1,0.1)))           

###inverse log-transform natural log-transformed (ln) HRV data######
    
    # nonID_mean <- exp(summary(m6_out)[[1]][6,1])
    # nonID_95_CI_low <- exp(summary(m6_out)[[2]][6,1])
    # nonID_95_CI_high <- exp(summary(m6_out)[[2]][6,5])
    # nonID_mean
    # nonID_95_CI_low
    # nonID_95_CI_high
    # 
    # ID_mean <- exp(summary(m6_out)[[1]][2,1])
    # ID_95_CI_low <- exp(summary(m6_out)[[2]][2,1])
    # ID_95_CI_high <- exp(summary(m6_out)[[2]][2,5])
    # ID_mean
    # ID_95_CI_low
    # ID_95_CI_high
    # 
    # doesn't work with negative numbers!!!
    # diff_mean <- exp(summary(m6_out)[[1]][3,1])
    # diff_95_CI_low <- exp(summary(m6_out)[[2]][3,1])
    # diff_95_CI_high <- exp(summary(m6_out)[[2]][3,5])
    # diff_mean
    # diff_95_CI_low
    # diff_95_CI_high
    

#######################DIAGNOSTICS###############################
######Convergence check###########


  #trace plots (look for fuzzy caterpillar)
    color_scheme_set("mix-blue-red")
    mcmc_trace(m6_out,
               pars=c("ID","non_ID","std_ID_v_non_ID"),
               facet_args=list(ncol=1,strip.position="left"))

#####Posterior predictive check of model fit#########

  #want mean Bayesian p-value around 0.5
    summary(m6_out)
    head(m6_out)

  #Look at min and max values to set correct limits
    m6_out_matrix <- as.matrix(m6_out)
    m6_out_matrix
    upr = round(max(c(m6_out_matrix[,4],m6_out_matrix[,5])), -1) + 10
    lwr = round(min(c(m6_out_matrix[,4],m6_out_matrix[,5])), -1) - 10

  #plot
    plot(m6_out_matrix[,4],m6_out_matrix[,5],
         xlim=c(lwr,upr),ylim=c(lwr,upr),
         xlab="T",ylab="T.pred",pch=".")
    lines(c(lwr,upr),c(lwr,upr))

#####Residuals#######

  #Subset data by group
    df_group <- function(data, treatment){
      df <- data %>% filter(Treatment==treatment)
      return(df)
    }
    df_clinical <- df_group(data_obj6, "clinical")
    df_healthy <- df_group(data_obj6, "none")

  #select mean inferred values
    mu <- function(row){
      mu <- summary(m6_out)[[1]][row,1]
      print(mu)
    }
    mu_clinical <- mu(2)
    mu_clinical
    mu_healthy <- mu(6)
    mu_healthy

    sd_clinical <- summary(m6_out)[[1]][7,1]
    sd_clinical
    sd_healthy <- summary(m6_out)[[1]][8,1]
    sd_healthy

  #Check assumption of constant variance of residuals
    res <- function(df, mu, group){
      res <- df[,3][[1]] - mu
      plot <- plot(seq(length(res)),
                   res,
                   ylim=c(-20, 20),
                   main=group,
                   ylab="residuals",
                   xlab="n") + abline(h=0,col="red")
      return(res)
    }
    res_clinical <- res(df_clinical, mu_clinical, "")
    res_healthy <- res(df_healthy, mu_healthy,"")

  #Check assumption of constant variance of standardized residuals
    res_std <- function(df, mu, sd, group){
      res_std <- (df[,3][[1]] - mu)/sd
      plot <- plot(seq(length(res_std)),
                   res_std,
                   ylim=c(-5, 5),
                   main=group,
                   ylab="standardized residuals",
                   xlab="n") + abline(h=0,col="red")
      return(res_std)
    }
    res_std_clinical <- res_std(df_clinical, mu_clinical, sd_clinical, "")
    res_std_healthy <- res_std(df_healthy, mu_healthy, sd_healthy,"")

  #Check assumption of normally distributed residuals
    qqnorm(res_std_clinical, main="");
    qqline(res_std_clinical)
    qqnorm(res_std_healthy, main="");
    qqline(res_std_healthy)

  #Compare std residuals plots to ensure assumption of equal variances met
    res_vert <- function(df, data){
      res <- ggplot(df,aes(x=Treatment,y=data)) +
        
        geom_point(alpha=.5) +
        ylab("standardized residuals") + ylim(-5,5)
    }
    res_vert_clinical <- res_vert(df_clinical, res_std_clinical)
    res_vert_healthy  <- res_vert(df_healthy, res_std_healthy)

    compare.plot <- res_vert_clinical + res_vert_healthy
    print(compare.plot)







