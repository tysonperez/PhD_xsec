#Objective 6b: 
#To assess whether greater changes in EEG metrics 
#are associated with greater changes in PROs.  

#Load libraries
library(rjags)  
library(readxl)   
library(tidyverse)
library(bayesplot)
library(bayestestR)
library(see)
library(patchwork)
library(rstanarm)

##CLINICAL##

# #Load & visualise clinical data
# clinical_data <- read_excel("Clinical.xlsx", sheet="Obj6")
# clinical_data
# 
# #Gather relevant data (baseline 2=t1) & add variables (tx)
# get_data_obj6 <- function(data){
#   #update last number in data[c(...)] to reflect column of interest
#   tmp = data[c(1,2,7,178)]
#   tmp$dv = tmp[,3][[1]]
#   tmp$iv = tmp[,4][[1]]
#   tmp$tx <- 1*(tmp$Group=="ID")
#   return(tmp)
# }
# data_obj6 <- get_data_obj6(clinical_data)
# data_obj6

##HEALTHY##
 
#Load & visualise healthy data
healthy_data <- read_excel("Healthy.xlsx",sheet="Obj6")
healthy_data

#Gather relevant data (baseline 2=t1) & add variables (tx)
get_data_obj6 <- function(data){
   #update last number in data[c(...)] to reflect column of interest
   tmp = data[c(1,2,7,163)]
   tmp$dv = tmp[,3][[1]]
   tmp$iv = tmp[,4][[1]]
   tmp$tx = 2*(tmp$Group=="non-ID")
   return(tmp)
 }
 data_obj6 <- get_data_obj6(healthy_data)
 data_obj6


######Regression plot of 500 random draws from posterior distribution using non-standardized variables######

#regress model
reg6 <- stan_glm(dv ~ iv,
                     data=data_obj6,
                     chains=3, iter=35000, warmup=10000)
summary(reg6)

#coerce model to data-frame & remove sigma
reg6_df <- reg6 %>% as_tibble() %>% rename(intercept=`(Intercept)`) %>% select(-sigma)
head(reg6_df)

#aesthetics
n_draws <- 500
alpha_level <- .15
color_draw <- "grey60"
color_mean <- "#3366FF"  

#plot
ggplot(data_obj6) +
  aes(x=iv,y=dv) +
  coord_cartesian(ylim=c(-2,5)) +                 #restrict y-axis to focus on centre of data
  geom_abline(aes(intercept=intercept,slope=iv),  #plot random sample rows from sim df as grey,semi-transparent lines
              data=sample_n(reg6_df,n_draws),
              color=color_draw,
              alpha=alpha_level) +
  geom_abline(intercept=mean(reg6_df$intercept),  #plot mean of intercept & slope coefficients in blue
              slope=mean(reg6_df$iv),
              size=1,
              color=color_mean) +
  geom_point(aes(color=Group)) + 
  labs(x="right dlPFC beta (12.5-30 Hz) log-CSD",
       y="HADS-D")
       #title="Visualization of 500 Regression Lines from the Posterior Distribution")

######Standardized regression coefficient######

#standardize the variables
data_obj6b_std <- data.frame(scale(data_obj6[,5:6]))
head(data_obj6b_std)

#regress model
reg6_std <- stan_glm(dv ~ iv, data=data_obj6b_std,
                     chains=3, iter=35000, warmup=10000)
summary(reg6_std)

#extract std reg coefficient from model object
b1 <- as.matrix(reg6_std)
head(b1)

#plot posterior distribution of std slope
post_plot <- function(data,parameter){
  color_scheme_set("brightblue")
  mcmc_areas(data, 
             prob=0.95, 
             point_est="mean",
             pars=parameter) +
    labs(title="Posterior distributions", 
         subtitle="with means & 95% credible intervals")
}
post_plot(b1,c("iv"))

#describe the posterior of std slope
describe_posterior(reg6_std,
                   centrality = "mean",
                   ci_method="hdi",
                   rope_ci = 1,
                   rope_range = c(-0.05,0.05))

#plot rope of std slope
plot(rope(reg6_std, ci=1))

#####Diagnostics#######

summary(reg6_std)           #Rhat & effective sample size
pp_check(reg6_std, "stat")  #posterior predictive distribution


