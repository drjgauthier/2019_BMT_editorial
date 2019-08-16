library(dplyr)
library(rms)
theme_set(theme_bw())

sim <- function(){
  #Simulate biomarker values
  biomarker <- round(runif(200,0,1000),digits=0)
  biomarker_m <- cut2(biomarker,g=2) #cut across median
  biomarker_q <- cut2(biomarker,g=4) #quartiles
  quart<-data.frame(quantile(biomarker,probs = c(0.05,0.25,0.5,0.75,0.95))) #Get qnatiles of biomarker
  truelogit <- sin(biomarker/150)
  prob <- exp(truelogit)/(1 + exp(truelogit))
  runis <- runif(200,0,1)
  sim_y <- ifelse(runis < prob,1,0)
  df_logit<-data.frame(biomarker,truelogit)
  df = data.frame(sim_y=sim_y,biomarker=biomarker,biomarker_m=biomarker_m,
                biomarker_t=biomarker_t)
  m<-glm(sim_y ~ biomarker_m,data=df,family=binomial)
  q<-glm(sim_y ~ biomarker_q,data=df,family=binomial)
  s<-lrm(sim_y ~ rms::rcs(biomarker,c(quart[1,1],quart[2,1],quart[3,1],quart[4,1],quart[5,1])),data=df)
  r <- (quart[5,1]-quart[1,1])^2
  b2_s <- s$coefficients[3]/r
  b3_s <- s$coefficients[4]/r
  b4_s <- s$coefficients[5]/r
  b5_s <- (b2_s*(quart[1,1]-quart[5,1])+
             b3_s*(quart[2,1]-quart[5,1])+
             b4_s*(quart[3,1]-quart[5,1]))/(quart[5,1]-quart[4,1])
  b6_s <- (b2_s*(quart[1,1]-quart[4,1])+
             b3_s*(quart[2,1]-quart[4,1])+
             b4_s*(quart[3,1]-quart[4,1]))/(quart[4,1]-quart[5,1])
  
  l<-list(df_logit,m$coefficients,q$coefficients,s$coefficients,b2_s,b3_s,b4_s,b5_s,b6_s)
return(l)
}

set.seed(1950)
biomarker <- round(runif(200,0,1000),digits=0)
biomarker_m <- cut2(biomarker,g=2) #median
biomarker_t <- cut2(biomarker,g=3) #tertiles
biomarker_q <- cut2(biomarker,g=4) #quartiles
quart<-data.frame(quantile(biomarker,probs = c(0.05,0.25,0.5,0.75,0.95)))
truelogit <- sin(biomarker/150)
trueprob <- exp(truelogit)/(1 + exp(truelogit))
df <- data.frame(biomarker,biomarker_m,biomarker_q,trueprob)
sims <- replicate(1000,sim(),simplify = FALSE)
str(sims)

#Average predicted intercepts and betas using f_m models
purrr::map(sims,c(2,1)) %>% bind_cols() %>% mutate(mean_int_m=rowMeans(.)) %>% select(mean_int_m) -> mean_int_m
purrr::map(sims,c(2,2)) %>% bind_cols() %>% mutate(mean_beta_m=rowMeans(.)) %>% select(mean_beta_m) -> mean_beta_m
#Get pred logits from averaged betas from log model using median x
pred_logit0_m <- mean_int_m + 0*mean_beta_m
pred_prob0_m <- exp(pred_logit0_m)/(1 + exp(pred_logit0_m))
pred_logit1_m <- mean_int_m + 1*mean_beta_m
pred_prob1_m <- exp(pred_logit1_m)/(1 + exp(pred_logit1_m))

#df_m$pred_logit <- ifelse(biomarker>=median(biomarker))
df %>% mutate(pred_prob=ifelse(biomarker_m==biomarker_m[1],pred_prob0_m,pred_prob1_m)) -> df
df <- data.frame(df,unlist(df$pred_prob))
df <- df %>% select(-5) %>% rename(pred_prob_m=unlist.df.pred_prob.)

#Graph predicted against simulated probabilities
df %>% select(biomarker,trueprob,pred_prob_m) %>% tidyr::gather(key=prob_type,value=prediction,-biomarker) -> dfm_g
m_p <- ggplot(dfm_g,aes(x=biomarker,y=prediction,col=prob_type)) + 
  geom_line(aes(col=prob_type))+
  labs(x="Simulated fictinin-1 value",y="Probability of aGvHD")+
  scale_color_manual(values=c("red","blue"),labels=c("Predicted probability","Assumed-true probability"))+
  theme(legend.title = element_blank(),legend.text = element_text(size=15),axis.title.x = element_text(size=15),axis.text.x=element_text(size=12),axis.title.y = element_text(size=15),axis.text.y=element_text(size=12))+
  ggtitle("Fictinin-1 dichotomised across median")

#Average predicted intercepts and betas using f_m models
purrr::map(sims,c(3,1)) %>% bind_cols() %>% mutate(mean_int_q=rowMeans(.)) %>% select(mean_int_q) -> mean_int_q 
purrr::map(sims,c(3,2)) %>% bind_cols() %>% mutate(mean_beta1_q=rowMeans(.)) %>% select(mean_beta1_q) -> mean_beta1_q 
purrr::map(sims,c(3,3)) %>% bind_cols() %>% mutate(mean_beta2_q=rowMeans(.)) %>% select(mean_beta2_q) -> mean_beta2_q 
purrr::map(sims,c(3,4)) %>% bind_cols() %>% mutate(mean_beta3_q=rowMeans(.)) %>% select(mean_beta3_q) -> mean_beta3_q 

#Get pred logits from averaged intercepts and betas from log model using quartiles of x
logit1_q <- mean_int_q
prob1_q <- exp(logit1_q)/(1 + exp(logit1_q))
logit2_q <- mean_int_q+1*mean_beta1_q
prob2_q <- exp(logit2_q)/(1 + exp(logit2_q))
logit3_q <- mean_int_q+1*mean_beta2_q
prob3_q <- exp(logit3_q)/(1 + exp(logit3_q))
logit4_q <- mean_int_q+1*mean_beta3_q
prob4_q <- exp(logit4_q)/(1 + exp(logit4_q))

df %>% mutate(pred_prob_q=case_when(
  biomarker_q == "[  2,250)"  ~ 0.6600773,
  biomarker_q == "[250,489)" ~ 0.6238583,
  biomarker_q == "[489,741)" ~ 0.3216825,
  biomarker_q == "[741,999]" ~ 0.4045998,
)) -> df

#Graph predicted against simulated probabilities
df %>% select(biomarker,trueprob,pred_prob_q) %>% tidyr::gather(key=prob_type,value=prediction,-biomarker) -> dfq_g
q_p <- ggplot(dfq_g,aes(x=biomarker,y=prediction,col=prob_type)) + 
  geom_line(aes(col=prob_type))+
  labs(x="Simulated fictinin-1 value",y="Probability of aGvHD")+
  scale_color_manual(values=c("red","blue"),labels=c("Predicted probability","Predicted probability"))+
  theme(legend.title = element_blank(),legend.text = element_text(size=15),axis.title.x = element_text(size=15),axis.text.x=element_text(size=12),axis.title.y = element_text(size=15),axis.text.y=element_text(size=12))+
  ggtitle("Fictinin-1 categorised across quartiles")

#Average predicted intercepts and betas using f_s (splined) models
purrr::map(sims,c(4,1)) %>% bind_cols() %>% mutate(mean_int_s=rowMeans(.)) %>% select(mean_int_s) -> mean_int_s # 1.50
purrr::map(sims,c(4,2)) %>% bind_cols() %>% mutate(mean_b1_s=rowMeans(.)) %>% select(mean_b1_s) -> mean_b1_s # -0.00209
purrr::map(sims,c(5,1)) %>% bind_cols() %>% mutate(mean_b2_s=rowMeans(.)) %>% select(mean_b2_s) -> mean_b2_s # -0.0126
purrr::map(sims,c(6,1)) %>% bind_cols() %>% mutate(mean_b3_s=rowMeans(.)) %>% select(mean_b3_s) -> mean_b3_s # 0.0297
purrr::map(sims,c(7,1)) %>% bind_cols() %>% mutate(mean_b4_s=rowMeans(.)) %>% select(mean_b4_s) -> mean_b4_s # 0.0125
purrr::map(sims,c(8,1)) %>% bind_cols() %>% mutate(mean_b5_s=rowMeans(.)) %>% select(mean_b5_s) -> mean_b5_s # 0.0125
purrr::map(sims,c(9,1)) %>% bind_cols() %>% mutate(mean_b6_s=rowMeans(.)) %>% select(mean_b6_s) -> mean_b6_s # 0.0125

coef_s<-data.frame(mean_int_s,mean_b1_s,mean_b2_s,mean_b3_s,mean_b4_s,mean_b5_s,mean_b6_s)
df_s <- data.frame(df,coef_s)

#Compute predicted logits for each biomarker value
df_s$pred_logit_s <- df_s$mean_int_s+
  df_s$mean_b1_s*df_s$biomarker+
  df_s$mean_b2_s*pmax(df_s$biomarker-quart[1,1],0)^3+
  df_s$mean_b3_s*pmax(df_s$biomarker-quart[2,1],0)^3+
  df_s$mean_b4_s*pmax(df_s$biomarker-quart[3,1],0)^3+
  df_s$mean_b5_s*pmax(df_s$biomarker-quart[4,1],0)^3+
  df_s$mean_b6_s*pmax(df_s$biomarker-quart[5,1],0)^3 
df_s$pred_prob_s <- exp(df_s$pred_logit_s)/(1 + exp(df_s$pred_logit_s))
df_s %>% select(biomarker,trueprob,pred_prob_s) %>% tidyr::gather(key=prob_type,value=prediction,-biomarker) -> dfs_g

#Graph predicted against simulated probabilities
s_p <- ggplot(dfs_g,aes(x=biomarker,y=prediction,col=prob_type)) + 
  geom_line(aes(col=prob_type))+
  labs(x="Simulated fictinin-1 value",y="Probability of aGvHD")+
  scale_color_manual(values=c("red","blue"),labels=c("Predicted probability","Assumed-true probability"))+
  theme(legend.title = element_blank(),legend.text = element_text(size=15),axis.title.x = element_text(size=15),axis.text.x=element_text(size=12),axis.title.y = element_text(size=15),axis.text.y=element_text(size=12))+
  ggtitle("Modelling fictinin-1 with a restricted cubic spline")
s_p

fig2 <- ggpubr::ggarrange(m_p,q_p,s_p,
                  ncol = 3, nrow = 1,labels="AUTO",common.legend =TRUE) #1500x500
                  
annotate_figure(fig2,top = text_grob(""),fig.lab="Figure 2. Predictions of the probability of aGvHD using logistic regression (simulated dataset)",fig.lab.size=16,fig.lab.face = "bold",fig.lab.pos = c("top.left")) 
