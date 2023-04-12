#######################################
## Centered Partition Processes: Informative Priors for Clustering
## Sally Paganin
## January 2021
#######################################
library(CPLogit)


## Function to get the Maximum at Posteriori classification
getClassMAP = function(clustering_mcmc, burn, tab = FALSE){
  n_iter0 = dim(clustering_mcmc)[2]
  class_tmp = apply(clustering_mcmc[, (burn+1):n_iter0], 1, function(x) table(factor(x, levels = 0:(H-1))))
  class = apply(t(class_tmp), 1, function(x) names(x)[which.max(x)])
  if(tab){ return(list(class = class, class_tmp = class_tmp))} else { return(class) }
}
#######################################
## load data
#load("simulation_discrete.RData")
load("simulation_continuos.RData")

Z_0 = Z_0 -1

###################################################
## 0) Separate LM
coef_out = list(); mod = list();
coef_lm = list()
for(i in 1:length(Y_list))
{
  model = lm(Y_list[[i]] ~., data = as.data.frame(X_list[[i]]))
  mod[[i]] = model
  coef_out[[i]] = model$coef

  df_coef = data.frame(cbind(coef_out[[i]][-1], confint.default(mod[[i]])[-1,]))
  colnames(df_coef) = paste0(c("mean_", "low_", "up_"), i)

  coef_lm[[i]] = df_coef
}

###################################################
## 1) Estimate a linear regression using the true grouping 
## to use as comparison baseline
###################################################

coef_grouped = list(); 
g0 = Z_0 +1
for(i in 1:length(unique(g0)))
{
  sel = which(g0 == i)
  Y = do.call(c, Y_list[sel])
  X = do.call(rbind, X_list[sel])

  model = lm(Y ~., data = as.data.frame(X))

  df_coef = data.frame(cbind(model$coef, confint.default(model)))
  colnames(df_coef) = c("mean", "low", "up")
  coef_grouped[[i]] = df_coef[-1, ]
}

###################################################
## 2) Bayesian modeling - DP process (no prior information on clustering)
###################################################
## number of coefficients
p = dim(X_list[[1]])[2]

## Upper bound for the number of clusters
H_upper = 12
## number of clusters at initialization
H = 7
## Dirichlet Process concentration parameters
alpha = 1
## Priors for regressions coefficients 
## intercept ~ N(a_prior, tau_prior = sigma^(-2))
a_prior = 0;
tau_prior = 1;

## beta_coeff ~ MVN (b_prior, Q_prior)

b_prior = array(0, dim = p)
Q_prior = diag(array(2, dim = p))
##################################
# Run gibbs sampling using DP process
n_iter = 5000
psi_par = 0    ## Equivalent to DP process 

set.seed(1)
res = gibbsLinearCP(X_list, Y_list, Z_0, psi_par, H, H_upper ,alpha,  a_prior, tau_prior,  b_prior, Q_prior,  n_iter)
burn = 1000
# check classification
cl = getClassMAP(res$clustering, burn, tab = TRUE)
table(cl$class, Z_0)

### Coefficient estimates

data_coef = array(0, dim = c(length(Y_list), p,  n_iter))
for(t in 1:n_iter)
{
  data_coef[,, t] = res$beta[res$clustering[, t] + 1, ,t]
}

odds_logit = apply(data_coef[,,burn:n_iter], c(1,2), mean)
odds_logit_low = apply(data_coef[,,burn:n_iter], c(1,2), function(x) quantile(x, 0.025))
odds_logit_up = apply(data_coef[,,burn:n_iter], c(1,2), function(x) quantile(x, 0.975))

# Estimates
dp_res = list(mean = odds_logit, lower = odds_logit_low, upper = odds_logit_up)

# Pariwise allocation matrix
pair_mat_dp = mcclust::comp.psm(t(res$clustering[, burn:n_iter] + 1))

###################################################
## 3) Bayesian modeling - CP process, with baseline DP EPPF
## using  prior information on clustering with psi = 17
###################################################s
n_iter = 5000
psi_par = 17
set.seed(4)

res = gibbsLinearCP(X_list, Y_list, Z_0, psi_par, H , H_upper, alpha,  a_prior, tau_prior,  b_prior, Q_prior,  n_iter)

burn = 1000
# check classification
cl = getClassMAP(res$clustering, burn, tab = TRUE)

data_coef = array(0, dim = c(length(Y_list), p,  n_iter))
for(t in 1:n_iter)
{
  data_coef[,, t] = res$beta[res$clustering[, t] + 1, ,t]
}

odds_logit = apply(data_coef[,,burn:n_iter], c(1,2), mean)
odds_logit_low = apply(data_coef[,,burn:n_iter], c(1,2), function(x) quantile(x, 0.025))
odds_logit_up = apply(data_coef[,,burn:n_iter], c(1,2), function(x) quantile(x, 0.975))

cp_res = list(mean = odds_logit, lower = odds_logit_low, upper = odds_logit_up)

# Pairwise allocation matrix
pair_mat_cp = mcclust::comp.psm(t(res$clustering[, burn:n_iter] + 1))

############################
## using  prior information on clustering with psi = 17
############################
n_iter = 5000
psi_par = 15
set.seed(88)
res = gibbsLinearCP(X_list, Y_list, Z_0, psi_par, H , H_upper, alpha,  a_prior, tau_prior,  b_prior, Q_prior,  n_iter)

burn = 1000
# check classification
cl = getClassMAP(res$clustering, burn, tab = TRUE)
table(cl$class, Z_0)

data_coef = array(0, dim = c(length(Y_list), p,  n_iter))
for(t in 1:n_iter)
{
  data_coef[,, t] = res$beta[res$clustering[, t] + 1, ,t]
}

odds_logit = apply(data_coef[,,burn:n_iter], c(1,2), mean)
odds_logit_low = apply(data_coef[,,burn:n_iter], c(1,2), function(x) quantile(x, 0.025))
odds_logit_up = apply(data_coef[,,burn:n_iter], c(1,2), function(x) quantile(x, 0.975))

cp_res2 = list(mean = odds_logit, lower = odds_logit_low, upper = odds_logit_up)

# Pairwise allocation matrix
pair_mat_cp2 = mcclust::comp.psm(t(res$clustering[, burn:n_iter] + 1))

###################################################
## 4) Producing plots in the Supplementary Material
###################################################s
library(ggplot2)
df_clust = reshape2::melt(pair_mat_dp)

df_clust$Var1 <- factor(df_clust$Var1, levels = 1:12)
df_clust$Var2 <- factor(df_clust$Var2, levels = 1:12)

p_dp =  ggplot(df_clust)+ geom_tile(aes(Var1,Var2,fill= value)) + scale_fill_gradient2(low= "white", high = "grey50") +
 theme_minimal() + theme(axis.title.x=element_blank(),axis.title.y=element_blank()) + theme(legend.position="none") +
 labs(title = paste("Dirichlet process"))+ theme(plot.title = element_text(hjust = 0.5)) +
 theme(text = element_text(size = 15))

p_dp


df_clust = reshape2::melt(pair_mat_cp)

df_clust$Var1 <- factor(df_clust$Var1, levels = 1:12)
df_clust$Var2 <- factor(df_clust$Var2, levels = 1:12)

p_cp =  ggplot(df_clust)+ geom_tile(aes(Var1,Var2,fill= value)) + scale_fill_gradient2(low= "white", high = "grey50") +
 theme_minimal() + theme(axis.title.x=element_blank(),axis.title.y=element_blank()) + theme(legend.position="none") +
 labs(title = paste("Centered Partition (\u03C8 = 17)"))+ theme(plot.title = element_text(hjust = 0.5))+
 theme(text = element_text(size = 15))
p_cp


df_clust = reshape2::melt(pair_mat_cp2)

df_clust$Var1 <- factor(df_clust$Var1, levels = 1:12)
df_clust$Var2 <- factor(df_clust$Var2, levels = 1:12)

p_cp2 =  ggplot(df_clust)+ geom_tile(aes(Var1,Var2,fill= value)) + scale_fill_gradient2(low= "white", high = "grey50") +
 theme_minimal() + theme(axis.title.x=element_blank(),axis.title.y=element_blank()) + theme(legend.position="none") +
 labs(title = paste("Centered Partition (\u03C8 = 15)"))+ theme(plot.title = element_text(hjust = 0.5))+
 theme(text = element_text(size = 15))
p_cp2
##################################
tt = coef_grouped[Z_0 + 1]
lm_res = sapply(tt, function(x) x$mean )

# dim(lm_res)
# dim(dp_res$mean)

dp_diff = as.data.frame(dp_res$mean - t(lm_res))
cp_diff = as.data.frame(cp_res$mean - t(lm_res))
cp_diff2 = as.data.frame(cp_res2$mean - t(lm_res))

dp_diff$id = 1:12
cp_diff$id = 1:12
cp_diff2$id = 1:12

df = reshape2::melt(dp_diff, id.vars = "id")

summary(df$value)

p1 = ggplot(df, aes(x = id, y = value, group = id)) +
geom_boxplot(fill = "grey90") +
scale_x_continuous(breaks = 1:12, labels = 1:12) + scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2, 2)) +
xlab("") + ylab("") + theme_minimal() +
theme(text = element_text(size = 20), strip.text.x = element_text(size = 20 ))
p1

df_cp = reshape2::melt(cp_diff, id.vars = "id")

summary(df_cp$value)

p2 = ggplot(df_cp, aes(x = id, y = value, group = id)) +
geom_boxplot(fill = "grey90") +  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2, 2)) +
scale_x_continuous(breaks = 1:12, labels = 1:12)  + 
xlab("") + ylab("") + theme_minimal() +
theme(text = element_text(size = 20), strip.text.x = element_text(size = 20 ))

p2

df_cp2 = reshape2::melt(cp_diff2, id.vars = "id")

p3 = ggplot(df_cp2, aes(x = id, y = value, group = id)) +
geom_boxplot(fill = "grey90") +
scale_x_continuous(breaks = 1:12, labels = 1:12) +  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2, 2)) +
xlab("") + ylab("") + theme_minimal() +
theme(text = element_text(size = 20), strip.text.x = element_text(size = 20 ))

p3

library(cowplot)
# p_cp, p2, 
pr = plot_grid(p_dp, p1, p_cp2, p3, ncol = 2, nrow = 2, rel_widths = c(0.5, 1))
pr
save_plot("SP_Fig1_Simulation_truec0.pdf", pr, ncol = 2, nrow = 2, device = cairo_pdf, base_aspect_ratio = 1, base_width =8, base_height =5)


###################################################
## 5) Repeat with wrong prior elicitation
###################################################s

cent = c(0,1,2,3,0,1,2,3,0,1,2,3)
vi_distC(cent, Z_0)
vi_distC(cent, rep(0,12))
table(Z_0, cent)
psi_par = 15
set.seed(2)


res = gibbsLinearCP(X_list, Y_list, cent, psi_par, H, H_upper  ,alpha,  a_prior, tau_prior,  b_prior, Q_prior,  n_iter)


burn =1000
# check classification
cl = getClassMAP(res$clustering, burn, tab = TRUE)
table(cl$class)
table(cl$class, Z_0)
table(cl$class, cent)

dp_clust = c(0,0,0, rep(1, 9))
vi_distC(as.numeric(cl$class), Z_0)
vi_distC(as.numeric(cl$class), dp_clust)
vi_distC(as.numeric(cl$class), cent)

# Pairwise allocation matrix
pair_mat_cpw = mcclust::comp.psm(t(res$clustering[, burn:n_iter] + 1))


df_clust = reshape2::melt(pair_mat_cpw)

df_clust$Var1 <- factor(df_clust$Var1, levels = 1:12)
df_clust$Var2 <- factor(df_clust$Var2, levels = 1:12)

p_w =  ggplot(df_clust)+ geom_tile(aes(Var1,Var2,fill= value)) + scale_fill_gradient2(low= "white", high = "grey50") +
 theme_minimal() + theme(axis.title.x=element_blank(),axis.title.y=element_blank()) + theme(legend.position="none") +
 labs(title = paste("Centered Partition (wrong guess)"))+ theme(plot.title = element_text(hjust = 0.5))+
 theme(text = element_text(size = 20))

p_w

###################################
data_coef = array(0, dim = c(length(Y_list), p,  n_iter))
for(t in 1:n_iter)
{
  data_coef[,, t] = res$beta[res$clustering[, t] + 1, ,t]
}

odds_logit = apply(data_coef[,,burn:n_iter], c(1,2), mean)
odds_logit_low = apply(data_coef[,,burn:n_iter], c(1,2), function(x) quantile(x, 0.025))
odds_logit_up = apply(data_coef[,,burn:n_iter], c(1,2), function(x) quantile(x, 0.975))

cpw_res = list(mean = odds_logit, lower = odds_logit_low, upper = odds_logit_up)

cpw_diff = as.data.frame(cpw_res$mean - t(lm_res))

cpw_diff$id = 1:12

df = reshape2::melt(cpw_diff, id.vars = "id")

summary(df$value)

pw = ggplot(df, aes(x = id, y = value, group = id)) +
geom_boxplot(fill = "grey90") +
scale_x_continuous(breaks = 1:12, labels = 1:12) +scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2, 2)) +
xlab("") + ylab("") + theme_minimal() +
theme(text = element_text(size = 20), strip.text.x = element_text(size = 20 ))

pw
###########
pr = plot_grid(p_w, pw, ncol = 2, nrow = 1, rel_widths = c(0.5, 1))
pr
save_plot("SP_Fig2_Simulation_wrongc0.pdf", pr, ncol = 2, nrow = 1, device = cairo_pdf, base_aspect_ratio = 1.5, base_height = 6)