######################################
# SIMULATIONS FOR LINEAR REGRESSION
######################################
set.seed(1)

p = 40
## Simu 1
# split_n = list(c(80, 50, 70), c(50, 50, 50),  c(70, 100, 80), c(50, 100, 80))
## Simu 2
split_n = list(c(5, 20, 70), c(10, 40, 50),  c(10, 50, 70), c(40, 10, 80))

n_j = unlist(lapply(split_n, sum))

# generate dichotomous data
## change some of the variables
X = list()
for(i in 1:length(n_j))
{
  x = matrix(0, n_j[i], p)
  for(j in 1:p)
  {
    for(k in 1:n_j[i])
    {
       # x[, j] = rbinom(n_j[i], 1, 0.5)
      x[k, j] = rnorm(1) 
    }
  }
  X[[i]] = x
}

# ## coefficients for groups - hyperparamters -- works kinda
vars = list(var1 = 1:10, var2 = 6:13, var3 = 11:20, var4 = c(1:5, 15:19))


beta_1 = jitter(rep(0, p)); beta_1[vars[[1]]] = c(rep(2, 5), rep(-2, 5))
beta_2 = jitter(rep(0, p)); beta_2[vars[[2]]] = c(rep(-2, 5), 1.8, 1.8, 1.8)
beta_3 = jitter(rep(0, p)); beta_3[vars[[3]]] = c(rep(-1.8, 5), rep(1.8, 5))
beta_4 = jitter(rep(0, p)); beta_4[vars[[4]]] = c(rep(2, 5), rep(1.8, 5))


beta_list = list(beta_1, beta_2, beta_3, beta_4)

## generate logist scaled probabilities
pi_list = list()
Y = list()
for(i in 1:length(n_j))
{
  pi_j = matrix(0, nrow = n_j[i], ncol = p)
  X_j = X[[i]]

  beta_vec = as.numeric(t(beta_list[[i]]))
  linear_pred =  X_j%*%beta_vec
  Y[[i]] = rnorm(n_j[i], X_j%*%beta_vec, 1)
}



coef_out = list()
for(i in 1:length(n_j))
{
  prova = lm(Y[[i]] ~., data = as.data.frame(X[[i]]))
  coef_out[[i]] = cbind(prova$coef, c(sum(beta_list[[i]]), beta_list[[i]]))
}
coef_out
##################################
# Random allocation to groups
Z_0 = c(1, 1, 1, 2, 2,2, 3, 3, 3, 4, 4 ,4)
X_list = list()
Y_list = list()
index = list()
###

### Random split
# split_n = list(c(100, 50, 50), c(50, 50, 100),  c(70, 100, 80), c(50, 100, 100))
lapply(split_n, sum)

## 1
d = 1
index_data = 1:dim(X[[d]])[1]
for(j in 1:length(split_n[[d]]))
{
  index = sample(index_data, split_n[[d]][j], replace = F)
  X_list[[j]] = X[[d]][index, ]
  Y_list[[j]] = Y[[d]][index]
  index_data = setdiff(index_data, index)
}
## 2
d = 2
index_data = 1:dim(X[[d]])[1]
for(j in 1:length(split_n[[d]]))
{
  index = sample(index_data, split_n[[d]][j], replace = F)
  X_list[[3 + j]] = X[[d]][index, ]
  Y_list[[3 + j]] = Y[[d]][index ]
  index_data = setdiff(index_data, index)
}
## 3
d = 3
index_data = 1:dim(X[[d]])[1]
for(j in 1:length(split_n[[d]]))
{
  index = sample(index_data, split_n[[d]][j], replace = F)
  X_list[[6 + j]] = X[[d]][index, ]
  Y_list[[6 + j]] = Y[[d]][index ]
  index_data = setdiff(index_data, index)
}
## 3
d = 4
index_data = 1:dim(X[[d]])[1]
for(j in 1:length(split_n[[d]]))
{
  index = sample(index_data, split_n[[d]][j], replace = F)
  X_list[[9 + j]] = X[[d]][index, ]
  Y_list[[9 + j]] = Y[[d]][index ]
  index_data = setdiff(index_data, index)
}

###################
# CHECK
coef_out = list()
for(i in 1:12)
{
  prova = lm(Y_list[[i]] ~., data = as.data.frame(X_list[[i]]))
  coef_out[[i]] = prova$coef
}
coef_out
beta_list

coef_pr = list()
for(i in 1:length(n_j))
{
  Y = do.call("c", Y_list[which(Z_0 == i)])
  X = do.call("rbind", X_list[which(Z_0 == i)])
  prova = lm(Y ~., data = as.data.frame(X))
  coef_pr[[i]] = cbind(prova$coef, c(sum(beta_list[[i]]), beta_list[[i]]))
}
coef_pr

save(Y_list, X_list, Z_0, beta_list, file = "simulation_continuos.RData")