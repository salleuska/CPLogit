######################################
# SIMULATIONS FOR LOGISTIC BORROWING
######################################
# MORE CHALLENGING SETTING
set.seed(45)

p = 10
# n_j = c(1000, 1200, 1100, 1200)
n_j = c(1000, 1200, 1100, 1200)

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
      x[k, j] = ifelse(runif(1) > 0.5, 0, 1)
      # x[k, j] = ifelse(rnorm(1) >0, 0, 1)
    }
  }
  X[[i]] = x
}

# ## coefficients for groups - hyperparamters -- works kinda
vars = list(var1 = c(1,2,3,4), var2 = c(4,5,6), var3 = 9:10, var4 = c(1,2,9,10))
beta_1 = rep(0, p); beta_1[vars[[1]]] = log(c(1.5, 0.5,1.3, 1.3))
beta_2 = rep(0, p); beta_2[vars[[2]]] = log(c(1.5, 0.7, 1.5))
beta_3 = rep(0, p); beta_3[vars[[3]]] = log(c(1.5, 0.5))
beta_4 = rep(0, p); beta_4[vars[[4]]] = log(c(1.5,0.5, 1.5, 0.5))
#beta_4 = rep(0, p); beta_4[vars[[4]]] = log(c(0.8, 1.2, 1.5, 1.2))

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
  pi_list[[i]] = exp(linear_pred)/(1 + exp(linear_pred))
  Y[[i]] = rbinom(n_j[i], 1, pi_list[[i]])
}
lapply(pi_list, mean)

X[[1]][, c(3,9)] = X[[1]][,c(9,3)]
X[[1]][, 9] = gtools::permute(X[[1]][, 5])
##
X[[2]][, 8] = gtools::permute(X[[2]][, 9])


coef_out = list()
for(i in 1:length(n_j))
{
  prova = glm(Y[[i]] ~., data = as.data.frame(X[[i]]), family = "binomial" )
  coef_out[[i]] = cbind(prova$coef, c(sum(beta_list[[i]]), beta_list[[i]]))
}
lapply(coef_out, exp)

##################################
# Random allocation to groups
Z_0 = c(1, 1, 1, 2, 2,2, 3, 3, 3, 4, 4 ,4)
X_list = list()
Y_list = list()
index = list()
###
### Random split
split_n = list(c(100, 600, 200), c(300, 100, 100),  c(500,100, 200), c(200, 200, 200))

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
  prova = glm(Y_list[[i]] ~., data = as.data.frame(X_list[[i]]), family = "binomial" )
  coef_out[[i]] = prova$coef
}
lapply(coef_out, function(x) exp(x[-1]))
lapply(beta_list, function(x) exp(x))

coef_pr = list()
for(i in 1:length(n_j))
{
  Y = do.call("c", Y_list[which(Z_0 == i)])
  X = do.call("rbind", X_list[which(Z_0 == i)])
  prova = glm(Y ~., data = as.data.frame(X), family = "binomial" )
  coef_pr[[i]] = cbind(prova$coef, c(sum(beta_list[[i]]), beta_list[[i]]))
}
lapply(coef_pr, exp)

save(Y_list, X_list, Z_0, beta_list, file = "simulation_discrete.RData")
