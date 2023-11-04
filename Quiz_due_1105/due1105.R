#### packages
library(tidyverse)
library(GGally)
library(xtable)

#### Simulation settings
set.seed(88)
N <- 100
e <- rnorm(N) 
e2 <- rnorm(N)
x1 <- rnorm(N, sd=2)
x2 <- rnorm(N, sd=4)
x3 <- rnorm(N, sd=2) #corr w/ x4
x4 <- .8*x3 + rnorm(N) # corr w/ x3
x5 <- runif(N, -3, 3)
x6 <- runif(N, -3, 3)
x7 <- .8*x6 + rnorm(N) # corr w/ x6

### Case1
y <- 5 + 1.5*x1 + 1.5*x2+ 1.5*x3 + 0.9*x4 + 1.5*x5 + 1.5*x6 + 
  0.9*x7  + e

### Case2
y2 <- 3 +  1.5*x1 + 1.5*x3 + 1.5*x4  + 
  1.5*x6 + e2

dta <- data.frame(e, e2, x1, x2, x3, x4, x5, x6, x7, y, y2)


#### Descriptive plotting: (Fig 1)
dta %>%
  select(c(starts_with("x"),"y")) %>%
  ggpairs() +
  ggtitle("Case1: Full model is true model")

dta %>%
  select(c(starts_with("x"),"y2")) %>%
  ggpairs() + 
  ggtitle("Case2: x2, x5, x7 are irrelevant vars")


## Case 1

##### Functions
# column names
cnames <- c("full_model", paste0("drop_x", 1:7)) 
rnames <- c("Intercept", paste0("x", 1:7))

# true beta values
beta <- c(5, rep(1.5,3),.9, rep(1.5,2),.9) 
names(beta) <-  rnames


# dataframe to design matrix
dmtx <-  function(df, drop_y = TRUE){
  if(drop_y){
    .d <- df[-1]
  } else{.d <- df}
  mtx <- .d %>% 
    mutate(intercept = rep(1, nrow(df))) %>%
    relocate(intercept) %>% 
    as.matrix()
  return(mtx)
}
# inv(XTX)
xtxinv <- function(df, drop_y =TRUE){
  mtx <- dmtx(df, drop_y = drop_y)
  return(
    solve(t(mtx)%*%mtx)
  )
}
# XTX
xtx <- function(df, drop_y =TRUE){
  mtx <- dmtx(df, drop_y = drop_y)
  return(
    t(mtx)%*%mtx
  )
}


f <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7
fm <- lm(f, dta) #full model

# Create a list of full model and submodels that dropped one variable
mlist <- list(full_model = fm)
for (i in 1:7) {
  f.drop.i <- as.formula(
    paste("y ~", paste(paste0("x", 1:7)[-i], collapse = " + "))
  )
  mlist[[paste0("drop_x", i)]] <- lm(f.drop.i, dta)
}

### bias and variance
df.coeff <- tibble(coef = rnames)
df.bias <- tibble(coef = rnames)
df.var <- tibble(coef = rnames)
.ivs <- paste0("x", 1:7)
for (i in 1:8) {
  .mddf <- mlist[[i]]$model 
  xtxI <- .mddf  %>% xtxinv()
  #   print(i)
  ## variance
  if (i==1){
    mdf <- dta[.ivs] 
  }else {mdf <- dta[.ivs[-i]]}
  .var <- xtxinv(mdf, drop_y = F) %>%
    diag()

  ## bias
  if (i ==1){ #full model
    .xtx <- .mddf  %>% xtx()
    .coeff <- xtxI %*% .xtx %*% beta
    .bias <- .coeff - beta

    
  }else if (i!=1 & i< length(beta)){ #1-7
    x1t.mtx <- .mddf %>%
      dmtx %>% t() #X1t
    x2.mtx <- dta[paste0("x",i-1)] %>%
      as.matrix() # missing vars(X2)
    .beta<- beta[i+1] %>% #
      as.matrix() 
    .b <- xtxI %*% x1t.mtx %*% x2.mtx %*% .beta
    # Inserting NA    
    .bias <- c(.b[1:(i-1)], NA, .b[i:length(.b)])  
    .coeff <-  beta + .bias
    .var <- c(.var[1:(i-1)], NA, .var[i:length(.var)])
  }else{ # 8
    x1t.mtx <- .mddf %>% dmtx %>% t() #X1t
    x2.mtx <- dta[paste0("x",i-1)] %>% 
      as.matrix() # missing vars
    .beta<- beta[i] %>% #missing beta
      as.matrix() 
    .b <- xtxI %*% x1t.mtx %*% x2.mtx %*% .beta
    # Inseting NA
    .bias <- c(.b[1:i-1],NA)
    .coeff <-  beta + .bias
    .var <- ifelse(i==1, .var ,c(.var, NA))
  }
  .cn <- cnames[i]
  df.coeff[.cn] <- .coeff
  df.bias[.cn] <- .bias
  df.var[.cn] <-.var
}

#### Tables
df.bias %>% 
  relocate(coef) %>%
  xtable(
          digits =3, label = "bias_a")

true_m <- matrix(rep(beta, 8), ncol =8)
.m <- df.bias[-1] %>% as.matrix()/true_m

.m %>% 
  as.data.frame() %>% 
  mutate(coef = rnames) %>% 
  relocate(coef)%>%
  xtable(caption = "",
         digits =2, label = "bias_f")


df.var %>% 
  relocate(coef) %>%
  xtable(caption = "",
          digits =3, label = "var_a")

### MSE
mtb <- df.bias[-1] %>% as.matrix()
.mtv <- df.var %>% as.matrix()
.m <- round(((.mtb^2)+.mtv), 3) #MSE

df.mse <- data.frame(.m) %>%
  mutate(coef = rnames) %>% 
  relocate(coef) 
#### Table
df.mse %>%
  xtable(
    caption = "MSE",
    digits = 3,
    label = "MSE_a"
  )

### EstVariance
EstV <- c()
for (i in 1:8){
  md <- mlist[[i]]
  .var <- md$residuals^2 %>% sum() / md$df.residual
  EstV <- c(EstV, .var) %>% round(3)
}
names(EstV) <- cnames
EstV %>% as.matrix() %>% 
  t() %>% 
  as.data.frame() %>% 
  xtable(caption = "Estimated Variance",
         digits = 3,
         label = "evar_a")

# Case 2
## true beta
beta_ir <-c(3, 1.5, 0, 1.5, 0, 0, 1.5, 0.9) 
names(beta_ir) <- rnames
##　Generate all models

irm_list <- list()
vec <- c(2, 5, 7)
# Generate all possible subsets
all_subsets <- unlist(lapply(1:length(vec),
                             function(n) combn(vec, n, simplify = FALSE)),
                              recursive = FALSE)

# Add the empty set
all_subsets <- c(all_subsets, list(c()))  #of irrivat

# add irri. vars in model
for (i in 1:length(all_subsets)) {
  .c <- all_subsets[[i]]
  .c <- c(c(1, 3, 4, 6) ,.c) %>% sort()
  .f <- as.formula(
    paste("y2 ~", paste(paste0("x", .c), collapse = " + "))
  )
  irm_list[[paste0("model", i)]] <- lm(.f, dta)
}



# function
fillna <- function(.v, .pos){
  #orignal vector, position
  .pos <-  vec[!(vec %in% .pos)] # complement, vec = c(2,5,7)
  if(length(.pos)!=0){
    for (i in 1:length(.pos)){
      # print(paste("len:",length(.v)))
      if (.pos[i]!= 7){
        .v <- c(.v[1:.pos[i] ], NA, .v[ (.pos[i]+1) : length(.v)])
      }else{.v <- c(.v[1:.pos[i]], NA )}
    }
  }
  return(.v)
}
## bias and variance
df.coef.irr <- tibble(coef=rnames)
df.bias.irr <- tibble(coef=rnames)
df.var.irr <- tibble(coef=rnames)
for (i in 1:length(irm_list)){
  .mddf <- irm_list[[i]]$model
  ## bias & beta_hat
  xtxI <- .mddf  %>% xtxinv() 
  .xtx <-  .mddf  %>% xtx() 
  .drop <- vec[!(vec %in% all_subsets[[i]])] + 1 # complement & its position

  if (length(.drop)>0){ #computing beta_hat
    .coef <- xtxI %*% .xtx %*% beta_ir[-.drop] 
  }else{ #all covariate
    .coef <- xtxI %*% .xtx %*% beta_ir
  }
  
  for (j in .drop){ # fill na
    if ( j == length(beta_ir) ){ # x7
      .coef <-  c(.coef[1:(j-1)], NA)
    }else{
      .coef <-  c(.coef[1:(j-1)], NA, .coef[j:length(.coef)])
    }
  }
  
  #variance
  var <- xtxinv(.mddf) %>%
    diag()
  .pos <- all_subsets[[i]]
  .vec <- fillna(var, .pos) # fill na

  df.var.irr[paste0("model", i)] = .vec 
  df.coef.irr[paste0("model", i)] = .coef 
  df.bias.irr[paste0("model", i)] = .coef - beta_ir
}

### Tables
df.bias.irr <- df.bias.irr %>% 
  relocate(coef, model8) 
colnames(df.bias.irr) <- c("Coef","TrueModel","addx2","addx5","addx7",
                   "addx2x5","addx2x7","addx5x7", "addx2x5x7")

df.bias.irr %>% xtable(
           label = "bias_b")

df.var.irr <- df.var.irr %>% 
  relocate(model8, .after = coef)
colnames(df.var.irr) <- c("Coef","TrueModel","addx2","addx5","addx7",
                          "addx2x5","addx2x7","addx5x7", "addx2x5x7")
df.var.irr %>% xtable(
                 digits = 3,
                 label = "var_b")

### EstVariance
EstV.irr <- c()
for (i in 1:8){
  md <- irm_list[[i]]
  .var <- md$residuals^2 %>% sum() / md$df.residual
  EstV.irr <- c(EstV.irr, .var) %>% round(3)
}
names(EstV.irr) <- c("TrueModel","addx2","addx5","addx7",
                     "addx2x5","addx2x7","addx5x7", "addx2x5x7")
EstV.irr %>% as.matrix() %>% 
  t() %>% 
  as.data.frame() %>% 
  xtable(digits = 3,
         label = "evar_b")


