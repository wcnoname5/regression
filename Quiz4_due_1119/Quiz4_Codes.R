#Packages
library(tidyverse)
library(stargazer)
library(GGally)
library(glmnet)

# Q1. (3.14, 5.5, 12.8)
### plot
E3_7 <- read_table("E3.7.txt",n_max = 20)

E3_7[-1] %>% ggpairs()

## 3.14
full_mod <- lm(data = E3_7, 
             y ~ x_1 + x_2 + x_3 + x_4 + x_5)
reduce_mod <- lm(data = E3_7, y ~ x_3 + x_5)  

stargazer(full_mod, reduce_mod, report = "vc*t",
          star.cutoffs = c(.05,.01,.001))

## 5.5

###qqplot
par(mfrow=c(1,2))
qqnorm(full_mod$residuals,
       main = "Normal Q-Q plot of residuals\n\ in Full Model"
       );qqline(full_mod$residuals, col=2)

qqnorm(reduce_mod$residuals,
       main = "Normal Q-Q plot of residuals\n in Reduced Model"
       );qqline(reduce_mod$residuals, col=2)

## 12.8 
## filter out outliers
mask <- E3_7$y<1.5 & E3_7$x_4>60

### PCA
# Following Codes were assist by ChatGPT
# Standardize the predictor variables (optional but recommended for PCR)
scaled_data <- scale(E3_7[-c(1,7)])
.scaled_data <- scale(E3_7[mask,-c(1,7)])
# Perform principal component analysis (PCA)
pca_result <- prcomp(scaled_data)
.pca_result <- prcomp(.scaled_data)
# Use the selected principal components as predictors
pcr_data <- predict(pca_result, newdata = scaled_data)[, 1:2]
.pcr_data <- predict(.pca_result, newdata = .scaled_data)[, 1:2]
# Combine with the response variable (mpg in this example)
response_variable <- E3_7$y
pcr_data <- cbind(response_variable, pcr_data)
# Fit a linear regression model using the principal components
pcr_model <- lm(response_variable ~ ., data = as.data.frame(pcr_data))
.pcr_model <- lm(V1 ~ .,
                 data = as.data.frame(cbind(E3_7$y[mask], .pcr_data)))
# Display the summary of the PCR model

stargazer(pcr_model,.pcr_model,
          star.cutoffs = c(.05, .01,.001), type = "text")
stargazer(pcr_model,.pcr_model,
          star.cutoffs = c(.05, .01,.001))

### Ridge
.N <- nrow(E3_7)
# Standardized the covariates
std_iv <- scale(as.matrix(E3_7[-c(1,7)]))

# Fit Ridge regression model using cross-validation to choose lambda
ridge_model <- cv.glmnet(std_iv, E3_7$y, alpha = 0, #alpha=0 ridge
                         nfolds = .N) #LOO cv

# Print the cross-validated results

par(mfrow= c(1,2))
plot(ridge_model, main="Full Sample(n=20)")


.std_iv <- scale(as.matrix(E3_7[mask,-c(1,7)]))

# Fit Ridge regression model using cross-validation to choose lambda
.ridge_model <- cv.glmnet(.std_iv, E3_7$y[mask], alpha = 0, #alpha=0 ridge
                         nfolds = 18) #LOO cv

# Print the cross-validated results

plot(.ridge_model, main = "Drop Outliers(n=18)")

## compare effect of ouliers
.cffull <- coef(ridge_model, s= ridge_model$lambda.1se) %>% as.matrix()
.cfdrop<- coef(.ridge_model, s= ridge_model$lambda.1se) %>% as.matrix()

ss <- cbind(.cffull,.cfdrop) %>%
  as.data.frame() 

colnames(ss) <- c("Full Sample", "Drop Outlier")
ss %>%   xtable::xtable(digits = 3,
                         caption = "coefficients of ridge regression",
                         label = "ridge_coef")

# Q2. (3.11, 3.18)
### read data
E3_6 <-
  read_table("E3.6.txt",
    n_max =50) %>%
  mutate(sal_rat = Y84/Y83, sal_diff = Y84 -Y83)

## 3.11
### plotting
E3_6 %>%
  mutate(across(c(3:5,7),log))%>%
  select(-c(1)) %>% 
  ggpairs() +ggtitle("SHARES, REV, INV log transformed")

### regressions
OLS_rat <-  lm(sal_rat ~ Y83 + SHARES + REV + INC + AGE,
               E3_6)
OLS_diff <-  lm(sal_diff ~ Y83 + SHARES + REV + INC + AGE,
               E3_6)
# add 1/Y83 term
OLS_diffl <- lm(sal_diff ~ Y83 + I(1/Y83) + SHARES + REV + INC + AGE,
               E3_6)
OLS_ratl <-  lm(sal_rat ~ Y83+ I(1/Y83) + SHARES + REV + INC + AGE,
               E3_6)
stargazer(OLS_rat, OLS_ratl, OLS_diff, OLS_diffl,# type = "text",
          star.cutoffs = c(.05,.01,.001),
          digits=2, dep.var.labels = c("salary ratio", "salary difference"))

## 3.18
### dealing with data
E3_9 <- read_table("E3.9.txt", col_types = list(y_1 = col_double()),
    skip = 1, n_max = 9)
# w/o aminophylline
E3_9na <- E3_9[1:3] %>%  
  mutate(ID = as.factor(1:nrow(E3_9)),
        ap = 0) 
# with aminophylline
E3_9a <-E3_9[4:6] %>% 
  rename(all_of(c(x_1 = "x_1_1", x_2 = "x_2_1", y = "y_1"))
    ) %>% 
  mutate(ID = as.factor(1:nrow(E3_9)),
        ap = 1) 
dta3_9 <- bind_rows(E3_9na, E3_9a)

#regression

m318 <- lm(y~ ap + x_1 + x_2 , dta3_9)
# add participant effect (\pi_s)
m318_rm <- lm(y ~ ap + x_1 + x_2 + ID, dta3_9) 
stargazer(m318, m318_rm,
          star.cutoffs = c(.05,.01,.001))

# Q4. (7.7, 7.8)
## 7.7
#data
E7_7<- tibble(
  year = seq(1790, 1970, by=10),
  pop = c(3929, 5308, 7239, 9638, 12866, 17069, 23191,
          31443, 39818, 50155, 62947, 75994, 91972, 105710,
          122775, 131669, 151325, 179323, 203211)
)
# time series object
e77ts <- ts(data = c(3929, 5308, 7239, 9638, 12866, 17069, 23191,
                     31443, 39818, 50155, 62947, 75994, 91972, 105710,
                     122775, 131669, 151325, 179323, 203211),
            start = 1790, end = 1970, frequency = 1/10)
#Durbin-Watson Test
library(car)
car::durbinWatsonTest(
  lm(pop~year, E7_7)
)
# ar(1), poly to 3rd degree
ar1m <- arima(sqrt(e77ts), order = c(1, 0, 0), 
              xreg = poly(time(e77ts), 3)
)
ar2m <- arima(e77ts, order = c(1, 0, 0), 
              xreg = poly(time(e77ts), 3)
)
stargazer(ar1m, ar2m, # type ="text",
          star.cutoffs = c(.05, .01, .001),
          dep.var.labels = c("sqrt(y)","y")
)
## 7.8

#first order difference
E7_7$dpop = c(NA, diff(E7_7$pop))

# For visulization
# plot(E7_7$year, E7_7$dpop,xlab = "Year",
#      ylab = "Increased populaiton",
#      main = "Population change per decade in US")
# plot(E7_7$year, log10(E7_7$dpop), xlab = "Year",
#      ylab = "Increased populaiton(log10)",
#      main = "Population change per decade in US")

pop_diff_model <- arima(diff(e77ts),
                        order = c(1, 0, 0), 
                        xreg = time(diff(e77ts))
)
predict(pop_diff_model, newxreg = 1980)$pred
