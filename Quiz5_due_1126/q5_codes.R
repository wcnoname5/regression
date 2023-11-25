# Quiz5
### Packages
library(tidyverse)
library(stargazer)
### Data
dta <- matrix(
  c(
    c(2.0340, 49, 16), c(2.7183, 24, 9), c(3.7062, 48, 34),
    c(4.3929, 49, 47), c(5.0028, 50, 47), c(5.4739, 8, 8)
  ),
  ncol = 3,
  byrow = TRUE,
  dimnames = list(NULL,
                  c("dose", "group_size", "positive_response")
  )
) %>% 
  as.data.frame() %>%
  mutate(
    #neg_response  = group_size - positive_response,
    positive_prop  = positive_response/group_size
  )

# Q1
## all models
model_list <- list()
.link <-  c("logit","probit","cloglog")
for (i in 1:3){
  .lf <- .link[i]
  .m <- glm(data=dta,
            positive_prop ~ dose, family = binomial(.lf),
            weights = group_size
  )  
  model_list <- c(model_list, list(.m))
  
} 
names(model_list) <- .link
### quasi-lik
for (i in 1:3){
  .lf <- .link[i]
  .m <- glm(data=dta,
            positive_prop ~ dose, family = quasibinomial(.lf),
            weights = group_size
  )  
  model_list <- c(model_list, list(.m))
} 
names(model_list)[4:6] <- paste0("quasibinom_",.link)
## Regression table
stargazer(model_list, type= "text",
          star.cutoffs = c(.05, .01,.001),
          title = "Models",
          label = "glms")


###Visulization

dta %>%
  mutate(se_per = sqrt(positive_prop * (1 - positive_prop) / group_size)) %>% 
  # Create the bar chart with error bars
  ggplot(aes(x = dose, y = positive_prop),
         fill = "grey35")+ #, fill = as.factor(dose))) +
  geom_point(size = 3) +
  stat_smooth(data = dta, aes(weight = group_size),method="glm",
              method.args=list(family=binomial("logit")),
              se = FALSE, color = "red", lwd=1.1, lty=1) +
  stat_smooth(data = dta, aes(weight = group_size),method="glm",
              method.args=list(family=binomial("probit")),
              se = FALSE, color = "darkcyan", lwd=1.1, lty=2)+
  stat_smooth(data = dta, aes(weight = group_size),method="glm",
              method.args=list(family=binomial("cloglog")),
              se = FALSE, color = "darkorange", lwd=1.1, lty=3)+
  geom_errorbar(aes(ymin = positive_prop - se_per,
                    ymax = positive_prop + se_per),
                position = position_dodge(0.9), width = 0.2) +
  scale_color_manual(values = c("red", "darkcyan", "darkorange"),
                     name = c("Logit", "Probit", "Cloglog")) +
  labs(x = "Dose",
       y = "Proportion of Positive Response") +
  theme_minimal()

####Function LLOCV

SPSE<- function(data, link = "logit" , quasi = F){
  .SPSE <- 0
  for (i in 1:nrow(data)){
    if (quasi){
      .m <- glm(data = data[-i,], #LLO
                positive_prop ~ dose, family = quasibinomial(link),
                weights = group_size
      )  
    }else{
      .m <- glm(data = data[-i,], #LLO
                positive_prop ~ dose, family = binomial(link),
                weights = group_size
      )
    }
    .p <- predict(.m , data[i,], type = "response")
    .true <-  data[i,"positive_prop"]
    .SPSE <-  .SPSE + (.p-.true)^2
    #print(paste(i,.SPSE))
  }
  return(as.numeric(.SPSE))
}

## Model fits
.name <- c()
.a <- c()
.rr <- c()
vSPSE <- c()
#res and deviance
for (i in model_list[1:6]){
  .r <-  residuals(i, type = "pearson")
  .rr <- c(.rr, .r)
  .p <- sum(
    (residuals(i, type = "pearson"))^2
  )
  #print("Pearson's chi-square:")
  .d <- i$deviance
  .a <- c(.a, .p, .d)
}
# LLO
for (j in c(T,F)){
  for (linkfunc in .link){
    if (j) .name <- c(.name, paste0("quasi_",linkfunc))
    else .name <- c(.name, linkfunc)
    #print(paste("quasi"))
    .name <- c(.name, linkfunc)
    vSPSE <- c(vSPSE, SPSE(dta, linkfunc, T))
  }
}

.rr <- .rr %>% matrix(nrow = 6) %>% 
  as.data.frame()
colnames(.rr) <- names(model_list)
.rr <- .rr %>% 
  mutate(group = 1:6) %>% 
  relocate(group) %>%
  pivot_longer(cols = logit:cloglog, 
               names_to = "model",
               values_to = "residuals")  


.a <- .a %>% matrix(nrow = 2) %>%
  rbind(vSPSE)%>% 
  as.data.frame(row.names = c("Pearson's chi-square","Deviance", "SPSE"))

colnames(.a) <- names(model_list)
.a  %>% xtable::xtable(digits = 3)

##Visualization - residuals

ggplot(.rr, aes(x = as.factor(group), y = residuals,
                pch = model,color = model))+
  geom_point(size = 3 ) +
  geom_line(aes(group= model), lwd = .5, lty =4) +
  geom_abline(intercept = 0, slope = 0, color = "grey30",lty=1 ) + 
  scale_shape_manual(values=c(18, 19, 17)) + 
  scale_color_manual(values = c("red", "darkcyan", "darkorange"))+
  labs(x = "Group",
       title = "Pearson's Residuals")+
  theme_classic()
## refooting
matrix(c( 
  model_list$logit$coefficients /( pi/sqrt(3)), 
  model_list$probit$coefficients ,
  c(model_list$cloglog$coefficients[1]+0.5772,
    model_list$cloglog$coefficients[2]
  )/ (pi/sqrt(6))),
  byrow = F,
  nrow = 2, dimnames = list( c("Intercept", "dose"),
                             c("Logit", "Probit", "C_log-log"))) %>% 
  xtable::xtable(digits = 3)


#Q2
### Visulizaiton: functions of g

g <- function(mu, a){
  .g <- log(((1-mu)^(-a)-1)/a)
  return(.g)
}
alp<- c(.1, .5, 1, 2, 5)
.grey <- seq(80,15,-15)[1:5]
for (i in 1:5) {
  curve(g(x, a = alp[i]), add = (i!=1),ylim = c(-4, 4),
        lwd = 2.5, col = paste0("grey", .grey[i]),
        lty = 6,
        main = "Link Function g",
        xlab = "mu (pi)", ylab = "eta")
}
curve(qlogis(x), add = T ,ylim = c(-4, 4),
      lwd = 2.5, col = "red", lty=2)
curve(qnorm(x), add = T ,ylim = c(-4, 4),
      lwd = 2.5, col = "darkcyan", lty=3)
curve(log(-log(1-x)), add = T ,ylim = c(-4, 4),
      lwd = 2.5, col = "darkorange", lty=4)
legend("bottomright",
       legend = c(paste("alpha",alp, sep="="), "Logit", "Probit", "C log-log"),
       col = c(paste0("grey", .grey), "red", "darkcyan", "darkorange"),
       lty = c(rep(6 ,5), 2, 3, 4), lwd = 2.5, border = "black",
       cex = .70)


### Algorithm
# Self defined link function g
gf <- function(alpha=1){
  if (alpha < 0) {
    cat("Error, alpha <0 is not allowed")
  }
  linkfun <- function(mu){
    if (alpha != 0) return(log(((1-mu)^(-alpha)-1)/alpha))  
    else return(log(-log(1-mu)))
  }
  
  linkinv <- function(eta){
    if (alpha != 0) return(1-(alpha*exp(eta)+1)^(-1/alpha))
    else return((1-exp(-exp(eta))))
  }
  ## derivative of invlink wrt eta
  mu.eta <- function(eta){
    if (alpha!=0) return((alpha*exp(eta)+1)^((-1/alpha)-1)*exp(eta))
    else return( exp(-exp(eta))*exp(eta) )
  } 
  valideta <- function(eta) TRUE
  link <- paste0("g(", alpha, ")")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link),
            class = "link-glm")
}

# This allow a shortcut for calling glm
working_model <- function(alpha){
  return (
    glm(data = dta,
        positive_prop ~ dose, family = binomial(link = gf(alpha=alpha)),
        weights = group_size)
  )
} 
# Main Dish
find_alpha <- function(u=5, l=0, crit=10^(-5), N=100){
  eps <- 10
  itinterv <-  c(l, u)
  i = 0 #iteration counter
  while (eps > crit){ #stopping criterion
    md <- mean(itinterv)
    .range <- itinterv[2] - itinterv[1]
    #calculate deviance
    dev <- c(working_model(itinterv[1])$deviance,
             working_model(md)$deviance,
             working_model(itinterv[2])$deviance)
    if (which.min(dev) ==1){ 
      if (itinterv[1] - (md/(i*10))<0) itinterv <- c(itinterv[1], md) #not allow <0
      else itinterv <- c(itinterv[1] - (md/(i*10)), md)
    }
    else if(which.min(dev) == 3) itinterv <- c(md, itinterv[2] + (md/(i*10)))
    else itinterv <- c(md-(.range/4), md+(.range/4))
    #cat(itinterv)
    eps <- abs(working_model(itinterv[2])$deviance - working_model(itinterv[1])$deviance)
    #print(eps)
    i <- i+1
    if(i>N){print("Over 100 iteration"); break}
  }
  #print(i)
  alphaset <- c(itinterv[1], mean(itinterv), itinterv[2])
  dev <- c(working_model(itinterv[1])$deviance,
           working_model(mean(itinterv))$deviance,
           working_model(itinterv[2])$deviance)
  return(alphaset[which.min(dev)])
}
find_alpha()
