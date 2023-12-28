setwd('/Users/zhangying/Desktop/哥大/Fall 2023/GR5221/project')
data=get(load("ar1.s.RData"))
#A,Plot the data and comment on any interesting features
#There are two interesting features of the graph of the residuals. The first is the absence of any discernible trend. 
#The second is the smoothness of the graph. (In particular, there are long stretches of residuals that have the same sign. 
#This would be very unlikely to occur if the residuals were observations of iid noise with zero mean.) 
#Smoothness of the graph of a time series is generally indicative of the existence of some form of dependence among the observations.
#Such dependence can be used to advantage in forecasting future values of the series.
plot(data,xlab="",ylab="",main="Data series"); 
fit.data=lm(data~1+time(data))
lines((1:100),fit.data$fitted.values,lty=1,col="red")
y=fit.data$residuals
plot(y,type="o",pch=22,lty=1,pty=2,xlab="x",ylab="",main="Residuals from a linear fit")
abline(h=0)
acf(y,lag.max=40,main="ACF of residuals from a linear fit")
Box.test(y,lag=20,type="Box-Pierce")# test independence

# try to  difference data at lag 1
data.d=diff(data)
plot(data.d,xlab="",ylab="",main="differenced Data series"); abline(h=0)
acf(data.d,lag.max=40,main="ACF of differenced data series")
acf(data.d,lag.max=40,type="partial",main="PACF of differenced data series")

#B,Based on the sample ACF and PACF plots, explain why an AR(1) model is plausible to fit the data.
# The sample PACF shown in Figure strongly suggests fitting an AR(1) model to the 
#mean-corrected data Xt = Yt − 45.7743.
par(mfrow=c(2,1))
acf(data,lag.max=40,main="ACF")
acf(data,lag.max=40,type="partial",main="PACF")
#by test,we detect the series is not stationary, and not have trned and sesaonal component,noise is not iid
adf.test(data)
kpss.test(data, null="Trend")
seasonal.test(data,seasonal = c("ocsb","ch","hegy"),alpha = 0.05)

#C,Find the Yule–Walker estimates of 𝜇, 𝜙, and 𝜎 ଶ.
#In this example we shall consider the problem of fitting an AR process 
#directly to the data without first removing any component.
#(Xt-45.7743)-0.9613(X(t-1)-45.7743)=Zt and Zt follow WN(0,8.270659)
acf2(data.mc,48)
mean(data)
data.mc=data-mean(data)
fit.data.1=ar(data.mc,order.max=1) # default method is "yule-walker"
fit.data.1$ar; fit.data.1$var.pred
#D Find a 95% confidence interval for 𝜙 using the result 
#the confidence interval (0.90698,1.01557)
n=length(data)
fit.data.1$ar+(1.96*(fit.data.1$var.pred)^(1/2))/((var(data))^(1/2)*sqrt(n))
fit.data.1$ar-(1.96*(fit.data.1$var.pred)^(1/2))/((var(data))^(1/2)*sqrt(n))
#Data preprogressing
model <- arima(data, order = c(1,0,0), method = "ML")
phi_hat <- model$coef[1] # phi_hat
mu_hat <- model$coef[2]  # mu_hat
z_hat <- data - mu_hat - phi_hat * (c(mu_hat, head(data, -1)) - mu_hat)#calculate residual
B <- 500   # sample number
n <- length(data)
phi_hat_star <- numeric(B)  # ans phi_hat
#simulation
for (b in 1:B) {
  success <- FALSE
  while (!success) {
    z_star <- sample(z_hat[-1], n - 1, replace = TRUE)#select from residul
    
    x_star <- numeric(n)
    x_star[1] <- data[1] 
    
    for (t in 2:n) {
      x_star[t] <- mu_hat + phi_hat * (x_star[t - 1] - mu_hat) + z_star[t - 1]
    }#create the dataset
    
    
    fit_result <- tryCatch({
      arima(x_star, order = c(1,0,0), include.mean = TRUE, method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(fit_result) && !is.na(fit_result$coef[1])) {
      phi_hat_star[b] <- fit_result$coef[1]
      success <- TRUE
    }#to get 500 available parameter value
  }
}
phi_hat_star
hist(phi_hat_star, breaks = 30, main = "Histogram of 500 phi_hat*", xlab = "phi_hat*", col = "blue")
mean <- mean(phi_hat_star, na.rm = TRUE)
sd <- sd(phi_hat_star, na.rm = TRUE)
curve(dnorm(data, mean = mean, sd = sd), add = TRUE, col = "red", lwd = 2)#normal distribution curve
diffs <- phi_hat_star - phi_hat
ci_lower <- quantile(diffs, 0.025)
ci_upper <- quantile(diffs, 0.975)
ci_phi <- c(phi_hat - ci_upper, phi_hat - ci_lower)
ci_phi
