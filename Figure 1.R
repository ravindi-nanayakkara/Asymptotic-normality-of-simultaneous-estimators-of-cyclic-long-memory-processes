pck <-
  c("waveslim",
    "wmtsa",
    "latex2exp")
lapply(pck, library, character.only = TRUE)

set.seed(654321)

Ts <- 1000
n0 <- min(floor(Ts/10), 100)

n <- Ts + 2*n0

#The function  dwpt.sim simulates a seasonal persistent process(cyclic long-memory time series) using the discrete wavelet packet transform(DWPT).
ts1 <- dwpt.sim(n, "mb16", .4, 0.1,epsilon=.001)

######## Figure 1(a) - This figure gives the plot of a realization of a cyclic long-memory time series.
plot(ts1, type="l", xlab="Time", ylab="Value")

######## Figure 1(b) - This figure gives the plot of the periodogram corresponding to the simulated cyclic long-memory time series.
plot(0:(n/2)/n, per(ts1), type="l", xlab="Frequency", 
     ylab="Value")

######## Figure 1(c) - This figure gives the plot of the sample covariance function.
acf(ts1, lag.max=100, ylim=c(-.6,1), main =" ")

#The function wavCWT() computes the continuous wavelet transform of a time series.
ts1.cwt <- wavCWT(ts1)

######## Figure 1(d) - This figure gives the plot of the wavelet coefficients.
plot(ts1.cwt)

