
##########################
#### FUNCTION GIRTEST ####
##########################

## Paper by Ganics, Inoue, Rossi (2018)
## Code by Oriol Gonzalez and Marko Irisarri


# Symmetric range CI (inspired by John T Kent (2005))
girtest=function(Fstat,K_2) {
  # Fstat is the F-statistic of the instruments in the first stage
  # K_2 is the number of instruments


  # Create the chi-squared function and then invert the chi (not the chi-squared!)
  F=function(u,p,lambda) pchisq(u^2,p,lambda^2)                 # cdf of noncentrality parameter
  Finv=function(prob,p,lambda) sqrt(qchisq(prob,p,lambda^2))    # quantile


  options(warn=-1)                                            # remove warnings in output

  alpha <- 0.05                                               #Tables in Ganics et al. (2018) only for 5% nomial level
  p <- K_2
  y <- sqrt(Fstat*K_2)

  #Function to obtain the CI for the concentration parameter
  u0=Finv(1-alpha,p,0)
  if(y<=u0) ll=0
  else {
    f=function(b) F(y,p,y-b)-F(max(y-2*b,0),p,y-b)-(1-alpha)
    bb=uniroot(f,c(0,y))$root
    ll=y-bb
  }
  f=function(b) F(y+2*b,p,y+b)-F(y,p,y+b)-(1-alpha)
  big=1; obj=-1
  while(obj<0) {big=2*big; obj=f(big)}
  bbb=uniroot(f,c(0,big))$root
  lu=y+bbb
  ci_mu=c(ll^2/p,lu^2/p)                                      #CI concentration parameter chi^2 (note the correction!)


  #Tables with critical values:
  bias <- read.csv("https://raw.githubusercontent.com/OriolGC/girtest/master/Table_cv_b.csv", header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  bias[,1] <- NULL
  x <- as.numeric(gsub("X", "", colnames(bias)))
  x[1] = 0

  size <- read.csv("https://raw.githubusercontent.com/OriolGC/girtest/master/Table_cv_s.csv", header = TRUE, sep = ";")
  size[,1] <- NULL
  y <- as.numeric(gsub("X", "", colnames(size)))
  y[1] = 0


  #Take the corresponding values from the table corresponding to the previous CI
  low_mu <- round(ci_mu[1],digits = 2)
  high_mu <- round(ci_mu[2],digits = 2)


  column_low_bias =  which(abs(x - low_mu) == min(abs(x - low_mu)))
  column_high_bias = which(abs(x - high_mu) == min(abs(x - high_mu)))

  column_low_size =  which(abs(y - low_mu) == min(abs(y - low_mu)))
  column_high_size = which(abs(y - high_mu) == min(abs(y - high_mu)))

  bias_low <- as.numeric(unlist(bias[p-1, column_high_bias]))
  bias_high <- as.numeric(unlist(bias[p-1, column_low_bias]))

  size_low <- as.numeric(as.character(unlist(size[p, column_high_size])))
  size_high <- as.numeric(as.character(unlist(size[p, min(column_low_size)])))


  #Create the final table and display the results
  cat(underline$bold("\nConfidence intervals by Ganics, Inoue and Rossi (2018)\n\n"))
  cat(c("Concentration Parameter", "[",low_mu,";", high_mu,"]"),"\n")
  cat(c("Bias                   ", "[",bias_low,";", bias_high,"]"),"\n")
  cat(c("Size distortion        ", "[",size_low,";", size_high,"]"),"\n")

}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("read_excel"))

## Example of use: Replication Angrist and Krueger (1991)/Bound et al. (1995)
# Angrist and Krueger (1991): 28 instruments and F-stat = 1.61 at the 5% level
# Bound et al. (1995): 3 instruments and F-stat = 13.49 at the 5% level
# girtest(1.61,28)
# girtest(13.49,3)


