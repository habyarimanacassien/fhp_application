rm(list=ls())

# Useful libraries
library(MittagLeffleR)
library(pracma)
library("rlist")
require(parallel)
library(ggplot2)
library(LaplacesDemon)
library(arm)
library(lubridate)
library(dplyr)


################################################################################
# Function to simulate observation events for FHP
simulate_events <- function(lambda, alpha, beta, Time){
  #It simulates and returns fractional Hawkes process number of observations
  #and time realization points
  # Hawkes process time realization
  library(MittagLeffleR)
  library(LaplacesDemon)
  library(hawkesbow)
  
  x <- hawkes(Time,fun=lambda,repr=alpha,family=function(n) {rml(n,beta)})
  SimPts <- c(0, x$p)
  
  return (length(SimPts))
  #end function
}

############ Function to simulate observation points for FHP    ###############
simulate_points <- function(lambda, alpha, beta, Time){
  #It simulates and returns fractional Hawkes process number of observations
  #and time realization points
  # Hawkes process time realization
  library(MittagLeffleR)
  library(LaplacesDemon)
  library(hawkesbow)
  
  x <- hawkes(Time,fun=lambda,repr=alpha,family=function(n) {rml(n,beta)})
  SimPts <- c(0, x$p)
  
  return (SimPts)
  #end function 
}

################################################################################
#                SUMMARY STATISTICS & DISTANCE FUNCTIONS TO USE                #
################################################################################

# 1. Function to compute_diff_log_nbr_distance (LN)
compute_diff_log_nbr_distance <- function(distr1, distr2){
  # It compute and return the wasserstein distance between two distributions
  n2 <- length(distr2)
  n1 <- length(distr1)
  distance = abs(log(n2)-log(n1))
  return (distance)
  #End function
}

# 2. Function to compute_wasserstein_distance (WS)
compute_wasserstein_distance <- function(distr1, distr2, Time){
  # It compute and return the wasserstein distance between two distributions
  n = min(length(distr1),length(distr2))
  m = max(length(distr1),length(distr2))
  
  if (length(distr1) > length(distr2)) {
    if (n==0){
      distance = (m-n)*Time - sum(na.omit(distr1[n+1:m]))
    } else {
      distance = sum(abs(distr1[1:n]-distr2[1:n])) + 
        (m-n)*Time - sum(na.omit(distr1[n+1:m]))
    }
  } else if (length(distr1) < length(distr2)) {
    if (n==0){
      distance = (m-n)*Time - sum(na.omit(distr2[n+1:m]))
    } else {
      distance = sum(abs(distr1[1:n]-distr2[1:n])) + 
        (m-n)*Time - sum(na.omit(distr2[n+1:m]))
    }
  } else {
    distance = sum(abs(distr1-distr2))
  }
  return (distance)
  #End function
}

# 3. Function to compute_meandiff_grt_q90_distance (MD90)
compute_meandiff_grt_q90_distance <- function(distr1, distr2){
  # It compute and return the diff between mean of event time diff
  # that great than their q90 for two distributions 
  q1 <- unname(quantile(diff(distr1), probs = .90))
  s1 <- mean(diff(distr1)[diff(distr1) > q1])
  if (length(distr2) > 2){
    q2 <- unname(quantile(diff(distr2), probs = .90))
    s2 <- mean(diff(distr2)[diff(distr2) > q2])
    distance <- abs(s2 - s1)
  } else {distance <- Inf}
  return (distance)
  #End function
}

# 4. Function to compute_meandiff_less_q50_distance (MD50)
compute_meandiff_less_q50_distance <- function(distr1, distr2){
  # It compute and return the diff between mean of event time diff
  # that great than their q90 for two distributions 
  q1 <- unname(quantile(diff(distr1), probs = .50))
  s1 <- mean(diff(distr1)[diff(distr1) < q1])
  
  if (length(distr2) > 2){
    q2 <- unname(quantile(diff(distr2), probs = .50))
    s2 <- mean(diff(distr2)[diff(distr2) < q2])
    distance <- abs(s2 - s1)
  } else {distance <- Inf}
  return (distance)
  #End function
}

# 5. Function to digamma_distance (DIGAD)
digamma_distance <- function(distr1, distr2){
  if (length(distr2) > 1){
    s1 = sum(digamma(diff(distr1)), na.rm = TRUE)
    s2 = sum(digamma(diff(distr2)), na.rm = TRUE)
    if(is.finite(s1) & is.finite(s2)){
      distance = abs(s1-s2)
    }else{
      distance = Inf
    }
  } else {
    distance = Inf
  }
  return (distance)
  #End function
}

################################################################################
#                      FUNCTION TO COMPUTE THE LIKELIHOOD                      #
################################################################################

# Likelihood function  for FRACTIONAL HAWKES PROCESS
compute_logLikelihood <- function(distr, lambda, alpha, beta, Time){
  library(MittagLeffleR)
  k = length(distr)
  ll1 = 0
  for (i in 1:k){
    si = 0
    if (i > 1){
      for(j in 1:(i-1)){
        si = si + ((distr[i]-distr[j])^(beta-1))*mlf(-(distr[i]-distr[j])^beta,beta,beta,1)
      }}
    ll1 = ll1 + log(lambda + alpha*si) - alpha*(1 - mlf(-(distr[k]-distr[i])^beta,beta,1,1))
  }
  ll <- ll1 - lambda*Time
  #  likelihood <- exp(ll)
  return (ll)
  #end function
}


################################################################################
data_finance <- read.csv('GE.csv') #Read the file
price <- data_finance$Volume#Adj_Close #Extract the adjusted close price
log_ret <- diff(log(price)) #Compute the logarithmic returns
abs_log_ret <- abs(log_ret) #Compute the absolute value of logarithmic returns

# Add variable 'abslogret' to the data set 'data_finance'
data_finance$abs_log_ret <- c(0, abs_log_ret)

excess_tres <- quantile(abs_log_ret, probs = .9) #Compute the 90th percentile of logarithmic returns

# Select data at which the logarithmic returns exceed a treshold value (exceedances)
data <- data_finance %>% filter(abs_log_ret > excess_tres ) %>% select(Date, abs_log_ret)

# Convert 'Date' to the 'YYYY-MM-DD' format
data$time_real <- mdy(data$Date)  # Converts 'Date' to Date type in 'YYYY-MM-DD' format

# Set a start date
start_time <- data_finance$Date[1]
start_time_real <- mdy(data_finance$Date[1])

# Add a start date to the data frame
# Define the new row to add at the first index
new_row <- data.frame(
  Date = start_time,         
  abs_log_ret = 0,               # 0 value for 'abs_log_ret' start time!
  time_real = start_time_real       
)

# Bind the new row to the top of the existing dataframe
data <- rbind(new_row, data)

# Calculate the difference in days and create a new variable 'time_absolute_days'
data$time_absolute_days <- as.numeric(difftime(data$time_real, start_time_real, units = "days"))

# Add 'time_absolute_months' and 'count' variables
data$time_absolute_months <- data$time_absolute_days/30.44
data$count <- seq(0, (length(data$Date)-1))

# Arrange variables
data <- data %>% select(Date, (ncol(data)-3):(ncol(data)), 2:(ncol(data)-3))

################## Finance exceedances Time Realization ###############

# Save the figure in a PDF file
pdf(file = "Fig_Finance_exceedances_time_realization.pdf", width = 8.5, height = 4)

# Adjust margins to create extra space at the bottom
par(mar = c(6, 4, 4, 2) + 0.1)

# Define custom labels (ensure same length as breaks)
custom_labels <- seq(min(data$time_real), max(data$time_real), by = "5 years")  # Adjust based on your data

# Extract the start and end count values for the height of the rectangle
start_count <- data$count[data$time_real == as.Date("2006-10-30")]
end_count <- data$count[data$time_real == as.Date("2007-08-30")]

# Plot with jump function
ggplot(data, aes(x = time_real, y = count)) +
  geom_step(color = "blue", linewidth = 0.6) +
  labs(x = "Time", y = "Number of Events") +
  theme_bw() +
  scale_x_date(breaks = seq(min(data$time_real), max(data$time_real), by = "5 years"),
               labels = custom_labels) +
  annotate("rect", xmin = as.Date("2006-10-30"), xmax = as.Date("2007-08-30"),
           ymin = start_count, ymax = end_count, alpha = 0.2, fill = "red")

# Close the PDF device
dev.off()


#############################ZOOMING
# Save the figure in a PDF file
pdf(file = "Fig_Finance_exceedances_time_realization_zoomed.pdf", width = 8.5, height = 4)

# Adjust margins to create extra space at the bottom
par(mar = c(6, 4, 4, 2) + 0.1)
# Define custom labels (ensure same length as breaks)
custom_labels <- seq(as.Date("2006-10-30"), as.Date("2007-08-30"), by = "2 months")  # Adjust based on your data
library(dplyr)
data_zoom <- data %>% filter(time_real >= as.Date("2006-10-30") & time_real <= as.Date("2007-08-30"))
data_zoom$count <- seq(1, length(data_zoom$count))
ggplot(data_zoom, aes(x = time_real, y = count)) +
  geom_step(color = "blue", linewidth = 0.6) +
  labs(x = "Time", y = "Number of Events") +
  theme_bw() +
  scale_x_date(breaks = seq(as.Date("2006-10-305"), as.Date("2007-08-30"), by = "2 months"),
               labels = custom_labels)
# Close the PDF device
dev.off()



#################### Finance Distribution Using FHP2PP ######################

##### 3.  FHP2PP_LAMBDA
proc <- "FHP2PP"
par <- "Lambda"
dataset <- "Fig_Finance_exceedances"

SimPoints <- data$time_absolute_days   # Time realizations (in days)
Time <- ceiling(max(SimPoints))

Niter <- 100000 #Number of iterations
lambda <- length(SimPoints)/Time   # Number of Finance exceedances per day   
alpha <- 0.0001
beta <- .5

number_of_cores <- 24

# Check if the file exists in the current working directory
if (!file.exists(paste0("data_",proc,"_",par,"_",dataset,"_NSimPoints",".rdata"))) {
  
  # Number of Observed data (NSimPoints) with a true parameter value
  NSimPoints <- c()
  
  for(i in 1:Niter){
    n <- simulate_events(lambda=lambda, alpha=alpha, beta=beta, Time=Time)
    NSimPoints <- c(NSimPoints, n)
    cat("Case 1: Iteration ", i, " of ", Niter, " is completed!\n")
  }
  
  # Save data of 'NSimPoints'
  df1 <- data.frame(NSimPoints)
  save(df1, file=paste0("data_",proc,"_",par,"_",dataset,"_NSimPoints",".rdata"))
  cat("File for 'NSimPoints' has been created and saved.\n")
} else {
  # LOAD THE DATA SAVED
  load(file=paste0("data_",proc,"_",par,"_",dataset,"_NSimPoints",".rdata"))
  
  # Extract data
  NSimPoints <- df1$NSimPoints
  cat("File for 'NSimPoints' already exists & extracted.\n")
}

mode <- length(which(Mode(NSimPoints) == NSimPoints))
y2 <- mode/length(NSimPoints)
ymax <- y2*1.05

# Save FIGURE in a pdf file
pdf(file = "Fig_Finance_distribution_FHP2PP.pdf", width=8, height=4.2)
par(mfrow=c(1,1), oma=c(1,1,1,1), xaxs="i")  # Set xaxs to 'i' for exact range

# Create plot function
fig <- function(){
  discrete.histogram(NSimPoints, axes = FALSE, cex.lab = 1,
                     xlab = "Number of events", xlim = range(NSimPoints))  # Set xlim to range of data
  
  # Set axis labels with intervals to reduce overlap
  x_labels <- pretty(range(NSimPoints), n = 10)  # Adjust `n` for fewer/more labels
  axis(2, cex.axis=1)
}

fig()

# Close the pdf device
dev.off()

########################## PARAMETER ESTIMATION ###############################

## Exact Posterior
# parameter values interval: 
dx <- 1/2500  #bin size
LL <- (lambda-lambda) + dx  # lower limit of the breaks interval
UL <- (lambda+lambda) + dx # upper limit of the breaks interval
breaks_interval <- seq(LL-dx,UL+dx,by=dx)
mids_interval <- seq(LL-dx/2, UL+dx/2, by=dx)

## Compute 'logLikelihood'

# Check if the file exists in the current working directory
if (!file.exists(paste0("data_",proc,"_",par,"_",dataset,"_exact_posterior",".rdata"))) {
  
  # Loglikelihood for FHP
  logLikelihood <- c()
  for (m in 1:length(mids_interval)){
    cat("Case of FHP: Iteration ", m, " of ", length(mids_interval), " is completed!\n")
    logL <- compute_logLikelihood(lambda = mids_interval[m], distr = SimPoints, alpha = alpha, beta = beta, Time = Time)
    logLikelihood <- c(logLikelihood, logL)
  }
  
  likelihood <- exp(logLikelihood - mean(logLikelihood))  # translation of the logLikelihood to remove impossibilities
  posterior <- likelihood # For uniform prior
  exact_posterior <- posterior
  exact_posterior <- exact_posterior/sum(exact_posterior)/(dx) #Normalization
  
  # Save data of 'likelihood'
  df2 <- data.frame(exact_posterior)
  save(df2, file=paste0("data_",proc,"_",par,"_",dataset,"_exact_posterior",".rdata"))
  
  cat("Files for 'likelihoods' have been created and saved.\n")
} else {
  # LOAD THE DATA SAVED
  load(file=paste0("data_",proc,"_",par,"_",dataset,"_exact_posterior",".rdata"))
  
  # Extract data
  exact_posterior <- df2$exact_posterior
  cat("Files for 'exact_posterior' already exist & extracted.\n")
}

## ABC Posterior

# Check if the file exists in the current working directory
if (!file.exists(paste0("data_",proc,"_",par,"_",dataset,"_SimPointsth_list",".rdata"))) {
  
  # We generate a uniformly distributed prior for the parameter LAMBDA
  ParameterSample_values <- runif(Niter,min=(lambda-lambda), max=(lambda+lambda))
  
  # Generate data for each parameter value
  cl <- makeCluster(number_of_cores)  # number_of_cores is the number of clusters (cores)
  SimPointsth_list <- parSapply(cl , ParameterSample_values, simulate_points, alpha=alpha, beta=beta, Time=Time)
  stopCluster(cl)
  
  # Save data of 'ParameterSample_values'
  df4 <- data.frame(ParameterSample_values)
  save(df4, file=paste0("data_",proc,"_",par,"_",dataset,"_ParameterSample_values",".rdata"))
  
  # Save data of 'SimPointsth_list'
  save(SimPointsth_list, file=paste0("data_",proc,"_",par,"_",dataset,"_SimPointsth_list",".rdata")) 
  
  cat("Files for 'ParameterSample_values' and 'SimPointsth_list' have been created and saved.\n")
} else {
  # LOAD THE DATA SAVED
  load(file=paste0("data_",proc,"_",par,"_",dataset,"_ParameterSample_values",".rdata"))
  load(file=paste0("data_",proc,"_",par,"_",dataset,"_SimPointsth_list",".rdata"))
  
  # Extract data
  ParameterSample_values <- df4$ParameterSample_values
  cat("Files for 'ParameterSample_values' and 'SimPointsth_list' already exist & extracted.\n")
}


# Compute distance

# Check if the file exists in the current working directory
if (!file.exists(paste0("data_",proc,"_",par,"_",dataset,"_data_frame1",".rdata"))) {
  
  cl <- makeCluster(number_of_cores)  # Open the cluster 
  dist1 <- "LN"
  out_distance1 <- parSapply(cl, SimPointsth_list, compute_diff_log_nbr_distance, distr1 = SimPoints)
  stopCluster(cl)  # Stop the cluster at the end
  
  cl <- makeCluster(number_of_cores)  # Open the cluster 
  dist2 <- "WS"
  out_distance2 <- parSapply(cl, SimPointsth_list, compute_wasserstein_distance, distr1 = SimPoints, Time = Time)
  stopCluster(cl)  # Stop the cluster at the end
  
  cl <- makeCluster(number_of_cores)  # Open the cluster 
  dist3 <- "MDG90"
  out_distance3 <- parSapply(cl, SimPointsth_list, compute_meandiff_grt_q90_distance, distr1 = SimPoints)
  stopCluster(cl)  # Stop the cluster at the end
  
  cl <- makeCluster(number_of_cores)  # Open the cluster 
  dist4 <- "MDL50"
  out_distance4 <- parSapply(cl, SimPointsth_list, compute_meandiff_less_q50_distance, distr1 = SimPoints)
  stopCluster(cl)  # Stop the cluster at the end
  
  cl <- makeCluster(number_of_cores)  # Open the cluster 
  dist5 <- "DIGAD"
  out_distance5 <- parSapply(cl, SimPointsth_list, digamma_distance, distr1 = SimPoints)
  stopCluster(cl)  # Stop the cluster at the end
  
  # DATAFRAMES FOR DISTANCES & SUMMARY STATISTICS, AND THEIR NAMES 
  data_frame1 <- data.frame(out_distance1, out_distance2, out_distance3, out_distance4, out_distance5)
  data_frame2 <- data.frame(dist1, dist2, dist3, dist4, dist5)
  
  # Save data 
  save(data_frame1, file=paste0("data_",proc,"_",par,"_",dataset,"_data_frame1",".rdata"))
  save(data_frame2, file=paste0("data_",proc,"_",par,"_",dataset,"_data_frame2",".rdata"))
  
  cat("Files for 'data_frame1' and 'data_frame2' of distances have been created and saved.\n")
} else {
  # LOAD THE DATA SAVED
  load(file=paste0("data_",proc,"_",par,"_",dataset,"_data_frame1",".rdata"))
  load(file=paste0("data_",proc,"_",par,"_",dataset,"_data_frame2",".rdata"))
  
  # Extract data
  data_frame1 <- data_frame1
  data_frame2 <- data_frame2
  cat("Files for 'data_frame1' and 'data_frame2' already exist & extracted.\n")
}


### ESTIMATION

# Choice on selection rate
rate_v <- c(.01, .025, .05)
rate_v <- c(.01, .05, .10)

################################################################################
# Save FIGURES in a pdf or png file
pdf(file = paste0("Fig", proc, "", par, "", dataset ,"_N", Niter, ".pdf"), width=12.5, height=9.5)
par(mar = c(4, 4, 3, 1))        # Reduce space around plots
par(mfrow=c(3,5),oma=c(3,3,3,3))

for (j in 1:length(rate_v)){
  
  for (idf in 1:length(data_frame1)){ # 
    
    out_distance <- as.numeric(data_frame1[[idf]])
    dist <- as.character(data_frame2[[idf]])
    
    rate <- rate_v[j]
    d <- quantile(out_distance, probs = rate)
    
    # Select parameter values corresponding to distances 
    # less than threshold distance d
    abc_posterior <- ParameterSample_values[which(out_distance < d)]
    
    ## Compute kullback–Leibler divergence
    
    density_LIK <- exact_posterior
    
    h <- hist(abc_posterior, breaks = breaks_interval, plot = FALSE)
    
    density_ABC <- h$density
    ymax <- max(density_LIK, density_ABC)
    
    #create plot
    if(dist != "LN"){
      h1 <- hist(abc_posterior, plot = FALSE)#probability=TRUE, plot = FALSE) 
      h2 <- hist(abc_posterior + (lambda-mean(h1$mids)), probability=TRUE, ylim=c(0,1.40*ceiling(ymax)), breaks = 25, 
                 col="darkolivegreen", border="#333333", xlab="", ylab="", main="")
      abc_posterior <- abc_posterior + (lambda-mean(h1$mids))
    } else {
      h2 <- hist(abc_posterior, probability=TRUE, ylim=c(0,1.40*ceiling(ymax)), breaks = 25,
                 col="darkolivegreen", border="#333333", xlab="", ylab="", main="")      
    }
    ###########
    
    # Generate the histogram with custom breaks
    hh <- hist(abc_posterior, breaks = breaks_interval, plot = FALSE)
    density_ABC <- hh$density
    
    ### Penalise zero (0) values of the densities
    eps <- 1e-10
    for (i in 1:length(density_ABC)){
      if (density_ABC[i]==0){
        density_ABC[i]=eps
      }
      if (density_LIK[i]==0){
        density_LIK[i]=eps
      }
    }
    
    ## We can find the KLD using hist densities
    KLD_value <- KLD(density_LIK,density_ABC)  ## KLD(Likelihood, ABC)
    
    dist_used = paste0("Dist : ",dist)
    threshold_dist = bquote("Threshold "*epsilon == .(round(d,3)))
    sel_rate = paste0("Sel Rate = ",round(rate*100,1),"%")
    dkl_ABC_LIK <- paste0("KLD(ABC, LIK_FHP) = ",round(KLD_value[[5]], 3))
    average <- paste0("Post Mean = ",round(mean(abc_posterior),3))
    stand_dev <- paste0("Post STD = ",round(std(abc_posterior),3))
    
    ###########
    lines(h$mids,exact_posterior,col='red',type="l", lwd=2)
    seg_max <- max(max(h2$density), max(exact_posterior))
    segments(x0=lambda,y0=0,x1=lambda,y1=seg_max,col='green',lwd=2)
    
    # if (dist %in% c("DIGAD")) {
    #   text(-0.005+min(abc_posterior), ymax*1.30, dkl_ABC_LIK , pos = 4, col = "#330000")
    #   text(-0.005+min(abc_posterior), ymax*1.15, average , pos = 4, col = "#330000")
    #   text(-0.005+min(abc_posterior), ymax*1.00, stand_dev , pos = 4, col = "#330000")
    # } else {
      text(min(abc_posterior), ymax*1.35, dkl_ABC_LIK , pos = 4, col = "#330000")
      text(min(abc_posterior), ymax*1.20, average , pos = 4, col = "#330000")
      text(min(abc_posterior), ymax*1.05, stand_dev , pos = 4, col = "#330000")
    # }
    ###return(hist_gr) 
    mtext(paste0(dist),side=3,line=0,outer=FALSE,las=1, col="darkblue", lwd = 0.5) 
  }
  mtext(paste0("Sel Rate = ",round(rate*100,1),"%"),side=3,line=2,outer=FALSE,las=1, col="darkred")
  mtext(expression(Lambda[0]),side=1,line=0,outer=TRUE,las=0)
  mtext("Density",side=2,line=0,outer=TRUE,las=0)
}

dev.off()
