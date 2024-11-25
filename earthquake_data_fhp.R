rm(list=ls())

# Load necessary libraries
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

# Likelihood function  for FRACTIONAL HAWKES PROCESS
compute_logLikelihood_etas <- function(distr, lambda, alpha, gamma, beta, c, distr2, Time){
  library(MittagLeffleR)
  k = length(distr)
  m0 = min(distr2, na.rm = TRUE)
  ll1 = 0
  for (i in 1:k){
    si = 0
    if (i > 1){
      for(j in 1:(i-1)){
        si = si + exp(gamma*(distr2[j]-m0))*c^beta*((distr[i]-distr[j])^(beta-1))*mlf(-c^beta*(distr[i]-distr[j])^beta,beta,beta,1)
      }}
    ll1 = ll1 + log(lambda + alpha*si) - alpha*exp(gamma*(distr2[k]-m0))*(1 - mlf(-c^beta*(distr[k]-distr[i])^beta,beta,1,1))
  }
  ll <- ll1 - lambda*Time
  #  likelihood <- exp(ll)
  return (ll)
  #end function
}

################################################################################
#Read the data
earthquake_data <- read.csv("earthquake_Rwanda.csv", header = TRUE, sep = ",")

# Parse 'time_absolute' as date format and save as 'time_real'
#earthquake_data$time_real <- as.Date(earthquake_data$time_absolute, format = "%m/%d/%Y %H:%M")
#earthquake_data$time_real <- as.Date(earthquake_data$time_absolute, format = "%m/%d/%y %I:%M %p")
earthquake_data$time_real <- as.Date(dmy_hm(earthquake_data$time_absolute))

# Filter
earthquake_data <- earthquake_data %>% filter(mag >= 4.5)

count <- earthquake_data$count
magnitude <- earthquake_data$mag

# Check that 'time_real' is parsed as Date format and no NA values remain
earthquake_data$time_real <- as.Date(earthquake_data$time_real)

#####################
# MAPPING
library(sf)
library(dplyr)

# List of shapefile paths
country_files <- c(
  "../Shapefile_Africa/gadm41_BDI_shp/gadm41_BDI_0.shp",
  "../Shapefile_Africa/gadm41_COD_shp/gadm41_COD_0.shp",
  "../Shapefile_Africa/gadm41_COG_shp/gadm41_COG_0.shp",
  "../Shapefile_Africa/gadm41_RWA_shp/gadm41_RWA_0.shp",
  "../Shapefile_Africa/gadm41_TZA_shp/gadm41_TZA_0.shp",
  "../Shapefile_Africa/gadm41_UGA_shp/gadm41_UGA_0.shp"
)

# Read all shapefiles and combine them & Use quiet = TRUE to suppress messages
region_of_interest <- lapply(country_files, function(file) st_read(file, quiet = TRUE)) %>%
  bind_rows()

#•	[-7.885, 3.382] Latitude
#•	[24.434, 33.794] Longitude

# Define the bounding box
bounding_box <- st_bbox(c(
  xmin = 24.434,  # Minimum longitude
  xmax = 33.794,  # Maximum longitude
  ymin = -7.885,  # Minimum latitude
  ymax = 3.382    # Maximum latitude
), crs = st_crs(4326))  # Use WGS84 coordinate system (EPSG:4326)

# Validate and fix geometry
region_of_interest <- st_make_valid(region_of_interest)

# region_of_interest <- countries
region_of_interest_zoomed <- st_crop(region_of_interest, bounding_box)

# Convert the magnitude data to an sf object
earthquake_sf <- st_as_sf(earthquake_data, coords = c("longitude", "latitude"), crs = 4326)

# Create a rectangle polygon from the bounding_box
bounding_rectangle <- st_as_sfc(bounding_box)

######################
# Replace the long country name with "DRC" for better visualization
region_of_interest_zoomed <- region_of_interest_zoomed %>%
  mutate(COUNTRY = ifelse(COUNTRY == "Democratic Republic of the Congo", "DRC", COUNTRY))

# Calculate centroids for each country in the region_of_interest_zoomed dataset
region_of_interest_zoomed <- region_of_interest_zoomed %>%
  mutate(centroid = st_centroid(geometry),
         longitude = st_coordinates(centroid)[, 1],
         latitude = st_coordinates(centroid)[, 2])

############
# Save the figure in a PDF file
pdf(file = "Fig_earthquake_Rwanda_map.pdf", width = 7.5, height = 7.0)

# Adjust margins to create extra space at the bottom
par(mar = c(6, 4, 4, 2) + 0.1)

# Plot the map with circles sized and colored by earthquake magnitude
ggplot() +
  geom_sf(data = region_of_interest_zoomed, fill = "white", color = "blue") +
  geom_sf(data = bounding_rectangle, color = "black", fill = NA, size = 4) +
  geom_point(data = earthquake_data, 
             aes(x = longitude, y = latitude, size = mag, color = mag), 
             alpha = 0.7) +
  geom_text(data = region_of_interest_zoomed, 
            aes(x = longitude, y = latitude, label = COUNTRY), 
            size = 3, color = "black")+#, fontface = "bold") +  # Add country names
  scale_size_continuous(name = "Magnitude size", range = c(0.1, 2.9)) +  # Adjust circle size range
  scale_color_gradient(name = "Magnitude color", low = "green", high = "red") +  # Custom color gradient
  labs(x = "Longitude",
       y = "Latitude") +
  theme_minimal()

dev.off()


############

# Save the figure in a PDF file
pdf(file = "Fig_earthquake_Rwanda_time_realization.pdf", width = 8.5, height = 4)

# Adjust margins to create extra space at the bottom
par(mar = c(6, 4, 4, 2) + 0.1)

# Define custom labels (ensure same length as breaks)
custom_labels <- seq(as.Date("1974-11-10"), as.Date("2024-11-10"), by = "10 years")  # Adjust based on your data

# Extract the start and end count values for the height of the rectangle
start_count <- earthquake_data$count[earthquake_data$time_real == as.Date("2000-03-03")]
end_count <- earthquake_data$count[earthquake_data$time_real == as.Date("2002-12-13")]

# Plot with jump function
ggplot(earthquake_data, aes(x = time_real, y = count)) +
  geom_step(color = "blue", linewidth = 0.6) + 
  labs(x = "Time", y = "Number of Events") +
  theme_bw() +
  scale_x_date(breaks = seq(as.Date("1974-11-10"), as.Date("2024-11-10"), by = "10 years"),
               labels = custom_labels) +
  annotate("rect", xmin = as.Date("2000-03-03"), xmax = as.Date("2002-12-13"),
           ymin = start_count, ymax = end_count, alpha = 0.2, fill = rgb(1, 0, 0, alpha = 1)) #rgb(1, 0, 0, alpha = 1) Pure red with full opacity


# Close the PDF device
dev.off()

#############################ZOOMING
# Save the figure in a PDF file
pdf(file = "Fig_earthquake_Rwanda_time_realization_zoomed.pdf", width = 8.5, height = 4)

# Adjust margins to create extra space at the bottom
par(mar = c(6, 4, 4, 2) + 0.1)
# Define custom labels (ensure same length as breaks)
custom_labels <- seq(as.Date("2000-03-03"), as.Date("2002-12-13"), by = "6 months")  # Adjust based on your data
library(dplyr)
earthquake_data_zoom <- earthquake_data %>% filter(time_real >= as.Date("2000-03-03") & time_real <= as.Date("2002-12-13"))
earthquake_data_zoom$count <- seq(1, length(earthquake_data_zoom$count))

ggplot(earthquake_data_zoom, aes(x = time_real, y = count)) +
  geom_step(color = "blue", linewidth = 0.6) +
  labs(x = "Time", y = "Number of Events") +
  theme_bw() +
  scale_x_date(breaks = seq(as.Date("2000-03-03"), as.Date("2002-12-13"), by = "6 months"),
               labels = custom_labels)
# Close the PDF device
dev.off()

#################### Earthquake Distribution Using FHP2PP ######################
##### 3.  FHP2PP_LAMBDA
proc <- "FHP2PP"
par <- "Lambda"
dataset <- "RealDataEarthquakeRwanda"

SimPoints <- earthquake_data$time_absolute_months   # Time realizations (in months)
Time <- ceiling(max(SimPoints))

Niter <- 100000 #Number of iterations
lambda <- length(SimPoints)/Time   # Number of earthquakes per month
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
#pdf(file = "Fig_earthquake_Rwanda_distribution_FHP2PP.pdf", width=12, height=4.2)
pdf(file = "Fig_earthquake_Rwanda_distribution_FHP2PP.pdf", width=8, height=4.2)
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
dx <- 1/250  #bin size
LL <- (lambda-lambda) + dx  # lower limit of the breaks interval
UL <- (lambda+lambda) + dx # upper limit of the breaks interval
breaks_interval <- seq(LL-dx,UL+dx,by=dx)
mids_interval <- seq(LL-dx/2, UL+dx/2, by=dx)

## Compute 'logLikelihood'

# Check if the file exists in the current working directory
if (!file.exists(paste0("data_",proc,"",par,"",dataset,"_logLik_ex_post",".rdata"))) {
  
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
  
  # Loglikelihood for SFHP
  logLikelihood2 <- c()
  for (m in 1:length(mids_interval)){
    cat("Case of SFHP: Iteration ", m, " of ", length(mids_interval), " is completed!\n")
    logL <- compute_logLikelihood_etas(lambda = mids_interval[m], distr = SimPoints, distr2 = magnitude, alpha = alpha, gamma = 0.5, beta = beta, c = 1, Time = Time)
    logLikelihood2 <- c(logLikelihood2, logL)
  }
  
  likelihood2 <- exp(logLikelihood2 - mean(logLikelihood2))  # translation of the logLikelihood2 to remove impossibilities
  posterior2 <- likelihood2 # For uniform prior 
  exact_posterior2 <- posterior2
  exact_posterior2 <- exact_posterior2/sum(exact_posterior2)/(dx) #Normalization
  
  # Save data of 'likelihood'
  df2 <- data.frame(logLikelihood, logLikelihood2, exact_posterior, exact_posterior2)
  save(df2, file=paste0("data_",proc,"",par,"",dataset,"_logLik_ex_post",".rdata"))
  cat("Files for 'logsLiks & likelihoods' have been created and saved.\n")
} else {
  # LOAD THE DATA SAVED
  load(file=paste0("data_",proc,"",par,"",dataset,"_logLik_ex_post",".rdata"))
  
  # Extract data
  logLikelihood <- df2$logLikelihood
  logLikelihood2 <- df2$logLikelihood2
  exact_posterior <- df2$exact_posterior
  exact_posterior2 <- df2$exact_posterior2
  cat("Files for 'logsLiks & likelihoods' already exist & extracted.\n")
}

############################
# Compute 'AIC' and 'BIC'
logLik1 <- max(logLikelihood)         # Log-likelihood for FHP
logLik2 <- max(logLikelihood2)        # Log-likelihood for SFHP
n <- length(mids_interval)      # Number of observations
k1 <- 1                         # Number of parameters in FHP
k2 <- 1                         # Number of parameters in SFHP

# Compute AIC and BIC for Model 1
AIC1 <- 2 * k1 - 2 * logLik1
BIC1 <- log(n) * k1 - 2 * logLik1

# Compute AIC and BIC for Model 2
AIC2 <- 2 * k2 - 2 * logLik2
BIC2 <- log(n) * k2 - 2 * logLik2

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
    density_LIK2 <- exact_posterior2
    
    h <- hist(abc_posterior, breaks = breaks_interval, plot = FALSE)
    
    density_ABC <- h$density
    ymax <- max(density_LIK, density_ABC, density_LIK2)
    
    #create plot
    if(dist != "LN"){
      h1 <- hist(abc_posterior, plot = FALSE)#probability=TRUE, plot = FALSE) 
      h2 <- hist(abc_posterior + (lambda-mean(h1$mids)), probability=TRUE,  ylim=c(0,1.40*ceiling(ymax)), breaks = 25, 
                 col="darkolivegreen", border="#333333", xlab="", ylab="", main="")
      abc_posterior <- abc_posterior + (lambda-mean(h1$mids))
    } else {
      h2 <- hist(abc_posterior, probability=TRUE,  ylim=c(0,1.40*ceiling(ymax)), breaks = 25,  
                 col="darkolivegreen", border="#333333", xlab="", ylab="", main="")        
    }
    ###########
    
    if (dist %in% c("MDL50", "DIGAD")) {
      # Create a sequence of breaks with a fixed number of intervals
      num_breaks <- length(exact_posterior)  # Number of desired breaks
      custom_breaks <- seq(min(abc_posterior), max(abc_posterior), length.out = num_breaks + 1) 
      
      # Generate the histogram with custom breaks
      hh <- hist(abc_posterior, breaks = custom_breaks, plot = FALSE)
      density_ABC <- hh$density
    } else {
      # Generate the histogram with custom breaks
      hh <- hist(abc_posterior, breaks = breaks_interval, plot = FALSE)
      density_ABC <- hh$density
    }
    
    ### Penalise zero (0) values of the densities
    eps <- 1e-10
    for (i in 1:length(density_ABC)){
      if (density_ABC[i]==0){
        density_ABC[i]=eps
      }
      if (density_LIK[i]==0){
        density_LIK[i]=eps
      }
      if (density_LIK2[i]==0){
        density_LIK2[i]=eps
      }
    }
    
    ## We can find the KLD using hist densities
    KLD_value <- KLD(density_LIK,density_ABC)  ## KLD(Likelihood, ABC)
    KLD_value2 <- KLD(density_LIK2,density_ABC)  ## KLD(Likelihood2, ABC)
    
    dist_used = paste0("Dist : ",dist)
    threshold_dist = bquote("Threshold "*epsilon == .(round(d,3)))
    sel_rate = paste0("Sel Rate = ",round(rate*100,1),"%")
    dkl_ABC_LIK <- paste0("KLD(ABC, LIK_FHP) = ",round(KLD_value[[5]], 3))
    dkl_ABC_LIK2 <- paste0("KLD(ABC, LIK_SFHP) = ",round(KLD_value2[[5]], 3))
    average <- paste0("Post Mean = ",round(mean(abc_posterior),3))
    stand_dev <- paste0("Post STD = ",round(std(abc_posterior),3))
    
    ###########
   lines(h$mids,exact_posterior2,col='blue3',type="l", lwd=4)
    lines(h$mids,exact_posterior,col='red',type="l", lwd=2)
    seg_max <- max(max(h2$density), max(exact_posterior), max(exact_posterior2))
    segments(x0=lambda,y0=0,x1=lambda,y1=seg_max,col='green',lwd=2)
    
    if (dist %in% c("MDL50", "DIGAD")) {
      text(-0.0015, ymax*1.32, dkl_ABC_LIK , pos = 4, col = "#330000")
      text(-0.0015, ymax*1.19, dkl_ABC_LIK2 , pos = 4, col = "#330000")
      text(-0.0015, ymax*1.06, average , pos = 4, col = "#330000")
      text(-0.0015, ymax*0.93, stand_dev , pos = 4, col = "#330000")
    } else {
      text(-0.0005+min(abc_posterior), ymax*1.36, dkl_ABC_LIK , pos = 4, col = "#330000")
      text(-0.0005+min(abc_posterior), ymax*1.25, dkl_ABC_LIK2 , pos = 4, col = "#330000")
      text(-0.0005+min(abc_posterior), ymax*1.14, average , pos = 4, col = "#330000")
      text(-0.0005+min(abc_posterior), ymax*1.03, stand_dev , pos = 4, col = "#330000")
    }
    ###return(hist_gr) 
    mtext(paste0(dist),side=3,line=0,outer=FALSE,las=1, col="darkblue", lwd = 0.5) 
  }
  mtext(paste0("Sel Rate = ",round(rate*100,1),"%"),side=3,line=2,outer=FALSE,las=1, col="darkred")
  mtext(expression(Lambda[0]),side=1,line=0,outer=TRUE,las=0)
  mtext("Density",side=2,line=0,outer=TRUE,las=0)
}

dev.off()

# Print results
cat("AIC_FHP = ", AIC1, ", AIC_SFHP = ", AIC2, ", Delta_AIC = ", AIC1-AIC2, "\n")
cat("BIC_FHP = ", BIC1, ", BIC_SFHP = ", BIC2, ", Delta_BIC = ", BIC1-BIC2, "\n")
