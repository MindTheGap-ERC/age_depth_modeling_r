######################################################################################################
## This R script is supplied in context of the GLOMAR short course "Age Models and Geochronology: An Introduction to Different Age-depth Modelling Approaches"
## It provides numerical recipes for some widely-adopted techniques for basic age-depth modeling 
## by dr. David De Vleeschouwer and dr. Christian Zeeden
######################################################################################################

######################################################################################################
# The lines until 493 load the R scrips in Supplement to Zeeden et al. 2018, https://doi.org/10.1016/j.quageo.2017.10.001
# exercises/examples start thereafter
# This is not the most common way of getting the functions into R, but for this course with non-R-experts the easiest way
# The function 'recalculate.ages' was adjusted to deal with a zero-age at the top
######################################################################################################

#' Function to create random error part PDF.
#'
#' This function creates a probability density function of the random
#' error part from a set of luminescence data sets
#'
#' @param data \code{List} or \code{data.frame} object, data set to be
#' processed. If more than one data set is used, all data frames must be
#' organised as a list object. Each data set must be of two columns: age
#' and age uncertainty.
#'
#' @param n.MC \code{\link{numeric}} vector, number of Monte Carlo
#' simulations.
#'
#' @param uncertainty.max \code{\link{numeric}} vector, maximum random
#' uncertainty part (%).
#'
#' @param n \code{\link{numeric}} vector, number of output values,
#' i.e. resolution.
#'
#' @param plot \code{\link{logical}} vector, optional plot output.
#'
#' @param ... Optional further parameters that may passed to the plot option.
#'
#' @return \code{Data frame} with probability density function
#'
#' @section Function version: 0.5.1
#'
#' @author Christian Zeeden, RWTH Aachen (Germany),\cr
#' Michael Dietze, GFZ Potsdam (Germany),\cr
#' Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @examples
#'
#' ## create input file list
#' files <- list.files(pattern = "*.txt")
#'
#' ## read input data
#' data <- lapply(X = files, FUN = read.table)
#'
#' ## create and plot pdf
#' pdf <- create.pdf(data = data)
#'
#' ## create pdf only for one data set
#' pdf <- create.pdf(data = data[[1]])
#'
#' ## create pdf assuming the data errors are ten times underestimated
#' pdf <- create.pdf(data = data[[1]],
#'                   uncertainty.max = 1000)
#'
create.pdf <- function(
  
  data,
  n.MC = 10000,
  uncertainty.max = 100,
  n = 100,
  plot = TRUE,
  ...
) {
  
  ## check input data
  if(is.data.frame(data) == TRUE) {
    data <- list(data)
  }
  
  for(i in 1:length(data)) {
    if(is.data.frame(data[[i]]) == TRUE) {
      if(is.numeric(data[[i]][1,1]) == FALSE) {
        stop("Datasets must be numeric throughout!")
      }
    } else if(is.matrix(data[[i]] == TRUE)) {
      if(is.numeric(data[[i]][1,1]) == FALSE) {
        stop("Datasets must be numeric throughout!")
      }
    } else {
      stop("Dataset contains neither data frames nor matrices!")
    }
  }
  
  ## create random error vector
  error.random <- seq(from = 0,
                      to = uncertainty.max,
                      length.out = n)
  
  inversions.observed <- numeric(length = length(x = data))
  
  for(i in 1:length(data)) {
    inversions.observed[i] <- sum(diff(data[[i]][,1]) < 0)
    
  }
  
  ## create result variables
  Result <- matrix(nrow = n,
                   ncol = length(data) + 1)
  
  ## fill first column with percentages of random error
  Result[,1] <- (0:(n - 1)) * uncertainty.max / 100
  
  ## loop through all files
  for(i in 1:length(data)) {
    
    ## create model data array (3-dimensional array)
    M <- array(data = 0,
               dim = c(nrow(x = data[[i]]),
                       n.MC,
                       length(x = error.random)))
    
    ## sort data by age
    data.sorted <- data[[i]][order(data[[i]][,1]),]
    
    ## fill model data array with random samples
    for(j in 1:length(error.random)) {
      M[,,j] <- rnorm(n = nrow(data.sorted) * n.MC,
                      mean = data.sorted[,1],
                      sd = data.sorted[,2]/2 * error.random[j] / 100)
    }
    
    ## find age inversions
    inversions.modelled <-
      M[2:nrow(data[[i]]),,] - M[1:(nrow(data[[i]]) - 1),,] < 0
    
    ## calculate differences of inversions to empirical data set
    inversions.difference <- apply(X = inversions.modelled,
                                   MARGIN = 3,
                                   FUN = colSums) - inversions.observed[i]
    
    ## replace differences with NA  or 1
    inversions.difference <- ifelse(test = inversions.difference == 0,
                                    yes = 1,
                                    no = NA)
    
    ## count number of inversions equal to inversions in model for error.random
    Result[,i+1]  <- colSums(x = inversions.difference,  na.rm = TRUE) /
      sum(inversions.difference, na.rm = TRUE)
  }
  
  ## create output data set
  probability <- Result[,2]
  
  if( i > 1) {
    
    probability <- rowSums(Result[,-1])
  }
  
  pdf <- data.frame(error.random = error.random,
                    probability = probability)
  
  ## optionally plot pdf
  if (plot == TRUE) {
    plot(pdf,
         type = "l",
         xlab = "% Random uncertainty",
         ylab = "Probability")
  }
  
  ## return function output
  return(list(pdf = pdf))
}
###########################################################################
###########################################################################
###########################################################################

#' Function to create random error part PDF.
#'
#' This function recalculates a set of input data (ages and associated
#' errors) by drawing random samples based on the input data but applying
#' the pdf-defined random error fraction. Only stratigraphically
#' consistent cases are kept and used to derive a robust model with mean
#' ages and their user-defined confidence intervals.
#'
#' @param data \code{Numeric} matrix or data frame, ages with errors
#' (2 sigma) to be recalculated.
#'
#' @param pdf \code{Numeric} matrix or data frame, probability density
#' function of the random error part of data. Output part of the function
#' \code{create.pdf()} that contains the pdf values (e.g., \code{pdf$pdf}).
#'
#' @param n.MC.min \code{Numeric} scalar, minimum number of Monte Carlo
#' simulations.
#'
#' @param n.MC.max \code{Numeric} scalar, maximum number of Monte Carlo
#' simulations.
#'
#' @param limits \code{Numeric} vector of length two, upper and lower limit of
#' assumed random error part. This limit will truncate the \code{pdf}
#' argument and might be constrained by laboratory experiences on the random
#' error range. Default is \code{c(0, 1)}, i.e., 0 to 100 %.
#'
#' @param probability = c(0.025, 0.975), \code{Numeric} vector of length two,
#' error probability interval for the modelled age uncertainties.
#'
#' @param verbose \code{Logical} scalar, optional screen messages about
#' model progress.
#'
#' @param plot \code{Logical} scalar, optional plot output.
#'
#' @param ... Optional further parameters that may passed to the plot option.
#'
#' @return \code{Numeric} data frame with recalculated mean ages along with
#' confidence interval boundaries.
#'
#' @section Function version: 0.6.0
#'
#' @author Christian Zeeden, RWTH Aachen (Germany),\cr
#' Michael Dietze, GFZ Potsdam (Germany),\cr
#' Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @examples
#'
#' ## create input file list
#' files <- list.files(pattern = "*.txt")
#'
#' ## read input data
#' data <- lapply(X = files, FUN = read.table)
#'
#' ## create random error pdf
#' pdf <- create.pdf(data = data)
#'
#' ## recalculate second input data set
#' data.ages.recal <- recalculate.ages(data = data[[2]],
#'                                     pdf = pdf$pdf)
#'
#' ## recalculate ages assuming a random error between 90 and 100 %
#' data.ages.recal <- recalculate.ages(data = data[[2]],
#'                                     pdf = pdf$pdf,
#'                                     limits = c(0.9, 1))
#'
recalculate.ages <- function(
  data,
  pdf,
  n.MC.min = 1000,
  n.MC.max = 10000,
  limits = c(0, 1),
  probability = c(0.025, 0.975),
  verbose = TRUE,
  plot = TRUE,
  ...
) {
  
  ## check if dependent library is installed
  if(!require("Hmisc")) {
    stop("Package Hmisc not installed!")
  }
  
  ## check input data
  if(is.matrix(data) == FALSE & is.data.frame(data) == FALSE) {
    stop("Data set must be matrix or data frame!")
  }
  
  if(is.matrix(pdf) == FALSE & is.data.frame(pdf) == FALSE) {
    stop("Data pdf must be matrix or data frame!")
  }
  
  ## remove NA and Inf-values from pdf
  pdf <- pdf[!is.na(pdf[,2]),]
  
  ## optionally truncate pdf by user-defined limits
  limits_pdf <- quantile(x = pdf$error.random, probs = limits)
  pdf <- pdf[pdf$error.random >= limits_pdf[1] &
               pdf$error.random <= limits_pdf[2],]
  
  ## check/set age uncertainty interval
  if(length(probability) != 2 | is.numeric(probability) == FALSE) {
    stop("Error probability must be a numeric vector of length two!")
  }
  
  ## check error probability interval for integrity
  if(probability[1] >= probability[2]) {
    stop("Error probability interval must be positive!")
  }
  
  ## check correct upper limit of random uncertainty part
  if(max(pdf[,1] >100)) {
    stop("Random part of uncertainty must be <=100%")
  }
  
  ## add boundaries to data set in the distance of
  ## 5 * (maximum age - minumim age) / amount data
  ## this is somewhat arbitrarily chosen to be not too close or distant
  ## from the real data. Introducing boundaries reduces the effect of
  ## different dataset lengths, and were introduced in OxCal by C.
  ## Bronk Ramsey
  
  ## calculate mean age difference between samples
  delta <- (max(data[,1]) - min(data[,1])) / nrow(data)
  delta_unc <- (max(data[,2]) - min(data[,2])) / nrow(data)
  
  ## calculate lower boundary age and error
  data.lower <- c(data[1,1] - 5 * delta,
                  (data[1,1] - 5 * delta) * mean(data[,2] / data[,1]))
    
  # set lower boundary to zero if otherwise in the future
  if(data.lower[1] <0)
  {data.lower = c(0,0)}
  
  ## calculate upper boundary age and error
  ### fixed for cases where a date or uncertainty is zero
  data.upper <- c(data[nrow(data),1] + 5 * delta,
                  (data[nrow(data),1] + 5 * delta) *
                    mean(Dates[which(Dates[,1]>0 & Dates[,2]>0),2] / Dates[which(Dates[,1]>0 & Dates[,2]>0),1]))
  
  
  ## append boundary values to data set
  data <- rbind(data.lower,
                data,
                data.upper)
  
  if(data[1,1] == data[2,1]) {data <- data[-1,]}
  
  ## define output data sets
  age.output <- numeric(0)
  weighting.factor <- numeric(0)
  e.random.output <- numeric(0)
  
  ## define number of samples, will be modified in loops
  n.sample <- n.MC.min
  
  ## set counter value
  i <- 1
  
  while(i <= (n.MC.max - 1)) {
    
    ## sample random error part from PDF Loess
    e.random <- sample(x = pdf[,1],
                       prob = pdf[,2],
                       size = n.sample,
                       replace = TRUE) / 100
    
    ## fill matrix with resampled ages
    age.resampled <- matrix(
      data = rnorm(n = nrow(data) * n.sample,
                   mean = rep(x = data[,1],
                              times = n.sample),
                   sd = rep(x = data[,2],
                            times = n.sample) / 2 *
                     e.random),
      nrow = nrow(data),
      ncol = n.sample)
    
    ## infer ordered data sets
    values.ordered <- !apply(X = age.resampled,
                             MARGIN = 2,
                             FUN = is.unsorted,
                             age.resampled)
    ID.ordered <- seq(from = 1,
                      to = length(values.ordered))[values.ordered]
    
    ## calculate weighting factor for ordered data sets
    weighting.factor <-
      c(weighting.factor,
        1 / ((age.resampled[nrow(age.resampled),ID.ordered] -
                age.resampled[1,ID.ordered])^(nrow(data))))
    
    ## write ordered data sets and random errors to output variable
    age.output <- cbind(age.output,
                        age.resampled[,ID.ordered])
    e.random.output <- c(e.random.output,
                         e.random[ID.ordered])
    
    ## update n.MC.min/n.sample
    n.sample <- round(x = n.sample * (n.MC.min - ncol(age.output)) /
                        length(ID.ordered),
                      digits = 0)
    n.sample <- ifelse(n.sample > n.MC.max, n.MC.max, n.sample)
    
    ## optionally print modell progress
    if(verbose == TRUE) {
      print(paste("Ordered data sets: ",
                  ncol(age.output),
                  " (",
                  round(x = ncol(age.output) / n.MC.min * 100,
                        digits = 0),
                  " % of n = ",
                  n.MC.min,
                  "), new sample size: ",
                  n.sample, sep = ""))
    }
    
    ## test if sucessful n.min is reached and update counter
    if(ncol(age.output) >= n.MC.min) {
      
      i <- n.MC.max
    } else {
      
      i <- i + 1
    }
  }
  
  ## get minimum data set size
  l_min <- ifelse(test = ncol(age.output) < n.MC.min,
                  yes = ncol(age.output),
                  no = n.MC.min)
  
  ## truncate output data set to defined length
  age.output <- age.output[,1:l_min]
  e.random.output <- e.random.output[1:l_min]
  e.syst.out <- (1 - e.random.output)
  weighting.factor <- weighting.factor[1:l_min] / min(weighting.factor)
  
  ## calculate weighted confidence envelopes
  ages.improved <- cbind(
    apply(X = age.output,
          MARGIN = 1,
          FUN = wtd.quantile,
          w = weighting.factor,
          probs=0.5),
    apply(X = age.output,
          MARGIN = 1,
          FUN = wtd.quantile,
          w = weighting.factor,
          probs = probability[1]) -
      (100 - 100 * mean(e.random.output)) * data[,2] / 100,
    apply(X = age.output,
          MARGIN = 1,
          FUN = wtd.quantile,
          w = weighting.factor,
          probs = probability[2]) +
      (100 - 100 * mean(e.random.output)) * data[,2] / 100
  )
  
  
  ## define row and column names for output data set
  rownames(ages.improved) <- seq(from = 1,
                                 to = nrow(ages.improved))
  
  colnames(ages.improved) <- c("mean",
                               as.character(probability[1]),
                               as.character(probability[2]))
  
  ## optionally plot output
  if(plot == TRUE) {
    
    ## generate empty plot
    plot(NA,
         xlim = c(1, nrow(ages.improved)),
         ylim = range(c(ages.improved[,2],
                        rev(ages.improved[,3]),
                        data[,1] + data[,2],
                        data[,1] - data[,2])),
         pch = NA,
         main = "Age recalculation",
         xlab = "Sample number",
         ylab = "Age")
    
    ## add model error polygon
    polygon(x = c(1:nrow(data), nrow(data):1),
            y = c(ages.improved[,2], rev(ages.improved[,3])),
            col = "grey",
            lty = 0)
    
    ## add average ages
    lines(x = 1:nrow(data),
          y = ages.improved[,1],
          col = "grey40",
          lwd = 2)
    
    ## add lines between age points
    lines(x = 1:nrow(data),
          y = data[,1],
          lty = 2)
    
    ## add age points
    points(x = 1:nrow(data),
           y = data[,1])
    
    ## add measured error bars
    segments(x0 = 1:nrow(data),
             y0 = data[,1] - data[,2],
             x1 = 1:nrow(data),
             y1 = data[,1] + data[,2])
    
    ## add legend
    legend(x = "topleft",
           legend = c("original", "recalculated"),
           fill = c("NA", "grey"),
           lty = c(1, 0),
           border = c(NA, NA))
  }
  
  ## return function output
  return(as.data.frame(ages.improved))
}

######################################################################################################
# Let's take a look at the GLOMAR dataset
######################################################################################################
dates = matrix(c(0,9,34,46,66,72,0,44.8,116,129,523,797.3,1,10,20,40,120,140), nrow = 6, ncol = 3)
colnames(dates) <- c("Depth (m)", "Age (ka)", "Uncertainty (2sigma, kyr)")
dates = as.data.frame(dates)
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim=c(0,1000),xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}

######################################################################################################
# assess if uncertainty is probably systematic or random
######################################################################################################

PDF_GLOMAR <- create.pdf(as.data.frame(cbind(dates[,2], dates[,3]/2)))[[1]]

PDF_unif <- PDF_GLOMAR
PDF_unif[[2]] <- rep(0.1, 100)

PDF_random <- PDF_GLOMAR
PDF_random[[2]] <- rep(0, 100)
PDF_random[[2]][100] <- 1

######################################################################################################
# run ADMin age-model
######################################################################################################
Dates <- as.data.frame(dates[,2:3])

ADMIn_GLOMAR <- recalculate.ages(Dates, pdf=PDF_GLOMAR, plot=TRUE)
ADMIn_GLOMAR2 <- recalculate.ages(Dates, pdf=PDF_unif, plot=TRUE)
ADMIn_GLOMAR3 <- recalculate.ages(Dates, pdf=PDF_random, plot=TRUE)


# plot a bit nicer, and together with depth scale
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(110,0), xlim=c(0,2000), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ADMIn_GLOMAR[1:6,1], dates$`Depth (m)`, col = "red", cex=2)
points(ADMIn_GLOMAR[1:6,2], dates$`Depth (m)`, col = "red", cex=1)
points(ADMIn_GLOMAR[1:6,3], dates$`Depth (m)`, col = "red", cex=1)

points(ADMIn_GLOMAR[7,1], 100, col = "red", cex=2)
points(ADMIn_GLOMAR[7,2], 100, col = "red", cex=1)
points(ADMIn_GLOMAR[7,3], 100, col = "red", cex=1)


