etasim1 <- function(bvalue, nd, threshold = 0.0, reference = 0.0, param)
{
# bvlalue                #  b-value of G-R law
# nd                     #  the number of the simulated events
# threshold              #  threshold magnitude
# reference              #  reference magnitude
# param                  #  initial estimates of parameters
  np <- length(param)
  if (np != 5)
    stop("Number of parameters is worng.")
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  p <- param[5]

  ic <- 0
  tstart <- 0
  xm <- rep(0, nd)
  zz <- rep(0, nd)

  z <- .Fortran(C_etasimf,
                as.integer(ic),
                as.double(bvalue),
                as.double(tstart),
                as.integer(nd),
                as.double(threshold),
                as.double(reference),
                as.double(a),
                as.double(b),
                as.double(c),
                as.double(d),
                as.double(p),
                as.double(xm),
                as.double(zz),
                mag = double(nd),
                time = double(nd),
                probx = double(1))

  sim <- matrix(, ncol = 9, nrow = nd)
  sim[,1] <- 1:nd
  sim[,2] <- rep(0, nd)
  sim[,3] <- rep(0, nd)
  sim[,4] <- z$mag
  sim[,5] <- z$time
  sim[,6] <- rep(0, nd)
  sim[,7] <- rep(0, nd)
  sim[,8] <- rep(0, nd)
  sim[,9] <- rep(0, nd)
  sim <- data.frame(sim)
  names(sim) <- c("no.", "longitude", "latitude", "magnitude", "time", "depth",
                  "year", "month", "days")

  return(etasim = sim)
}


etasim2 <- function(etas, tstart, threshold = 0.0, reference = 0.0, param)
{
  n <- dim(etas)[1]
  mag <- etas[, 4]       #  magnitude
  time <- etas[, 5]      #  time from the main shock in days
# tstart                #  
# threshold             #  threshold magnitude
# reference             #  reference magnitude
# param                 #  initial estimates of parameters
  np <- length(param)
  if (np != 5)
    stop("Number of parameters is worng.")
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  p <- param[5]

  nd <- 0
  mag1 <- NULL
  time1 <- NULL
  for (i in 1:n)
      if (mag[i] >= threshold) {
        nd <- nd + 1
        mag1 <- c(mag1, mag[i])
        time1 <- c(time1, time[i])
      }

  ic <- 1
  bvalue <- 0
  xx <- rep(0, nd)
  probx <- 0

  z <- .Fortran(C_etasimf,
                as.integer(ic),
                as.double(bvalue),
                as.double(tstart),
                as.integer(nd),
                as.double(threshold),
                as.double(reference),
                as.double(a),
                as.double(b),
                as.double(c),
                as.double(d),
                as.double(p),
                as.double(mag1),
                as.double(time1),
                mag2 = double(nd),
                time2 = double(nd),
                probx = double(1))

  if (z$probx > 1)
    stop(gettextf("prob = %f > 1", z[[3L]]))

  sim <- matrix(, ncol = 9, nrow = nd)
  sim[,1] <- 1:nd
  sim[,2] <- rep(0, nd)
  sim[,3] <- rep(0, nd)
  sim[,4] <- z$mag2
  sim[,5] <- z$time2
  sim[,6] <- rep(0, nd)
  sim[,7] <- rep(0, nd)
  sim[,8] <- rep(0, nd)
  sim[,9] <- rep(0, nd)
  sim <- data.frame(sim)
  names(sim) <- c("no.", "longitude", "latitude", "magnitude", "time",
                  "depth", "year", "month", "days")

  return(etasim2 = sim)
}

