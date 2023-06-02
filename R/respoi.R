respoi <- function(time, mag, param, zts, tstart, zte, threshold = 0.0,
                   plot = TRUE)
{
# time                   #  time from the main shock in days
  nd <- length(time)
# mag                    #  magnitude

# param                  #  estimates of parameters
  np <- length(param) + 1
  if (np != 5)
    stop("Number of parameters is worng.")
  param <- c(param[1:3],0,param[4])

# zts                    #  the start of the precursory period
# zstart                 #  the start of the target period
# zte                    #  the end of the target period

  dep <- rep(0, nd)
  xp <- rep(0, nd)
  yp <- rep(0, nd)

  z <- .Fortran(C_respoif,
                as.double(time),
                as.double(mag),
                as.double(dep),
                as.double(xp),
                as.double(yp),
                as.integer(nd),
                as.double(param),
                as.integer(np),
                as.double(zts),
                as.double(zte),
                as.double(tstart),
                as.double(threshold),
                mag1 = double(nd),
                dep1 = double(nd),
                xp1 = double(nd),
                yp1 = double(nd),
                ntstar = integer(1),
                xx = double(nd),
                x = double(nd),
                nn = integer(1))

  n <- z$nn
  cn <- c(1:n) - z$ntstar
  ti <- z$x[1:n]
  mag1 <- z$mag1[1:n]

  if (plot == TRUE) {
    t <- threshold
    level <- min(0, cn[1])
    xrange <- c(min(ti), max(ti))
    mgrange <- max(cn) / 4
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(max(cn), max(ti)))
    plot(xrange, yrange, type = "n", main = "Omori-Utsu Residual",
         xlab = "Transformed Time", ylab = "Cumulative Number of Events",
         lwd = 1, pty = "s", xaxs = 'r')
    points(ti, 1:length(ti)+cn[1]+1, type = 's') 
    mgmax <- max(mag1 - t + 1)
    mag1 <- mag1 - t + 0.5
    segments(ti, bottom, ti, mag1/mgmax*mgrange+bottom)	
    abline(h = bottom)
    abline(h = 0)
    abline(v = 0, lty = 2)
    abline(0,1, lty = 1, col = 'red')
    s <- tstart
  }

  respoi.out <- list(trans.time = ti, cnum = cn)
  return(respoi.out)
}


respoi2 <- function(etas, param, zts, tstart, zte, threshold = 0.0, plot = TRUE)
{
  if (dim(etas)[2] != 9)
    stop("etas is a sequential data with 9 columns-format.")

# param                  #  estimates of parameters
  np <- length(param) + 1
  if (np != 5)
    stop("Number of parameters is worng.")
  eparam <- c(param[1:3], 0, param[4])

  nd <- dim(etas)[1]
#  xp <- etas[,2]         #  longitude
#  yp <- etas[,3]         #  latitude
#  mag <- etas[,4]        #  magnitude
#  time <- etas[,5]       #  time from the main shock in days
#  dep <- etas[,6]        #  depth
  xp <- rep(0, nd)
  yp <- rep(0, nd)
  mag <- rep(0, nd)
  time <- rep(0, nd)
  dep <- rep(0, nd)
  for (i in 1: nd) {
    xp[i] <- etas[i, 2]
    yp[i] <- etas[i, 3]
    mag[i] <- etas[i, 4]
    time[i] <- etas[i, 5]
    dep[i] <- etas[i, 6]
  }

  z <- .Fortran(C_respoif,
                as.double(time),
                as.double(mag),
                as.double(dep),
                as.double(xp),
                as.double(yp),
                as.integer(nd),
                as.double(eparam),
                as.integer(np),
                as.double(zts),
                as.double(zte),
                as.double(tstart),
                as.double(threshold),
                mag1 = double(nd),
                dep1 = double(nd),
                xp1 = double(nd),
                yp1 = double(nd),
                ntstar = integer(1),
                xx = double(nd),
                x = double(nd),
                nn = integer(1))

  n <- z$nn
  cn <- c(1:n) - z$ntstar
  ti <- z$x[1:n]
  mag1 <- z$mag1[1:n]

  if (plot == TRUE) {
    t <- threshold
    level <- min(0, cn[1])
    xrange <- c(min(ti), max(ti))
    mgrange <- max(cn) / 4
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(max(cn), max(ti)))
    plot(xrange, yrange, type = "n", main = "Omori-Utsu Residual",
         xlab = "Transformed Time", ylab = "Cumulative Number of Events",
         lwd = 1, pty = "s", xaxs = 'r')
    points(ti, 1:length(ti)+cn[1]+1, type = 's') 
    mgmax <- max(mag1 - t + 1)
    mag1 <- mag1 - t + 0.5
    segments(ti, bottom, ti, mag1/mgmax*mgrange+bottom)	
    abline(h = bottom)
    abline(h = 0)
    abline(v = 0, lty = 2)
    abline(0, 1, lty = 1, col = 'red')
    s <- tstart
  }

  res <- matrix(, ncol = 7, nrow = n)
  res[,1] <- cn
  res[,2] <- z$xp1[1:n]
  res[,3] <- z$yp1[1:n]
  res[,4] <- z$mag1[1:n]
  res[,5] <- z$xx[1:n]
  res[,6] <- z$dep1[1:n]
  res[,7] <- ti
  res <- data.frame(res)
  names(res) <- c("no.", "longitude", "latitude", "magnitude", "time", "depth",
                  "trans.time")

  respoi2.out <- list(resData = res)
  return(respoi2.out)

}
