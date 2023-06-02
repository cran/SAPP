etarpp <- function(time, mag, threshold = 0.0, reference = 0.0, parami,
                   zts = 0.0, tstart, zte, ztend = NULL, plot = TRUE)
{
# time                   #  time from the main shock in days
  n <- length(time)
# mag                    #  magnitude
# threshold              #  threshold magnitude
# reference              #  reference magnitude
# parami                 #  initial estimates of parameters
  np <- length(parami)
  if (np != 5)
    stop("Number of parameters is worng.")
# zts                    #  the start of the precursory period
# tstart                 #  the start of the target period
# zte                    #  the end of the target period
# ztend                  #  the end of the predictionperiod
  if (is.null(ztend))
    ztend <- time[n]

  mm <- 0
  nn <- 0
  time0 <- NULL
  mag0 <- NULL
  time1 <- NULL
  mag1 <- NULL
  for (i in 1:n) 
    if (mag[i] >= threshold) {
      mm <- mm + 1
      time0 <- c(time0, time[i])
      mag0 <- c(mag0, mag[i])
      if (time[i] >= zts && time[i] <= ztend) {
        nn <- nn + 1
        time1 <- c(time1, time[i])
        mag1 <- c(mag1, mag[i])
      }
    }

  z <- .Fortran(C_etarppf,
                as.double(time1),
                as.double(mag1),
                as.double(reference),
                as.integer(nn),
                as.double(parami),
                as.integer(np),
                as.double(zts),
                as.double(ztend),
                as.double(tstart),
                x = double(nn),
                ntst = integer(1))

  x <- z$x
  ntstar <- z$ntst
  cnum <- 1:nn - ntstar

  if (plot == TRUE) {
# r.seisetas
## scan("work.etas")
    mgmin <- min(mag0)
    zero <- rep(0, mm)
    cn <- 1:mm
    cna <- append(cn, 0, after = 0)
    cn1 <- cna[1:mm]
    tia <- append(time0, 0, after = 0)
    ti1 <- tia[1:mm]
    xrange <- c(min(time0), ztend)
    mgrange <- max(cn) / 4
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(max(cn), max(x-ntstar)))
    mtitle <- paste("ETAS Fit and Prediction\nM>=", threshold, "S=", zts, "T=",
                    zte, "Tend=", ztend)
    plot(xrange, yrange, type = "n", main = mtitle,
         xlab = "Ordinary Time (Days)", ylab = "Cumulative Number of Events",
         lwd = 1, xlim = c(xrange), pty = "s", xaxs = "r", yaxs = "r")
## scan("work.res")
    lines(time1, x + ntstar, type = 'l', col = 'red')
    segments(ti1, cn1, time0, cn1)
    segments(time0, cn1, time0, cn)
    mgmax <- max(mag0 - mgmin + 1)
    mag2 <- mag0 - mgmin + 0.5
    segments(time0, bottom, time0, mag2/mgmax*mgrange+bottom)
    abline(h = 0)
    abline(h = bottom)
    mark0 <- tstart; abline(v = mark0, lty = 2)
#   mark1 <- ngmle ; abline(v = mark1, lty = 2)
    mark2 <- ztend;  abline(v = mark2, lty = 2)
    abline(v = tstart, lty = 2)
    abline(v = zte, lty = 2)
#    text(max(time0)*0.3, max(cn)*0.95, paste('M>=',mgmin,'S=',zts,'T=',zte,'Tend=',ztend))

# r.retas
## scan("work.res")
#    X11()
    par(ask = TRUE)
    level <- min(0, cnum[1])
    cn <- 1:mm + cnum[1]
    ti <- x
    xrange <- c(min(ti), max(ti))
    mgrange <- max(cn / 4)
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(max(cn), max(ti)))
    mtitle <- paste("ETAS Residual\nM>=", threshold, "S=", zts, "T=", zte,
                    "Tend=", ztend)
    plot(xrange, yrange, type = "n", main = mtitle, xlab = "Transformed Time",
         ylab = "Cumulative Number of Events", lwd = 1, pty = "s", xaxs = "r")
    points(x, 1:nn+cnum[1]+1, type = 's')
    mgmax <- max(mag1 - threshold + 1)
    mag3 <- mag1 - threshold + 0.5
    segments(ti, bottom, ti, mag3/mgmax*mgrange+bottom)
    abline(h = bottom)
    abline(h = 0)
    abline(v = 0, lty = 2)
    abline(0, 1, lty = 1, col = 'red')
    timax <- max(ti[time1 <= zte])
    abline(v = timax, lty = 2)
#    text(max(ti)*0.4, max(cn)*0.9, paste('M>=',threshold,'S=',zts,'T=',zte,'Tend=',ztend))
    par(ask = FALSE)
  }

  etarpp.out <- list(trans.time = x, no.tstart=ntstar)
  return(etarpp.out)
}


etarpp2 <- function(etas, threshold = 0.0, reference = 0.0, parami, zts= 0.0,
                    tstart, zte, ztend = NULL, plot = TRUE)
{
  n <- dim(etas)[1]
# threshold              #  threshold magnitude
# reference              #  reference magnitude
# parami                 #  initial estimates of parameters
  np <- length(parami)
  if (np != 5)
    stop("Number of parameters is worng.")
# zts                    #  the start of the precursory period
# tstart                 #  the start of the target period
# zte                    #  the end of the target period
# ztend                  #  the end of the predictionperiod
  if (is.null(ztend))
    ztend <- time[n]

  xp <- etas[, 2]        #  longitude
  yp <- etas[, 3]        #  latitude
  mag <- etas[, 4]       #  magnitude
  time <- etas[, 5]      #  time from the main shock in days
  dep <- etas[, 6]       #  depth


  mm <- 0
  nn <- 0
  time0 <- NULL
  mag0 <- NULL
  xp1 <- NULL
  yp1 <- NULL
  mag1 <- NULL
  time1 <- NULL
  dep1 <- NULL
  for (i in 1:n ) 
    if (mag[i] >= threshold) {
      mm <- mm + 1
      time0 <- c(time0, time[i])
      mag0 <- c(mag0, mag[i])
      if (time[i] >= zts && time[i] <= ztend) {
        nn <- nn + 1
        xp1 <- c(xp1, xp[i])
        yp1 <- c(yp1, yp[i])
        mag1 <- c(mag1, mag[i])
        time1 <- c(time1, time[i])
        dep1 <- c(dep1, dep[i])
      }
    }

###  x <- rep(0, nn)
###  ntstar <- 0

  z <- .Fortran(C_etarppf,
                as.double(time1),
                as.double(mag1),
                as.double(reference),
                as.integer(nn),
                as.double(parami),
                as.integer(np),
                as.double(zts),
                as.double(ztend),
                as.double(tstart),
                x = double(nn),
                ntst = integer(1))

  x <- z$x
  ntstar <- z$ntst
  cnum <- 1:nn - ntstar

  if (plot == TRUE) {
# r.seisetas
## scan("work.etas")
    mgmin <- min(mag0)
    zero <- rep(0, mm)
    cn <- 1:mm
    cna <- append(cn, 0, after = 0)
    cn1 <- cna[1:mm]
    tia <- append(time0, 0, after = 0)
    ti1 <- tia[1:mm]
    par(pty = "s", xaxs = "r", yaxs = "r")
    xrange <- c(min(time0), ztend)
    mgrange <- max(cn) / 4
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(max(cn), max(x-ntstar)))
    mtitle <- paste("ETAS Fit and Prediction\nM>=", threshold, "S=", zts, "T=",
                    zte, "Tend=", ztend)
    plot(xrange, yrange, type = "n", main = mtitle,
         xlab = "Ordinary Time (Days)",ylab = "Cumulative Number of Events",
         lwd = 1, xlim = c(xrange))
## scan("work.res")
    lines(time1, x+ntstar, type = 'l', col = 'red')
    segments(ti1, cn1, time0, cn1)
    segments(time0, cn1, time0, cn)
    mgmax <- max(mag0 - mgmin + 1)
    mag2 <- mag0 - mgmin + 0.5
    segments(time0, bottom, time0, mag2/mgmax*mgrange+bottom)
    abline(h = 0)
    abline(h = bottom)
    mark0 <- tstart; abline(v = mark0, lty = 2)
#   mark1 <- ngmle ; abline(v = mark1, lty = 2)
    mark2 <- ztend;  abline(v = mark2, lty = 2)
    abline(v = tstart, lty = 2)
    abline(v = zte, lty = 2)
#    text(max(time0)*0.3, max(cn)*0.95, paste('M>=',mgmin,'S=',zts,'T=',zte,'Tend=',ztend))

# r.retas
## scan("work.res")
#    X11()
    par(ask = TRUE)
    level <- min(0, cnum[1])
    cn <- 1:mm+cnum[1]
    ti <- x
    xrange <- c(min(ti), max(ti))
    mgrange <- max(cn / 4)
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(max(cn), max(ti)))
    par(pty = "s", xaxs = "r")
    mtitle <- paste("ETAS Residual\nM>=", threshold, "S=", zts, "T=", zte,
                    "Tend=", ztend)
    plot(xrange, yrange, type = "n", main = mtitle, xlab = "Transformed Time",
         ylab = "Cumulative Number of Events", lwd = 1)
    points(x, 1:nn+cnum[1]+1, type = 's')
    mgmax <- max(mag1 - threshold + 1)
    mag3 <- mag1 - threshold + 0.5
    segments(ti, bottom, ti, mag3/mgmax*mgrange+bottom)
    abline(h = bottom)
    abline(h = 0)
    abline(v = 0, lty = 2)
    abline(0, 1, lty = 1, col = 'red')
    timax <- max(ti[time1 <= zte])
    abline(v = timax, lty = 2)
#    text(max(ti)*0.4,max(cn)*0.9,paste('M>=',threshold,'S=',zts,'T=',zte,'Tend=',ztend))
    par(ask = FALSE)
  }

  res <- matrix(, ncol = 7, nrow = nn)
  res[, 1] <- cnum
  res[, 2] <- xp1
  res[, 3] <- yp1
  res[, 4] <- mag1
  res[, 5] <- time1
  res[, 6] <- dep1
  res[, 7] <- x
  res <- data.frame(res)
  names(res) <- c("no.", "longitude", "latitude", "magnitude", "time", "depth",
                  "trans.time")

  etarpp2.out <- list(resData = res)
  return(etarpp2.out)
}
