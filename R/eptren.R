eptren <- function(data, mag = NULL, threshold = 0.0, nparam, nsub, cycle = 0,
                   tmpfile = NULL, nlmax = 1000, plot = TRUE)
{

# data                   #  point process data
  n <- length(data)      #  total number of variable
# mag                    #  magnitude
# threshold              #  threshold magnitude
# nparam                 #  maximum number of parameters
# nsub                   #  number of subdivisions in either (0,t) or (0,cycle)
                         #  for the numerical integration of intensity.
# cycle                  # > 0 : periodicity to be investigated in days
  if (cycle == 0)
    nfunct = 1     #  trend  (exponential polynomial trend fitting)
  if (cycle > 0)
    nfunct = 2     #  cycle  (exponential fourier series fitting)
 
  nn <- 0
  xx <- NULL
  if (is.null(mag)) {
    xx <- data
    magnitude <- data
    nn <- n
    t <- data[nn]
    if (nfunct == 2)
      for (i in 1:nn)
        magnitude[i] <- magnitude[i] - as.integer(data[i]/t) * t
  } else {
    magnitude <- NULL
    for (i in 1:n)
      if ((mag[i] >= threshold) && (data[i] >= 0)) {
        nn <- nn + 1
        xx <- c(xx, data[i])
        magnitude <- c(magnitude, mag[i])
      }
    t <- data[nn]
  }

  nmax <- nparam
  if (nfunct == 2)
    nmax <- 2 * nparam - 1
  np <- 101
  nlm <- nlmax
  if (is.null(tmpfile))
    nlm <- 0

  z <- .Fortran(C_eptrenf,
                as.double(xx),
                as.double(t),
                as.integer(nn),
                as.integer(nfunct),
                as.integer(nparam),
                as.integer(nsub),
                as.double(cycle),
                xa = double(nparam*nmax),
                aic = double(nparam),
                aicm = double(1),
                morder = integer(1),
                xval = double(np),
                fval = double(np),
                x = double(nparam*nmax),
                g = double(nparam*nmax),
                id = integer(nlmax),
                rmd = double(nlmax),
                eee = double(nlmax),
                nl = integer(1),
                as.integer(nmax),
                as.integer(np),
                as.integer(nlm))

  xa <- array(z$xa, dim = c(nmax, nparam))
  param <- list()
  for (i in 1:nparam) {
    m <- i
    if (nfunct == 2)
      m <- m * 2 - 1
    param[[i]] <- xa[1:m, i]
  }
 
  if (plot == TRUE) {
    par(mfrow = c(2, 1))
    data1 <-matrix(, nrow = nn, ncol = 2)  
    data1[, 1] <- xx
    data1[, 2] <- magnitude
    if (is.null(mag))
      data1[1:nn, 2] <- rep(1.0, nn)
    data2 <- matrix(, nrow = np, ncol = 2)
    data2[, 1] <- z$xval
    data2[, 2] <- z$fval
    if (cycle == 0) {
      plot(data2, type = "l", main = "EPTREN-Trend", xlab = "Time",
           ylab = "Intensity Rate")
      plot(data1, type = "h", main = "", xlab = "Time", ylab = "Magnitude")
    } else if (cycle > 0) {
      for (i in 1:nn)
        data1[i, 1] <- data1[i, 1] - as.integer(data1[i, 1] / cycle) * cycle
      plot(data2, type = 'l', main = "EPTREN > EXP(CYCLE)",
           xlab = "Superposed Occurrence Times", ylab = "Intensity Rate")
      plot(data1, type = "h", main = "", xlab = "Superposed Occurrence Times",
           ylab = "Magnitude") 
    }
    par(mfrow = c(1, 1))
  }

  if (nlm > 0) {
    nl <- z$nl
    x <- array(z$x, c(nmax, nparam))
    g <- array(z$g, c(nmax, nparam))
    id <- z$id[1:nl]
    ramda <- z$rmd[1:nl]
    ee <- z$eee[1:nl]
    print.process(nfunct, nparam, x, g, id, ramda, ee, tmpfile)
  }

  eptren.out <- list(aic = z$aic, param = param, aicmin = z$aicm,
                     maice.order = z$morder, time = z$xval, intensity = z$fval)
  class(eptren.out) <- "eptren"
  return(eptren.out)
}

print.eptren <- function(x, ...)
{

  nparam <- length(x$aic)
  for (i in 1:nparam) {
    cat(sprintf("\n AIC\t%f\n", x$aic[i]))
    cat(" parameters\n")
    ii <- length(x$param[[i]])
    pmsg <- NULL
    for (j in 1:ii) {
      pmsg <- paste(pmsg, gettextf("  %e", x$param[[i]][j]))
      if ((j%%5 == 0) || (j == ii)) {
        cat(pmsg, "\n", sep = "")
        pmsg <- NULL
      }
    }
  }

  cat(sprintf("\n\n minimum AIC\t%f\n", x$aicmin))
  cat(" parameters\n")
  i <- x$maice.order
  ii <- length(x$param[[i]])
  pmsg <- NULL
  for (j in 1:ii) {
    pmsg <- paste(pmsg, gettextf("  %e",x$param[[i]][j]))
    if ((j%%5 == 0) || (j == ii)) {
      cat(pmsg, "\n", sep = "")
      pmsg <- NULL
    }
  }
}

print.process <- function(nfunct, n, x, g, id, ramda, ee, tmpfile)
{
  nl <- length(ramda)
  if (nfunct == 0) {
    np <- n
    ist <- 0
  } else if (nfunct == 1) {
    np <- 1
    ist <- 1
  } else if (nfunct == 2) {
    np <- 1
    ist <- 2
  }
  j2 <- 1

  for (i in 1:nl) {
    k <- id[i]
    if (k > 0 & k < 7)
      write(sprintf("lambda = %e     e%i = %e", ramda[i], k, ee[i]), tmpfile,
            append = TRUE);
    if (k == 330)
      write(sprintf("lambda = %e     log likelihood  = %e", ramda[i], ee[i]),
                    tmpfile, append = TRUE);
    if (k == 340)
      write(sprintf("\n                          log-likelihood = %e",
                    ramda[i]), tmpfile, append = TRUE);
    if (k == -1) {
      write("\n ----- x -----", tmpfile, append = TRUE)
      out <- NULL
      for (j in 1:np)
        out <- paste(out, sprintf("  %e", x[j,j2]))
      write(out, tmpfile, append = TRUE)

      write("\n *** gradient ***", tmpfile, append = TRUE)
      out <- NULL
      for (j in 1:np)
        out <- paste(out, sprintf("  %e", g[j,j2]**2))
      out <- paste(out, "\n")
      write(out, tmpfile, append = TRUE)

      np <- np+ist
      j2 <- j2+1
    }
  }
  write(sprintf("................................ Process steps  1- %d ................................\n", nl),
        tmpfile, append = TRUE)
}
