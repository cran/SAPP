ptspec <- function(data, nfre, prdmin, prd, nsmooth = 1, pprd, interval,
                   plot = TRUE)
{
# data                #  data of events
  n <- length(data)   #  number of events of data set
# nfre                #  number of sampling frequencies of spectra
  nh1 <- nfre+1
# prdmin              #  the minimum periodicity of the sampling
# prd                 #  a periodicity for calculating the rayleigh probability  
# smooth              #  number for smoothing of periodgram
  is <- nsmooth
# pprd                #  particular periodcities to be investigated among others
  nt <- length(pprd)  #  number of particular periodcities to be investigated
# interval            #  length of observed time interval of events

  z <- .Fortran(C_ptspecf,
                as.double(data),
                as.integer(n),
                as.double(interval),
                as.double(pprd),
                as.double(prdmin),
                as.double(prd),
                as.integer(nfre),
                as.integer(nt),
                as.integer(is),
                prb = double(1),
                r1 = double(1),
                rwx = double(1),
                rwy = double(1),
                phs = double(1),
                wt = double(nt),
                ht = double(nt),
                w = double(nh1),
                h = double(nh1),
                g = double(nh1))

  wt <- z$wt
  ht <- z$ht
  w <- z$w
  h <- z$h
  g <- z$g

  if (plot == TRUE) {
    par(mfrow = c(2, 1))
    pt.spec <- matrix(, nrow = nh1-1, ncol = 2)
    pt.spec[, 1] <- w[2:nh1]
    pt.spec[, 2] <- h[2:nh1]
    ymin <- min(0.0, pt.spec[, 2])

    ram <- nh1 / interval
    sgm <- log(interval / prdmin / pi) - log(0.05)
    sgm1 <- log10(ram * sgm) * 10.e0
    sgm <- log(interval /prdmin / pi) - log(0.1)
    sgm2 <- log10(ram * sgm) * 10.e0
    sgm <- log(interval / prdmin / pi) - log(0.01)
    sgm3 <- log10(ram * sgm) * 10.e0
 
    confidence.bar <- c(sgm1, sgm2, sgm3)
    vertical.bar <- matrix(, nrow = nt, ncol = 2)
    vertical.bar[, 1] <- wt / (2 * pi)
    vertical.bar[, 2] <- ht
    pt.spec[, 1] <- pt.spec[, 1] / (2 * pi)
    plot(pt.spec, pch = 1, cex = 0.5, xlab = "Frequency", ylab = "DB")
    points(vertical.bar, pch = "+", cex = 2, col = 2)
    points(vertical.bar, lty = 2, type = "h")
    abline(h=confidence.bar, lty = 3)

    plot(pt.spec[, 1], 10^(pt.spec[, 2] / 10), pch = 1, cex = 0.5,
        xlab = "Frequency", ylab = "")
    points(vertical.bar[, 1], 10^(vertical.bar[, 2] / 10), pch = "+", cex = 2,
          col = 2)
    points(vertical.bar[, 1], 10^(vertical.bar[, 2] / 10), lty = 2, type = "h")
    abline(h = 10^(confidence.bar / 10), lty = 3)
    par(mfrow = c(1, 1))
  }

  ptspec.out <- list(f = w, db = h, power = g, rayleigh.prob=z$prb,
                   distance = z$r1, rwx = z$rwx, rwy = z$rwy,
                   phase = z$phs, n = n, nfre = nfre, prdmin = prdmin,
                   nsmooth = nsmooth, interval = interval)

  class(ptspec.out) <- "ptspec"
  return(ptspec.out)
}


print.ptspec <- function(x, ...)
{
  om <- 2 * pi / x$prdmin
  cat(sprintf("\n maximum frequency= %e divided into %d points\n", om, x$nfre))

  ave <- x$n / x$interval
  cat(sprintf(" average= %e\n\n", ave))

  cat(sprintf(" rayleigh probability = %f\n", x$rayleigh.prob))
  cat(sprintf(" distance = %f\n", x$distance))
  cat(sprintf(" rwx = %f\trwy = %f\n", x$rwx, x$rwy))
  cat(sprintf(" phase = %f\n", x$phase))

  nh1 <- length(x$f)
  cat("\n frequency\tD.B.\t\tpower\n")
  for(i in 2:nh1)
    cat(sprintf(" %f\t%f\t%f\n", x$f[i], x$db[i], x$power[i]))

  if (x$nsmooth == 1) {
    xi <- x$f / (2 * pi)
    lt4 <- 4 / x$interval
    cat("\n\n list of high powers(>4) with period < t/4\n")
    cat(" frequency\tperiods\tpowers\n")
    for (i in 1:nh1)
      if (xi[i] >= lt4 && x$power[i] >= 4) {
        prd4 <- 1 / xi[i]
        cat(sprintf(" %f\t%f\t%f\n", xi[i],prd4,x$power[i]))
      }
    cat("\n\n")
  }
}
