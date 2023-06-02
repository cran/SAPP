etasap <- function(time, mag, threshold = 0.0, reference = 0.0, parami,
                   zts = 0.0, tstart, zte, approx = 2, tmpfile = NULL,
                   nlmax = 1000, plot = TRUE)
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
# approx                 #  the level for approximation version (1, 2, 4, 8, 16)
                         #  =0 : the exact version

  nn <- 0
  xx <- NULL
  magnitude <- NULL
  amx1 <- threshold
  amx2 <- 10.0

  for (i in 1:n ) 
    if (mag[i] >= amx1 && mag[i] <= amx2)
      if (time[i] >= zts && time[i] <= zte) {
        nn <- nn + 1
        xx <- c(xx, time[i])
        magnitude <- c(magnitude, mag[i])
      }

  nfunct <- 9
  if (approx == 0)
    nfunct <- 4
  nlm <- nlmax
  if (is.null(tmpfile))
    nlm <- 0

  z <- .Fortran(C_etasapf,
                as.double(xx),
                as.double(magnitude),
                as.integer(nn),
                as.double(reference),
                as.double(threshold),
                as.double(parami),
                as.integer(np),
                as.double(zts),
                as.double(zte),
                as.double(tstart),
                as.integer(nfunct),
                as.integer(approx),
                f = double(1),
                x = double(np),
                g = double(np),
                aic2 = double(1),
                id = integer(nlm), 
                ee = double(nlm),
                x1 = double(nlm*np),
                nl = integer(1), 
                as.integer(nlm))

  x <- z$x
  pa <- list(mu = x[1], K = x[2], c = x[3], alpha = x[4], p = x[5])

  if (nlm > 0) {
    nl <- z$nl
    id <- z$id[1:nl]
    g <- z$g
    ee <- z$ee[1:nl]
    x0 <- array(z$x1, c(np, nlmax))
    x1 <- x0[1:np, 1:nl]
    print.process9(id, x, g, ee, x1, tmpfile)
  }

  if (plot == TRUE) {
    mag1 <- min(magnitude)
    ti <- xx
    zero <- rep(0, nn)
    cn <- 1:nn
    xrange <- c(min(xx), max(xx))
    mgrange <- max(cn)/4
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(cn))
    plot(xrange, yrange, type = "n", main = "Seismicity CUM & M-T plot", 
         xlab = "Ordinary Time (Days)", ylab = "Cumulative Number of Events",
         lwd = 1, pty = "s", xaxs = 'r', yaxs = 'i')
    lines(ti, cn, type = 'S')
    mgmax <- max(magnitude - mag1 + 1)
    magnitude1 <- magnitude - mag1 + 0.5
    segments(xx, bottom, xx, magnitude1/mgmax*mgrange+bottom)	
    abline(h = 0)
    abline(h = bottom)
  }

  etasap.out <- list(ngmle = z$f, aic2 = z$aic2, param = pa)
  return(etasap.out)
}


print.process9 <- function(id, x, g, ee, x1, tmpfile)
{
  np <- dim(x1)[1]
  nl <- length(id)

  for (i in 1:nl) {
    k <- id[i] 
      if (k == 8) {
        out <- sprintf(" -ll = %e", ee[i])
        for (j in 1:np)
          out <- paste(out, sprintf("\t%e", x1[j,i]**2))
        write(out, tmpfile, append = TRUE)
      }
    if (k == 330) 
      write(sprintf(" lambda = %e     -LL =%e   %e   %e", ee[i], x1[1,i],
                    x1[2,i], x1[3,i]), tmpfile, append = TRUE) 
    if (k == 340)
      write(sprintf("\n\tinitial (-1)*Log-Likelihood = %e", x1[1,i]), tmpfile, append = TRUE)
    if (k == 600) {
      write("\n ----- x -----", tmpfile, append = TRUE)
      out <- NULL
      for (j in 1:np)
        out <- paste(out, sprintf("  %e", x[j]))
      write(out, tmpfile, append = TRUE)
      write("\n\n *** gradient ***", tmpfile, append = TRUE)
      out <- NULL
      for (j in 1:np)
        out <- paste(out, sprintf("  %e", g[j]**2))
      out <- paste(out, "\n")
      write(out, tmpfile, append = TRUE)
    }
  }
  write(sprintf("................................ Process steps  1- %d ................................\n", nl),
        tmpfile, append = TRUE)
}
