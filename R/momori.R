
momori <- function(data, mag = NULL, threshold = 0.0, tstart, tend, parami,
                   tmpfile = NULL, nlmax = 1000)
{

# data                   #  point process data
  n <- length(data)      #  total number of variable
# mag                    #  magnitude
# threshold              #  threshold magnitude
# tstart                 #  the start of the target perio
  zts = tstart
  zte=tend
# tend                   #  the end of the target period
# parami                 #  initial estimates of parameters
  np <- length(parami)
  np1 <- np + 1
  if (np1 != 5)
    stop("Number of parameters is worng.")
  parami <- c(parami[1:3], 0, parami[4])

  nn <- 0
  xx <- NULL
  magnitude <- NULL

  if (is.null(mag)) {
    xx <- data
    nn <- n
  } else {
    for (i in 1:n) 
      if (mag[i] >= threshold) 
        if ((data[i] >= zts) && (data[i] <= zte)) {
          nn <- nn + 1
          xx <- c(xx, data[i])
          magnitude <- c(magnitude, mag[i])
        }
  }

  nfunct <- 6
  ncount <- 1
  nlm <- nlmax
  if (is.null(tmpfile))
    nlm <- 1
  kn <- (np-1)/3

  z <- .Fortran(C_momorif,
                as.double(xx),
                as.integer(nn),
                as.double(parami),
                as.integer(np),
                as.double(zts),
                as.double(zte),
                as.integer(ncount),
                as.integer(nfunct),
                ff = double(1),
                x = double(2*np),
                g = double(2*np),
                pa = double(np),
                ahaic = double(ncount),
                t0 = double(1),
                ti = double(kn),
                ak = double(kn),
                c = double(kn),
                p = double(kn),
                cls = double(kn),
                id = integer(nlm),
                rmd = double(nlm),
                x1 = double(np*nlm),
                h = double(2*np*np),
                hf = double(4*np*np),
                nl = integer(1),
                as.integer(nlm))

#  param <- c(t=z$x[1], k=z$x[2], c=z$x[3], p=x0, cls=z$x[4])
  pa1 <- z$pa
  pa2 <- list(t_i = z$t0, K = z$ak, c = z$c, p = z$p, cls = z$cls)

  if (nlm > 1) {
    x <- array(z$x, c(np, 2))
    g <- array(z$g, c(np, 2))
    nl <- z$nl
    id <- z$id[1:nl]
    ramda <- z$rmd[1:nl]
    xx0 <- array(z$x1, c(np, nlm))
    xx1 <- xx0[1:np, 1:nl]
    h <- array(z$h, c(np, np, 2))
    hf <- array(z$hf, c(np, np, 2, 2))
    print.process6(id, x, g, h, hf, ramda, xx1, tmpfile)
  }

  momori.out <- list(param = pa1, ngmle = z$ff, aic = z$ahaic, plist = pa2)
  return(momori.out)
}

print.process6 <- function(id, x, g, h, hf, ramda, xx1, tmpfile)
{
  np <- dim(x)[1]
  nl <- length(id)
  j2 <- 1

  for (i in 1:nl) {
    k <- id[i]
    if (k == 8) {
      out <- sprintf(" -ll = %e", ramda[i])
      for (j in 1:np)
        out <- paste(out, sprintf("\t%e", xx1[j, i]**2))
      write(out, tmpfile, append = TRUE)
    } 
    if (k == 330) 
      write(sprintf(" lambda = %e     -LL =%e   %e   %e", ramda[i],xx1[1,i],xx1[2,i],xx1[3,i]),
                    tmpfile, append = TRUE) 
    if (k == 340)
      write(sprintf("\n\tinitial (-1)*Log-Likelihood = %e", xx1[1,i]), tmpfile, append = TRUE)
    if (k == 600) {
      write("\n ----- x -----", tmpfile, append = TRUE)
      out <- NULL
      for (j in 1:np)
        out <- paste(out, sprintf("  %e", x[j,j2]))
      write(out, tmpfile, append = TRUE)

      write("\n *** gradient ***", tmpfile, append = TRUE)
      out <- NULL
      for (j in 1:np)
        out <- paste(out, sprintf("  %e", g[j,j2]))
      write(out, tmpfile, append = TRUE)

      write("\n *** estimated inverse hessian gradient ***", tmpfile, append = TRUE)
      for (j0 in 1:np) {
         out <- NULL
         for (j1 in 1:np)
           out <- paste(out, sprintf("  %e", h[j0,j1,j2]))
         write(out, tmpfile, append = TRUE)
      }

      write("\n *** fisher matrix ***", tmpfile, append = TRUE)
      for (j0 in 1:np)  {
        out <- NULL
        for (j1 in 1:np)
          out <- paste(out, sprintf("  %e", hf[j0,j1,j2,1]))
        write(out, tmpfile, append = TRUE)
      }

      write("\n *** inverse fisher ***", tmpfile, append = TRUE)
      for (j0 in 1:np)  {
        out <- NULL
        for (j1 in 1:np)
          out <- paste(out, sprintf("  %e", hf[j0,j1,j2,2]))
        out <- paste(out, "\n")
        write(out, tmpfile, append = TRUE)
      }
      j2 <- j2 + 1
    }
  }
  write(sprintf("................................ Process steps  1- %d ................................\n", nl),
        tmpfile, append = TRUE)
}
