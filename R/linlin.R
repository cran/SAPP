linlin <- function(external, self.excit, interval, c, d, ax = NULL, ay = NULL,
                   ac = NULL, at = NULL, opt = 0, tmpfile = NULL, nlmax = 1000)
{
  yy <- external       # input data set
  xx <- self.excit     # self-exciting data set
  t <- interval        # length of time interval in which events take place

  x <- c(c,d)
  nn <- length(xx)
  mm <- length(yy)
  if (is.null(ax)) {
    kkx <- 0
    ax <- 0.0
  }  else {
    x <- c(x, ax)
    kkx <- length(ax)  # kxx-1 is the order of lgp trf; xx --> xx
  }
  if (is.null(ay)) {
    kky <- 0
    ay <- 0.0
  } else {
    x <- c(x, ay)
    kky <- length(ay)  # kxy-1 is the oeder of lgp trf; yy --> xx
  }
  if (is.null(ac)) {
    kkc <- 0
    ac <- 0.0
  } else {
    x <- c(x, ac)
    kkc <- length(ac)  # kxz-1 is the order of polynomial for xx data
  }
  if (is.null(at)) {
    kkt <- 0
    at <- 0.0
  } else {
    x <- c(x, at)
    kkt <- length(at)  # kxz-1 is the order of polynomial for xx data
  }
  n <- length(x)

# opt        # =0 : minimize the likelihood with fixed exponential coefficient c
             # =1 :  not fixed d

  kmax <- max(kkx, kky) + 1
  kmax <- max(kmax, 3)
  nlm <- nlmax
  if (is.null(tmpfile))
    nlm <- 0

  z <- .Fortran(C_linlinf,
                as.integer(n),
                as.double(x),
                as.integer(opt),
                as.double(t),
                as.integer(nn),
                as.integer(mm),
                as.double(xx),
                as.double(yy),
                as.integer(kkx),
                as.integer(kky),
                as.integer(kmax),
                as.integer(kkc),
                as.integer(kkt),
                as.integer(nlm),
                x1 = double(n),
                x2 = double(n),
                aic = double(1),
                f = double(1),
                prb = double(1),
                r1 = double(1),
                rwx = double(1),
                rwy = double(1),
                phs = double(1),
                px = double(5*n),
                pg = double(5*n),
                id = integer(nlm),
                rmd = double(nlm),
                eee = double(nlm),
                nl = integer(1),
                ier = integer(1))
	
  ier <- z$ier
  if (ier == -1)
    stop(" subroutine funct : n or kkx or kky kkc kkt error")

# initial estimates
  x1 <- z$x1
  c1 <- x1[1]
  d1 <- x1[2]
  ax1 <- NULL
  ay1 <- NULL
  ac1 <- NULL
  at1 <- NULL
  if (kkx > 0)
    ax1 <- x1[3:(kkx+2)]
  if (kky > 0)
    ay1 <- x1[(kkx+3):(kkx+kky+2)]
  if (kkc > 0)
    ac1 <- x1[(kkx+kky+3):(kkx+kky+kkc+2)]
  if (kkt > 0)
    at1 <- x1[(kkx+kky+kkc+3):n]

# final estimates
  x2 <- z$x2
  c2 <- x2[1]
  d2 <- x2[2]
  ax2 <- NULL
  ay2 <- NULL
  ac2 <- NULL
  at2 <- NULL
  if (kkx > 0)
    ax2 <- x2[3:(kkx+2)]
  if (kky > 0)
    ay2 <- x2[(kkx+3):(kkx+kky+2)]
  if (kkc > 0)
    ac2 <- x2[(kkx+kky+3):(kkx+kky+kkc+2)]
  if (kkt > 0)
    at2 <- x2[(kkx+kky+kkc+3):n]

  if (nlm > 0) {
    nfunct <- 0
    nl <- z$nl
    x <- array(z$px, c(n, 5))
    g <-array(z$pg, c(n, 5))
    id <- z$id[1:nl]
    ramda <- z$rmd[1:nl]
    ee <- z$eee[1:nl]
    print.process(nfunct, n, x, g, id, ramda, ee, tmpfile)
  }

  linlin.out <- list(c1 = c1, d1 = d1, ax1 = ax1, ay1 = ay1, ac1 = ac1,
                    at1 = at1, c2 = c2, d2 = d2, ax2 = ax2, ay2 = ay2,
                    ac2 = ac2, at2=at2, aic2 = z$aic, ngmle = z$f,
                    rayleigh.prob = z$prb, distance = z$r1, rwx = z$rwx,
                    rwy = z$rwy, phase = z$phs)

  class(linlin.out) <- "linlin"
  return(linlin.out)
}

print.linlin <- function(x, ...)
{
  cat(sprintf("\n rayleigh probability = %f\n", x$rayleigh.prob))
  cat(sprintf(" distance = %f\n", x$distance))
  cat(sprintf(" rwx = %f\t\trwy = %f\n", x$rwx, x$rwy))
  cat(sprintf(" phase = %f\n", x$phase))

  kkx <- length(x$ax1)
  kky <- length(x$ay1)
  kkc <- length(x$ac1)
  kkt <- length(x$at1)

  cat("\n initial estimates\n")
  cat(sprintf(" c, d\t%f\t%f\n", x$c1, x$d1))
  if (kkx > 0) {
    pmsg <- " ax"
    for (i in 1:kkx)
      pmsg <- paste(pmsg, gettextf("\t%e", x$ax1[i]))
    cat(pmsg, "\n", sep = "")
  }
  if (kky > 0) {
    pmsg <- " ay"
    for (i in 1:kky)
      pmsg <- paste(pmsg, gettextf("\t%e", x$ay1[i]))
    cat(pmsg, "\n", sep = "")
  }
  if (kkc > 0) {
    pmsg <- " ac"
    for (i in 1:kkc)
      pmsg <- paste(pmsg, gettextf("\t%e", x$ac1[i]))
    cat(pmsg, "\n", sep = "")
  }
  if (kkt > 0) {
    pmsg <- " at"
    for (i in 1:kkt)
      pmsg <- paste(pmsg, gettextf("\t%e", x$at1[i]))
    cat(pmsg, "\n", sep = "")
  }

  cat("\n final outputs\n")
  cat(sprintf(" c, d\t%f\t%f\n", x$c2, x$d2))
  if (kkx > 0) {
    pmsg <- " ax"
    for (i in 1:kkx)
      pmsg <- paste(pmsg, gettextf("\t%e", x$ax2[i]))
    cat(pmsg, "\n", sep = "")
  }
  if (kky > 0) {
    pmsg <- " ay"
    for (i in 1:kky)
      pmsg <- paste(pmsg, gettextf("\t%e", x$ay2[i]))
    cat(pmsg, "\n", sep = "")
  }
  if (kkc > 0) {
    pmsg <- " ac"
    for (i in 1:kkc)
      pmsg <- paste(pmsg, gettextf("\t%e", x$ac2[i]))
    cat(pmsg, "\n", sep = "")
  }
  if (kkt > 0) {
    pmsg <- " at"
    for (i in 1:kkt)
      pmsg <- paste(pmsg, gettextf("\t%e", x$at2[i]))
    cat(pmsg, "\n", sep = "")
  }

  cat(sprintf("\n negative max likelihood = %f\n", x$ngmle))
  cat(sprintf(" AIC/2 = %f\n\n", x$aic2))

}
