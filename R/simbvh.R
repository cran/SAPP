simbvh <- function(interval, axx = NULL, axy = NULL, axz = NULL, ayx = NULL,
                   ayy = NULL, ayz = NULL, c, d, c2, d2, ptxmax, ptymax)
{
  t <- interval          # length of time interval in which events take place

# axx(i),axy(i),ayx(i),ayy(i),axz(i),ayz(i)
#    coefficients of the corresponding polynomials

  if (is.null(axx)) {
    axx <- 0
    kxx <- 0
  } else {
    kxx <- length(axx)     # kxx-1 is the order of lgp trf; xx --> xx
  }
  if (is.null(axy)) {
    axy <- 0
    kxy <- 0
  } else {
    kxy <- length(axy)     # kxy-1 is the oeder of lgp trf; yy --> xx
  }
  if (is.null(ayx)) {
    ayx <- 0
    kyx <- 0
  } else {
    kyx <- length(ayx)     # kyx-1 is the oeder of lgp trf; xx --> yy
  }
  if (is.null(ayy)) {
    ayy <- 0
    kyy <- 0
  } else {
    kyy <- length(ayy)     # kyy-1 is the oeder of lgp trf; yy --> yy
  }
  if (is.null(axz)) {
    axz <- 0
    kxz <- 0
  } else {
    kxz <- length(axz)     # kxz-1 is the order of polynomial for xx data
  }
  if (is.null(ayz)) {
    ayz <- 0
    kyz <- 0
  } else {
    kyz <- length(ayz)     # kyz-1 is the order of polynomial for yy data
  }

# c,d,c2,d2              # exponential coefficients in lgp corresponding to
                         # xx,xy,yx and yy, respectively

# ptxmax                 # upper bounds of trend polynomials corresponding to xz
# ptymax                 # upper bounds of trend polynomials corresponding to yz
  kmax <- max(kxx, kxy, kyx, kyy)
  kmax <- max(kmax, 3)

  nnmax <- 10000
  mmmax <- 10000

  z <- .Fortran(C_simbvhf,
                as.integer(kxx),
                as.integer(kxy),
                as.integer(kxz),
                as.integer(kyx),
                as.integer(kyy),
                as.integer(kyz),
                as.double(t),
                as.double(c),
                as.double(d),
                as.double(c2),
                as.double(d2),
                as.double(axx),
                as.double(axy),
                as.double(axz),
                as.double(ayx),
                as.double(ayy),
                as.double(ayz),
                as.double(ptxmax),
                as.double(ptymax),
                as.integer(kmax),
                x = double(nnmax),
                y = double(mmmax),
                ii1 = integer(1),
                jj1 = integer(1),
                err = double(1),
                ier = integer(1),
                as.integer(nnmax),
                as.integer(mmmax))

  if (z$ier != 0) 
    stop("Simulated data length is greater than 10000.")

  err <- z$err
  if (err != 0.)
    warning(gettextf("Are ptxmax & ptymax correct? prob=%f", err), domain = NA)

  nx <- z$ii1
  ny <- z$jj1
  simbvh.out <- list(x=z$x[1:nx], y=z$y[1:ny] )
  return(simbvh.out)
}









