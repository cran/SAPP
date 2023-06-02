
linsim <- function(data, interval, c, d, ax, ay, at, ptmax)
{
  kxx <- length(ax)      # kxx-1 is the order of lgp trf; xx --> xx
  kxy <- length(ay)      # kxy-1 is the oeder of lgp trf; yy --> xx
  kxz <- length(at)      # kxz-1 is the order of polynomial for xx data
  t <- interval          # length of time interval in which events take place
# c,d                    # exponential coefficients in lgp corresponding to xx
                         # and xy, respectively
# ax,ay,at               # coefficients of the corresponding polynomials
  mm <- length(data)
  yy <- c(data, t)
# ptmax                  # ptxmax: an upper bound of trend polynomial
  kmax <- max(kxx, kxy, kxz)
  kmax <- max(kmax, 3)

  z <- .Fortran(C_linsimf,
                as.integer(kxx),
                as.integer(kxy),
                as.integer(kxz),
                as.double(t),
                as.double(c),
                as.double(d),
                as.double(ax),
                as.double(ay),
                as.double(at),
                as.double(yy),
                as.integer(mm),
                as.double(ptmax),
                as.integer(kmax),
                xx = double(2*mm),
                ii1 = integer(1),
                jj1 = integer(1),
                err = double(1),
                ier = integer(1))

  ier <- z$ier
  if (ier != 0)
    stop("Simulated data length is greater than 2*(original data length)")
  err <- z$err
  if (err != 0.)
    stop(gettextf("'ptmax' is incorrect. prob = %f", err), domain = NA)

  simul <- z$xx[1:z$ii1]
  input <- data[1:z$jj1]

  linsim.out <- list(in.data=input, sim.data=simul)
  return(linsim.out)
}
