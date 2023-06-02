pgraph <- function(data, mag, threshold = 0.0, h, npoint, days, delta = 0.0,
                   dmax = 0.0, separate.graphics = FALSE)
{
  nfunct <- 0
  isi <- 0
  n1 <- length(data)     #  total number of variable

  if (delta == 0.0 || dmax == 0.0) {
    dmax <- data[n1] / 4
    delta <- dmax / 100
  }
  td <- data[n1] / delta
  kmax <- as.integer(td/16.0 + 2)

  zd <- rep(0, n1)
  xmg <- rep(0, n1)
  xmg1 <- 0.0
  xmg2 <- 10.0
  nn <- 0
  for (i in 1:n1)
    if (mag[i]==threshold || mag[i]>threshold)
    if (mag[i]==xmg2 || mag[i]<xmg2)
    if (mag[i]==xmg1 || mag[i]>xmg1)
    if (data[i]== 0.0 || data[i]>0.0) {
        nn <- nn + 1
        zd[nn] <- data[i]
        xmg[nn] <- mag[i]
#        if (xmg[nn]== 0.0 ) xmg[nn]=6.0
    }
  zd <- zd[1:nn]
  xmg <- xmg[1:nn]
  nn1 <- nn - 1

  z <- .Fortran(C_pgraphf,
                as.integer(nfunct),
                as.integer(isi),
                as.double(zd),
                as.integer(nn),
                as.integer(npoint),
                as.double(days),
                as.double(h),
                as.double(delta),
                as.double(dmax),
                as.integer(kmax),
                xtau = double(2*nn),
                y = double(2*nn),
                kn = integer(1),
                xl = double(nn1),
                xx = double(nn1*6),
                ydev = double(nn1),
                ui = double(nn1),
                ncn = double(nn1),
                sui = double(nn1),
                xp = double(4),
                xrate = double(npoint+1),
                dlt = double(1),
                xtime = double(kmax),
                yvar = double(kmax*5),
                sigma = double(kmax),
                k = integer(1),
                ier = integer(1))

  kn <- z$kn
  k <- z$k

  tn <- data[n1] / as.double(nn)
  xtau <- z$xtau[1:kn]
  y <- z$y[1:kn]
  xl <- z$xl[1:nn1] * tn
  xx <- array(z$xx, dim = c(nn1, 6)) * tn
  xx <- xx[1:nn1, 1:6]
  ydev <- z$ydev[1:nn1]
  ui <- z$ui[1:nn1]
  ncn <- z$ncn[1:nn1]
  sui <- z$sui[1:nn1]
  xp <- z$xp
  xrate <- z$xrate[1:npoint]
  dlt <- z$dlt
  xtime <- z$xtime[1:k]
  yvar <- array(z$yvar, dim = c(5, kmax))
  yvar <- yvar[1:5, 1:k]
  sigma <- z$sigma[1:k]

  outputPtGraph(zd, xmg, h, kn, xtau, y, xl, xx, ydev, ui, ncn, sui, xp, npoint,
             dlt, xrate, k, xtime, sigma, yvar, separate.graphics)

  pgraph.out <- list(cnum = zd, lintv = xl, tau = xtau, nevent = y,
                     survivor = xx, deviation = ydev, normal.cnum = ncn,
                     normal.lintv = ui, success.intv = sui, occur = xrate,
                     time = xtime, variance = sigma, error = yvar)
  return(pgraph.out)
}


outputPtGraph <- function(zd, xmg, h, kn, xtau, y, xl, xx, ydev, ui, ncn, sui,
                          xp, npoint, dlt, xrate, k, xtime, sigma, yvar,
                          separate.graphics)
{
  nn <- length(zd)
  nn1 <- nn-1
  rnn <- as.double(nn)
  ymax <- nn
  xmax <- zd[nn]
  nband1 <- 1.35810 * sqrt(rnn)
  nband2 <- 1.62762 * sqrt(rnn)

# r.pgCumMT
  par(mfrow=c(2, 1))
  plot(x = zd, y = c(1:nn), type = "s", xlim = c(0, xmax), ylim = c(0, ymax),
       main = "Cumulative Curve & M-T plot",
       xlab = "Time", ylab = "Cumulative Number")
  points(c(0, xmax), c(0, ymax), lty = 2, type = "l")
  points(c(0, xmax), c(0, ymax)+nband1, col = 2, lty = 2, type = "l")
  points(c(0, xmax), c(0, ymax)+nband2, col = 2, lty = 2, type = "l")
  points(c(0, xmax), c(0, ymax)-nband1, col = 2, lty = 2, type = "l")
  points(c(0, xmax), c(0, ymax)-nband2, col = 2, lty = 2, type = "l")
  
  plot(x = zd, y = xmg, type = "h", xlim = c(0, xmax), xlab = "Time",
       ylab = "Magnitude")

# r.pgPTnum
  par(ask = TRUE)
  if (separate.graphics == TRUE)
    dev.new()
  par(mfrow = c(2, 1))
  plot(x = xtau, y = y, type = "l", xlim = c(0, nn),
       xlab = "tau = Time x (Total Number of Events) / (Time End)",
       ylab = "Number of Events in [tau, tau+h]")
  abline(h = seq(-3, 3, 1), lty = 2, col = 2)
  title(main = paste("Normalized number of point in [t, t+h], where h=", h, sep = ""))
  plot(x = zd, y = xmg, type = "h", xlim = c(0, xmax),
       xlab = "Ordinary or Transformed Time", ylab = "Magnitude")

# r.pgSurviv
  if (separate.graphics == TRUE)
    dev.new()
  par(mfrow = c(2, 1))
  yy <- log10(c(nn1:1))
  plot(x = xl, y = yy, type = "p", pch = "+", cex = 0.5,
       main="Survivor Function", xlab = "Interval Length",
       ylab = "Cumulative Number", axes = F)
  for (j in 1:6)
    points(xx[, j], yy, type = "l", col = 2)
  box()
  axis(1)
  n9 <- c(1:9)
  axis(2, at = log10(c(n9, n9*10, n9*10^2, n9*10^3, n9*10^4, 10^5)),
       labels = rep("", 46))
  axis(2, at = log10(c(1, 10, 100, 1000, 10000, 10000)),
       labels = c(expression(10^0), expression(10^1), expression(10^2),
       expression(10^3), expression(10^4), expression(10^5)))

  xx <- c(1:nn1)
  plot(x = xx, y = ydev, type = "n",
       main = "Deviation of Survivor Function from the Poisson",
       xlab = "Interval Length", ylab = "Deviations") 
  points(xx, ydev, type = "p", pch = "+", cex = 0.5)
  abline(h = seq(-3, 3, 1), lty = 2, col = 2)
  abline(h = 0, lty = 1)

# r.pgInterP
# Inter1
  if (separate.graphics == TRUE)
    dev.new()
  par(mfrow = c(2, 1))
  xmax <- ui[1]
  ymax <- ncn[1]
  nband1 <- nband1 / nn1
  nband2 <- nband2 / nn1
  plot(x = ui, y = ncn, type = "l", xlim = c(0, xmax), ylim = c(0, ymax),
       xaxs = "i", yaxs = "i", main = "Interval-Length Distribution",
       xlab = "U(i) = exp{-(Normalized Interval Length)}",
       ylab = "Normalized Cumulative Number")
  points(c(0, xmax), c(0, ymax)+nband1, type = "l", lty=2, col = 2)
  points(c(0, xmax), c(0, ymax)+nband2, type = "l", lty=2, col = 2)
  points(c(0, xmax), c(0, ymax)-nband1, type = "l", lty=2, col = 2)
  points(c(0, xmax), c(0, ymax)-nband2, type = "l", lty=2, col = 2)
  points(c(0, xmax), c(0, ymax), type = "l", lty=2, col = 2)

# Inter2
  data<- matrix(, ncol = 2, nrow = nn1-1)
  data[,1] <- sui[1:(nn1-1)]
  data[,2] <- sui[2:nn1]
  plot(data, type = "p", pch = "+", xaxs = "i", yaxs = "i",
       main="Successive Pair of Intervals", xlab = "U(i)", ylab = "U(i+1)")

# r.pgPalm
  if (separate.graphics == TRUE)
    dev.new()
  par(mfrow = c(1, 1))
  band <- xp
  data <- matrix(, ncol = 2, nrow = npoint)
  data[, 1] <- c(1:npoint) * dlt
  data[, 2] <- xrate
  plot(data, type = "p", pch = "*", main = "Palm Intensity",
       xlab = "Elapsed Time After Each Event", ylab = "Occurrence Rate")
  points(data, type = "l")
  abline(h = band, lty = 2)
  abline(h = mean(band))

# r.pgVTC
  if (separate.graphics == TRUE)
    dev.new()
  plot(x = xtime, y = sigma, type = "p", pch = "+", cex = 0.5,
       ylim = range(pretty(c(sigma, yvar))), main = "Variance-Time Curve",
       xlab = "Time", ylab = "Var{N(0,Time)}")
  points(xtime, yvar[1, ], type = "l", col = 2)
  for (i in 2:5)
    points(xtime, yvar[i, ], type = "l", lty=2, col = 2)
  abline(h = 0)
  abline(v = 0)
  par(ask = FALSE)
}
