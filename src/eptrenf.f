cc      program eptren
      subroutine eptrenf(xx,t,nn,nfunct,nb,ni,cycle,xa,aic,aicmin,
     & imin,xval,fval,px,g,id,rmd,eee,nl,nmax,np,nlmax)
c
      include 'sapp_f.h'
c
c     this program carries out maximum likelihood estimates of
c   intensity rates of either exponential plynomial or exponential
c   fourier series of non-stationary poisson process models.
c   the optimal oreders are selected automatically by the use of
c   akaike information criterion(aic). a subroutine for plotting
c   a graph of the estimated lines are attached.
c
c     structure of the program
c
c          eptren
c             i---inputs
c             i---reduc1
c             i---reduc2
c             i---davidn----------funct
c             i      i---hesian
c             i      +---linear---funct
c             i---fincal
c             +---output---printr---trenfn
c                             +-----cyclfn
c
c     this program is designed by y. ogata, and programed by y. ogata
c   and k. katsura, inst. statist. math., tokyo. (31/01/85)
c
c     reference of the model, see
c   maclean,c.j.(1974) "estimation and testing og an exponential
c   polynomial rate function within the non-stationary poisson
c   process." biometrika 61, pp.81-86.
c
cx      implicit real * 8 (a-h,o-z)
cc      dimension xx(9000),x(100),aic(20),xa(100,20)
cc      dimension amg(9000)
cc      external funct1,funct2
cx      dimension xx(nn),x(nmax),aic(nb),xa(nmax,nb)
cx      dimension amg(nn)
cx      dimension rx(nmax,nfunct),sx(ni+1,nmax,nfunct)
cx      dimension xval(np),fval(np)
c
cx      dimension px(nmax,nb),g(nmax,nb)
cx      dimension id(nlmax),rmd(nlmax),eee(nlmax)
c
      integer :: nn, nfunct, nb, ni, imin, nlmax, id(nlmax), nl,
     1           nmax, np
      real(8) :: xx(nn), t, cycle, xa(nmax,nb), aic(nb), aicmin,
     1           xval(np), fval(np), px(nmax,nb), g(nmax,nb),
     2           rmd(nlmax), eee(nlmax)
      real(8) :: x(nmax), amg(nn), rx(nmax,nfunct), sx(ni+1,nmax,nfunct)

      nl = 0
cx      do 5 i = 1,nlmax
cx    5 id(i) = 0
      id(1:nlmax) = 0
c
cc      call inputs(xx,amg,t,nn,nfunct,ipl,x,nb,ni,cycle)
cc      if(nfunct.eq.1) call reduc1(t,xx,nn,nb,ni)
cc      if(nfunct.eq.2) call reduc2(t,xx,nn,nb,ni,cycle)
      if(nfunct.eq.1) call reduc1(t,xx,nn,nb,ni,rx,sx,ns)
      if(nfunct.eq.2) call reduc2
     &      (t,xx,nn,nb,ni,cycle,rx(1,1),sx(1,1,1),rx(1,2),sx(1,1,2),ns)
      do 10 i=1,nb
      ii=i
      if(nfunct.eq.2) ii=2*i-1
      do 20 j=1,ii
      x(j)=0.0
   20 continue
cc      if(nfunct.eq.1) call davidn(x,ii,ihes,funct1)
cc      if(nfunct.eq.2) call davidn(x,ii,ihes,funct2)
cx      call davidn1(x,ii,ihes,nfunct,rx,sx,ns,nmax,px(1,i),g(1,i),
      call davidn1(x,ii,nfunct,rx,sx,ns,nmax,px(1,i),g(1,i),
     &               id,rmd,eee,nl,nlmax)
      call fincal(ii,x,aic(i),xa(1,i),t,nfunct)
   10 continue
cc      call output(xx,amg,nn,t,xa,aic,nb,ipl,nfunct,cycle)
      call output(xx,amg,nn,t,xa,aic,nb,nfunct,cycle,aicmin,imin,
     & xval,fval,nmax,np)
      return
      end
cc      subroutine reduc1(t,xx,nn,nb,ni)
      subroutine reduc1(t,xx,nn,nb,ni,rxz,sxz,ns)
c
c     reduction of the data for the (exp)polynomial trend analysis
c
c     input
c           t: length of the observed time interval
c           xx: data of events
c           nn: number of data
c           nb: the largest number of parameters
c           ni: number of subdivisions in (0,t)
c     output
c           rxz:  values of polynomial elements at each events
c           sxz:  values of polynomial at each selected points
c           delta:  length of subdivision by the selected points
CX      implicit real * 8(a-h,o-z)
cc      common /rd1fn1/delta,rxz(20),sxz(4001,20),ns
cc      dimension rxz(20),sxz(4001,20)
      integer :: nn, nb, ni, ns
      real(8) :: t, xx(nn), rxz(nb), sxz(ni+1,nb)
      real(8) :: delta
      common /rd1fn1/delta
cx      dimension rxz(nb),sxz(ni+1,nb)
cx      dimension xx(1)
cx      dimension xx(nn)
      ns=ni
cx      do 10 j=1,nb
      do 11 j=1,nb
      rxz(j)=0.0
      do 10 i=1,nn
      rxz(j)=rxz(j)+(xx(i)/t)**(j-1)
   10 continue
   11 continue
      delta=1.0d0/ni
      sxz(1,1)=1.d0
      do 30 j=2,nb
cx   30 sxz(1,j)=0.0
      sxz(1,j)=0.0
   30 continue
cx      do 20 i=2,ni+1
      do 21 i=2,ni+1
      do 20 j=1,nb
      sxz(i,j)=((i-1)*delta)**(j-1)
   20 continue
   21 continue
      return
      end
cc      subroutine funct1(n,a,f,g,ifg)
      subroutine funct1(n,a,f,g,ifg,rxz,sxz,ns,nmax)
c
c     negative log likelihood function and gradients for log linearly
c     parametrized intensity (exponential polynomial trend)
c
c   input
c        n: number of parameters
c        ns: number of subdivision in (0,1)
c        a: parameters
c        delta: length of subdivisions in (0,1)
c        rxz: values of elements at each events
c        sxz: values of elements at each selected points
c   output
c        f: value of the function
c        g: gradient vecter
c        ifg: index for restrictions
c
cx      implicit real * 8 (a-h,o-z)
cc      common /rd1fn1/delta,rxz(20),sxz(4001,20),ns
      integer :: n, ifg, ns, nmax
      real(8) :: a(n), f, g(n), rxz, sxz
      real(8) :: delta, r, ff, aic, sd
      common /rd1fn1/delta
      common     / ddd /  r , ff , aic , sd
      dimension rxz(nmax),sxz(ns+1,nmax)
cc      dimension g(1),gs(30),a(1)
cx      dimension g(1),gs(n),a(1)
cx      dimension g(n),gs(n),a(n)
      real(8) :: gs(n), fxx, ssxx, rmd, exprm, exprmd
      ifg=0
      fxx=0.0
      do 20 j=1,n
cx   20 fxx=fxx+a(j)*rxz(j)
      fxx=fxx+a(j)*rxz(j)
   20 continue
c
      ssxx=0.0
cx      do 10 j=1,n
cx   10 gs(j)=0.0
      gs(1:n)=0.0
      do 30 i=1,ns+1
      rmd=0.0
      do 40 j=1,n
cx   40 rmd=rmd+a(j)*sxz(i,j)
      rmd=rmd+a(j)*sxz(i,j)
   40 continue
      if(rmd.gt.100.d0) go to 80
      exprm=exp(rmd)
      exprmd=exprm
      if(i.eq.1.or.i.eq.ns+1) exprmd=exprm/2
      ssxx=ssxx+exprmd
      do 50 j=1,n
cx   50 gs(j)=gs(j)+sxz(i,j)*exprmd
      gs(j)=gs(j)+sxz(i,j)*exprmd
   50 continue
   30 continue
      f=-fxx+ssxx*delta
      ff=f
      do 60 j=1,n
cx   60 g(j)=-rxz(j)+gs(j)*delta
      g(j)=-rxz(j)+gs(j)*delta
   60 continue
      return
   80 continue
      f=1.d30
      ifg=1
      return
      end
cc      subroutine reduc2(t,xx,nn,nb,ni,cycle)
      subroutine reduc2(t,xx,nn,nb,ni,cycle,rxc,sxc,rxs,sxs,ns)
c
c   reduction of the data for the exponential fourier trend analysis
c
c     input
c           t: length of the observed time interval
c           cycle: periodicity
c           xx: data of events
c           nn: number of data
c           nb: the largest number of parameters
c           ni: number of subdivisions in (0,cycle)
c     output
c           rxz:  values of fourier elements at each events
c           sxz:  values of fourier at each selected points
c           delta:  length of subdivision by the selected points
c
cx      implicit real * 8(a-h,o-z)
cc      common /rd2fn2/delta,rxc(20),sxc(4001,20),rxs(20),sxs(4001,20),tr,
cc     &               it,ns,nnd
      integer :: nn, nb, ni, ns
      real(8) :: t, xx(nn), cycle, rxc(nb), sxc(ni+1,nb), rxs(nb),
     1           sxs(ni+1,nb)
      real(8) :: delta, tr, pi
      common /rd2fn2/delta,tr,it,nnd
cx      dimension rxc(1),sxc(ni+1,nb*2-1),rxs(1),sxs(ni+1,nb*2-1)
cx      dimension xx(1)
cx      dimension rxc(nb),sxc(ni+1,nb),rxs(nb),sxs(ni+1,nb)
cx      dimension xx(nn)
      data pi/3.14159265358979d0/
      ns=ni
      nnd=nn
c     cycle=365.d0
      it=int(t/cycle)
      tr=t-it*cycle
cx      do 10 k=1,nb
      do 11 k=1,nb
      rxc(k)=0.0
      rxs(k)=0.0
      do 10 i=1,nn
      rxc(k)=rxc(k)+cos(2*pi*k*xx(i)/cycle)
      rxs(k)=rxs(k)+sin(2*pi*k*xx(i)/cycle)
   10 continue
   11 continue
      delta=cycle/ni
cx      do 30 k=2,nb
cx      sxc(1,k)=0.0
cx   30 sxs(1,k)=0.0
      sxc(1,2:nb)=0.0
      sxs(1,2:nb)=0.0
cx      do 20 i=1,ni+1
      do 21 i=1,ni+1
      do 20 k=1,nb
c     sxc(i,k)=cos(2*pi*k*(i-0.5d0)*delta/cycle)
      sxc(i,k)=cos(2*pi*k*(i-1)*delta/cycle)
c     sxs(i,k)=sin(2*pi*k*(i-0.5d0)*delta/cycle)
      sxs(i,k)=sin(2*pi*k*(i-1)*delta/cycle)
   20 continue
   21 continue
      return
      end
cc      subroutine funct2(n,a,f,g,ifg)
      subroutine funct2(n,a,f,g,ifg,rxc,sxc,rxs,sxs,ns,nmax)
c
c     negative log likelihood function and gradients for log linearly
c     parametrized intensity (periodicity)
c
c     input
c        n: number of parameters
c        nnd: number of events
c        ns: number of subdivision in (0,cycle)
c        a: parameters
c        delta: length of subdivisions in (0,cycle)
c        rxz: values of elements at each events
c        sxz: values of elements at each selected points
c     output
c        f: value of the function
c        g: gradient vecter
c        ifg: index for restrictions
c
cx      implicit real * 8 (a-h,o-z)
cc      common /rd2fn2/delta,rxc(20),sxc(4001,20),rxs(20),sxs(4001,20),tr,
cc     &               it,ns,nnd
      integer :: n, ifg, ns, nmax
      real(8) :: a(n), f, g(n), rxc(nmax), sxc(ns+1,nmax),
     1           rxs(nmax), sxs(ns+1,nmax)
      real(8) :: delta, tr, r, ff, aic, sd
      common /rd2fn2/delta,tr,it,nnd
      common     / ddd /  r , ff , aic , sd
cx      dimension rxc(nmax),sxc(ns+1,nmax),rxs(nmax),sxs(ns+1,nmax)
cc      dimension g(20),a(1)
cx      dimension g(n),a(1)
cc      dimension gs(20),gc(20),gcp(20),gsp(20)
cx      dimension g(n),a(n)
cx      dimension gs(n/2),gc(n/2),gcp(n/2),gsp(n/2)
      real(8) :: gs(n/2), gc(n/2), gcp(n/2), gsp(n/2), fxx, ssxx, ssxxp,
     1           rmd, exprm, exprmd
      ifg=0
      fxx=a(1)*nnd
      n2=(n-1)/2
      if(n2.eq.0) go to 25
      do 20 j=1,n2
cx   20 fxx=fxx+a(2*j)*rxc(j)+a(2*j+1)*rxs(j)
      fxx=fxx+a(2*j)*rxc(j)+a(2*j+1)*rxs(j)
   20 continue
   25 continue
c
      ssxx=0.0
      ssxxp=0.0
      g(1)=1.0d0
      if(n2.eq.0) go to 15
cx      do 10 j=1,n2
cx      gc(j)=0.0
cx      gcp(j)=0.0
cx      gs(j)=0.0
cx   10 gsp(j)=0.0
      gc(1:n2)=0.0
      gcp(1:n2)=0.0
      gs(1:n2)=0.0
      gsp(1:n2)=0.0
   15 continue
      do 30 i=1,ns+1
      rmd=a(1)
      if(n2.eq.0) go to 45
      do 40 j=1,n2
cx   40 rmd=rmd+a(2*j)*sxc(i,j)+a(2*j+1)*sxs(i,j)
      rmd=rmd+a(2*j)*sxc(i,j)+a(2*j+1)*sxs(i,j)
   40 continue
   45 continue
      if(rmd.gt.100.d0) go to 80
      exprm=exp(rmd)
      exprmd=exprm
      if(i.eq.1.or.i.eq.ns+1) exprmd=exprm/2
      ssxx=ssxx+exprmd
      if(i*delta.le.tr) ssxxp=ssxx
      if(n2.eq.0) go to 55
      do 50 j=1,n2
      gc(j)=gc(j)+sxc(i,j)*exprmd
      if(i*delta.le.tr) gcp(j)=gc(j)
      gs(j)=gs(j)+sxs(i,j)*exprmd
      if(i*delta.le.tr) gsp(j)=gs(j)
   50 continue
   55 continue
   30 continue
c
      f=-fxx+(it*ssxx+ssxxp)*delta
      g(1)=-nnd+(it*ssxx+ssxxp)*delta
      if(n2.eq.0) go to 65
      do 60 j=1,n2
      g(2*j)=-rxc(j)+(it*gc(j)+gcp(j))*delta
cx   60 g(2*j+1)=-rxs(j)+(it*gs(j)+gsp(j))*delta
      g(2*j+1)=-rxs(j)+(it*gs(j)+gsp(j))*delta
   60 continue
   65 continue
      ff=f
      return
   80 continue
      f=1.d30
      ifg=1
      return
      end
cc      subroutine  davidn( x,n,ihes,funct )
cx      subroutine  davidn1( x,n,ihes,nfunct,rx,sx,ns,nmax,px,g,
      subroutine  davidn1( x,n,nfunct,rx,sx,ns,nmax,px,g,
     &                          id,rmd,eee,nl,nlmax )
c
c     minimization by davidon-fletcher-powell procedure
c
c     this subroutine was copied from timsac 78.
c
c     ----------------------------------------------------------------
c     the following subroutines are directly called by this subroutine
c          funct
c          hesian
c          linear
c     ----------------------------------------------------------------
c     inputs:
c            x: vector of initial values
c            k: dimension of the vector x
c            ihes: =0   inverse of hessian matrix is not available
c                  =1   inverse of hessian matrix is available
c
c     output:
c            x: vector of minimizing solution
c
cx      implicit  real * 8  ( a-h , o-z )
cc      dimension  x(82) , dx(82) , g(82) , g0(82) , y(82)
cc      dimension  h(82,82) , wrk(82) , s(82)
      integer :: n, nfunct, ns, nmax, nlmax, id(nlmax), nl
      real(8) :: x(n), rx(nmax,nfunct), sx(ns+1,nmax,nfunct),
     1           px(n), g(n), rmd(nlmax), eee(nlmax)
      real(8) :: r , f , aic , sd
cx      dimension  x(n) , dx(n) , g(n) , g0(n) , y(n)
cx      dimension  h(n,n) , wrk(n) , s(n)
cx      dimension  rx(nmax,nfunct), sx(ns+1,nmax,nfunct)
      common     / ccc /  isw,ipr
      common     / ddd /  r , f , aic , sd
c
cx      dimension  px(n)
cx      dimension  id(nlmax), rmd(nlmax), eee(nlmax)
      real(8) :: dx(n), g0(n), y(n), h(n,n), wrk(n), s(n),
     1           tau1, tau2, eps1, eps2, ramda, const1, xm,
     2           sum, s1, s2, ss, stem, ds2, gtem, ed, xmb
c
cc      external funct
      data  tau1 , tau2  /  1.0d-5 , 1.0d-5  /
      data  eps1 , eps2  / 1.0d-5 , 1.0d-5  /
      ramda = 0.5d0
      const1 = 1.0d-70
c
c          initial estimate of inverse of hessian
c
      h(1:n,1:n) = 0.0d00
      s(1:n) = 0.0d00
      dx(1:n) = 0.0d00
      do  20   i=1,n
cx      do  10   j=1,n
cx   10 h(i,j) = 0.0d00
cx      s(i) = 0.0d00
cx      dx(i) = 0.0d00
cx   20 h(i,i) = 1.0d00
      h(i,i) = 1.0d00
   20 continue
      isw = 0
c
cc      call  funct( n,x,xm,g,ig )
      if( nfunct.eq.1 )  call  funct1( n,x,xm,g,ig,rx,sx,ns,nmax )
      if( nfunct.eq.2 )  call  funct2( n,x,xm,g,ig,rx(1,1),sx(1,1,1),
     &                                       rx(1,2),sx(1,1,2),ns,nmax )
c
c     write( 6,340 )     xm , sd , aic
c
c          inverse of hessian computation (if available)
c
c     if( ihes .eq. 1 )   call  hesian( x,n,h )
c
      icc = 0
c      iteration
11110 continue
      icc = icc + 1
      do  11111   ic=1,n
      if( ic .eq. 1 .and. icc .eq. 1 )     go to 120
c
      do  40   i=1,n
cx   40 y(i) = g(i) - g0(i)
      y(i) = g(i) - g0(i)
   40 continue
      do  60   i=1,n
      sum = 0.0d00
      do  50   j=1,n
cx   50 sum = sum + y(j) * h(i,j)
      sum = sum + y(j) * h(i,j)
   50 continue
cx   60 wrk(i) = sum
      wrk(i) = sum
   60 continue
      s1 = 0.0d00
      s2 = 0.0d00
      do  70   i=1,n
      s1 = s1 + wrk(i) * y(i)
cx   70 s2 = s2 + dx(i) * y(i)
      s2 = s2 + dx(i) * y(i)
   70 continue
      if( s1.le.const1 .or. s2.le.const1 )  go to 900
      if( s1 .le. s2 )     go to 100
c
c          update the inverse of hessian matrix
c
c               ---  davidon-fletcher-powell type correction  ---
c
cx      do  90   i=1,n
      do  91   i=1,n
      do  90   j=i,n
      h(i,j) = h(i,j) + dx(i)*dx(j)/s2 - wrk(i)*wrk(j)/s1
cx   90 h(j,i) = h(i,j)
      h(j,i) = h(i,j)
   90 continue
   91 continue
      go to  120
c
c               ---  fletcher type correction  ---
c
  100 continue
      stem = s1 / s2 + 1.0d00
cx      do  110   i=1,n
      do  111   i=1,n
      do  110   j=i,n
      h(i,j) = h(i,j)- (dx(i)*wrk(j)+wrk(i)*dx(j)-dx(i)*dx(j)*stem)/s2
cx  110 h(j,i) = h(i,j)
      h(j,i) = h(i,j)
  110 continue
  111 continue
c
c
c
  120 continue
      ss = 0.0d00
      do  150   i=1,n
      sum = 0.0d00
      do  140   j=1,n
cx  140 sum = sum + h(i,j)*g(j)
      sum = sum + h(i,j)*g(j)
  140 continue
      ss = ss + sum * sum
cx  150 s(i) = -sum
      s(i) = -sum
  150 continue
c
c
      s1 = 0.0d00
      s2 = 0.0d00
      do  170   i=1,n
      s1 = s1 + s(i)*g(i)
cx  170 s2 = s2 + g(i)*g(i)
      s2 = s2 + g(i)*g(i)
  170 continue
cc      ds2 = sqrt(s2)
cc      gtem = abs(s1) / ds2
      ds2 = dsqrt(s2)
      gtem = dabs(s1) / ds2
      if( gtem .le. tau1  .and.  ds2 .le. tau2 )     go to  900
      if( s1 .lt. 0.0d00 )     go to  200
      h(1:n,1:n) = 0.0d00
      do  190   i=1,n
cx      do  180   j=1,n
cx  180 h(i,j) = 0.0d00
      h(i,i) = 1.0d00
cx  190 s(i) = -s(i)
      s(i) = -s(i)
  190 continue
  200 continue
c
      ed = xm
c
c          linear  search
c
cc      call  linear( x,s,ramda,ed,n,ig,funct )
      call  linear1( x,s,ramda,ed,n,ig,nfunct,rx,sx,ns,nmax,
     &                id,rmd,eee,nl,nlmax )
c
cc      write( 6,330 )     ramda , f
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 330
         rmd(nl) = ramda
         eee(nl) = f
      end if
c
      s1 = 0.0d00
      do  210   i=1,n
      dx(i) = s(i) * ramda
      s1 = s1 + dx(i) * dx(i)
      g0(i) = g(i)
cx  210 x(i) = x(i) + dx(i)
      x(i) = x(i) + dx(i)
  210 continue
      xmb = xm
      isw = 0
c
cc      call  funct( n,x,xm,g,ig )
      if( nfunct.eq.1)  call  funct1( n,x,xm,g,ig,rx,sx,ns,nmax )
      if( nfunct.eq.2)  call  funct2( n,x,xm,g,ig,rx(1,1),sx(1,1,1),
     &                                       rx(1,2),sx(1,1,2),ns,nmax )
c
      s2 = 0.d0
      do  220     i=1,n
cx  220 s2 = s2 + g(i)*g(i)
      s2 = s2 + g(i)*g(i)
  220 continue
      if( sqrt(s2) .gt. tau2 )   go to  11111
      if( xmb/xm-1.d0 .lt. eps1  .and.  sqrt(s1) .lt. eps2 )  go to 900
11111 continue
      if( icc .ge. 5 )     go to 900
      go to 11110
  900 continue
cc      write( 6,600 )
cc      write( 6,610 )     (x(i),i=1,n)
cc      write( 6,601 )
cc      write( 6,610 )     (g(i),i=1,n)
cx      if( ifg.eq.1 )  then
      do 910 i=1,n
cx  910 px(i) = x(i)
      px(i) = x(i)
  910 continue
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = -1
      end if
      return
cx  330 format( 'lambda =',d15.7,5x,'log likelihood =',d23.15)
cx  340 format( 28x,'log-likelihood =',d23.15,5x,'sd =',d22.15,5x,
cx     *  'aic =',d23.15 )
cx  600 format( /'-----  x  -----' )
cx  601 format( /'***  gradient  ***' )
cx  610 format( 5d13.5 )
      end
c
c
cc      subroutine  linear( x,h,ram,ee,k,ig,funct )
      subroutine  linear1( x,h,ram,ee,k,ig,nfunct,rx,sx,ns,nmax,
     &                     id,rmd,eee,nl,nlmax )
c
c     this subroutine was copied from timsac 78.
c
c     this subroutine performs the linear search along the direction
c     specified by the vector h
c   -----------------------------------------------------------------
c     the following subroutine is directly called by this subroutine:
c         funct
c   -----------------------------------------------------------------
c
c     inputs:
c        x:       vector of position
c        h:       search direction
c        k:       dimension of vector x
c
c     outputs:
c        ram:     optimal step width
c        e2:      minimum function value
c        ig:      error code
c
cx      implicit  real  *8 ( a-h,o-z )
cx      integer  return,sub
cc      dimension  x(100) , h(100) , x1(82)
cc      dimension  g(82)
      integer :: k, ig, nfunct, ns, nmax, nlmax, id(nlmax), nl
      real(8) :: x(k), h(k), ram, ee, rx(nmax,nfunct),
     1           sx(ns+1,nmax,nfunct), rmd(nlmax), eee(nlmax)
cx      dimension  x(k) , h(k) , x1(k)
cx      dimension  g(k)
cx      dimension  rx(nmax,nfunct), sx(ns+1,nmax,nfunct)
      common     / ccc /  isw , ipr
c
cx      dimension  rmd(nlmax), eee(nlmax), id(nlmax)
      integer :: return, sub
      real(8) :: x1(k), g(k), const2, hnorm, e1, e2, e3, ram1,
     1           ram2, ram3, a1, a2, a3, b1, b2 
c
      isw = 1
      ipr = 7
      if( ram .le. 1.0d-30 )  ram = 0.01d0
      const2 = 1.0d-60
      hnorm = 0.d0
      do 10  i=1,k
cx   10 hnorm = hnorm + h(i)**2
      hnorm = hnorm + h(i)**2
   10 continue
      hnorm = sqrt( hnorm )
c
      ram2 = ram
      e1 =ee
      ram1 = 0.d0
c
      do 20  i=1,k
cx   20 x1(i) = x(i) + ram2*h(i)
      x1(i) = x(i) + ram2*h(i)
   20 continue
cc      call  funct( k,x1,e2,g,ig )
      if( nfunct.eq.1 )  call  funct1( k,x1,e2,g,ig,rx,sx,ns,nmax )
      if( nfunct.eq.2 )  call  funct2( k,x1,e2,g,ig,rx(1,1),sx(1,1,1),
     &                                       rx(1,2),sx(1,1,2),ns,nmax )
cc      if(ipr.ge.7)  write(6,2)  ram2,e2
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 2
         rmd(nl) = ram2
         eee(nl) = e2
      end if
c
      if( ig .eq. 1 )  go to  50
      if( e2 .gt. e1 )  go to 50
   30 ram3 = ram2*2.d0
      do 40  i=1,k
cx   40 x1(i) = x(i) + ram3*h(i)
      x1(i) = x(i) + ram3*h(i)
   40 continue
cc      call  funct( k,x1,e3,g,ig )
      if( nfunct.eq.1 )  call  funct1( k,x1,e3,g,ig,rx,sx,ns,nmax )
      if( nfunct.eq.2 )  call  funct2( k,x1,e3,g,ig,rx(1,1),sx(1,1,1),
     &                                       rx(1,2),sx(1,1,2),ns,nmax )
      if( ig.eq.1 )  go to  500
cc      if( ipr.ge.7 )  write(6,3)  ram3,e3
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 3
         rmd(nl) = ram3
         eee(nl) = e3
      end if
      if( e3 .gt. e2 )  go to 70
      ram1 = ram2
      ram2 = ram3
      e1 = e2
      e2 = e3
      go to 30
c
   50 ram3 = ram2
      e3 = e2
      ram2 = ram3*0.1d0
      if( ram2*hnorm .lt. const2 )  go to  400
      do 60  i=1,k
cx   60 x1(i) = x(i) + ram2*h(i)
      x1(i) = x(i) + ram2*h(i)
   60 continue
cc      call  funct( k,x1,e2,g,ig )
      if( nfunct.eq.1 )  call  funct1( k,x1,e2,g,ig,rx,sx,ns,nmax )
      if( nfunct.eq.2 )  call  funct2( k,x1,e2,g,ig,rx(1,1),sx(1,1,1),
     &                                       rx(1,2),sx(1,1,2),ns,nmax )
cc      if(ipr.ge.7)  write(6,4)  ram2,e2
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 4
         rmd(nl) = ram2
         eee(nl) = e2
      end if
c     write(6,8) e1,e2,e3,ee
cx    8 format(1h ,5d20.10)
      if( e2.gt.e1 )  go to 50
c
cc   70 assign 80 to return
   70 return = 80
      go to 200
c
   80 do 90  i=1,k
cx   90 x1(i) = x(i) + ram*h(i)
      x1(i) = x(i) + ram*h(i)
   90 continue
cc      call  funct( k,x1,ee,g,ig )
      if( nfunct.eq.1 )  call  funct1( k,x1,ee,g,ig,rx,sx,ns,nmax )
      if( nfunct.eq.2 )  call  funct2( k,x1,ee,g,ig,rx(1,1),sx(1,1,1),
     &                                       rx(1,2),sx(1,1,2),ns,nmax )
cc      if(ipr.ge.7)  write(6,5)  ram,ee
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 5
         rmd(nl) = ram
         eee(nl) = ee
      end if
c
      ifg = 0
cc      assign  300 to  sub
cc      assign 200 to sub
cc   95 assign 130 to return
      sub = 200
   95 return = 130
      if( ram .gt. ram2 )  go to 110
      if( ee .ge. e2 )  go to 100
      ram3 = ram2
      ram2 = ram
      e3 =e2
      e2 =ee
cc      go to  sub,( 200,300 )
      if( sub .eq. 200 ) go to 200
      if( sub .eq. 300 ) go to 300
c
  100 ram1 = ram
      e1 = ee
cc      go to  sub,( 200,300 )
      if( sub .eq. 200 ) go to 200
      if( sub .eq. 300 ) go to 300
c
  110 if( ee .le. e2 )  go to 120
      ram3 = ram
      e3 = ee
cc      go to  sub,( 200,300 )
      if( sub .eq. 200 ) go to 200
      if( sub .eq. 300 ) go to 300
c
  120 ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = ee
cc      go to  sub,( 200,300 )
      if( sub .eq. 200 ) go to 200
      if( sub .eq. 300 ) go to 300
c
  130 do 140  i=1,k
cx  140 x1(i) = x(i) + ram*h(i)
      x1(i) = x(i) + ram*h(i)
  140 continue
cc      call  funct( k,x1,ee,g,ig )
      if( nfunct.eq.1 )  call  funct1( k,x1,ee,g,ig,rx,sx,ns,nmax )
      if( nfunct.eq.2 )  call  funct2( k,x1,ee,g,ig,rx(1,1),sx(1,1,1),
     &                                       rx(1,2),sx(1,1,2),ns,nmax )
cc      if( ipr.ge.7 )  write(6,6)  ram,ee
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 6
         rmd(nl) = ram
         eee(nl) = ee
      end if
cc      assign 200 to sub
      sub = 200
      ifg = ifg+1
      ifg = 0
      if( ifg .eq. 1 )  go to 95
c
      if( e2 .lt. ee )  ram = ram2
      return
c
c      -------  internal subroutine sub1  -------
  200 a1 = (ram3-ram2)*e1
      a2 = (ram1-ram3)*e2
      a3 = (ram2-ram1)*e3
      b2 = (a1+a2+a3)*2.d0
      b1 = a1*(ram3+ram2) + a2*(ram1+ram3) + a3*(ram2+ram1)
      if( b2 .eq. 0.d0 )  go to 210
      ram = b1 /b2
cc      go to return ,( 80,130 )
      if( return .eq. 80 ) go to 80
      if( return .eq. 130 ) go to 130
c
  210 ig = 1
      ram = ram2
      return
c
c      -------  internal subroutine sub2  -------
c
  300 if( ram3-ram2 .gt. ram2-ram1 )  go to 310
      ram = (ram1+ram2)*0.5d0
cc      go to return ,( 80,130 )
      if( return .eq. 80 ) go to 80
      if( return .eq. 130 ) go to 130
c
  310 ram = (ram2+ram3)*0.5d0
cc      go to return ,( 80,130 )
      if( return .eq. 80 ) go to 80
      if( return .eq. 130 ) go to 130
c ------------------------------------------------------------
c
  400 ram = 0.d0
      return
c ------------------------------------------------------------
c
  500 ram = (ram2+ram3)*0.5d0
  510 do 520  i=1,k
cx  520 x1(i) = x(i) + ram*h(i)
      x1(i) = x(i) + ram*h(i)
  520 continue
cc      call  funct( k,x1,e3,g,ig )
      if( nfunct.eq.1)  call  funct1( k,x1,e3,g,ig,rx,sx,ns,nmax )
      if( nfunct.eq.2)  call  funct2( k,x1,e3,g,ig,rx(1,1),sx(1,1,1),
     &                                       rx(1,2),sx(1,1,2),ns,nmax )
cc      if( ipr.ge.7 )  write(6,7)  ram,e3
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 7
         rmd(nl) = ram
         eee(nl) = e3
      end if
      if( ig.eq.1 )  go to 540
      if( e3.gt.e2 )  go to 530
      ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = e3
      go to 500
c
  530 ram3 = ram
      go to 70
c
  540 ram = (ram2+ram)*0.5d0
      go to 510
c
c ------------------------------------------------------------
cx    1 format( 1h ,'lambda =',d18.10, 10x,'e1 =',d25.17 )
cx    2 format( 1h ,'lambda =',d18.10, 10x,'e2 =',d25.17 )
cx    3 format( 1h ,'lambda =',d18.10, 10x,'e3 =',d25.17 )
cx    4 format( 1h ,'lambda =',d18.10, 10x,'e4 =',d25.17 )
cx    5 format( 1h ,'lambda =',d18.10, 10x,'e5 =',d25.17 )
cx    6 format( 1h ,'lambda =',d18.10, 10x,'e6 =',d25.17 )
cx    7 format( 1h ,'lambda =',d18.10, 10x,'e7 =',d25.17 )
      e n d
      subroutine fincal(n,x,aic,xa,t,nfunct)
cx      implicit real * 8 (a-h,o-z)
      integer :: n, nfunct
      real(8) :: x(n), aic, xa(n), t
      real(8) :: r, xm, aicc, sd
      common     / ddd /  r , xm , aicc , sd
cx      dimension x(1),xa(1)
cx      dimension x(n),xa(n)
      if(nfunct.eq.2) go to 20
      xa(1)=x(1)-log(t)
      if(n.eq.1) go to 40
      do 10 i=2,n
      xa(i)=x(i)/t**(i-1)
   10 continue
   40 continue
      aic=2*xm+2*n
      return
   20 do 30 i=1,n
      xa(i)=x(i)
   30 continue
      aic=2*xm+2*n
      return
      end
cc      subroutine output( xx,amg,nn,t,xa,aic,n,ipl,nfunct,cycle )
      subroutine output( xx,amg,nn,t,xa,aic,n,nfunct,cycle,acmin,
     & imin,xval,fval,nmax,np )
cx      implicit real * 8 (a-h,o-z)
cc      dimension xx(1),amg(1),aic(1),xa(100,20)
cx      dimension xx(1),amg(1),aic(1),xa(nmax,n)
cx      dimension xval(1),fval(1)
      integer :: nn, n, nfunct, imin, nmax, np 
      real(8) :: xx(nn), amg(nn), t, xa(nmax,n), aic(n), cycle, acmin,
     1           xval(np),fval(np)
      real(8) :: tt, xd, x1 
cx      dimension xx(nn),amg(nn),aic(n),xa(nmax,n)
cx      dimension xval(np),fval(np)

      acmin=1.d10
      do 10 i=1,n
cc      write(6,1) aic(i)
cx    1 format(////' a i c ',f15.2)
      ii=i
      if(nfunct.eq.2) ii=2*i-1
cc      write(6,2) (xa(j,i),j=1,ii)
cx    2 format(//' parameters '/(2x,5d13.5))
      if(acmin.lt.aic(i)) go to 10
      acmin=aic(i)
      imin=i
   10 continue
      iimin=imin
      tt=t
      if(nfunct.eq.2) iimin=2*imin-1
      if(nfunct.eq.2) tt=cycle
cc      if(nfunct.eq.1) open(7,file='out.eptrend1')
cc      if(nfunct.eq.1) write(7,*) 'eptren1'
cc      if(nfunct.eq.2) open(7,file='out.epcycle1')
cc      if(nfunct.eq.2) write(7,*) 'epcycle1'
c     x1=1.0
      do 20 i=1,nn
      xd=xx(i)
      if(nfunct.eq.2) xd=xx(i)-int(xx(i)/tt)*tt
      x1=amg(i)
cc      write(7,*) xd,x1
   20 continue
cc      close(7)
cc      call printr(tt,xa(1,imin),iimin,nfunct)
      call printr(tt,xa(1,imin),iimin,xval,fval,nfunct,np)
cc      write(6,3) iimin,(xa(i,imin),i=1,iimin)
cx    3 format(i10/(3d21.13))
      return
      end
      subroutine trenfn(xa,x,y,n)
cx      real * 8 xa(1),yy,x,y
cx      real * 8 xa(n),yy,x,y
      real(8) :: xa(n), yy, x, y
      yy=xa(1)
      if(n.eq.1) go to 20
      do 10 i=2,n
       yy=yy+xa(i)*x**(i-1)
   10 continue
   20 continue
      y=exp(yy)
      return
      end
      subroutine cyclfn(xa,x,y,n)
cx      real * 8 xa(1),yy,x,y,pi
cx      real * 8 xa(n),yy,x,y,pi
      real(8) :: xa(n), yy, x, y, pi
      data pi/3.14159265358979d0/
      yy=xa(1)
      if(n.eq.1) go to 20
      n2=(n-1)/2
      do 10 i=1,n2
       yy=yy+xa(2*i)*cos(2*pi*i*x)+xa(2*i+1)*sin(2*pi*i*x)
   10 continue
   20 continue
cc      y=exp(yy)
      y=dexp(yy)
      return
      end
cc      subroutine printr(t,xa,n,nfunct)
      subroutine printr(t,xa,n,x,y,nfunct,nn)
cc      real * 8 t,xa,xx,yy
cc      dimension x(2000),y(2000),xa(1)
cx      implicit real * 8 (a-h,o-z)
cx      dimension x(1),y(1),xa(1)

cx      dimension x(nn),y(nn),xa(nn)
      integer :: n, nfunct, nn
      real(8) :: t, xa(nn), x(nn), y(nn)
      real(8) :: ymin, ymax, xx, yy
cx      character*1 xo
cx      data xo/'o'/
cc      nn=101
      ymin=0.0
      ymax=0.0
      do 10 i=1,nn
      x(i)=t*(i-1)/(nn-1)
      xx=x(i)
      if(nfunct.eq.2) xx=1.d0*(i-1)/nn
      if(nfunct.eq.1) call trenfn(xa,xx,yy,n)
      if(nfunct.eq.2) call cyclfn(xa,xx,yy,n)
      y(i)=yy
      if(ymin.gt.y(i)) ymin=y(i)
      if(ymax.lt.y(i)) ymax=y(i)
   10 continue
cx      x1=ymin
cx      x2=ymin+(ymax-ymin)*2/4
cx      x3=ymax
cc      write(6,1) x1,x2,x3
cx    1 format(//'x-values  f-values  ',f8.5,2(10x,e11.4)/
cx     &       21x,'+',4('---------+'))
cc      if(nfunct.eq.1) open(7,file='out.eptrend2')
cc      if(nfunct.eq.1) write(7,*) 'eptrend2'
cc      if(nfunct.eq.2) open(7,file='out.epcycle2')
cc      if(nfunct.eq.2) write(7,*) 'epcycle2'
cx      do 20 i=1,nn
cx      l=(y(i)-ymin)*40/(ymax-ymin)+0.5
cx      l=int((y(i)-ymin)*40/(ymax-ymin)+0.5)
cc      if(l.gt.0) write(6,2) x(i),y(i),(xo,j=1,l)
cc      if(l.le.0) write(6,2) x(i),y(i)
cx    2 format(2(1x,d9.4),' i',40a1)
cc      write(7,*) x(i),y(i)
cx   20 continue
cc      close(7)
      return
      end

