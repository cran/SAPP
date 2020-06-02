cc      program etasap
      subroutine etasapf(xx,xmg,nd,xmag0,amx1,xini,n,zts,zte,
     &  tstart0,nfunct0,iappr0,f,x,g,aic2,id,ee,xx1,nl,nlmax)
c
      include 'sapp_f.h'
c-----------------------------------------------------------------------
c Maximum likelihood procedure for the ETAS point process model.
c Subroutine FUNC4 corresponds to the exact likelihood and FUNC9 to the
c approximate version: see one of the references below.
c 
c     References:
c     Ogata, Y. (1988). J. Amer. Statist. Soc., 83, pp. 9-27.
c       ---     (1989). Tectonophysics, 169, pp. 159-174.
c       ---     (1992). J. Geophys. Res. 97, pp. 19845-19871.
c     Ogata, Y., Matsu'ura S.R., Katsura, K. (1993). submitted to 
c                Geophys. Res. Letters.
c-----------------------------------------------------------------------
cx      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=17777, npara=5)
cc      common/xyod/xx(ldata),xmg(ldata),xmag0
cc      common/param/xini(npara),n
cx      dimension xx(nd),xmg(nd)
cx      dimension xini(n)
      integer :: nd, n, nfunct0, iappr0, nlmax, id(nlmax), nl
      real(8) :: xx(nd), xmg(nd), xmag0, amx1, xini(n), zts, zte,
     1           tstart0, f, x(n), g(n), aic2, ee(nlmax), xx1(n,nlmax)
      real(8) :: t, tstart
cxx      common /range/tstart,ntstar
      common /range3/ tstart,ntstar
      common /kkxy/kkx,kky,kkt
cxx      common t,nn,mm,iappr,nfunct
      common /etasap/ t,nn,mm,iappr,nfunct
c
cx      dimension  x(n), g(n)
cx      dimension id(nlmax),ee(nlmax),xx1(n,nlmax)
c
      nl = 0
cx      do 30 i = 1,nlmax
cx   30 id(i) = 0
      id(1:nlmax) = 0
c
cc      call input
cx      call input(xx,xmg,nd,xmag0,amx1,xini,n,zts,zte,tstart0,
cx     &           nfunct0,iappr0)
      call input(xx,xmg,nd,amx1,zts,zte,tstart0,nfunct0,iappr0)
c
cc      do 10 i=2,nn
cc      if(xx(i).ge.xx(i-1)) go to 10
cc      write(6,*) 'reverse occurrence time'
cc      write(6,*) i,xx(i),xx(i-1),xmg(i),xmg(i-1)
cc   10 continue
cc      write(6,8) nfunct
cc    8 format(1h ,' funct = ',i5)
cc      write(6,9) t,nn,mm
cc      write(6,5) xmag0
cx    5 format(1h ,'reference magnitudes; xmag0',5x,f10.4)
cx    3 format(1h ,10f12.4/(1h ,10f12.4))
cx    9 format(1h ,'t,nn,mm',5x,f10.4,2i6)
cx    2 format(f10.2,2i10)
cx    1 format(8f10.2)
c
cc      call finout
      call finout(xx,xmg,xmag0,nn,xini,n,f,x,g,aic2,id,ee,xx1,nl,nlmax)
c
cx   20 continue
      return
      end
c***********************************************************************
cc      subroutine input
cx      subroutine input(xx,xmg,nd,xmag0,amx1,xini,n,zts,zte,tstart0,
      subroutine input(xx,xmg,nd,amx1,zts,zte,tstart0,
     &                 nfunct0,iappr0)
c
c       Reading parameters
c
cc      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=17777, npara=5)
cc      common/param/xini(npara),n
cx      dimension  xini(n)
cc      dimension z(ldata),amg(ldata)
cc      character*80 hypodata
cc      common /xyod/xx(ldata),xmg(ldata),xmag0
cx      dimension  xx(nd),xmg(nd)
      integer :: nd, nfunct0, iappr0
      real(8) :: xx(nd), xmg(nd), amx1, zts, zte, tstart0
      real(8) :: t, tstart, smg, amct, bvl
      common /kkxy/kkx,kky,kkt
cc      common /fukasa/dep(ldata)
cxx      common t,nn,mm,iappr,nfunct
cxx      common /range/tstart,ntstar
      common /etasap/ t,nn,mm,iappr,nfunct
      common /range3/ tstart,ntstar
cc      open(unit=1,file='etas.open')
c     open(unit=1,file='./etaspc.open')
c     open(unit=1,file='siminit.dat')
*     read(1,112) hypodata
cx  112 format(a)
cc      amx2=10.0
cc      read(1,*) nfunct,iappr
cc      read(1,*) zts,zte,tstart
      nfunct=nfunct0
      iappr=iappr0
      tstart=tstart0
cc      read(1,*) amx1,xmag0
cc      write(6,*) nfunct,iappr
cc      write(6,*) zts,zte,tstart
cc      write(6,*) amx1,xmag0
cc      n=5
cc      read(1,*) (xini(i),i=1,n)
cc      close(unit=1)
cc      write(6,*) '      minmag, dep1, ts, tend,tstart'
cc      write(6,5) amx1,zts,zte,tstart
cx    1 format(1h ,20a4)
cx    4 format(8f10.3)
cx    5 format(1h ,8f10.3)
cx    9 format(3i10)
cx    2 format(f10.0,i10)
cx    3 format(20a4)
c
cc      call zisin(t,nd,z,amg,dep,hypodata)
c
      t=zte-zts
      tstart=tstart-zts
      nnn=0
      nn=0
      ntstar=0
      smg=0.0
      do 10 i=1,nd
cc      if(amg(i).lt.amx1.or.amg(i).gt.amx2) go to 10
c     if(dep(i).lt.depx1.or.dep(i).gt.depx2) go to 10
cc      if(z(i).lt.zts.or.z(i).gt.zte) go to 10
      nn=nn+1
cc      write(6,1002) nn,i,amg(i),z(i),dep(i)
cx 1002 format(2i7,4f12.5,f5.0)
cc      if(z(i).lt.tstart) ntstar=nn
cc      xx(nn)=z(i)-zts
      if(xx(i).lt.tstart) ntstar=nn
      xx(nn)=xx(i)-zts
cc      xmg(nn)=amg(i)
cc      dep(nn)=dep(i)
      amct=amx1
      if(xmg(nn).ge.amct) smg=smg+(xmg(nn)-amct+.05)
      if(xmg(nn).ge.amct) nnn=nnn+1
   10 continue
cx  111 format(3(i5,f12.5,f4.1))
      mm=nd
      bvl=nnn/log(10.)/smg
cc      write(6,*) 'read #data; selected #data; #tstart; b-value; M_0'
cc      write(6,*)  nd, nn, ntstar, bvl, amct
      if(zte.ne.0.0) t=zte
      t=zte-zts
      return
      end
c***********************************************************************
cc      subroutine finout
      subroutine finout(xx,xmg,xmag0,ldata,xini,n,ff,x,g,aic20,
     &   id,ee,x1,nl,nlmax)
cx      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=17777,npara=5)
cc      external func4,func9
cc      common/xyod/xx(ldata),xmg(ldata),xmag0
      external func4,func91
cx      dimension xx(ldata),xmg(ldata)
      integer :: ldata, n, nlmax, id(nlmax), nl
      real(8) :: xx(ldata), xmg(ldata), xmag0, xini(n), ff, x(n),
     1           g(n), aic20, ee(nlmax), x1(n,nlmax)
      real(8) :: t, f, aic2
cxx      common t,nn,mm,iappr,nfunct
      common /etasap/ t,nn,mm,iappr,nfunct
cc      common/param/xini(npara),n
cx      dimension xini(n)
cxx      common /ddd/ f,aic2
      common /ddd3/ f,aic2
      common /kkxy/kkx,kky,kkt
cc      dimension x(npara)
cx      dimension x(n)
c
cx      dimension  g(n), id(nlmax), ee(nlmax), x1(n,nlmax)
c
      do 10 i=1,nn
cx   10 xmg(i)=xmg(i)-xmag0
      xmg(i)=xmg(i)-xmag0
   10 continue
c
      do 60 i=1,n
cx   60 x(i)=xini(i)
      x(i)=xini(i)
   60 continue
cc      write(6,1020)   n
cc      write(6,1030)  (x(i),i=1,n)
c
      do 70 i=1,n
cx   70 x(i)=sqrt(x(i))
      x(i)=sqrt(x(i))
   70 continue
c
      do 30 ii=1,1 
c
cc      if(nfunct.eq.4) call davidn(x,n,4,func4)
cc      if(nfunct.eq.9) call davidn(x,n,9,func9)
      if(nfunct.eq.4)
cx     &       call davidn9(xx,xmg,ldata,x,n,4,func4,g,id,ee,x1,nl,nlmax)
     &       call davidn9(xx,xmg,ldata,x,n,func4,g,id,ee,x1,nl,nlmax)
      if(nfunct.eq.9)
cx     &       call davidn9(xx,xmg,ldata,x,n,9,func91,g,id,ee,x1,nl,nlmax)
     &       call davidn9(xx,xmg,ldata,x,n,func91,g,id,ee,x1,nl,nlmax)
c
   30 continue
c
      do 80 i=1,n
cx   80 x(i)=x(i)**2
      x(i)=x(i)**2
   80 continue
c
cc      write(6,1040) f,(x(i),i=1,n)
      ff=f
      aic2=f+n
cx 1005 format(e20.10,3i10)
cc      write(6,1001) aic2
      aic20=aic2
cx 1001 format(1h ,'aic/2 =',e20.10)
cx   20 continue
      return
cx 1000 format(3i10,2f15.6)
cx 1010 format(8f10.4)
cx 1020 format(1h ,'n=',i3)
cx 1030 format(1h ,'x=',6e12.5)
cx 1040 format(1h , 'neg max lklhd=',1 e16.7
cx     3      /' max lklhd est.=',9e12.5/('                 ',9e12.5))
cx 1050 format(4d20.13)
cx 1060 format(e25.15)
      end
c***********************************************************************
cc      subroutine  davidn( x,n,ihes,funct )
cx      subroutine  davidn9( xx,xmg,nn,x,n,ihes,funct,
      subroutine  davidn9( xx,xmg,nn,x,n,funct,
     &                          g,id,ee,xx1,nl,nlmax)
c
c          minimization by davidon-fletcher-powell procedure
c
c       the following subroutines are directly called by this subroutine
c             funct
c             hesian
c             linear
c          inputs:
c             x:       vector of initial values
c             k:       dimension of the vector x
c
c          output:
c             x:       vector of minimizing solution
c
cx      implicit  real * 8  ( a-h , o-z )
      parameter(npara=5)
      external funct
cc      dimension  x(npara) , dx(npara) , g(npara) , g0(npara) , y(npara)
cc      dimension  h(npara,npara) , wrk(npara) , s(npara)
cx      dimension  x(n) , dx(n) , g(n) , g0(n) , y(n)
cx      dimension  h(n,n) , wrk(n) , s(n)
cx      dimension  xx(nn) , xmg(nn)
cx      dimension  id(nlmax), ee(nlmax), xx1(n,nlmax)
      integer :: nn, n, id(nlmax), nl, nlmax
      real(8) :: xx(nn), xmg(nn), x(n), g(n), ee(nlmax), xx1(n,nlmax)
      real(8) :: f, aic2, tau1, tau2, eps1, eps2
cxx      common /ccc/ isw,ipr
cxx      common /ddd/ f,aic2
      common /ddd3/ f,aic2
      data  tau1 , tau2  /  1.0d-5 , 1.0d-5  /
      data  eps1 , eps2  / 1.0d-5 , 1.0d-5  /
      real(8) :: dx(n), g0(n), y(n), h(n,n), wrk(n), s(n), ramda,
     1           const1, xm, sum, s1, s2, stem, ss, ds2, gtem, ed,
     2           xmb
      ramda = 0.5d0
      const1 = 1.0d-70
c
c          initial estimate of inverse of hessian
c
cgk
      iccc=0
cx22221 continue
      iccc=iccc+1
cgk
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
cxx      isw = 0
c
cc      call  funct( n,x,xm,g,ig )
      call funct( xx,xmg,nn,n,x,xm,g,ig )
c
cc      write( 6,340 )     xm
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 340
         xx1(1,nl) = xm
      end if
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
      ds2 = dsqrt(s2)
      gtem = dabs(s1) / ds2
c     write(6,610)gtem,ds2
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
cx      call  linear9( xx,xmg,nn,x,s,ramda,ed,n,ig,funct,lu,ifg )
      call  linear9( xx,xmg,nn,x,s,ramda,ed,n,ig,funct,
     &                id,ee,xx1,nl,nlmax )
c
cc      write( 6,330 )     ramda , ed , s1 , s2
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 330
         ee(nl) = ramda
         xx1(1,nl) = ed
         xx1(2,nl) = s1
         xx1(3,nl) = s2
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
cxx      isw = 0
c
cc      call  funct( n,x,xm,g,ig )
      call funct( xx,xmg,nn,n,x,xm,g,ig )
c
      s2 = 0.d0
      do  220     i=1,n
cx  220 s2 = s2 + g(i)*g(i)
      s2 = s2 + g(i)*g(i)
  220 continue
      if( dsqrt(s2) .gt. tau2 )   go to  11111
      if( xmb/xm-1.d0 .lt. eps1  .and.  dsqrt(s1) .lt. eps2 )  go to 900
11111 continue
      if( icc .ge. 5 )     go to 900
      go to 11110
  900 continue
cc      write( 6,600 )
cc      write( 6,610 )     (x(i)**2,i=1,n)
cc      write( 6,601 )
cc      write( 6,610 )     (g(i)**2,i=1,n)
         if( nl.lt.nlmax ) then
            nl = nl+1
            id(nl) = 600
         end if
      return
cx  330 format( 1h ,'lambda =',d15.7,5x,'-LL =',d23.15,2x,d9.2,2x,d9.2)
cx  340 format( 1h ,4x,'initial (-1)*Log-Likelihood =',d23.15)
cx  600 format( 1h ,'-----  x  -----' )
cx  601 format( 1h ,'***  gradient  ***' )
cx  610 format( 1h ,10d13.5 )
      end
c***********************************************************************
cc      subroutine  linear( x,h,ram,ee,k,ig,funct )
      subroutine  linear9( xx,xmg,nn,x,h,ram,ee,k,ig,funct,
     &                         id,eee,xx1,nl,nlmax )
c
c     this subroutine performs the linear search along the direction spe
c     by the vector h
c----------------------------------------------------------------------
c       the following subroutine is directly called by this subroutine:
c             funct
c----------------------------------------------------------------------
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
cc      parameter(npara=5)
      external funct
cx      integer  return,sub
cc      dimension  x(npara) , h(npara) , x1(npara)
cc      dimension  g(npara)
      integer :: nn, k, ig, id(nlmax), nl, nlmax
      real(8) :: xx(nn), xmg(nn), x(k), h(k), ram, ee, eee(nlmax),
     1           xx1(k,nlmax)
cx      dimension  x(k) , h(k) , x1(k)
cx      dimension  g(k)
cxx      common /ccc/ isw,ipr
cx      dimension  xx(nn) , xmg(nn)
cx      dimension  id(nlmax), eee(nlmax), xx1(k,nlmax)
      integer :: return, sub
      real(8) :: x1(k), g(k), const2, hnorm, ram1, ram2, ram3,
     1           e1, e2, e3, a1, a2, a3, b1, b2
c
cxx      isw = 1
cxx      ipr=7
      if( ram .le. 1.0d-30 )  ram = 0.01d0
      const2 = 1.0d-60
      hnorm = 0.d0
      do 10  i=1,k
cx   10 hnorm = hnorm + h(i)**2
      hnorm = hnorm + h(i)**2
   10 continue
      hnorm = dsqrt( hnorm )
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
      call  funct( xx,xmg,nn,k,x1,e2,g,ig )
c     if(ipr.ge.7)  write(6,2)  ram2,e2
cc      if(ipr.ge.7)  write(6,8)  e2,(x1(i)**2,i=1,k)
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         eee(nl) = e2
         do 21 i=1,k
cx   21    xx1(i,nl) = x1(i)
         xx1(i,nl) = x1(i)
   21    continue
      end if

cx    8 format(' -ll=',d13.5,1x,5d12.5)
c
      if( ig .eq. 1 )  go to  50
      if( e2 .gt. e1 )  go to 50
   30 ram3 = ram2*2.d0
      do 40  i=1,k
cx   40 x1(i) = x(i) + ram3*h(i)
      x1(i) = x(i) + ram3*h(i)
   40 continue
cc      call  funct( k,x1,e3,g,ig )
      call  funct( xx,xmg,nn,k,x1,e3,g,ig )
      if( ig.eq.1 )  go to  500
c     if( ipr.ge.7 )  write(6,3)  ram3,e3
cc      if(ipr.ge.7)  write(6,8)  e3,(x1(i)**2,i=1,k)
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         eee(nl) = e3
         do 41 i=1,k
cx   41    xx1(i,nl) = x1(i)
         xx1(i,nl) = x1(i)
   41    continue
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
      call  funct( xx,xmg,nn,k,x1,e2,g,ig )
c     if(ipr.ge.7)  write(6,4)  ram2,e2
cc      if(ipr.ge.7)  write(6,8)  e2,(x1(i)**2,i=1,k)
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         eee(nl) = e2
         do 61 i=1,k
cx   61    xx1(i,nl) = x1(i)
         xx1(i,nl) = x1(i)
   61    continue
      end if

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
      call  funct( xx,xmg,nn,k,x1,ee,g,ig )
c     if(ipr.ge.7)  write(6,5)  ram,ee
cc      if(ipr.ge.7)  write(6,8)  ee,(x1(i)**2,i=1,k)
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         eee(nl) = ee
         do 91 i=1,k
cx   91    xx1(i,nl) = x1(i)
         xx1(i,nl) = x1(i)
   91    continue
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
      if( sub.eq.200 ) go to 200
      if( sub.eq.300 ) go to 300
c
  100 ram1 = ram
      e1 = ee
cc      go to  sub,( 200,300 )
      if( sub.eq.200 ) go to 200
      if( sub.eq.300 ) go to 300
c
  110 if( ee .le. e2 )  go to 120
      ram3 = ram
      e3 = ee
cc      go to  sub,( 200,300 )
      if( sub.eq.200 ) go to 200
      if( sub.eq.300 ) go to 300
c
  120 ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = ee
cc      go to  sub,( 200,300 )
      if( sub.eq.200 ) go to 200
      if( sub.eq.300 ) go to 300
c
  130 do 140  i=1,k
cx  140 x1(i) = x(i) + ram*h(i)
      x1(i) = x(i) + ram*h(i)
  140 continue
cc      call  funct( k,x1,ee,g,ig )
      call  funct( xx,xmg,nn,k,x1,ee,g,ig )
c     if( ipr.ge.7 )  write(6,6)  ram,ee
cc      if(ipr.ge.7)  write(6,8)  ee,(x1(i)**2,i=1,k)
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         eee(nl) = ee
         do 141 i=1,k
cx  141    xx1(i,nl) = x1(i)
         xx1(i,nl) = x1(i)
  141    continue
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
      if( return.eq.80 ) go to 80
      if( return.eq.130 ) go to 130
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
      if( return.eq.80 ) go to 80
      if( return.eq.130 ) go to 130
c
  310 ram = (ram2+ram3)*0.5d0
cc      go to return ,( 80,130 )
      if( return.eq.80 ) go to 80
      if( return.eq.130 ) go to 130
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
      call  funct( xx,xmg,nn,k,x1,e3,g,ig )
c     if( ipr.ge.7 )  write(6,7)  ram,e3
cc      if(ipr.ge.7)  write(6,8)  e3,(x1(i)**2,i=1,k)
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         eee(nl) = e3
         do 521 i=1,k
cx  521    xx1(i,nl) = x1(i)
         xx1(i,nl) = x1(i)
  521    continue
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
      end
c***********************************************************************
cc      subroutine func4(n,b,f,h,ifg)
      subroutine func4(xx,xmg,ldata,n,b,f,h,ifg)
c----------------------------------------------------------------------
c     likelihood function of modified oomori type
c     hawkes' point process
c     the optimization w.r.t. a5 is not yet succeeded
c     at the date of 26th oct. 1981.
c     -----this is succeeded at 23 dec.1983-----
c----------------------------------------------------------------------
cx      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=17777, npara=5)
cc      common/xyod/xx(ldata),xmg(ldata),xmag0
cx      dimension xx(ldata),xmg(ldata)
      integer :: ldata, n, ifg
      real(8) :: xx(ldata), xmg(ldata), b(n), f, h(n)
      real(8) :: t, fff, aic2, tstart
cxx      common t,nn,mm,iappr
cxx      common /ddd/fff,aic2
cxx      common /range/tstart,ntstar
      common /etasap/ t,nn,mm,iappr
      common /ddd3/ fff,aic2
      common /range3/ tstart,ntstar
      common /kkxy/kkx,kky,kkt
cc      dimension b(npara),h(npara),ggt(npara),gt(npara),at(npara)
cx      dimension b(n),h(n),ggt(n),gt(n),at(n)
      real(8) :: ggt(n), gt(n), at(n), a1, a2, a3, a4, a5, ff,
     1           d1ff, d2ff, d3ff, d4ff, d5ff, rmdti1, ramdai, rmdi,
     2           rmd1i, rmdmi, rmdli, ft, d3ft, d4ft, d5ft, fs, d1fs,
     3           d2fs, d3fs, d4fs, d5fs, stsum, g1, g2, g3, g4, g5
      ifg=0
      a1=b(1)**2
      a2=b(2)**2
      a3=b(3)**2
      a4=b(4)**2
      a5=b(5)**2
c
      if(a4.gt.10.d0) go to 50
      if(a5.gt.3.d0) go to 50
c
        ff=0.0
      d1ff=0.0
      d2ff=0.0
      d3ff=0.0
      d4ff=0.0
      d5ff=0.0
         d1fs=0.0
      do 110 k=1,kkt
cx  110 at(k)=b(k+5)
      at(k)=b(k+5)
  110 continue
cx      do 120 k=1,kkt
cx  120 ggt(k)=0.0
      ggt(1:kkt)=0.0
      rmdti1=0.0
c
      if(xx(1).lt.tstart) go to 130
c
      do 140 k=1,kkt
cx  140 rmdti1=rmdti1+at(k)*(xx(1)/t)**k
      rmdti1=rmdti1+at(k)*(xx(1)/t)**k
  140 continue
      ramdai=a1+rmdti1
      if(ramdai.le.0.0) go to 50
      ff=log(ramdai)
      d1ff=1.d0/ramdai
      do 150 k=1,kkt
      ggt(k)=ggt(k)+(xx(1)/t)**k/ramdai
  150 continue
c
  130 continue
c
      do 10 i=2,nn
c
      rmdi=0.0
      rmd1i=0.0
      rmdmi=0.0
      rmdli=0.0
c
      do 20 j=1,i-1
      rmdi=rmdi+exp(a4*xmg(j))/(xx(i)-xx(j)+a3)**a5
      rmd1i=rmd1i+exp(a4*xmg(j))/(xx(i)-xx(j)+a3)**(a5+1)
      rmdmi=rmdmi+xmg(j)*exp(a4*xmg(j))/(xx(i)-xx(j)+a3)**a5
      rmdli=rmdli-exp(a4*xmg(j))/(xx(i)-xx(j)+a3)**a5
     &                          *log(xx(i)-xx(j)+a3)
   20 continue
c
      if(xx(i).lt.tstart) go to 10
c
      rmdti1=0.0
      if(kkt.eq.0) go to 160
      do 170 k=1,kkt
cx  170 rmdti1=rmdti1+at(k)*(xx(i)/t)**k
      rmdti1=rmdti1+at(k)*(xx(i)/t)**k
  170 continue
  160 continue
c
      ramdai=a1+a2*rmdi+rmdti1
      if(ramdai.le.0.0) go to 50
      ff=ff+log(ramdai)
      d1ff=d1ff+1.d0/ramdai
      d2ff=d2ff+rmdi/ramdai
      d3ff=d3ff-a2*a5*rmd1i/ramdai
      d4ff=d4ff+a2*rmdmi/ramdai
      if(a5.ne.1.d0) d5ff=d5ff+a2*rmdli/ramdai
      do 180 k=1,kkt
      ggt(k)=ggt(k)+(xx(i)/t)**k/ramdai
  180 continue
c
   10 continue
c
      d3ft=0.0
      d4ft=0.0
      d5ft=0.0
      ft=0.0
c
      if(a5.eq.1.d0) then
c
      do 30 i=ntstar+1,nn
      ft=ft+1.d0*(log(t-xx(i)+a3)-log(a3))*exp(a4*xmg(i))
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d4ft=d4ft+(log(t-xx(i)+a3)-log(a3))*xmg(i)*exp(a4*xmg(i))
   30 continue
c
      do 40 i=1,ntstar
      ft=ft+1.d0*(log(t-xx(i)+a3)-log(a3))*exp(a4*xmg(i))
      ft=ft-1.d0*(log(tstart-xx(i)+a3)-log(a3))*exp(a4*xmg(i))
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d3ft=d3ft-((tstart-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d4ft=d4ft+(log(t-xx(i)+a3)-log(a3))*xmg(i)*exp(a4*xmg(i))
      d4ft=d4ft-(log(tstart-xx(i)+a3)-log(a3))*xmg(i)*exp(a4*xmg(i))
   40 continue
c
      else
c
      do 60 i=ntstar+1,nn
      ft=ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *exp(a4*xmg(i))
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d4ft=d4ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*exp(a4*xmg(i))
      d5ft=d5ft+1.d0/(1.d0-a5)*(-(t-xx(i)+a3)**(1.d0-a5)
     &      *log(t-xx(i)+a3)+a3**(1.d0-a5)*log(a3))*exp(a4*xmg(i))
   60 continue
c
      do 70 i=1,ntstar
      ft=ft+1.d0/(1.d0-a5)*((     t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *exp(a4*xmg(i))
      ft=ft-1.d0/(1.d0-a5)*((tstart-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *exp(a4*xmg(i))
      d3ft=d3ft+((     t-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d3ft=d3ft-((tstart-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d4ft=d4ft+1.d0/(1.d0-a5)*((     t-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*exp(a4*xmg(i))
      d4ft=d4ft-1.d0/(1.d0-a5)*((tstart-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*exp(a4*xmg(i))
      d5ft=d5ft+1.d0/(1.d0-a5)*(-(     t-xx(i)+a3)**(1.d0-a5)
     &      *log(     t-xx(i)+a3)+a3**(1.d0-a5)*log(a3))*exp(a4*xmg(i))
      d5ft=d5ft-1.d0/(1.d0-a5)*(-(tstart-xx(i)+a3)**(1.d0-a5)
     &      *log(tstart-xx(i)+a3)+a3**(1.d0-a5)*log(a3))*exp(a4*xmg(i))
c
   70 continue
c
      endif
c
      fs=a1*(t-tstart)+a2*ft
      if(a1.gt.0.0) d1fs=(t-tstart)
      d2fs=ft
      d3fs=a2*d3ft
      d4fs=a2*d4ft
      d5fs=0.0
      if(a5.ne.1.d0) d5fs=a2*(ft/(1.d0-a5)+d5ft)
c
      stsum=0.0
      do 80 k=1,kkt
      stsum=stsum+at(k)*t*(     t/t)**(k+1)/(k+1)
      stsum=stsum-at(k)*t*(tstart/t)**(k+1)/(k+1)
   80 continue
c
      f=ff-fs-stsum
      g1=d1ff-d1fs
      g2=d2ff-d2fs
      g3=d3ff-d3fs
      g4=d4ff-d4fs
      g5=d5ff-d5fs
c
      do 90 k=1,kkt
      gt(k)=ggt(k)-t*(     t/t)**(k+1)/(k+1)
      gt(k)=ggt(k)+t*(tstart/t)**(k+1)/(k+1)
   90 continue
c
      f=-f
      h(1)=-g1*2.0d0*b(1)
      h(2)=-g2*2.0d0*b(2)
      h(3)=-g3*2.0d0*b(3)
      h(4)=-g4*2.0d0*b(4)
      h(5)=-g5*2.0d0*b(5)
      do 100 k=1,kkt
cx  100 h(k+5)=-gt(k)
      h(k+5)=-gt(k)
  100 continue
      if(a5.eq.1.d0) h(5)=0.0
      fff=f
c     write(6,1030) fff,a1,a2,a3,a4,a5
cx 1030 format(1h ,'f=',e12.5,'; x=',6e12.5)
cx    3 format(1h ,110x,d18.10)
cx    1 format(1h ,7d18.10)
      return
   50 continue
      ifg=1
      f=1.0d30
      return
      end
c***********************************************************************
cc      subroutine func9(n,b,f,h,ifg)
      subroutine func91(xx,xmg,ldata,n,b,f,h,ifg)
c-----------------------------------------------------------------------
c     likelihood function of "etas" model with the approximation
c     method: this is copied fron func2 and modified (1992 october)
c     improved for the effective use of 'iap' processor.
c     when is=16, the gradient were not converge because approximation
c     for a3 did not fit well. this subrourine has overcomed this
c     by making use of the direct diffrerential of lamdai: see pi#(.).
c-----------------------------------------------------------------------
cx      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=17777, npara=5)
cc      common/xyod/xx(ldata),xmg(ldata),xmag0
      integer :: ldata, n, ifg
      real(8) :: xx(ldata), xmg(ldata), b(n), f, h(n)
      real(8) :: t, fff, aic2, tstart
cx      dimension xx(ldata),xmg(ldata)
cxx      common t,nn,mm,iappr
cxx      common /ddd/fff,aic2
cxx      common /range/tstart,ntstar
      common /etasap/ t,nn,mm,iappr
      common /ddd3/ fff,aic2
      common /range3/ tstart,ntstar
cc      dimension b(npara),h(npara)
cx      dimension b(n),h(n)
cx      dimension xi1(144),xi2(144),wx1(144),wx2(144)
cx      dimension fi1(144),fi2(144),alf1(144),alf2(144),ci1(144),ci2(144)
cx      dimension rmd(ldata),rmdc(ldata),rmdm(ldata),rmdl(ldata)
      integer :: ixhiab
      real(8) :: xi1(144), xi2(144), wx1(144), wx2(144), fi1(144),
     1           fi2(144), alf1(144), alf2(144), ci1(144), ci2(144),
     2           rmd(ldata), rmdc(ldata), rmdm(ldata), rmdl(ldata),
     3           delta0, xi0, wx0, pi2, delta, a1, a2, a3, a4, a5,
     4           ff, d1ff, d2ff, d3ff, d4ff, d5ff, rmdi, rmdci, rmdmi,
     5           rmdli, fi0, ci0, alf0, cpg, cpg3, cpg5, qi0,
     6           gi0, hi0, blf0, qi1, qi2, gi1, hi1, gi2, hi2, gam,
     7           blf1, blf2, ramdai, d3ft, d4ft, d5ft, ft, effmag, fs,
     8           d1fs, d2fs, d3fs, d4fs, d5fs, g1, g2, g3, g4, g5
      data ixhiab /0/
      save ixhiab,delta0,xi0,xi1,xi2,wx0,wx1,wx2
c
      if(ixhiab.eq.0) call hiab(delta0,xi0,xi1,xi2,wx0,wx1,wx2)
c
      pi2=1.570796326794397d0
c take 'is' below as one of the 1,2,4,8 or 16:
c the larger is the faster but the less accurate.
      is=iappr
      im=is
c
      delta=delta0*im
cc      if(ixhiab.eq.0) write(6,*) 'is=',is
      ixhiab=1
      ifg=0
      a1=b(1)**2
      a2=b(2)**2
      a3=b(3)**2
      a4=b(4)**2
      a5=b(5)**2
c
      if(a5.gt.10.d0) go to 50
c
      ff=0.0
      d1ff=0.0
      d2ff=0.0
      d3ff=0.0
      d4ff=0.0
      d5ff=0.0
c
      if(ixhiab.eq.0) call hiab(delta,xi0,xi1,xi2,wx0,wx1,wx2)
      ixhiab=1
c
c the folloings are not redundant: consider the case where i=1.
      rmdi=0.0
      rmdci=0.0
      rmdmi=0.0
      rmdli=0.0
      fi0=0.0
      ci0=0.0
      alf0=0.0
      do 90 j=is,144,im
      fi1(j)=0.0
      fi2(j)=0.0
      ci1(j)=0.0
      ci2(j)=0.0
      alf1(j)=0.0
      alf2(j)=0.0
   90 continue
      cpg=a3**(-a5)/gam(0,a5)
      cpg3=-a5*a3**(-a5-1)/gam(0,a5)
      cpg5=-cpg*(log(a3)+gam(1,a5)/gam(0,a5))
      do 10 i=1,nn
      if(xx(i).eq.0.0.and.i.eq.1) go to 10
      if(i.eq.1) go to 10
      if (a4*xmg(i-1).gt.170) go to 50
cx      if (a4*xmg(i-1).gt.170) write(6,*) ' go to 50'
cx      if (xx(i)-xx(i-1).lt.0.0) write(6,*) 'reverse',i,xx(i),xmg(i)
c
      ci0=exp(-xi0/a3*(xx(i)-xx(i-1)))*
     &    (ci0+xi0/a3**2*(xx(i)-xx(i-1))*(fi0+exp(a4*xmg(i-1))))
      qi0=ci0*xi0**(a5-1)*exp(-xi0)
      fi0=(fi0+exp(a4*xmg(i-1)))*exp(-xi0/a3*(xx(i)-xx(i-1)))
      gi0=fi0*xi0**(a5-1)*exp(-xi0)
      hi0=gi0*log(xi0)
      alf0=(alf0+xmg(i-1)*exp(a4*xmg(i-1)))*exp(-xi0/a3*(xx(i)-xx(i-1)))
      blf0=alf0*xi0**(a5-1)*exp(-xi0)
      rmd(i)=gi0*wx0
      rmdc(i)=qi0*wx0
      rmdm(i)=blf0*wx0
      rmdl(i)=hi0*wx0
   10 continue
c
      do 20 j=is,144,im
      do 25 i=2,nn
c
      ci1(j)=exp(-xi1(j)/a3*(xx(i)-xx(i-1)))*
     &   (ci1(j)+xi1(j)/a3**2*(xx(i)-xx(i-1))*(fi1(j)+exp(a4*xmg(i-1))))
      qi1=ci1(j)*xi1(j)**(a5-1)*exp(-xi1(j))
      ci2(j)=exp(-xi2(j)/a3*(xx(i)-xx(i-1)))*
     &   (ci2(j)+xi2(j)/a3**2*(xx(i)-xx(i-1))*(fi2(j)+exp(a4*xmg(i-1))))
      qi2=ci2(j)*xi2(j)**(a5-1)*exp(-xi2(j))
      rmdc(i)=rmdc(i)+qi1*wx1(j)+qi2*wx2(j)
      fi1(j)=(fi1(j)+exp(a4*xmg(i-1)))*exp(-xi1(j)/a3*(xx(i)-xx(i-1)))
      gi1=fi1(j)*xi1(j)**(a5-1)*exp(-xi1(j))
      hi1=gi1*log(xi1(j))
      fi2(j)=(fi2(j)+exp(a4*xmg(i-1)))*exp(-xi2(j)/a3*(xx(i)-xx(i-1)))
      gi2=fi2(j)*xi2(j)**(a5-1)*exp(-xi2(j))
      hi2=gi2*log(xi2(j))
      rmd(i)=rmd(i)+gi1*wx1(j)+gi2*wx2(j)
      rmdl(i)=rmdl(i)+hi1*wx1(j)+hi2*wx2(j)
c
      alf1(j)=(alf1(j)+xmg(i-1)*exp(a4*xmg(i-1))) *
     &      exp(-xi1(j)/a3*(xx(i)-xx(i-1)))
      blf1=alf1(j)*xi1(j)**(a5-1)*exp(-xi1(j))
      alf2(j)=(alf2(j)+xmg(i-1)*exp(a4*xmg(i-1))) *
     &      exp(-xi2(j)/a3*(xx(i)-xx(i-1)))
      blf2=alf2(j)*xi2(j)**(a5-1)*exp(-xi2(j))
      rmdm(i)=rmdm(i)+blf1*wx1(j)+blf2*wx2(j)
   25 continue
   20 continue
c
c     do 15 i=1,nn
      do 15 i=ntstar+1,nn
      if(i.eq.1) go to 80
      rmdi=rmd(i)*pi2*delta
      rmdci=rmdc(i)*pi2*delta
      rmdmi=rmdm(i)*pi2*delta
      rmdli=rmdl(i)*pi2*delta
   80 continue
      ramdai=a1+a2*cpg*rmdi
      ff=ff+log(ramdai)
      d1ff=d1ff+1.d0 /ramdai
      d2ff=d2ff+1.d0 /ramdai* cpg *rmdi
      d3ff=d3ff+a2   /ramdai*(cpg3*rmdi+cpg*rmdci)
      d4ff=d4ff+a2   /ramdai* cpg *rmdmi
      if(a5.ne.1.d0)
     &d5ff=d5ff+a2   /ramdai*(cpg5*rmdi+cpg*rmdli)
   15 continue
c
      d3ft=0.0
      d4ft=0.0
      d5ft=0.0
      ft=0.0
c
      if(a5.eq.1.d0) then
c
      do 30 i=ntstar+1,nn
      effmag=exp(a4*xmg(i))
      ft=ft+1.d0*(log(t-xx(i)+a3)-log(a3))*effmag
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d4ft=d4ft+(log(t-xx(i)+a3)-log(a3))*xmg(i)*effmag
   30 continue
c
      do 40 i=1,ntstar
      effmag=exp(a4*xmg(i))
      ft=ft+1.d0*(log(t-xx(i)+a3)-log(a3))*effmag
      ft=ft-1.d0*(log(tstart-xx(i)+a3)-log(a3))*effmag
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d3ft=d3ft-((tstart-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d4ft=d4ft+(log(t-xx(i)+a3)-log(a3))*xmg(i)*effmag
      d4ft=d4ft-(log(tstart-xx(i)+a3)-log(a3))*xmg(i)*effmag
   40 continue
c
      else
c
      do 60 i=ntstar+1,nn
      effmag=exp(a4*xmg(i))
      ft=ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *effmag
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d4ft=d4ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*effmag
      d5ft=d5ft+1.d0/(1.d0-a5)**2 *
     &           ((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5)) * effmag
     &         +1.d0/(1.d0-a5)    * ( -(t-xx(i)+a3)**(1.d0-a5)
     &      *log(t-xx(i)+a3)+a3**(1.d0-a5)*log(a3) ) * effmag
   60 continue
c
      do 70 i=1,ntstar
      effmag=exp(a4*xmg(i))
      ft=ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *effmag
      ft=ft-1.d0/(1.d0-a5)*((tstart-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *effmag
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d3ft=d3ft-((tstart-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d4ft=d4ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*effmag
      d4ft=d4ft-1.d0/(1.d0-a5)*((tstart-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*effmag
      d5ft=d5ft+1.d0/(1.d0-a5)**2 *
     &           ((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5)) * effmag
     &         +1.d0/(1.d0-a5)    * ( -(t-xx(i)+a3)**(1.d0-a5)
     &      *log(t-xx(i)+a3)+a3**(1.d0-a5)*log(a3) ) * effmag
      d5ft=d5ft-1.d0/(1.d0-a5)**2 *
     &           ((tstart-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5)) * effmag
     &         -1.d0/(1.d0-a5)    * ( -(tstart-xx(i)+a3)**(1.d0-a5)
     &      *log(tstart-xx(i)+a3)+a3**(1.d0-a5)*log(a3) ) * effmag
   70 continue
c
      endif
c
      fs=a1*(t-tstart)+a2*ft
      d1fs=0.0
      if(a1.gt.0.0) d1fs=t-tstart
      d2fs=ft
      d3fs=a2*d3ft
      d4fs=a2*d4ft
      d5fs=0.0
      if(a5.ne.1.d0) d5fs=a2*d5ft
c
      f=ff-fs
      g1=d1ff-d1fs
      g2=d2ff-d2fs
      g3=d3ff-d3fs
      g4=d4ff-d4fs
      g5=d5ff-d5fs
c
      f=-f
      h(1)=-g1*2.0d0*b(1)
      h(2)=-g2*2.0d0*b(2)
      h(3)=-g3*2.0d0*b(3)
      h(4)=-g4*2.0d0*b(4)
      h(5)=-g5*2.0d0*b(5)
      if(a5.eq.1.d0) h(5)=0.0
      fff=f
c     write(6,1030) fff,a1,a2,a3,a4,a5
cx 1030 format(1h ,'f=',e12.5,'; x=',6e12.5)
cx    3 format(1h ,110x,d18.10)
cx    1 format(1h ,7d13.6)
      return
   50 continue
      ifg=1
      f=1.0d30
      return
      end
c***********************************************************************
      subroutine hiab(h,a0,a1,a2,b0,b1,b2)
cx      implicit real*8(a-h,o-z)
cx      dimension a1(144),a2(144),b1(144),b2(144)
      real(8) :: h, a0, a1(144), a2(144), b0, b1(144), b2(144)
      real(8) :: pi2, eh, ehi, eni, en, s1
      pi2=1.570796326794397d0
      h=1.d0/32
      eh=exp(h)
      ehi=1/eh
      eni=1.d0
      in=144
      inc=in
      a0=1.d0/exp(pi2)
      b0=2.d0*a0
      do 20 i=1,in
      eni=ehi*eni
      en=1.d0/eni
      s1=h*i
      a1(i)=1.d0/exp(pi2*(s1+en))
      b1(i)=a1(i)*(1.d0+en)
      a2(i)=exp(pi2*(s1-eni))
      b2(i)=a2(i)*(1.d0+eni)
   20 continue
      return
      end
c***********************************************************************
c     real*8 function gam(id,q)
cx      real*8 function dbgam(id,q)
      double precision function dbgam(id,q)
c hitac    real function dbgam*8(id,q)
cx      implicit real*8(a-h, o-z)
cx      dimension a(10)
      integer :: id
      real(8) :: q
      real(8) :: a(10)
      data a/ 0.99999 99998 71452d0, 0.42278 43615 29813d0,
     1        0.41183 94326 49605d0, 0.08158 87915 49927d0,
     2        0.07416 87114 09713d0, 0.00004 89152 06125d0,
     3        0.01038 02945 70428d0, -.00162 85524 78086d0,
     4        0.00082 46883 39196d0, -.00000 66427 76723d0 /
      real(8) :: fact, dfac, d2fac, p, x, gamm, gam1, dgam, d2gam,
     1           eps, gam
      fact=1
      dfac=0.0
      d2fac=0.0
      p=q
   30 if(p.ge.1.d0 .and. p.le.2.d0) then
        d2fac=d2fac*p+2*dfac
        dfac=dfac*p+fact
        fact=fact*p
        x=p-1
      else if (p.lt.1) then
        d2fac=d2fac*p+2*dfac
        dfac=dfac*p+fact
        fact=fact*p
        p=p+1
        go to 30
      else
        d2fac=(d2fac*(p-1))/(p-1)**2
     &        -2*(dfac*(p-1)-fact)/(p-1)**3
        dfac=(dfac*(p-1)-fact)/(p-1)**2
        dfac=(dfac*(p-1)-fact)/(p-1)**2
        fact=fact/(p-1)
        p=p-1
        go to 30
      end if
c
      gamm=0.0
      gam1=0.0
      dgam=0.0
      d2gam=0.0
      eps=0.1d-6
      do 10 i=1,10
      gamm=gamm+a(i)*x**(i-1)
      dgam=dgam+(i-1)*a(i)*x**(i-2)
      d2gam=d2gam+(i-1)*(i-2)*a(i)*x**(i-3)
   10 continue
cxx      if(id.eq.0) gam=gamm/fact
cxx      if(id.eq.1) gam=(dgam*fact-gamm*dfac)/fact**2
cxx      if(id.eq.2) gam=(d2gam*fact-gamm*d2fac)/fact**2
cxx     &             -2*(dgam*fact-gamm*dfac)/fact**3 *dfac
      if(id.eq.1) then
         gam=(dgam*fact-gamm*dfac)/fact**2
      else if(id.eq.2) then
         gam=(d2gam*fact-gamm*d2fac)/fact**2
     &             -2*(dgam*fact-gamm*dfac)/fact**3 *dfac
      else
         gam=gamm/fact
      end if
      dbgam=gam
      return
      end
c***********************************************************************
c     real*8 function qbgam(id,qq)
cx      real*8 function gam(id,qq)
      double precision function gam(id,qq)
c hitac     real function gam*8(id,q)
c     implicit real*16(a-h, o-z)
cx      implicit real* 8(a-h, o-z)
cx      real*8 qq
cx      dimension a(11),b(11)
      integer :: id
      real(8) :: qq
      real(8) :: a(11), b(11)
      data a/ -2 98354.32785 74342 13883 04376 59         d0,
     1        -2 38495.39700 18198 87246 87344 23         d0,
     2        -1 17049.47601 21780 68840 38544 45         d0,
     3        -  39494.45048 30157 19364 21824 091        d0,
     4        -  10466.99423 82752 14053 30650 531        d0,
     5        -   2188.21811 00718 16369 39479 5998       d0,
     6        -    380.51122 08641 73465 75849 22631      d0,
     7        -     52.83123 75563 58453 83718 97838 2    d0,
     8        -      6.12857 17637 04498 30688 94282 12   d0,
     9        -       .50280 18054 41681 24673 64198 75   d0,
     t        -       .03343 06032 23305 95274 51566 0112 d0 /
      data b/ -2 98354.32785 74342 13883 04385 24         d0,
     1        -1 12355.86087 48644 91134 23064 08         d0,
     2           53327.16689 11814 21574 85686 311        d0,
     3            8571.16049 89070 43851 96114 7763       d0,
     4        -   4734.86597 70282 11706 55681 977        d0,
     5             196.04976 12885 58583 89970 39621      d0,
     6             125.77333 67869 88864 59666 47426      d0,
     7        -     20.53126 15310 06727 64513 92906 7    d0,
     8               1.                                   d0,
     9               0.0                                  d0,
     t               0.0                                  d0 /
      real(8) :: fact, dfac, d2fac, p, x, gam1, gam2, dga1,
     1           dga2, d2ga1, d2ga2, eps
      fact=1
      dfac=0.0
      d2fac=0.0
      p=qq
   30 if(p.ge.1.d0 .and. p.le.2.d0) then
        d2fac=d2fac*p+2*dfac
        dfac=dfac*p+fact
        fact=fact*p
        x=p-1
      else if (p.lt.1.d0) then
        d2fac=d2fac*p+2*dfac
        dfac=dfac*p+fact
        fact=fact*p
        p=p+1
        go to 30
      else
        d2fac=d2fac/(p-1)+2*dfac/(p-1)**2*(-1)
     &        +fact/(p-1)**3*(-1)*(-2)
        dfac=dfac/(p-1)+fact/(p-1)**2*(-1)
        fact=fact/(p-1)
        p=p-1
        go to 30
      end if
c
      gam1=a(1)
      gam2=b(1)
      do 10 i=2,11
      gam1=gam1+a(i)*x**(i-1)
cx   10 gam2=gam2+b(i)*x**(i-1)
      gam2=gam2+b(i)*x**(i-1)
   10 continue
      dga1=a(2)
      dga2=b(2)
      do 40 i=3,11
      dga1=dga1+(i-1)*a(i)*x**(i-2)
      dga2=dga2+(i-1)*b(i)*x**(i-2)
   40 continue
      d2ga1=2*a(3)
      d2ga2=2*b(3)
      do 50 i=4,11
      d2ga1=d2ga1+(i-1)*(i-2)*a(i)*x**(i-3)
      d2ga2=d2ga2+(i-1)*(i-2)*b(i)*x**(i-3)
   50 continue
c
cxx      if(id.eq.0) gam=gam1/gam2/fact
      gam=gam1/gam2/fact
      if(id.eq.1) gam=(dga1*gam2*fact-gam1*dga2*fact-gam1*gam2*dfac)
     & / (gam2*fact)**2
      if(id.eq.2) gam=(d2ga1*gam2*fact+dga1*dga2*fact+dga1*gam2*dfac
     &                -dga1*dga2*fact-gam1*d2ga2*fact-gam1*dga2*dfac
     &                -dga1*gam2*dfac-gam1*dga2*dfac-gam1*gam2*d2fac)
     & / (gam2*fact)**2
     &             -2*(dga1*gam2*fact-gam1*dga2*fact-gam1*gam2*dfac)
     & / (gam2*fact)**3 * (dga2*fact+gam2*dfac)
      eps=0.1d-4
      return
      end
