c--------------------------------------------------------------------
c83-11-11-15:06:19/87-03-25-10:21 nhpossn fort 
cc      program momori
cx      subroutine momorif(xx,nni,xini,np,zts,zte,tstart,nc,nfuncti,
      subroutine momorif(xx,nni,xini,np,zts,zte,nc,nfuncti,
     & ff,x,g,pa,ahaic,t00,ti,ak,c,p,cls, id,rmd,xx1,h,hf,nl,nlmax)

c-----------------------------------------------------------------------
c     this program performs the maximum likelihood calculation and 
c     estimate of the parameters of modified omori formula assuming 
c     piecewise non-stationary poisson process. 
c 
c     this is made based on the 'nhpossn4'. 
c 
c     func5---------exponential decay poisson process 
c     func6---------the modified omori's model for a interval (t0,t1) 
c     func9---------model for detecting the 2nd aftershock 
c     func10--------model for a pair of 2nd aftershock sequences 
c 
c     this program is desined by y.ogata and programed by y.ogata and 
c     k.katsura.  see ogata (1983; j.p.e., vol.31,pp.115-124). 
c----------------------------------------------------------------------- 
c
      include 'sapp_f.h'
c
cx      implicit real * 8 (a-h,o-z) 
cc      real*4 time 
cc      common/xyod/xdumy,xx(19999) 
cc      common/y/y(50)
      integer :: nni, np, nc, nfuncti, nlmax, id(nlmax), nl
      real(8) :: xx(nni), xini(np+1), zts, zte, ff, x(np,2), g(np,2),
     1           pa(np), ahaic(nc), t00, ti((np-1)/3), ak((np-1)/3),
     2           c((np-1)/3),p((np-1)/3),cls((np-1)/3), rmd(nlmax),
     3           xx1(np,nlmax), h(np,np,2), hf(np,np,2,2)
      real(8) :: t, t0, t1, t2, t3
cxx      common t,nn,nfunct 
cxx      common/range/t0,t1,t2,t3
      common /momori/ t,nn,nfunct 
      common /range1/ t0,t1,t2,t3
cc      dimension xxxx(50),x(50),ahaic(30) 
cc      dimension axxx(19999),amag(19999)
cx      dimension x(np,2),pa(np),ahaic(nc),g(np,2)
cx      dimension xx(nni),xini(np+1),ti((np-1)/3)
cx      dimension ak((np-1)/3),c((np-1)/3),p((np-1)/3),cls((np-1)/3) 
cx      dimension id(nlmax), rmd(nlmax), xx1(np,nlmax)
cx      dimension h(np,np,2), hf(np,np,2,2)
c
      nn=nni
      nfunct=nfuncti
c
      nl = 0
cx      do 5 i = 1,nlmax
cx    5 id(i) = 0
      id(1:nlmax) = 0
c
cc      call input(nnnn,axxx,amag)
      t=zte-zts
      if(nfunct.lt.6) go to 155 
      t0=zts 
      t1=zte 
cx      bmag=amx1 
  155 continue 
c 
      kaisu=1 
      iend=0 
c     if(nfunct.eq.6) kaisu=100 
      do 9753 ijkl=1,kaisu 
cc      if(ijkl.ne.1) write(6,9751) 
cx 9751 format(1h ) 
      if(nn.eq.0) go to 9753
cc      call repara(ijkl,nnn,xxxx,x)
cc      call dav(nnn,x,ahaic)
cc      call output(nnn,x,ahaic)
cc      if(nfunct.ne.5) call sizes(nnn,x)
cx      call repara(xini,np+1,ijkl,nnn,x)
      call repara(xini,np+1,nnn,x)
      call dav6(nni,xx,nnn,x,g,nc,ahaic,pa,
     &  id,rmd,xx1,h,hf(1,1,1,1),hf(1,1,1,2),nl,nlmax)
cx      call output6(nnn,pa,nc,ahaic,ff)
      call output6(nnn,pa,ff)
      kn=(nnn-1)/3
      if(nfunct.ne.5) call sizes(nnn,pa,kn,t00,ti,ak,c,p,cls) 
c     call clock(time) 
c      write(6,9752) time 
cx 9752 format(1h ,'time= ',f10.3) 
 9753 continue 
cx 9754 continue 
      return 
      end 
cc      subroutine repara(ijkl,nnn,xxxx,x)
cx      subroutine repara(xini,n,ijkl,nnn,x) 
      subroutine repara(xini,n,nnn,x) 
c----------------------------------------------------------------------- 
card 7 nnn (i10) # of parameters 
card 8 (x(i),i=1,n) (8f10.2) initial estimates 
c----------------------------------------------------------------------- 
cx      implicit real * 8 (a-h,o-z) 
      integer :: n, nnn
      real(8) :: xini(n), x(n-1)
      real(8) :: t, t0, t1, t2, t3
cxx      common /range/t0,t1,t2,t3 
cxx      common t,nn,nfunct 
      common /range1/ t0,t1,t2,t3 
      common /momori/ t,nn,nfunct 
cc      dimension xxxx(50),x(50),xini(50)
cx      dimension xxxx(n-1),x(n-1),xini(n) 
      real(8) :: xxxx(n-1)
cc      n=5 
cc      read(1,*) (xini(i),i=1,n) 
cx 1010 format(7f10.4) 
      nnn=n-1 
      do 41 i=1,nnn 
      xxxx(i)=xini(i) 
      if(i.eq.nnn) xxxx(i)=xini(n) 
cx   41 x(i)=xxxx(i)
      x(i)=xxxx(i)
   41 continue
      do 45 i=1,nnn 
      if(nfunct.eq.5.or.nfunct.eq.6) x(i)=sqrt(x(i)) 
      if(nfunct.eq.9.and.x(i).ne.0.0) x(i)=log(x(i)) 
      if(nfunct.eq.10.and.x(i).ne.0.0) x(i)=log(x(i)) 
   45 continue 
      return 
      end 
c 
c 
c 
c 
c 
cc      subroutine dav(n,x,ahaic)
      subroutine dav6(nni,xx,n,x0,g,ncount,ahaic,x,
     &                     id,rmd,xx1,h,hf,hfi,nl,nlmax)

cx      implicit real * 8 (a-h,o-z) 
      external func5,func6,func9,func10
      integer :: nni, n, ncount, nlmax, id(nlmax), nl
      real(8) :: xx(nni), x0(n,2), g(n,2), ahaic(ncount), x(n),
     1           rmd(nlmax), xx1(n,nlmax), h(n,n,2), hf(n,n,2),
     2           hfi(n,n,2)
      real(8) :: t, t0, t1, t2, t3, f, aic
cxx      common t,nn,nfunct 
cxx      common/range/t0,t1,t2,t3 
      common /momori/ t,nn,nfunct 
      common /range1/ t0,t1,t2,t3 
cc      common/y/y(50)
cxx      common/ddd/f,aic 
      common /ddd1/ f,aic 
cc      dimension x(50),ahaic(30)
cx      dimension x(n),ahaic(ncount)
cx      dimension xx(nni), x0(n,2), g(n,2)
cx      dimension  id(nlmax), rmd(nlmax), xx1(n,nlmax)
cx      dimension  h(n,n,2), hf(n,n,2), hfi(n,n,2)
c
cc      write(6,1020) n 
cc      write(6,1030)  (x(i),i=1,n)
      do 20 ii=1,n
      x(ii)=x0(ii,1)
   20 continue
      if(nfunct.eq.5) mm=3 
      if(nfunct.eq.6) mm=4 
      if(nfunct.eq.9) mm=7 
      if(nfunct.eq.10) mm=10
      if(mm.lt.n) mm=n 
      do 30 ii=1,2
cc      if(nfunct.eq.5) call davidn(x,n,0,func5) 
cc      if(nfunct.eq.6) call davidn(x,n,0,func6) 
cc      if(nfunct.eq.9) call davidn(x,n,0,func9) 
cc      if(nfunct.eq.10) call davidn(x,n,0,func10) 
cx      if(nfunct.eq.5) call davidn6(xx,nni,x,n,mm,0,func5,g(1,ii),
      if(nfunct.eq.5) call davidn6(xx,nni,x,n,mm,func5,g(1,ii),
     &          id,rmd,xx1,h(1,1,ii),hf(1,1,ii),hfi(1,1,ii),nl,nlmax)
cx      if(nfunct.eq.6) call davidn6(xx,nni,x,n,mm,0,func6,g(1,ii),
      if(nfunct.eq.6) call davidn6(xx,nni,x,n,mm,func6,g(1,ii),
     &          id,rmd,xx1,h(1,1,ii),hf(1,1,ii),hfi(1,1,ii),nl,nlmax)
cx      if(nfunct.eq.9) call davidn6(xx,nni,x,n,mm,0,func9,g(1,ii),
      if(nfunct.eq.9) call davidn6(xx,nni,x,n,mm,func9,g(1,ii),
     &          id,rmd,xx1,h(1,1,ii),hf(1,1,ii),hfi(1,1,ii),nl,nlmax)
cx      if(nfunct.eq.10) call davidn6(xx,nni,x,n,mm,0,func10,g(1,ii),
      if(nfunct.eq.10) call davidn6(xx,nni,x,n,mm,func10,g(1,ii),
     &          id,rmd,xx1,h(1,1,ii),hf(1,1,ii),hfi(1,1,ii),nl,nlmax)
      do 25 jj=1,n
cx   25   x0(jj,ii) = x(jj)
        x0(jj,ii) = x(jj)
   25 continue
   30 continue 
cx   80 continue 
      ahaic(1)=aic 
      return 
cx 1020 format(1h ,3x,'input data'/1h ,5x,'n=',i3) 
cx 1030 format(1h ,                   5x,'x=',6e16.7) 
      end 
cc      subroutine output(n,x,ahaic)
cx      subroutine output6(n,x,ncount,ahaic,ff)
      subroutine output6(n,x,ff)
cx      implicit real * 8 (a-h,o-z)
      integer :: n
      real(8) :: x(n), ff
      real(8) :: t, t0, t1, t2, t3, f, aic 
      real(8) :: x0
cxx      common t,nn,nfunct 
cxx      common/range/t0,t1,t2,t3 
cxx      common/ddd/f,aic 
      common /momori/ t,nn,nfunct 
      common /range1/ t0,t1,t2,t3 
      common /ddd1/ f,aic 
cc      common/y/y(50) 
cc      dimension x(50),ahaic(30) 
cx      dimension x(n),ahaic(ncount)
cx      dimension x(n)
      do 70 i=1,n 
      if(nfunct.ne.9.and.nfunct.ne.10) x(i)=x(i)**2 
      if(x(i).eq.0.0) go to 70 
      if(nfunct.eq.9.or.nfunct.eq.10) x(i)=exp(x(i)) 
   70 continue 
ckk      write(6,1040) f,(x(i),i=1,n) 
      x0=0.0 
cc      write(6,1040) f,(x(i),i=1,3),x0,x(4)
c---  
      ff=f
c---
cc      ncount=1 
cc      do 110 iii=1,ncount 
cc      write(6,1080) iii,ahaic(iii) 
cc  110 continue 
cx 1080 format(1h ,i10,d20.10) 
      return 
cx 1000 format(3i10,2f15.6) 
cx 1010 format(8f10.4) 
cx 1020 format(1h ,3x,'input data'/1h ,5x,'n=',i3) 
cx 1030 format(1h ,                   5x,'x=',6e16.7) 
cx 1040 format(
cx     2      /1h ,'neg max lklhd=',1 e16.7
cx     3    /1h ,'max lklhd est.=',10e12.5/('                 ',10e12.5))
cx 1050 format(4d20.13) 
cx 1100 format(i10) 
cx 1060 format(e25.15) 
cx 1070 format(1h ,'  c = ',e25.15) 
      end 
cc      subroutine  linear( x,h,ram,ee,k,ig,funct ) 
      subroutine  linear6( xx,nn,x,h,ram,ee,k,ig,funct,
     &                         id,rmd,xx1,nl,nlmax)
c 
c     this subroutine performs the linear search along the direction spe 
c     by the vector h 
c       ---------------------------------------------------------------- 
c       the following subroutine is directly called by this subroutine: 
c             funct 
c       ---------------------------------------------------------------- 
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
      integer :: return, sub 
cc      dimension  x(1) , h(1) , x1(50) 
cc      dimension  g(50)
cx      dimension  x(1) , h(1) , x1(k)
      integer :: nn, k, ig, nlmax, id(nlmax), nl
      real(8) :: xx(nn), x(k), h(k), ram, ee, rmd(nlmax), xx1(k,nlmax)
      real(8) :: x1(k), g(k), const2, hnorm, ram1, ram2, ram3, 
     1           e1, e2, e3, a1, a2, a3, b1, b2
cx      dimension  x(k) , h(k) , x1(k)
cx      dimension  g(k)
cx      dimension  xx(nn) 
c
cx      dimension id(nlmax),rmd(nlmax),xx1(k,nlmax)
c
      external funct 
cxx      common     / ccc /  isw , ipr
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
      call  funct( nn,xx,k,x1,e2,g,ig ) 
c     if(ipr.ge.7)  write(6,2)  ram2,e2 
cc      if(ipr.ge.7)  write(6,8)  ram2,(x1(i)**2,i=1,k) 
cx    8 format(1h ,'-ll=',d13.5,1x,4d12.5) 
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         rmd(nl) = ram2
         do 21 i=1,k
cx   21    xx1(i,nl) = x1(i)
         xx1(i,nl) = x1(i)
   21    continue
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
      call  funct( nn,xx,k,x1,e3,g,ig ) 
      if( ig.eq.1 )  go to  500 
c     if( ipr.ge.7 )  write(6,3)  ram3,e3 
cc      if(ipr.ge.7)  write(6,8)  ram3,(x1(i)**2,i=1,k) 
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         rmd(nl) = ram3
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
      call  funct( nn,xx,k,x1,e2,g,ig ) 
c     if(ipr.ge.7)  write(6,4)  ram2,e2 
cc      if(ipr.ge.7)  write(6,8)  ram2,(x1(i)**2,i=1,k) 
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         rmd(nl) = ram2
         do 61 i=1,k
cx   61    xx1(i,nl) = x1(i)
         xx1(i,nl) = x1(i)
   61    continue
      end if
      if( e2.gt.e1 )  go to 50 
c 
cc   70 assign 80 to return 
   70 continue
      return = 80
      go to 200 
c 
   80 do 90  i=1,k 
cx   90 x1(i) = x(i) + ram*h(i) 
      x1(i) = x(i) + ram*h(i)
   90 continue
cc      call  funct( k,x1,ee,g,ig ) 
      call  funct( nn,xx,k,x1,ee,g,ig ) 
c     if(ipr.ge.7)  write(6,5)  ram,ee 
cc      if(ipr.ge.7)  write(6,8)  ram,(x1(i)**2,i=1,k) 
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         rmd(nl) = ram
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
      sub = 300
      sub = 200 
   95  continue
      return = 130
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
      call  funct( nn,xx,k,x1,ee,g,ig ) 
c     if( ipr.ge.7 )  write(6,6)  ram,ee 
cc      if(ipr.ge.7)  write(6,8)  ram,(x1(i)**2,i=1,k) 
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         rmd(nl) = ram
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
cc     go to return ,( 80,130 ) 
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
      call  funct( nn,xx,k,x1,e3,g,ig ) 
c     if( ipr.ge.7 )  write(6,7)  ram,e3 
cc      if(ipr.ge.7)  write(6,8)  ram,(x1(i)**2,i=1,k) 
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 8
         rmd(nl) = ram2
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
      e n d 
cc      subroutine  davidn( x,n,ihes,funct )
cx      subroutine  davidn6( xx,nni,x,n,m,ihes,funct,g,
      subroutine  davidn6( xx,nni,x,n,m,funct,g,
     &                          id,rmd,xx1,h,hf,hfi,nl,nlmax)
c 
c          minimization by davidon-fletcher-powell procedure 
c 
c       ---------------------------------------------------------------- 
c       the following subroutines are directly called by this subroutine 
c             funct 
c             hesian 
c             linear 
c       ---------------------------------------------------------------- 
c          inputs: 
c             x:       vector of initial values 
c             k:       dimension of the vector x 
c             ihes:    =0   inverse of hessian matrix is not available 
c                      =1   inverse of hessian matrix is available 
c 
c          output: 
c             x:       vector of minimizing solution 
c 
cx      implicit  real * 8  ( a-h , o-z ) 
cc      dimension  x(50) , dx(50) , g(50) , g0(50) , y(50) 
cc      dimension  h(50,50) , wrk(50) , s(50) 
cc      dimension  ht(50,50),hf(50,50)
cx      dimension  x(n) , dx(n) , g(m) , g0(n) , y(n)
cx      dimension  h(n,n) , wrk(n) , s(n) 
cx      dimension  hf(n,n)
c---
cx      dimension  xx(nni), hfi(n,n)
cx      dimension  id(nlmax), rmd(nlmax), xx1(n,nlmax) 
c---
      integer :: nni, n, m, nlmax, id(nlmax), nl
      REAL(8) :: xx(nni), x(n), g(m), rmd(nlmax), xx1(n,nlmax), h(n,n),
     1           hf(n,n), hfi(n,n)
      real(8) :: f, aic, t, tau1, tau2, eps1, eps2
      external funct
cxx      common     / ccc /  isw , ipr 
cxx      common     / ddd /   f , aic 
cxx      common t,nn,nfunct 
      common /ddd1/ f, aic 
      common /momori/ t,nn,nfunct 
      data  tau1 , tau2  /  1.0d-5 , 1.0d-5  / 
      data  eps1 , eps2  / 1.0d-5 , 1.0d-5  / 
      real(8) :: dx(n), g0(n), y(n), wrk(n), s(n), ramda, const1, sum,
     1           s1, s2, ss, ds2, stem, gtem, ed, xm, xmb, hdet
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
cxx      isw = 0
c
cc      call  funct( n,x,xm,g,ig ) 
      call  funct( nni,xx,m,x,xm,g,ig ) 
c 
cc      write( 6,340 )     xm 
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 340
         xx1(1,nl) = xm
      end if
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
      if( ic .eq. 1 .and. icc .eq. 1 )  go to 120 
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
      call  linear6( xx,nn,x,s,ramda,ed,m,ig,funct,id,rmd,xx1,
     & nl,nlmax)
cc      write( 6,330 )     ramda , f , s1 , s2 
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 330
         rmd(nl) = ramda
         xx1(1,nl) = f
         xx1(2,nl) = s1
         xx1(3,nl) = s2
      end if
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
      call  funct( nni,xx,m,x,xm,g,ig ) 
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
cc      write( 6,610 )     (x(i),i=1,n) 
cc      write( 6,601 ) 
cc      write( 6,610 )     (g(i),i=1,n) 
c     the followings are inserted by y.o.(13/10/82) 
cc      write(6,602) 
cc      do 651 i=1,n
cc  651 write(6,610) (h(i,j),j=1,n) 
cc      if(nfunct.ne.6) return 
cc      call fisher(x,n,hf)
cc      write(6,605) 
cx  605 format(/1h ,'***  fisher matrix  ***') 
cc      do 654 i=1,n 
cc  654 write(6,604) (hf(i,j),j=1,n)
cc      call invdet(hf,hdet,n,50)
cc      write(6,606) 
cx  606 format(/1h ,'***  inverse fisher  ***') 
cc      do 655 i=1,n 
cc  655 write(6,604) (hf(i,j),j=1,n) 
c
      if( nl.lt.nlmax ) then
         nl = nl+1
         id(nl) = 600
      end if
      if(nfunct.eq.6) then
         call fisher(x,m,hf)
         do 654 i=1,n 
            do 653 j=1,n
               hfi(i,j) = hf(i,j)
  653       continue
  654    continue
         call invdet(hfi,hdet,n,n)
      end if 
c
cx  602 format(/1h ,'***  estimated inverse hessian  ***') 
cx  604 format(1h ,8d15.5/(1h ,8d15.5)) 
cx  603 format(/1h ,'***  inverse hessian  ***') 
c-----------------------------------------------------------------------
      return 
cx  330 format( 1h ,'lambda =',d15.7,5x,'-LL =',d23.15,2x,d9.2,2x,d9.2) 
cx  340 format( 1h ,4x,'initial (-1)*Log-Likelihood =',d23.15) 
cx  600 format(/1h ,'-----  x  -----' ) 
cx  601 format(/1h ,'***  gradient  ***' ) 
cx  610 format( 1h ,10d13.5 ) 
      end 
c 
c 
c 
      subroutine invdet(x,xdet,mm,mj) 
c----------------------------------------------------------------------- 
c       the inverse and determinant of x computation 
c----------------------------------------------------------------l 11890 
c       inputs: 
c          x:     mm*mm square matrix 
c          mm:    dimension of x 
c          mj:    absolute dimension of x in the main program 
c 
c       outputs: 
c          x:     inverse of x 
c          xdet:  determinant of x 
c 
cx      implicit  real * 8  ( a-h , o-z ) 
cx      dimension x(mj,mj) 
cc      dimension  ids(100)
cx      dimension  ids(mm)
      integer :: mm, mj
      real(8) :: x(mj,mj), xdet 
      integer :: ids(mm)
      real(8) :: xmaxp, xc
      xdet = 1.0d00 
      do 10 l=1,mm 
c     pivoting at l-th stage 
      xmaxp=0.10000d-10 
      maxi=0 
      do 110 i=l,mm 
cx   1 if( abs(xmaxp) .ge. abs(x(i,l)) )     go to 110 
      if( abs(xmaxp) .ge. abs(x(i,l)) )     go to 110 
      xmaxp=x(i,l) 
      maxi=i 
  110 continue 
      ids(l)=maxi 
      if(maxi.eq.l) go to 120 
      if(maxi.gt.0) go to 121 
      xdet = 0.0d00 
      go to 140 
c     row interchange 
  121 do 14 j=1,mm 
      xc=x(maxi,j) 
      x(maxi,j)=x(l,j) 
cx   14 x(l,j)=xc
      x(l,j)=xc
   14 continue
      xdet=-xdet 
c 120 xdet=xdet*xmaxp 
  120 continue 
      xc = 1.0d00 / xmaxp 
      x(l,l)=1.0d00 
      do 11 j=1,mm 
cx   11 x(l,j)=x(l,j)*xc 
      x(l,j)=x(l,j)*xc
   11 continue
      do 12 i=1,mm 
      if(i.eq.l) go to 12 
      xc=x(i,l) 
      x(i,l) = 0.0d00 
      do 13 j=1,mm 
cx   13 x(i,j)=x(i,j)-xc*x(l,j) 
      x(i,j)=x(i,j)-xc*x(l,j)
   13 continue
   12 continue 
   10 continue 
      if(mm.gt.1) go to 123 
      go to 140 
c     column interchange 
  123 mm1=mm-1 
      do 130 j=1,mm1 
      mmj=mm-j 
      jj=ids(mmj) 
      if(jj.eq.mmj) go to 130 
      do 131 i=1,mm 
      xc=x(i,jj) 
      x(i,jj)=x(i,mmj) 
cx  131 x(i,mmj)=xc 
      x(i,mmj)=xc
  131 continue
  130 continue 
  140 return 
      end 
c 
c 
cc      subroutine func5(n,b,f,h,ifg) 
      subroutine func5(nni,xx,n,b,f,h,ifg)
c----------------------------------------------------------------------- 
c     likelihood function of exp-decay poisson process 
c     lammbda = a1 + a2*exp(-a3*t) 
c----------------------------------------------------------------------- 
cx      implicit real * 8 (a-h,o-z)
      integer :: nni, n, ifg
      real(8) :: xx(nni), b(n), f, h(n)
      real(8) :: t, ff, aic
cc      common/xyod/xdumy,xx(19999) 
cc      common/y/y(50)
cxx      common t,nn,nfunct 
cxx      common/ddd/ff,aic 
      common /momori/ t,nn,nfunct 
      common /ddd1/ ff,aic 
cc      dimension b(50),h(50),g(50)
cx      dimension b(n),h(n),g(n)
cx      dimension xx(nni)
      real(8) :: g(n), a1, a2, a3, uni, f1, gg1, gg2, gg3,
     1          ramdai, sasump, gs1, gs2, gs3
      ifg=0 
      a1=b(1)**2 
      a2=b(2)**2 
      a3=b(3)**2 
      f1=0.0 
      gg1=0.0 
      gg2=0.0 
      gg3=0.0 
      do 20 i=1,nn 
      ramdai=a1+a2*exp(-a3*xx(i)) 
      if(ramdai.le.0.0) go to 50 
      f1=f1+log(ramdai) 
      uni=1.0d00 
      gg1=gg1+uni/ramdai 
      gg2=gg2+exp(-a3*xx(i))/ramdai 
      gg3=gg3-a2*xx(i)*exp(-a3*xx(i))/ramdai 
   20 continue 
c 
c 
      sasump=(1.0d0-exp(-a3*t))/a3 
      gs1=t 
      gs2=sasump 
      gs3=-sasump*a2/a3+a2/a3*t*exp(-a3*t) 
      go to 240 
c 
c 
   50 continue 
      ifg=1 
      f=1.0d30 
      return 
c 
c 
  240 continue 
      f=f1-a1*t-a2*sasump 
      g(1)=gg1-gs1 
      g(2)=gg2-gs2 
      g(3)=gg3-gs3 
      f=-f 
      h(1)=-g(1) 
      h(1)=h(1)*2.0d00*b(1) 
      h(2)=-g(2) 
      h(2)=h(2)*2.0d00*b(2) 
      h(3)=-g(3) 
      h(3)=h(3)*2.0d00*b(3) 
      ff=f 
      na=0 
      do 800 i=1,n 
cx  800 if(b(i).ne.0.0) na=na+1 
      if(b(i).ne.0.0) na=na+1
  800 continue
      aic=ff+na 
cx    3 format(1h ,110x,d18.10) 
cx    1 format(1h ,7d18.10) 
      return 
      end 
c 
c 
cc      subroutine func6(n,b,f,h,ifg)
      subroutine func6(nni,xx,n,b,f,h,ifg) 
c----------------------------------------------------------------------- 
c     likelihood function of the omori's poisson process 
c     lammbda = a1 + a2/(a3+t)**a4 
c----------------------------------------------------------------------- 
cx      implicit real * 8 (a-h,o-z) 
      integer :: nni, n, ifg
      real(8) :: xx(nni), b(n), f, h(n)
      real(8) :: ff, aic, t, t0, t1, t2, t3
cc      common/xyod/xdumy,xx(19999) 
cxx      common/ddd/ff,aic 
      common /ddd1/ ff,aic 
cc      common/y/y(50)
cxx      common/range/t0,t1,t2,t3
      common /range1/ t0,t1,t2,t3
cc      common/grad/g(50) 
cxx      common t,nn,nfunct
      common /momori/ t,nn,nfunct 

cc      dimension b(50),h(50)
cx      dimension xx(nni),b(n),h(n),g(n)
      real(8) :: g(n), a1, a2, a3, a4, uni, f1, gg1, gg2,
     1           gg3, gg4, ramdai, sasump, gs1, gs2, gs3, gs4
      ifg=0
      a1=b(1)**2
      a2=b(2)**2
      a3=b(3)**2
      a4=b(4)**2
      if(a4.gt.5.0d00) go to 119
      if(a3.gt.10000.) go to 119
      f1=0.0
      gg1=0.0
      gg2=0.0
      gg3=0.0
      gg4=0.0
         gs4=0.0
CCC         sasump=0.0
      do 20 i=1,nn
      ramdai=a1+a2/(a3+xx(i))**a4
      if(ramdai.le.0.0) go to 50
      f1=f1+log(ramdai)
      uni=1.0d00
      gg1=gg1+uni/ramdai
      gg2=gg2+uni/ramdai/(a3+xx(i))**a4
      gg3=gg3-a2*a4/ramdai/(a3+xx(i))**(a4+uni)
      gg4=gg4-a2*log(a3+xx(i))/ramdai/(a3+xx(i))**a4
   20 continue
c
c
cxx      if(b(4).eq.uni) sasump=log(t1+a3)-log(t0+a3)
cxx      if(b(4).gt.uni) sasump=
cxx     & (uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4)
cxx      if(b(4).lt.uni) sasump=
cxx     & ((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      if(b(4).gt.uni) then
         sasump=(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4)
      else if (b(4).lt.uni) then
         sasump=((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      else
         sasump=log(t1+a3)-log(t0+a3)
      end if
      gs1=t1-t0
      gs2=sasump
      if(b(4).ne.uni) gs3=a2*(uni/(t1+a3)**a4-uni/(t0+a3)**a4)
      if(b(4).eq.uni) gs3=a2*(uni/(t1+a3)-uni/(t0+a3))
      if(b(4).gt.uni) gs4=
     & a2/(uni-a4)**2
     &  *(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))
     & +a2/(uni-a4)
     &  *(-log(t1+a3)/(t1+a3)**(a4-uni)+log(t0+a3)/(t0+a3)**(a4-uni))
      if(b(4).lt.uni) gs4=
     & a2/(uni-a4)**2
     &  *((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))
     & +a2/(uni-a4)
     &  *(-(t1+a3)**(uni-a4)*log(t1+a3)+(t0+a3)**(uni-a4)*log(t0+a3))
      go to 240
c
c
   50 continue
      ifg=1
      f=1.0d30
      return
c
c
  240 continue
      f=f1-a1*(t1-t0)-a2*sasump
      g(1)=gg1-gs1
      g(2)=gg2-gs2
      g(3)=gg3-gs3
      g(4)=gg4-gs4
      if(b(4).eq.1.0d00) g(4)=0.0
      f=-f
      h(1)=-g(1)
      h(1)=h(1)*2.0d00*b(1)
      h(2)=-g(2)
      h(2)=h(2)*2.0d00*b(2)
      h(3)=-g(3)
      h(3)=h(3)*2.0d00*b(3)
      h(4)=-g(4)
      h(4)=h(4)*2.0d00*b(4)
      ff=f
      na=0
      do 800 i=1,n
cx  800 if(b(i).ne.0.0) na=na+1
      if(b(i).ne.0.0) na=na+1
  800 continue
      aic=ff+na
cx    3 format(1h ,110x,d18.10)
cx    1 format(1h ,7d18.10)
      return
  119 continue
      f=1.0d50
      ifg=1
      return
      end
c
cc      subroutine func9(n,b,f,h,ifg)
      subroutine func9(nni,xx,n,b,f,h,ifg) 
c----------------------------------------------------------------------- 
c     likelihood function of the omori's poisson process 
c     lammbda = a1 + a2/(a3+t)**a4 + a5/(a6+t-t2)**a7 * i(t.gt.t2) 
c     a model for detecting the second aftershock 
c     this is not succeeded at the moment (82/11/22) 
c----------------------------------------------------------------------- 
cx      implicit real * 8 (a-h,o-z)
      integer :: nni, n, ifg
      real(8) :: xx(nni), b(n), f, h(n)
      real(8) :: ff, aic, t, t0, t1, t2, t3
cc      common/xyod/xdumy,xx(19999)
cxx      common/ddd/ff,aic 
      common /ddd1/ ff,aic 
cc      common/y/y(50)
cxx      common/range/t0,t1,t2,t3 
      common /range1/ t0,t1,t2,t3 
cc      common/grad/g(50)
cxx      common t,nn,nfunct 
      common /momori/ t,nn,nfunct 
cc      dimension b(50),h(50)
cx      dimension xx(nni),b(n),h(n),g(n)
      real(8) :: g(n), a1, a2, a3, a4, a5, a6, a7, uni, f1, t4,
     1           gg1, gg2, gg3, gg4, gg5, gg6, gg7, ramdai,
     2           sasump, sasumq, gs1, gs2, gs3, gs4, gs5, gs6, gs7
      ifg=0 
      a1=0.0 
      a1=b(1)**2 
      if(b(1).ne.0.0) a1=exp(b(1)) 
      a2=0.0 
      a2=b(2)**2 
      if(b(2).ne.0.0) a2=exp(b(2)) 
      a3=0.0 
      a3=b(3)**2 
      if(b(3).ne.0.0) a3=exp(b(3)) 
      a4=0.0 
      a4=b(4)**2 
      if(b(4).ne.0.0) a4=exp(b(4))
      a5=0.0 
      a5=b(5)**2 
      if(b(5).ne.0.0) a5=exp(b(5)) 
      a6=0.0 
      a6=b(6)**2 
      if(b(6).ne.0.0) a6=exp(b(6)) 
      a7=0.0 
      a7=b(7)**2 
      if(b(7).ne.0.0) a7=exp(b(7)) 
      if(b(7).eq.0.0) a7=exp(b(4)) 
      if(b(6).eq.0.0) a6=exp(b(3)) 
c     if(a4.gt.5.0d00) go to 119 
      if(a3.gt.10000.) go to 119 
      uni=1.0d00 
      f1=0.0 
      gg1=0.0 
      gg2=0.0 
      gg3=0.0 
      gg4=0.0 
      gg5=0.0 
      gg6=0.0 
      gg7=0.0 
         gs4=0.0
         gs7=0.0
      if(a7*log(a6+t1-t2).gt.150.) go to 50 
      if(a4*log(a3).lt.-150.) go to 50 
      if(a4*log(a3+t1).gt.150.) go to 50 
      do 20 i=1,nn 
      ramdai=a1+a2/(a3+xx(i))**a4 
      if(xx(i).gt.t2) 
     &ramdai=a1+a2/(a3+xx(i))**a4+a5/(a6+xx(i)-t2)**a7 
      if(ramdai.le.0.0) go to 50 
      gg1=gg1+uni/ramdai 
      gg2=gg2+uni/ramdai/(a3+xx(i))**a4 
      gg3=gg3-a2*a4/ramdai/(a3+xx(i))**(a4+uni) 
      gg4=gg4-a2*log(a3+xx(i))/ramdai/(a3+xx(i))**a4 
c 
      if(xx(i).le.t2) go to 10 
      if(ramdai.le.0.0) go to 50 
      gg5=gg5+uni/ramdai/(a6+xx(i)-t2)**a7 
      gg6=gg6-a5*a7/ramdai/(a6+xx(i)-t2)**(a7+uni) 
      gg7=gg7-a5*log(a6+xx(i)-t2)/ramdai/(a6+xx(i)-t2)**a7 
c 
   10 continue 
      f1=f1+log(ramdai) 
   20 continue 
c 
c 
cxx      if(b(4).eq.uni) sasump=log(t1+a3)-log(t0+a3) 
cxx      if(b(4).gt.uni) sasump= 
cxx     & (uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4) 
cxx      if(b(4).lt.uni) sasump= 
cxx     & ((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      if(b(4).gt.uni) then
         sasump=(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4) 
      else if(b(4).lt.uni) then
         sasump=((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      else
         sasump=log(t1+a3)-log(t0+a3) 
      end if
      gs1=t1-t0 
      gs2=sasump 
      if(b(4).ne.uni) gs3=a2*(uni/(t1+a3)**a4-uni/(t0+a3)**a4) 
      if(b(4).eq.uni) gs3=a2*(uni/(t1+a3)-uni/(t0+a3)) 
      if(b(4).gt.uni) gs4= 
     & a2/(uni-a4)**2 
     &  *(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni)) 
     & +a2/(uni-a4) 
     &  *(-log(t1+a3)/(t1+a3)**(a4-uni)+log(t0+a3)/(t0+a3)**(a4-uni)) 
      if(b(4).lt.uni) gs4= 
     & a2/(uni-a4)**2 
     &  *((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4)) 
     & +a2/(uni-a4) 
     &  *(-(t1+a3)**(uni-a4)*log(t1+a3)+(t0+a3)**(uni-a4)*log(t0+a3)) 
c 
      t3=t1-t2 
      t4=0.0
cxx      if(b(7).eq.uni) sasumq=log(t3+a6)-log(t4+a6) 
cxx      if(b(7).gt.uni) sasumq= 
cxx     & (uni/(t3+a6)**(a7-uni)-uni/(t4+a6)**(a7-uni))/(uni-a7) 
cxx      if(b(7).lt.uni) sasumq= 
cxx     & ((t3+a6)**(uni-a7)-(t4+a6)**(uni-a7))/(uni-a7)
      if(b(7).gt.uni) then
         sasumq=(uni/(t3+a6)**(a7-uni)-uni/(t4+a6)**(a7-uni))/(uni-a7)
      else if(b(7).lt.uni) then
         sasumq=((t3+a6)**(uni-a7)-(t4+a6)**(uni-a7))/(uni-a7)
      else
         sasumq=log(t3+a6)-log(t4+a6) 
      end if
      gs5=sasumq 
      if(b(7).ne.uni) gs6=a5*(uni/(t3+a6)**a7-uni/(t4+a6)**a7) 
      if(b(7).eq.uni) gs6=a5*(uni/(t3+a6)-uni/(t4+a6)) 
      if(b(7).gt.uni) gs7= 
     & a5/(uni-a7)**2 
     &  *(uni/(t3+a6)**(a7-uni)-uni/(t4+a6)**(a7-uni)) 
     & +a5/(uni-a7) 
     &  *(-log(t3+a6)/(t3+a6)**(a7-uni)+log(t4+a6)/(t4+a6)**(a7-uni)) 
      if(b(7).lt.uni) gs7= 
     & a5/(uni-a7)**2 
     &  *((t3+a6)**(uni-a7)-(t4+a6)**(uni-a7)) 
     & +a5/(uni-a7) 
     &  *(-(t3+a6)**(uni-a7)*log(t3+a6)+(t4+a6)**(uni-a7)*log(t4+a6)) 
c 
      go to 240 
c 
c 
   50 continue 
      ifg=1 
      f=1.0d30 
      return 
c 
c 
  240 continue 
      f=f1-a1*(t1-t0)-a2*sasump-a5*sasumq 
      g(1)=gg1-gs1 
      g(2)=gg2-gs2 
      g(3)=gg3-gs3 
      g(4)=gg4-gs4 
      g(5)=gg5-gs5 
      g(6)=gg6-gs6 
      g(7)=gg7-gs7 
      f=-f 
      h(1)=-g(1) 
      h(1)=h(1)*a1 
      if(b(1).eq.0.0) h(1)=0.0 
      h(2)=-g(2) 
      h(2)=h(2)*a2 
      if(b(2).eq.0.0) h(2)=0.0 
      h(3)=-g(3) 
      if(b(6).eq.0.0) h(3)=-(g(3)+g(6)) 
      h(3)=h(3)*a3 
      if(b(3).eq.0.0) h(3)=0.0 
      h(4)=-g(4) 
      if(b(7).eq.0.0) h(4)=-(g(4)+g(7)) 
      h(4)=h(4)*a4 
      if(b(4).eq.0.0) h(4)=0.0 
      h(5)=-g(5) 
      h(5)=h(5)*a5 
      if(b(5).eq.0.0) h(5)=0.0 
      h(6)=-g(6) 
      h(6)=h(6)*a6 
      if(b(6).eq.0.0) h(6)=0.0 
      h(7)=-g(7) 
      h(7)=h(7)*a7 
      if(b(7).eq.0.0) h(7)=0.0 
      ff=f 
      na=0 
      do 800 i=1,n 
cx  800 if(b(i).ne.0.0) na=na+1 
      if(b(i).ne.0.0) na=na+1
  800 continue
      aic=ff+na 
cx    3 format(1h ,110x,d18.10) 
cx    1 format(1h ,7d18.10) 
      return 
  119 continue 
      f=1.0d50 
      ifg=1 
      return 
      end 
c 
c 
c 
cc      subroutine func10(n,b,f,h,ifg)
      subroutine func10(nni,xx,n,b,f,h,ifg) 
c----------------------------------------------------------------------- 
c     likelihood function of the omori's poisson process 
c     lammbda = a1 + a2/(a3+t)**a4 + a5/(a6+t-t2)**a7 * i(t.gt.t2) 
c                  + a8/(a9+t-t3)**a10 * i(t.gt.t3) 
c     a model for detecting the third aftershock
c     this is not succeeded at the moment (82/11/22) 
c----------------------------------------------------------------------- 
cx      implicit real * 8 (a-h,o-z) 
cc      common/xyod/xdumy,xx(19999) 
      integer :: nni, n, ifg
      real(8) :: xx(nni), b(n), f, h(n)
      real(8) :: ff, aic, t, t0, t1, t2, t3
cxx      common/ddd/ff,aic 
      common /ddd1/ ff,aic 
cc      common/y/y(50) 
cxx      common/range/t0,t1,t2,t3 
      common /range1/ t0,t1,t2,t3 
cc      common/grad/g(50)
cxx      common t,nn,nfunct 
      common /momori/ t,nn,nfunct 
cc      dimension b(50),h(50)
cx      dimension xx(nni),b(n),h(n),g(n)
      real(8) :: g(n), a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, uni,
     1           f1, gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, gg10,
     2           ramdai, sasump, sasumq, sasumr, t02, t03, t13, t12,
     3           gs1, gs2, gs3, gs4, gs5, gs6, gs7, gs8, gs9, gs10
      do 111 i=1,10 
      if(b(i).gt.270.d0) go to 50 
  111 continue 
      ifg=0 
      a1=0.0 
      a1=b(1)**2 
      if(b(1).ne.0.0) a1=exp(b(1)) 
      a2=0.0 
      a2=b(2)**2 
      if(b(2).ne.0.0) a2=exp(b(2)) 
      a3=0.0 
      a3=b(3)**2 
      if(b(3).ne.0.0) a3=exp(b(3)) 
      a4=0.0 
      a4=b(4)**2 
      if(b(4).eq.0.0) a4=exp(b(10))
      if(b(4).ne.0.0) a4=exp(b(4)) 
      a5=0.0 
      a5=b(5)**2 
      if(b(5).ne.0.0) a5=exp(b(5)) 
      a6=0.0 
      a6=b(6)**2 
      if(b(6).ne.0.0) a6=exp(b(6)) 
      a7=0.0 
      a7=b(7)**2 
      if(b(7).ne.0.0) a7=exp(b(7)) 
      a8=0.0 
      a8=b(8)**2 
      if(b(8).ne.0.0) a8=exp(b(8)) 
      a9=0.0 
      a9=b(9)**2 
      if(b(9).ne.0.0) a9=exp(b(9)) 
      a10=0.0 
      a10=b(10)**2 
      if(b(10).ne.0.0) a10=exp(b(10)) 
      if(b(7).eq.0.0) a7=exp(b(4)) 
      if(b(10).eq.0.0) a10=exp(b(4)) 
      if(b(6).eq.0.0) a6=exp(b(3)) 
      if(b(9).eq.0.0) a9=exp(b(3)) 
      if(a3.gt.10000.) go to 119 
      uni=1.0d00 
      f1=0.0 
      gg1=0.0 
      gg2=0.0 
      gg3=0.0 
      gg4=0.0 
      gg5=0.0 
      gg6=0.0 
      gg7=0.0 
      gg8=0.0 
      gg9=0.0 
      gg10=0.0
         gs4=0.0
         gs7=0.0
         gs10=0.0
      if(a7*log(a6+t1-t2).gt.150.) go to 50 
      if(a4*log(a3).lt.-150.) go to 50 
      if(a4*log(a3+t1).gt.150.) go to 50 
      do 20 i=1,nn 
      ramdai=a1+a2/(a3+xx(i))**a4 
      if(xx(i).gt.t2) 
     &ramdai=a1+a2/(a3+xx(i))**a4+a5/(a6+xx(i)-t2)**a7 
      if(xx(i).gt.t3) 
     &ramdai=a1+a2/(a3+xx(i))**a4+a5/(a6+xx(i)-t2)**a7 
     &         +a8/(a9+xx(i)-t3)**a10 
      if(ramdai.le.0.0) go to 50 
      gg1=gg1+uni/ramdai 
      gg2=gg2+uni/ramdai/(a3+xx(i))**a4 
      gg3=gg3-a2*a4/ramdai/(a3+xx(i))**(a4+uni) 
      gg4=gg4-a2*log(a3+xx(i))/ramdai/(a3+xx(i))**a4 
c 
      if(xx(i).le.t2) go to 10 
      if(ramdai.le.0.0) go to 50 
      gg5=gg5+uni/ramdai/(a6+xx(i)-t2)**a7 
      gg6=gg6-a5*a7/ramdai/(a6+xx(i)-t2)**(a7+uni) 
      gg7=gg7-a5*log(a6+xx(i)-t2)/ramdai/(a6+xx(i)-t2)**a7 
c 
      if(xx(i).le.t3) go to 10 
      if(ramdai.le.0.0) go to 50 
      gg8=gg8+uni/ramdai/(a9+xx(i)-t3)**a10 
      gg9=gg9-a8*a10/ramdai/(a9+xx(i)-t3)**(a10+uni) 
      gg10=gg10-a8*log(a9+xx(i)-t3)/ramdai/(a9+xx(i)-t3)**a10 
c 
   10 continue 
      f1=f1+log(ramdai) 
   20 continue 
c 
c
cxx      if(b(4).eq.uni) sasump=log(t1+a3)-log(t0+a3) 
cxx      if(b(4).gt.uni) sasump= 
cxx     & (uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4) 
cxx      if(b(4).lt.uni) sasump= 
cxx     & ((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      if(b(4).gt.uni) then
         sasump=(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4) 
      else if(b(4).lt.uni) then
         sasump=((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      else
         sasump=log(t1+a3)-log(t0+a3) 
      end if
      gs1=t1-t0 
      gs2=sasump 
      if(b(4).ne.uni) gs3=a2*(uni/(t1+a3)**a4-uni/(t0+a3)**a4) 
      if(b(4).eq.uni) gs3=a2*(uni/(t1+a3)-uni/(t0+a3)) 
      if(b(4).gt.uni) gs4= 
     & a2/(uni-a4)**2 
     &  *(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni)) 
     & +a2/(uni-a4) 
     &  *(-log(t1+a3)/(t1+a3)**(a4-uni)+log(t0+a3)/(t0+a3)**(a4-uni)) 
      if(b(4).lt.uni) gs4= 
     & a2/(uni-a4)**2 
     &  *((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4)) 
     & +a2/(uni-a4) 
     &  *(-(t1+a3)**(uni-a4)*log(t1+a3)+(t0+a3)**(uni-a4)*log(t0+a3)) 
c 
      t12=t1-t2 
      t02=0.0
cxx      if(b(7).eq.uni) sasumq=log(t12+a6)-log(t02+a6) 
cxx      if(b(7).gt.uni) sasumq= 
cxx     & (uni/(t12+a6)**(a7-uni)-uni/(t02+a6)**(a7-uni))/(uni-a7) 
cxx      if(b(7).lt.uni) sasumq= 
cxx     & ((t12+a6)**(uni-a7)-(t02+a6)**(uni-a7))/(uni-a7)
      if(b(7).gt.uni) then
         sasumq= 
     &      (uni/(t12+a6)**(a7-uni)-uni/(t02+a6)**(a7-uni))/(uni-a7) 
      else if(b(7).lt.uni) then
         sasumq=((t12+a6)**(uni-a7)-(t02+a6)**(uni-a7))/(uni-a7)
      else
         sasumq=log(t12+a6)-log(t02+a6) 
      end if

      gs5=sasumq 
      if(b(7).ne.uni) gs6=a5*(uni/(t12+a6)**a7-uni/(t02+a6)**a7) 
      if(b(7).eq.uni) gs6=a5*(uni/(t12+a6)-uni/(t02+a6)) 
      if(b(7).gt.uni) gs7= 
     & a5/(uni-a7)**2 
     &  *(uni/(t12+a6)**(a7-uni)-uni/(t02+a6)**(a7-uni)) 
     & +a5/(uni-a7) 
     & *(-log(t12+a6)/(t12+a6)**(a7-uni)+log(t02+a6)/(t02+a6)**(a7-uni)) 
      if(b(7).lt.uni) gs7= 
     & a5/(uni-a7)**2 
     &  *((t12+a6)**(uni-a7)-(t02+a6)**(uni-a7)) 
     & +a5/(uni-a7) 
     & *(-(t12+a6)**(uni-a7)*log(t12+a6)+(t02+a6)**(uni-a7)*log(t02+a6)) 
c 
      t13=t1-t3 
      t03=0.0 
cxx      if(b(10).eq.uni) sasumr=log(t13+a9)-log(t03+a9) 
cxx      if(b(10).gt.uni) sasumr= 
cxx     & (uni/(t13+a9)**(a10-uni)-uni/(t03+a9)**(a10-uni))/(uni-a10) 
cxx      if(b(10).lt.uni) sasumr= 
cxx     & ((t13+a9)**(uni-a10)-(t03+a9)**(uni-a10))/(uni-a10) 
      if(b(10).gt.uni) then
         sasumr=
     &      (uni/(t13+a9)**(a10-uni)-uni/(t03+a9)**(a10-uni))/(uni-a10)
      else if(b(10).lt.uni) then
         sasumr=((t13+a9)**(uni-a10)-(t03+a9)**(uni-a10))/(uni-a10) 
      else
         sasumr=log(t13+a9)-log(t03+a9)
      end if 
      gs8=sasumr 
      if(b(10).ne.uni) gs9=a8*(uni/(t13+a9)**a10-uni/(t03+a9)**a10) 
      if(b(10).eq.uni) gs9=a8*(uni/(t13+a9)-uni/(t03+a9)) 
      if(b(10).gt.uni) gs10= 
     & a8/(uni-a10)**2 
     &  *(uni/(t13+a9)**(a10-uni)-uni/(t03+a9)**(a10-uni)) 
     & +a8/(uni-a10) 
     & *(-log(t13+a9)/(t13+a9)**(a10-uni)+log(t03+a9)/(t03+a9)**(a10-uni 
     &)) 
      if(b(10).lt.uni) gs10= 
     & a8/(uni-a10)**2 
     &  *((t13+a9)**(uni-a10)-(t03+a9)**(uni-a10)) 
     & +a8/(uni-a10) 
     & *(-(t13+a9)**(uni-a10)*log(t13+a9)+(t03+a9)**(uni-a10)*log(t03+a9 
     &)) 
c 
      go to 240 
c 
c 
   50 continue 
      ifg=1 
      f=1.0d30 
      return 
c 
c 
  240 continue 
      f=f1-a1*(t1-t0)-a2*sasump-a5*sasumq-a8*sasumr 
      g(1)=gg1-gs1 
      g(2)=gg2-gs2 
      g(3)=gg3-gs3 
      g(4)=gg4-gs4 
      g(5)=gg5-gs5 
      g(6)=gg6-gs6 
      g(7)=gg7-gs7 
      g(8)=gg8-gs8 
      g(9)=gg9-gs9 
      g(10)=gg10-gs10 
      if(b(4).eq.1.0d00) g(4)=0.0 
      if(b(7).eq.1.0d00) g(7)=0.0 
      f=-f 
      h(1)=-g(1) 
      h(1)=h(1)*a1 
      if(b(1).eq.0.0) h(1)=0.0 
      h(2)=-g(2) 
      h(2)=h(2)*a2 
      if(b(2).eq.0.0) h(2)=0.0 
      h(3)=-g(3) 
      if(b(6).eq.0.0) h(3)=-(g(3)+g(6)) 
      if(b(6).eq.0.0.and.b(9).eq.0.0) h(3)=-(g(3)+g(6)+g(9)) 
      h(3)=h(3)*a3 
      if(b(3).eq.0.0) h(3)=0.0 
      h(4)=-g(4) 
      if(b(7).eq.0.0.and.b(10).ne.0.0) h(4)=-(g(4)+g(7)) 
      if(b(7).ne.0.0.and.b(10).eq.0.0) h(4)=-(g(4)+g(10)) 
      if(b(7).eq.0.0.and.b(10).eq.0.0) h(4)=-(g(4)+g(7)+g(10)) 
      h(4)=h(4)*a4 
      if(b(4).eq.0.0) h(4)=0.0 
      h(5)=-g(5) 
      h(5)=h(5)*a5 
      if(b(5).eq.0.0) h(5)=0.0 
      h(6)=-g(6) 
      h(6)=h(6)*a6 
      if(b(6).eq.0.0) h(6)=0.0 
      h(7)=-g(7) 
      if(b(10).eq.0.0) h(7)=-(g(7)+g(10)) 
      h(7)=h(7)*a7 
      if(b(7).eq.0.0) h(7)=0.0 
      h(8)=-g(8) 
      h(8)=h(8)*a8 
      if(b(8).eq.0.0) h(8)=0.0 
      h(9)=-g(9) 
      h(9)=h(9)*a9 
      if(b(9).eq.0.0) h(9)=0.0 
      h(10)=-g(10) 
      if(b(4).eq.0) h(10)=-(g(4)+g(10)) 
      h(10)=h(10)*a10 
      if(b(10).eq.0.0) h(10)=0.0 
      ff=f 
      na=0 
      do 800 i=1,n 
cx  800 if(b(i).ne.0.0) na=na+1 
      if(b(i).ne.0.0) na=na+1
  800 continue
      aic=ff+na 
cx    3 format(1h ,110x,d18.10) 
cx    1 format(1h ,7d18.10) 
      go to 300 
  119 continue 
      f=1.0d50 
      ifg=1 
  300 continue 
      return 
      end 
c 
c 
c 
cx      function sf1(x,q) 
      double precision function sf1(x,q)
cx      implicit real * 8 (a-h,o-z)
      real(8) :: x, q
      sf1=x**(1.d0-q)/(1.d0-q) 
      return 
      end 
cx      function sf2(x,q)
      double precision function sf2(x,q)
cx      implicit real * 8 (a-h,o-z) 
      real(8) :: x, q, sf1
      sf2=sf1(x,q)*(log(x)-1.d0/(1.d0-q)) 
      return 
      end 
cx      function sf3(x,q)
      double precision function sf3(x,q)
cx      implicit real * 8 (a-h,o-z) 
      real(8) :: x, q, sf1, sf2
      sf3=sf1(x,q)*(log(x))**2-2.d0/(1.d0-q)*sf2(x,q) 
      return 
      end 
c 
c     real function gm*8 (x,q,c) 
cx      real*8 function gm (x,q,c) 
      double precision function gm (x,q,c)
cx      implicit real * 8 (a-h,o-z) 
      real(8) :: x, q, c
      real(8) :: gmi
      gm=0.0 
      if(x.eq.c) go to 20 
      gmi=x**(-q) 
      do 10 i=1,50 
      i1=i-1 
      if(i1.eq.0) i1=1 
c     gmi=gmi*x/i1 
      gmi=gmi*(x-c)/i1 
      gm=gm+(-1)**(i-1)*gmi/(i-q) 
      if(gmi/gm.lt.1.d-13) go to 20 
   10 continue 
   20 return 
      end 
c 
c     real function dgm*8 (x,q,c) 
cx      real*8 function dgm (x,q,c)
      double precision function dgm (x,q,c)
cx      implicit real*8(a-h,o-z) 
      real(8) :: x, q, c
      real(8) :: gm, dgmi
      dgm=0.0 
      if(x.eq.c) go to 30 
      dgmi=x**(-q) 
      do 10 i=1,50 
      i1=i-1 
      if(i1.eq.0) i1=1 
c     dgmi=dgmi*x/i1 
      dgmi=dgmi*(x-c)/i1 
      dgm=dgm+(-1)**i*dgmi/(i-q)**2 
      if(dgmi/dgm.lt.1.d-13) go to 20 
   10 continue 
c  20 dgm=-dgm-gm(x,q)*log(x) 
   20 dgm=-dgm-gm(x,q,c)*log(x) 
   30 continue 
      return 
      end 
c 
c 
c
       subroutine fisher(b,n,h) 
c----------------------------------------------------------------------- 
c     fisher information matrix of the modified omori model 
c     lammbda = a1 + a2/(a3+t)**a4 
c----------------------------------------------------------------------- 
cx      implicit real * 8 (a-h,o-z)
      integer :: n
      real(8) :: b(n), h(n,n)
      real(8) :: t, t0, t1, t2, t3
cxx      common/range/t0,t1,t2,t3 
cxx      common t,nn,nfunct 
      common /range1/ t0,t1,t2,t3 
      common /momori/ t,nn,nfunct 
cc      dimension b(50),h(50,50) 
cx      dimension b(n),h(n,n) 
      real(8) :: a1, a2, a3, a4, sf1, sf2, sf3
      a1=b(1)**2 
      a2=b(2)**2 
      a3=b(3)**2 
      a4=b(4)**2 
c---- careful about q = 1 in the functions --- 
      if(a1.ne.0.0) h(1,1)=(t1-t0)/a1-log(t1-t0)/a1 
      if(a1.eq.0.0) h(1,1)=1.0d0 
      h(1,2)=0.0 
      h(1,3)=0.0 
      h(1,4)=0.0 
      h(2,2)=(sf1(t1+a3,a4)-sf1(t0+a3,a4))/a2 
      h(2,3)=-a4*(sf1(t1+a3,a4+1)-sf1(t0+a3,a4+1)) 
      h(2,4)=-(sf2(t1+a3,a4)-sf2(t0+a3,a4)) 
      h(3,3)=a2*a4**2*(sf1(t1+a3,a4+2)-sf1(t0+a3,a4+2)) 
      h(3,4)=a2*a4*(sf2(t1+a3,a4+1)-sf2(t0+a3,a4+1)) 
      h(4,4)=a2*(sf3(t1+a3,a4)-sf3(t0+a3,a4)) 
cx      do 10 i=1,4 
      do 11 i=1,4
      do 10 j=i,4 
cx   10 h(j,i)=h(i,j) 
      h(j,i)=h(i,j)
   10 continue
   11 continue
      return 
      end 
c 
c 
c 
cc      subroutine sizes(n,x)
      subroutine sizes(n,x,kn,t00,ti,ak,c,p,cls) 
cx      implicit real * 8 (a-h,o-z) 
cc      dimension x(50),cls(20),ti(20) 
cc      dimension ak(20),p(20),c(20) 
      integer :: n, kn
      real(8) :: x(n), t00, ti(kn), ak(kn), c(kn), p(kn), cls(kn)
      real(8) :: t0, t1, t2, t3
cx      dimension x(n),cls(kn),ti(kn)
cx      dimension ak(kn),p(kn),c(kn) 
cxx      common/range/t0,t1,t2,t3 
      common /range1/ t0,t1,t2,t3 
      ti(1)=t2 
cc      ti(2)=t3 
      if(kn.ge.2) ti(2)=t3
cc      kn=(n-1)/3 
      do 10 k=1,kn 
      ak(k)=x(3*k-1) 
      c(k)= x(3*k) 
      if(c(k).eq.0.0) c(k)=c(k-1) 
      p(k)= x(3*k+1) 
      if(p(k).eq.0.0) p(k)=p(k-1) 
   10 continue 
      cls(1)=ak(1)*((t1+c(1))**(1-p(1))-c(1)**(1-p(1)))/(1-p(1)) 
      if(p(1).eq.1.d0) cls(1)=ak(1)*(log(t1+c(1))-log(c(1))) 
      if(kn.eq.1) go to 25 
      do 20 k=2,kn 
      if(p(k).eq.1.d0) then 
      cls(k)=ak(k)*(log(t1-ti(k-1)+c(k))-log(c(k))) 
      else 
      cls(k)=ak(k)*((t1-ti(k-1)+c(k))**(1-p(k))-c(k)**(1-p(k)))/(1-p(k)) 
      end if 
   20 continue 
   25 continue 
cc      write(6,3) 
cc      write(6,4) t0,ak(1),c(1),p(1),cls(1)
      t00=t0
cc      if(kn.eq.1) go to 35 
cc      do 30 k=2,kn 
cc      write(6,4) ti(k-1),ak(k),c(k),p(k),cls(k) 
cc   30 continue 
cc   35 continue 
cx    3 format(1h ,'    ti         k          c         p         cls') 
cx    4 format(1h ,5e11.4) 
      return 
      end 
 
