cc      program respoi
      subroutine respoif(z,xmg1,dep1,xp1,yp1,nd,xini,npara,zts,zte,
     &                       tstart,amx1,xmg,dep,xp,yp,ntstar,xx,x,nn)
c
      include 'sapp_f.h'
c
c-----------------------------------------------------------------------
c     residual of modified Omori Poisson process
c     copied from rpp (for ETAS) only subroutine residual is changed.
c-----------------------------------------------------------------------
cx      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=19999, npara=5)
cc      common/xyod/xx(ldata),xmg(ldata),xmag0
cc      common/param/xini(npara),n
cc      common /range/tstart,ntstar
cc      common t,nn,mm,iappr,nfunct
cx      dimension z(nd), dep(nd), xmg(nd), xp(nd), yp(nd), xini(npara)
cx      dimension xx(nd), x(nd)
cx      dimension xmg1(nd), dep1(nd), xp1(nd), yp1(nd)
      integer :: nd, npara, ntstar, nn
      real(8) :: z(nd), xmg1(nd), dep1(nd), xp1(nd), yp1(nd), 
     1           xini(npara), zts, zte, tstart, amx1, xmg(nd), 
     2           dep(nd), xp(nd), yp(nd), xx(nd), x(nd), t
c
      do 9 i=1,nd
         xmg(i) = xmg1(i)
         dep(i) = dep1(i)
         xp(i) = xp1(i)
         yp(i) = yp1(i)
    9 continue
c
cc      call input
      call input1(z,xmg,dep,xp,yp,nd,zts,zte,tstart,ntstar,amx1,xx,nn,t)
c
cc      do 10 i=2,nn
cc      if(xx(i).ge.xx(i-1)) go to 10
cc      write(6,*) 'reverse occurrence time'
cc      write(6,*) i,xx(i),xx(i-1),xmg(i),xmg(i-1)
cc   10 continue
cc      write(6,6) nfunct
cc      write(6,4) t,nn,mm
cc      write(6,5) xmag0
cc      write(6,*)
cx    6 format(1h ,' funct = ',i5)
cx    5 format(1h ,'reference magnitudes; xmag0',5x,f10.4)
cx    4 format(1h ,'(T,nn,mm) =',5x,f10.4,2i6)
cx    3 format(1h ,10f12.4/(1h ,10f12.4))
cx    2 format(f10.2,2i10)
cx    1 format(8f10.2)
c
cc      call residual
cx      call presidual(xini,npara,tstart,ntstar,t,xx,x,nn)
      call presidual(xini,npara,tstart,xx,x,nn)
c
cx   20 continue
      return
      end
c***********************************************************************
cc      subroutine input
      subroutine input1(z,amg,dep,xp,yp,nd,zts,zte,tstart,ntstar,amx1,
     &                  xx,nn,t)
c
c       Reading parameters
c
cx      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=19999, npara=5)
cc      common/param/xini(npara),n
cc      dimension z(ldata),amg(ldata)
cc      common /xyod/xx(ldata),xmg(ldata),xmag0
cc      common/epi/ xp(ldata),yp(ldata)
cc      common /fukasa/dep(ldata)
cx      dimension z(nd),xx(nd),amg(nd)
cx      dimension xp(nd),yp(nd),dep(nd)
      integer :: nd, ntstar, nn
      real(8) :: z(nd), amg(nd), dep(nd), xp(nd), yp(nd), zts, zte,
     1           tstart, amx1, xx(nd), t
cc      common t,nn,mm,iappr,nfunct
cc      common /range/tstart,ntstar
cc      character*60 hypodata
cc      open(unit=1,file='./aftpoi.open')
c     read(1,112) hypodata
cx  112 format(a)
cc      read(1,*) nfunct,iappr
cc      read(1,*) zts,zte,tstart
cc      read(1,*) amx1,xmag0
cc      n=5
cc      read(1,*) (xini(i),i=1,n)
cc      close(unit=1)
cc      write(6,*) ' cutoff-mag/  ref-mag/    t_0  /     T   / dead-time.'
cc      write(6,5) amx1,xmag0,zts,zte,tstart
cx    1 format(1h ,20a4)
cx    4 format(8f10.3)
cx    5 format(1h ,8f10.3)
cx    6 format(3i10)
cx    2 format(f10.0,i10)
cx    3 format(20a4)
c
cc      call zisin(t,nd,z,amg,dep,xp,yp,hypodata)
c
      t=zte-zts
      tstart=tstart-zts
      nnn=0
      nn=0
      ntstar=0
      do 10 i=1,nd
      if(amg(i).lt.amx1) go to 10
      if(z(i).lt.zts.or.z(i).gt.zte) go to 10
      nn=nn+1
      if(z(i).lt.tstart) ntstar=nn
      xx(nn)=z(i)-zts
cc      xmg(nn)=amg(i)
      amg(nn)=amg(i)
      dep(nn)=dep(i)
      xp(nn)=xp(i)
      yp(nn)=yp(i)
   10 continue
cc      write(6,*) 'input hypocenter data'
cc      write(6,7) (i,xx(i),xmg(i),i=1,nn)
cx    7 format(3(i5,f12.5,f4.1))
      mm=nd
cc      write(6,*) 'read #data; selected #data; #tstart; M_0'
cc      write(6,8)  nd, nn, ntstar, amx1
cx    8 format(i10,i16,i9,f6.2)
      if(zte.ne.0.0) t=zte
      t=zte-zts
      return
      end
c***********************************************************************
cc      subroutine residual
cx      subroutine presidual(a,npara,tstart,ntstar,t,xx,x,nn)
      subroutine presidual(a,npara,tstart,xx,x,nn)
cx      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=19999, npara=5)
cc      common/param/a(npara),n
cc      common /range/tstart,ntstar
cc      common t,nn,mm,iappr,nfunct
cc      common/xyod/xx(ldata),xmg(ldata),xmag0
cc      common/epi/ xp(ldata),yp(ldata)
cc      common /fukasa/dep(ldata)
cc      dimension x(ldata)
cx      dimension a(npara)
cx      dimension xx(nn), x(nn)
      integer :: npara, nn
      real(8) :: a(npara), tstart, xx(nn), x(nn)
      real(8) :: chtsta, ft, func41, func4p
cx      func41(t,a3)=(log(t+a3)-log(a3))
cx      func4p(t,a3,a5)=1.d0/(1-a5)*((t+a3)**(1-a5)-a3**(1-a5))
c
      chtsta=a(1)*tstart
      ft=0.0
c     do 30 j=1,ntstar
      if(a(5).eq.1.d0) ft=ft+func41(tstart,a(3))
      if(a(5).ne.1.d0) ft=ft+func4p(tstart,a(3),a(5))
c  30 continue
      chtsta=chtsta+a(2)*ft
c
      x(1)=a(1)*xx(1)-chtsta
*     do 10 i=2,nn
      do 10 i=1,nn
      ft=0.0
      if(a(5).eq.1.d0) ft=ft+func41(xx(i),a(3))
      if(a(5).ne.1.d0) ft=ft+func4p(xx(i),a(3),a(5))
      x(i)=a(1)*xx(i)+a(2)*ft-chtsta
   10 continue
c
c     write(6,1001) (i,xmg(i),x(i),i=1,nn)
cc      write(6,1005) (i-ntstar,xp(i),yp(i),xmg(i),xx(i),
cc     &                                         dep(i),x(i),i=1,nn)
cc      open(unit=1,file='work.res')
c     write(1,1003) (i-ntstar,xp(i),yp(i),xmg(i),x(i),xx(i),i=1,nn)
cc      write(1,1005) (i-ntstar,xp(i),yp(i),xmg(i),xx(i),
cc     &                                         dep(i),x(i),i=1,nn)
cc      close(unit=1)
cx 1001 format(i5,24x,2f12.5,5x)
cx 1003 format(i5,2f12.5,f4.1,2f15.5)
cx 1004 format(i5,5f12.5,5x)
cx 1005 format(i5,2f12.5,f12.1,f12.5,f8.2,2x,f12.5)
      return
      end
      double precision function func41(t,a3)
      real(8) :: t, a3
      func41=(log(t+a3)-log(a3))
      return
      end
      double precision function func4p(t,a3,a5)
      real(8) :: t, a3, a5
      func4p=1.d0/(1-a5)*((t+a3)**(1-a5)-a3**(1-a5))
      return
      end

