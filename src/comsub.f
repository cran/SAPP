      INCLUDE 'sapp_f.h'
C------------------------------------------------- 
C     cycle   ---  ptspec, linlin
C     comfac  ---  linsim, linlin, simbvh
C     initl   ---  linsim, simbvh
C     unifor  ---  linsim, simbvh
C     ptrend  ---  linsim, simbvh
C------------------------------------------------- 
cc      subroutine cycle(xx,nn,prd)
      subroutine cycle(xx,nn,prd,prb,r1,rwx,rwy,phs)
cx      implicit real*8(a-h,o-z)
cc      dimension xx(10000)
cx      dimension xx(nn)
      integer :: nn
      real(8) :: xx(nn), prd, prb, r1, rwx, rwy, phs
      real(8) :: pai, r2
      rwx=0.0
      rwy=0.0
      pai=3.14159265358979d0
      do 10 i=1,nn
      rwy=rwy+sin(2.d0*pai*xx(i)/prd)
      rwx=rwx+cos(2.d0*pai*xx(i)/prd)
   10 continue
      r2=rwx**2+rwy**2
      r1=sqrt(r2)
      phs=acos(rwx/r1)
      phs=phs/2.d0/pai*360.d0
      prb=0.0
      if(r2/nn.lt.100.d0) prb=exp(-r2/nn)
cc      write(6,1) prb,r1,rwx,rwy,phs
cx    1 format(1h ,' rayleigh prob =',d13.5,'   distance =',d13.5,
cx     &           '    rwx = ',d13.5,'     rwy = ',d13.5,
cx     &           '    phs = ',d13.5)
      return
      end
c
c
cc      subroutine comfac(lf)
      subroutine comfac(kmax,lf)
c
c     calculation of the combination factorials
c     n & m should be interpreted to be n-1 & m-1, respectively.
c
cc      dimension lf(51,51)
cx      dimension lf(kmax,kmax)
      integer :: kmax, lf(kmax,kmax)

      lf(1,1)=1
      lf(2,1)=1
      lf(2,2)=1
cc      do 10 n=3,51
      do 10 n=3,kmax
      lf(n,1)=1
      lf(n,n)=1
      do 20 m=2,n-1
      lf(n,m)=lf(n-1,m-1)+lf(n-1,m)
   20 continue
   10 continue
      return
      end
c
c
      subroutine initl(kx,ax,c,fmax)
c
c     calculation of an upper bound of a response function
c
cx      implicit real*8 (a-h,o-z)
cx      dimension ax(1)
cx      dimension ax(kx)
      integer :: kx
      real(8) :: ax(kx), c, fmax
      real(8) :: x, fx 
      fmax=0.0
      if(kx.eq.0) go to 120
      do 100 i=1,1000
      x=(1000-i)*6.0d0/c/1000
      fx=ax(1)
      if(kx.eq.1) go to 130
      do 110 k=2,kx
cx  110 fx=fx+ax(k)*x**(k-1)
      fx=fx+ax(k)*x**(k-1)
  110 continue
  130 continue
      fx=fx*exp(-c*x)
      if(fmax.lt.fx) fmax=fx
  100 continue
  120 continue
      fmax=fmax*1.5
      return
      end
c
c
cc      subroutine unifor(r)
      subroutine unifor(r,ir)
c     generation of pseudo-random numbers
cc      data ir/584287/
      integer :: ir
      real(4) :: r
      ir=ir*48828125
cx      if(ir) 10,20,20
      if(ir.lt.0) go to 10
      if(ir.ge.0) go to 20
   10 ir=(ir+2147483647)+1
   20 r=float(ir)*0.4656613e-9
      return
      end
c
c
      subroutine ptrend(x,ptx,axz,kxz)
c
c     trend polynomial; it is possible to includes here another
c     deterministic functions for an intensity of non-stationary
c     poisson process such as a fourier series for cyclic component.
c
cx      implicit real*8(a-h,o-z)
cx      dimension axz(1)
cx      dimension axz(kxz)
      integer :: kxz
      real(8) :: x, ptx, axz(kxz)
      ptx=0.0
      do 10 i=1,kxz
      ptx=ptx+axz(i)*x**(i-1)
   10 continue
      return
      end
