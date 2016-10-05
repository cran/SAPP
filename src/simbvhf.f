cc      program simbvh
      subroutine simbvhf(kxx,kxy,kxz,kyx,kyy,kyz,t,c,d,c2,d2,axx,axy,
cx     & axz,ayx,ayy,ayz,ptxmax,ptymax,kmax,xx,yy,i1,j1,err,nnmax,mmmax)
     & axz,ayx,ayy,ayz,ptxmax,ptymax,kmax,xx,yy,i1,j1,err,nnmax,mmmax,
     & ier)
c
      include 'sapp_f.h'
c
c     this program performs the simulation of bi-variate hawkes'
c  mutually exciting point processes.  the response functions are
c  parametrized by the laguerre type polynomials.
c
c     structure of the program
c
c          simbvh
c             i---input
c             i---initl
c             i---simda----comfac
c             i     i------duf
c             i     i------unifor
c             i     +------fx-------ptrend
c             +---output
c
c     this program is designed by y. ogata, and programed by y. ogata
c  and k. katsura, inst. statist. math., tokyo. (31/01/85)
c
c     references
c
c  y. ogata (1981). "on lewis' simulation method for point processes."
c     ieee information theory, vol. it-27, pp. 23-31.
c  a. g. hawkes (1971). "spectra of some self-exciting andmutually
c     exciting point processes."  biometrika, vol. 58, pp. 83-90.
c
cx      implicit real*8(a-h,o-z)
cc      dimension xx(10000),yy(10000)
cc      dimension axx(100),axy(100),axz(100),ayx(100),ayy(100),ayz(100)
cc      dimension ei(100),ej(100),fi(100),fj(100),lf(51,51)
cc      dimension ei2(100),ej2(100),fi2(100),fj2(100)
cx      dimension xx(nnmax),yy(mmmax)
cx      dimension axx(kxx),axy(kxy),axz(kxz),ayx(kyx),ayy(kyy),ayz(kyz)
cx      dimension ei(kmax),ej(kmax),fi(kmax),fj(kmax),lf(kmax,kmax)
cx      dimension ei2(kmax),ej2(kmax),fi2(kmax),fj2(kmax)
      integer :: kxx, kxy, kxz, kyx, kyy, kyz, kmax,i1, j1, nnmax,
     1           mmmax, ier
      real(8) :: t, c, d, c2, d2, axx(kxx), axy(kxy), axz(kxz),
     1           ayx(kyx), ayy(kyy), ayz(kyz), ptxmax, ptymax,
     2           xx(nnmax), yy(mmmax), err
      integer :: lf(kmax,kmax)
      real(8) :: ei(kmax), ej(kmax), fi(kmax), fj(kmax), ei2(kmax),
     1           ej2(kmax), fi2(kmax), fj2(kmax), fxxmax, fxymax,
     2           fyxmax, fyymax
c
c
cc      call input(kxx,kxy,kxz,kyx,kyy,kyz,t,c,d,c2,d2,axx,axy,axz,ayx,ayy
cc     &           ,ayz,ptxmax,ptymax,kmax)
c
      call initl(kxx,axx,c,fxxmax)
      call initl(kxy,axy,d,fxymax)
      call initl(kyx,ayx,c2,fyxmax)
      call initl(kyy,ayy,d2,fyymax)
c
      call simda(kxx,kxy,kxz,kyx,kyy,kyz,t,c,d,c2,d2,axx,axy,axz,ayx,ayy
     &           ,ayz,fxxmax,fxymax,fyxmax,fyymax,xx,yy,ei,ej,fi,fj,
cc     &           ei2,ej2,fi2,fj2,kmax,ptxmax,ptymax,lf,i1,j1)
cx     &           ei2,ej2,fi2,fj2,kmax,ptxmax,ptymax,lf,i1,j1,err)
     & ei2,ej2,fi2,fj2,kmax,ptxmax,ptymax,nnmax,mmmax,lf,i1,j1,err,ier)
c
cc      call output(xx,yy,i1,j1,t)
c
      return
      end
      subroutine simda(kxx,kxy,kxz,kyx,kyy,kyz,t,c,d,c2,d2,axx,axy,axz,
     &                 ayx,ayy,ayz,fxxmax,fxymax,fyxmax,fyymax,xx,yy,
     &                 ei,ej,fi,fj,ei2,ej2,fi2,fj2,kmax,ptxmax,
cc     &                 ptymax,lf,i1,j1)
cx     &                 ptymax,lf,i1,j1,err)
     &                 ptymax,nmax,mmax,lf,i1,j1,err,ier)
cx      implicit real*8(a-h,o-z)
cx      dimension axx(1),axy(1),axz(1),xx(1),yy(1),ei(1),ej(1),fi(1),fj(1)
cx      dimension ayx(1),ayy(1),ayz(1),ei2(1),ej2(1),fi2(1),fj2(1)
cx      dimension axx(kxx),axy(kxy),axz(kxz),xx(nmax),yy(mmax)
cx      dimension ei(kmax),ej(kmax),fi(kmax),fj(kmax)
cx      dimension ayx(kyx),ayy(kyy),ayz(kyz)
cx      dimension ei2(kmax),ej2(kmax),fi2(kmax),fj2(kmax)
cx      dimension lf(kmax,kmax)
cx      real*4r
      integer :: kxx, kxy, kxz, kyx, kyy, kyz, kmax, nmax, mmax,
     1           lf(kmax,kmax), i1, j1, ier
      real(8) :: t, c, d, c2, d2, axx(kxx), axy(kxy), axz(kxz), 
     1           ayx(kyx), ayy(kyy), ayz(kyz), fxxmax, fxymax, fyxmax,
     2           fyymax, xx(nmax), yy(mmax), ei(kmax), ej(kmax), 
     3           fi(kmax), fj(kmax), ei2(kmax), ej2(kmax), fi2(kmax),
     4           fj2(kmax), ptxmax, ptymax, err
      real(4) :: r
      real(8) :: uity, duity, e, x, xity, yity, probx, prob
c---
      ier=0
      err=0.0
      ir=584287
c---
      do 100 k=1,kmax
      ei(k)=0.0
      ej(k)=0.0
      fi(k)=0.0
      fj(k)=0.0
      ei2(k)=0.0
      ej2(k)=0.0
      fi2(k)=0.0
      fj2(k)=0.0
  100 continue
c---------------------
cc     call comfac(lf)
      call comfac(kmax,lf)
c---------------------
      x=0.0
      uity=ptxmax+ptymax+fxxmax+fxymax+fyxmax+fyymax
      i=0
      j=0
   30 continue
c-----------------------------------------------------------------------
cc      call duf(i,j,x,duity,xx,yy,axx,axy,ayx,ayy,kxx,kxy,kyx,kyy,c,d,
      call dufs(i,j,x,duity,xx,yy,axx,axy,ayx,ayy,kxx,kxy,kyx,kyy,c,d,
     &         c2,d2,ei,ej,fi,fj,ei2,ej2,fi2,fj2,ptxmax,ptymax)
c-----------------------------------------------------------------------
      if(duity.lt.uity) uity=duity
   20 continue
c-----------------------
cc      call unifor(r)
      call unifor(r,ir)
c-----------------------
   40 continue
      e=-alog(r)/uity
      x=x+e
      if(x.gt.t)go to 80
c-----------------------------------------------------------------------
cc      call fx(i,j,x,xity,yity,axx,axy,axz,ayx,ayy,ayz,kxx,kxy,kxz,kyx,
cc     &        kyy,kyz,c,d,c2,d2,lf,ei,ej,fi,fj,ei2,ej2,fi2,fj2,xx,yy)
      call fxs(i,j,x,xity,yity,axx,axy,axz,ayx,ayy,ayz,kxx,kxy,kxz,
     & kyx,kyy,kyz,kmax,c,d,c2,d2,lf,ei,ej,fi,fj,ei2,ej2,fi2,fj2,xx,yy)
c-----------------------------------------------------------------------
      probx=xity/uity
      prob=(xity+yity)/uity
      if(prob.le.1.0d00)go to 50
cc      write(6,1) prob
      err=prob
      go to 40
   50 continue
cx   60 continue
c----------------------
cc      call unifor(r)
      call unifor(r,ir)
c----------------------
      if(probx.le.r)go to 70
      i=i+1
c--------------------
      if(i.gt.nmax) then
         ier=-1
         return
      end if
c--------------------
      xx(i)=x
      do 10 k=1,kmax
      ei(k)=fi(k)
      ei2(k)=fi2(k)
   10 continue
      uity=uity+fxxmax+fyxmax
      go to 20
   70 if(prob.le.r) go to 30
      j=j+1
c--------------------
      if(j.gt.mmax) then
         ier=-2
         return
      end if
c--------------------
      yy(j)=x
      do 90 k=1,kmax
      ej(k)=fj(k)
      ej2(k)=fj2(k)
   90 continue
      uity=uity+fxymax+fyymax
      go to 20
   80 continue
      i1=i
      j1=j
cx    1 format(1h ,'warning: are ptxmax & ptymax correct? prob=',f10.5)
      return
      end
cc      subroutine duf(i,j,x,duity,xx,yy,axx,axy,ayx,ayy,kxx,kxy,kyx,kyy,
      subroutine dufs(i,j,x,duity,xx,yy,axx,axy,ayx,ayy,kxx,kxy,kyx,kyy,
     &         c,d,c2,d2,ei,ej,fi,fj,ei2,ej2,fi2,fj2,ptxmax,ptymax)
c
c     decreasing process except jumps, which is always greater than
c     the intensity process in subroutine fx.
c
cx      implicit real * 8 (a-h,o-z)
cx      dimension axx(1),axy(1),ei(1),ej(1),fi(1),fj(1)
cx      dimension ayx(1),ayy(1),ei2(1),ej2(1),fi2(1),fj2(1)
cx      dimension xx(1),yy(1)
cx      dimension axx(kxx),axy(kxy),ei(1),ej(1),fi(1),fj(1)
cx      dimension ayx(kyx),ayy(kyy),ei2(1),ej2(1),fi2(1),fj2(1)
cx      dimension xx(i),yy(j)
      integer :: i, j, kxx, kxy, kyx, kyy
      real(8) :: x, duity, xx(i), yy(j), axx(kxx), axy(kxy), ayx(kyx),
     1           ayy(kyy), c, d, c2, d2, ei(1), ej(1), fi(1), fj(1),
     2           ei2(1), ej2(1), fi2(1), fj2(1), ptxmax, ptymax
      real(8) :: cxp, cyp, cxp2, cyp2, bxx, cx, bxy, cy, byx, byy, dxxi,
     1           ecdxxi, ecdyxi, dyyj, ecdxyj, ecdyyj, xity, yity
      ixf=1
      iyf=1
      ixf2=1
      iyf2=1
      cxp=0.0
      cyp=0.0
      cxp2=0.0
      cyp2=0.0
      if(kxx.eq.0) go to 35
      do 30 ix=1,kxx
      bxx=max(axx(ix),0.d0)
      cx=bxx/(c/2)**(ix-1) *ixf
      cxp=max(cx,cxp)
      ixf=ixf*ix
   30 continue
   35 continue
      if(kxy.eq.0) go to 55
      do 50 iy=1,kxy
      bxy=max(axy(iy),0.d0)
      cy=bxy/(d/2)**(iy-1) *iyf
      cyp=max(cy,cyp)
      iyf=iyf*iy
   50 continue
   55 continue
      if(kyx.eq.0) go to 45
      do 40 ix=1,kyx
      byx=max(ayx(ix),0.d0)
      cx=byx/(c2/2)**(ix-1) *ixf2
      cxp2=max(cx,cxp2)
      ixf2=ixf2*ix
   40 continue
   45 continue
      if(kyy.eq.0) go to 65
      do 60 iy=1,kyy
      byy=max(ayy(iy),0.d0)
      cy=byy/(d2/2)**(iy-1) *iyf2
      cyp2=max(cy,cyp2)
      iyf2=iyf2*iy
   60 continue
   65 continue
      if(i.eq.0) go to 10
      dxxi=x-xx(i)
      ecdxxi=0.0d00
      if(-c/2*dxxi.ge.-20.0) ecdxxi=exp(-c/2*dxxi)
      fi(1)=ecdxxi*(ei(1)+1.d0)
      ecdyxi=0.0d00
      if(-c2/2*dxxi.ge.-20.0) ecdyxi=exp(-c2/2*dxxi)
      fi2(1)=ecdyxi*(ei2(1)+1.d0)
   10 continue
      if(j.eq.0) go to 20
      dyyj=x-yy(j)
      ecdxyj=0.0d00
      if(-d/2*dyyj.ge.-20.0) ecdxyj=exp(-d/2*dyyj)
      fj(1)=ecdxyj*(ej(1)+1.d0)
      ecdyyj=0.0d00
      if(-d2/2*dyyj.ge.-20.0) ecdyyj=exp(-d2/2*dyyj)
      fj2(1)=ecdyyj*(ej2(1)+1.d0)
   20 continue
      xity=ptxmax+cxp*fi(1)+cyp*fj(1)
      yity=ptymax+cxp2*fi2(1)+cyp2*fj2(1)
      duity=xity+yity
      return
      end
cc      subroutine fx(i,j,x,xity,yity,axx,axy,axz,ayx,ayy,ayz,kxx,kxy,kxz,
cc     &       kyx,kyy,kyz,c,d,c2,d2,lf,ei,ej,fi,fj,ei2,ej2,fi2,fj2,xx,yy)
      subroutine fxs(i,j,x,xity,yity,axx,axy,axz,ayx,ayy,ayz,kxx,kxy,
     &   kxz,kyx,kyy,kyz,kmax,c,d,c2,d2,lf,ei,ej,fi,fj,ei2,ej2,fi2,fj2,
     &   xx,yy)
c
c     intensity processes
c
ccx      implicit real * 8 (a-h,o-z)
cc      dimension axx(1),axy(1),axz(1),ei(1),ej(1),fi(1),fj(1),lf(51,51)
cx      dimension axx(1),axy(1),axz(1),ei(1),ej(1),fi(1),fj(1)
cx      dimension lf(kmax,kmax)
cx      dimension ayx(1),ayy(1),ayz(1),ei2(1),ej2(1),fi2(1),fj2(1)
cx      dimension xx(1),yy(1)
ccx      dimension axx(kxx),axy(kxy),axz(kxz)
ccx      dimension ei(kmax),ej(kmax),fi(kmax),fj(kmax)
ccx      dimension ayx(kyx),ayy(kyy),ayz(kyz)
ccx      dimension ei2(kmax),ej2(kmax),fi2(kmax),fj2(kmax)
ccx      dimension xx(i),yy(j)
      integer :: i, j, kxx, kxy, kxz, kyx, kyy, kyz, kmax, lf(kmax,kmax)
      real(8) :: x, xity, yity, axx(kxx), axy(kxy), axz(kxz),
     1           ayx(kyx), ayy(kyy), ayz(kyz), c, d, c2, d2, ei(kmax),
     2           ej(kmax), fi(kmax), fj(kmax), ei2(kmax), ej2(kmax),
     3           fi2(kmax), fj2(kmax), xx(i), yy(j)
      REAL(8) :: dxxi, ecdxxi, ff, ecdyxi, dyyj, eddxyj, eddyyj, ptx,
     1           pty
c
cc      dxxi=x-xx(i)
      dxxi=x
      if(i.ne.0) dxxi=x-xx(i)
      ecdxxi=0.0d00
      if(kxx.eq.0) go to 75
      if(-c*dxxi.ge.-20.0) ecdxxi=exp(-c*dxxi)
      do 70 ix=1,kxx
      ff=ei(1)*lf(ix,1)+1.d0
      if(ix.eq.1) go to 85
      do 80 jx=2,ix
      ff=ff*dxxi+ei(jx)*lf(ix,jx)
   80 continue
   85 continue
      fi(ix)=ecdxxi*ff
   70 continue
   75 continue
      ecdyxi=0.0d00
      if(kyx.eq.0) go to 175
      if(-c2*dxxi.ge.-20.0) ecdyxi=exp(-c2*dxxi)
      do 170 ix=1,kyx
      ff=ei2(1)*lf(ix,1)+1.d0
      if(ix.eq.1) go to 185
      do 180 jx=2,ix
      ff=ff*dxxi+ei2(jx)*lf(ix,jx)
  180 continue
  185 continue
      fi2(ix)=ecdyxi*ff
  170 continue
  175 continue
c
cc      dyyj=x-yy(j)
      dyyj=x
      if(j.ne.0) dyyj=x-yy(j)
      eddxyj=0.0d00
      if(kxy.eq.0) go to 95
      if(-d*dyyj.ge.-20.0) eddxyj=exp(-d*dyyj)
      do 90 ix=1,kxy
      ff=ej(1)*lf(ix,1)+1.d0
      if(ix.eq.1) go to 105
      do 100 jx=2,ix
      ff=ff*dyyj+ej(jx)*lf(ix,jx)
  100 continue
  105 continue
      fj(ix)=eddxyj*ff
   90 continue
   95 continue
      eddyyj=0.0d00
      if(kyy.eq.0) go to 195
      if(-d2*dyyj.ge.-20.0) eddyyj=exp(-d2*dyyj)
      do 190 ix=1,kyy
      ff=ej2(1)*lf(ix,1)+1.d0
      if(ix.eq.1) go to 205
      do 200 jx=2,ix
      ff=ff*dyyj+ej2(jx)*lf(ix,jx)
  200 continue
  205 continue
      fj2(ix)=eddyyj*ff
  190 continue
  195 continue
c
      call ptrend(x,ptx,axz,kxz)
      xity=ptx
      if(kxx.eq.0) go to 15
      do 10 k=1,kxx
      xity=xity+axx(k)*fi(k)
   10 continue
   15 continue
c
      if(kxy.eq.0) go to 25
      do 20 k=1,kxy
      xity=xity+axy(k)*fj(k)
   20 continue
   25 continue
c
      call ptrend(x,pty,ayz,kyz)
      yity=pty
      if(kyx.eq.0) go to 115
      do 110 k=1,kyx
      yity=yity+ayx(k)*fi2(k)
  110 continue
  115 continue
c
      if(kyy.eq.0) go to 125
      do 120 k=1,kyy
      yity=yity+ayy(k)*fj2(k)
  120 continue
  125 continue
c
      return
      end
