c***********************************************************************
c
c                     function yjm0
c
c     this function returns the value of the spherical
c     harmonic y(j,m) when cos(theta)=x and phi=0.
c
c     input:
c       x      cos(theta)
c       j,m    order of the spherical harmonic
c
c     output:
c       yjm0   value of Y(j,m|cos(theta),0)
c
c***********************************************************************
      module P_lm
       implicit integer (i-n) 
        private
        public::yjm0
        contains
      function yjm0(j,m,x)
      implicit real*8 (a-h,o-z)
      pi = dacos(-1.d0)
c-----------------------------------------------------------------------
c
c      -- first calculate the associated legendre function
c
c-----------------------------------------------------------------------
      mabs = abs(m)
      pjm0 = plgndr(j,mabs,x)
c-----------------------------------------------------------------------
c
c     -- calculate the spherical harmonics
c
c-----------------------------------------------------------------------
c
c     -- positive m
c
      factor = (2.*j+1) * factrl(j-mabs) / (4.*pi*factrl(j+mabs) )
      sqfact = sqrt(factor)
      yjm0 = (-1.)**mabs * sqfact * pjm0
c
c     -- negative m
c
      if(m .lt. 0) then
         phase = (-1.)**m
         yjm0  = phase * yjm0
      end if
      yjm0=sqrt(2.d0*pi)*yjm0
      return
      end
c***********************************************************************
c
c                     function dyjm0
c
c     this function returns the value of the derivative of the
c     spherical harmonic y(j,m) wrt theta when cos(theta)=x
c     and phi=0.
c
c     input:
c       x      cos(theta)
c       j,m    order of the spherical harmonic
c
c     output:
c       dyjm0  value of dY(j,m|cos(theta),0)/dtheta
c
c***********************************************************************
      function dyjm0(j,m,x)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n) 
      pi = dacos(-1.d0)
      if(j .eq. 0) then
         dyjm0 = 0.
         return
      end if
c-----------------------------------------------------------------------
c
c      -- first calculate the associated legendre functions
c
c-----------------------------------------------------------------------
      if(m .eq. 0) then
         pjmp = plgndr(j,1,x)
         factor = -factrl(j-1) / factrl(j+1)
         pjmm = factor * pjmp
      else if(m .gt. 0) then
         pjmp = plgndr(j,m+1,x)
         pjmm = plgndr(j,m-1,x)
      else if(m .lt. 0) then
         mabs = iabs(m+1)
         pjmp = plgndr(j,mabs,x)
         factor = (-1.)**mabs * factrl(j-mabs) / factrl(j+mabs)
         pjmp = factor * pjmp
         mabs = iabs(m-1)
         pjmm = plgndr(j,mabs,x)
         factor = (-1.)**mabs * factrl(j-mabs) / factrl(j+mabs)
         pjmm = factor * pjmm
      end if
      dpjm0 = .5 * ( pjmp - (j+m)*(j-m+1)*pjmm )
c-----------------------------------------------------------------------
c
c     -- calculate the derivative of the spherical harmonic
c
c-----------------------------------------------------------------------
      factor = (2.*j+1) * factrl(j-m) / (4.*pi*factrl(j+m) )
      sqfact = sqrt(factor)
      dyjm0 = -(-1.)**m * sqfact * dpjm0
      return
      end
c***********************************************************************
c
c                       function plgndr
c
c     This function was taken from the book "Numerical Recipes -
c     the Art of Scientific Computing " by W.H. Press, B.P. Flannery,
c     S.A. Teukolsky and W.T. Vetterling, published by Cambridge Univ.
c     Press, New York, 1986. See p.182.
c
c     Computes the associated Legendre polynomial P_{l,m}(x). Here,
c     m and l are integers satisfying (0 .le. m .le. l), while x
c     lies in the range (-1 .le. x .le. 1).
c
c***********************************************************************
      function plgndr(l,m,x)
      implicit real*8(a-h,o-z)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) then
         write(6,*) ' *** fatal error in function plgndr '
         write(6,*) ' *** bad arguements (l,m,x): ', l, m, x
      end if
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      end
c***********************************************************************
c
c                       function factrl
c
c     This function was taken from the book "Numerical Recipes -
c     the Art of Scientific Computing " by W.H. Press, B.P. Flannery,
c     S.A. Teukolsky and W.T. Vetterling, published by Cambridge Univ.
c     Press, New York, 1986. See p.182.
c
c     Computes the value of N! as a floating point number.
c
c***********************************************************************
      function factrl(n)
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension a(33)
      data ntop,a(1)/0,1./
      if (n.lt.0) then
        write(6,*) ' *** fatal error in routine factrl '
        write(6,*) ' *** negative factorial n: ', n
        stop
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
	x=n+1.d0
        factrl=exp(gammln(x))
      endif
      return
      end
c***********************************************************************
c
c                       function gammln
c
c     This function was taken from the book "Numerical Recipes -
c     the Art of Scientific Computing " by W.H. Press, B.P. Flannery,
c     S.A. Teukolsky and W.T. Vetterling, published by Cambridge Univ.
c     Press, New York, 1986. See p.182.
c
c     Computes the value ln[gamma(xx)] for xx>0. Full accuracy is
c     obtained for xx>1. For 0<xx<1, the reflection formula (6.1.4)
c     can be used first.
c
c***********************************************************************
      function gammln(xx)
      implicit real*8(a-h,o-z)
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)
      return
      end
c===========================================================
c	implicit real*8 (a-h,o-z)
c	dimension fac(0:100),pl(100)

c	fac(0)=1.D0
c	do i=1,100
c	  fac(i)=factrl(i)
c	end do
c
c	read*,l,m,x
c	call yjm0M(l,m,x,pl,fac)
c
c	m0=iabs(m)
c	do i=1,l-m0+1
c	  j=m0+i-1
c	  write(*,*) pl(i),yjm0(j,m,x),j,m
c	end do
c
c	stop
c	end
c***********************************************************************


      end module P_lm
