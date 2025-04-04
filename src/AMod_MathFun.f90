module Math_Functions
   use amod_pot_KrRb, only: refpot=>KrStar_Rb
   integer,parameter::rk = kind(0.0D0)
   real(rk),parameter:: pi = acos(-1.d0)
   private
   public::rbesjy,Gamma,Solve_DVRm,DVR_0,PO_DVR
   contains
   
   subroutine PO_DVR(nn, rmin, rmax, mu, n, nvib, nrot, PotMax, r, psi, e1, U, T, UPO, POTREF)
      integer, intent(in):: nn, n, nvib, nrot
      real(rk), intent(in)::rmin, rmax, mu, PotMax
      real(rk), intent(out):: r(n), psi(n, n), U(n, n), T(n, n), UPO(nn, n), POTREF(n)
      real(rk):: length
      
      real(rk), allocatable, dimension(:):: r0, pot0, E0
      real(rk), allocatable, dimension(:, :):: WF0, U0, T0, H0, X0, wf01
      real(rk), allocatable, dimension(:, :) :: WF1, X1, H1
      real(rk), allocatable, dimension(:)::E1
      integer:: i, j,  ierr
      real(rk), allocatable:: fv1(:), fv2(:),  C_j_j0(:)
      length = (rmax - rmin)
 
      allocate (r0(nn), source=0.0_rk)
      allocate (pot0(nn), source=0.0_rk)
      do i = 1, nn
         r0(i) = x_a(i, length, nn, rmin)
         call refpot(r0(i),pot0(i))
      end do
      !pot0 = pot0+(0.0_rk - minval(pot0))
      DO I = 1, nn
         if (pot0(i) > PotMax) pot0(i) = PotMax
      END DO
    
      write (*, *) 'The range of initial potential:$', minval(pot0), ' arrow ', maxval(pot0), '.$'
      allocate (E0(nn), source=0.0_rk)
      allocate (WF0(nn, nn), source=0.0_rk)
      allocate (U0(nn, nn), source=0.0_rk)
      allocate (T0(nn, nn), source=0.0_rk)
      allocate (H0(nn, nn), source=0.0_rk)
      allocate (wf01(nn, nn), source=0.0_rk)
      allocate (C_j_j0(nn), source=0.0_rk)
      !call DVR_0(nn, nn, rmin, rmax, mu, pot0, nvib, nrot, r0, wf01, E0, U0, T0, H0)
      call DVR_0(nn, nn, rmin, rmax, mu, pot0, 0, 0, r0, WF0, E0, U0, T0)

      ! do i = 1, nn
      !    C_j_j0(i) = sum(wf0(:, i)*wf01(:, 1 + nvib))
      !    !print*,C_j_j0(i),'C_j_j0'
      ! end do
      ! print *, 'Nor', sum(C_j_j0**2)

      deallocate (H0, pot0)
      do j = 1, n
         do i = 1, nn
            upo(i, j) = WF0(i, j)
         end do
      end do
      allocate (X0(nn, nn), source=0.0_rk)
      do i = 1, nn
         X0(i, i) = r0(i)
      end do
      allocate (X1(n, n), source=0.0_rk)
      write (*, *) 'The range of initial r:$', minval(r0), ' arrow ', maxval(r0), '.$'
      X1 = 0.0_rk
      X1 = matmul(transpose(UPO), matmul(X0, UPO))
      deallocate (X0)
      write (*, *) 'U_po', rank(UPO), size(upo(:, 1)), size(UPO(1, :))
      allocate (fv1(n), source=0.0_rk)
      allocate (fv2(n), source=0.0_rk)
      call rs(n, n, X1, r, 1, U, fv1, fv2, ierr); if (ierr /= 0) stop 'Error at :PO_DVR:RS !'
      deallocate (fv1, fv2, X1)
      do j = 1, n
         if (U(1, j) < 0) then
            do i = 1, n
               U(i, j) = -U(i, j)
            end do
         end if
      end do
      allocate (H1(n, n), source=0.0_rk)
      T = 0.0_rk
      do i = 1, n
         T(i, i) = E0(i)!*C_j_j0(i)
      end do

      deallocate (E0)
      write (*, *) 'The range of r:$', minval(r), ' arrow ', maxval(r), '.$'
      H1 = matmul(transpose(U), matmul(T, U))

      POTREF = 0.0_rk
     
      do i = 1, n
         call refpot(r(i),POTREF(I))
      end do
      !POTREF = POTREF+(0.0_rk-minval(POTREF))
      do i = 1, n
         if (POTREF(i) > POtMax) POTREF(i) = PotMax
      end do

      write (*, *) 'The range of Reference potential:$', minval(POTREF), ' arrow ', maxval(POTREF), '.$'
   
      do i = 1, n
         H1(i, i) = H1(i, i) + nrot*(nrot + 1.0_rk)/2.0_rk/mu/r(i)**2
      end do

      allocate (fv1(n), source=0.0_rk)
      allocate (fv2(n), source=0.0_rk)
      allocate (WF1(n, n), source=0.0_rk)
      !allocate(E1(n),source=0.0_rk)
      call rs(n, n, H1, E1, 1, WF1, fv1, fv2, ierr)
      deallocate (fv1, fv2, H1)
      !deallocate(E1)
      psi = WF1!(:,1+nvib)
      deallocate (WF1)
     
      !allocate(wft(nn,nn))
      !wft = (UPO,matmul())
      write (*, *) 'Get initial WF'
      allocate (fv1(nn))
      fv1 = matmul(UPO, matmul(u, psi(:, 1 + nvib)))
      open (66, file='WF_PODVR0.dat')
      do i = 1, size(fv1)
         write (66, '(2(g0,x))') r0(i), abs(fv1(i))**2
      end do
      close (66)
      deallocate (fv1, r0)
      write (*, *) 'Initial WF OK'
   end subroutine

   
   subroutine DVR_0(nDVR, nfbr, x_min, x_max, mu, pot, v, l, xe, WFs, ee, Uja0, T)
     
      integer, intent(in)::nDVR, nfbr
      real(rk), intent(in)::x_min, x_max, pot(nDVR)
      real(rk), intent(in)::mu
      integer, intent(in)::v, l
      real(rk), intent(out)::T(nfbr, nfbr)
      real(rk), intent(out)::WFs(nDVR, nFBR), ee(nDVR), xe(nDVR), Uja0(nfbr, nDVR)

      integer, parameter::nn = 1000
      real(rk), allocatable, dimension(:, :):: TFBR
      real(rk), allocatable, dimension(:)::xx, fv1, fv2
      real(rk):: length, Hh(nDVR, nDVR)
      integer::i, j, ierr
      if (size(pot) /= nDVR) stop 'The size of pot neq grid (DVR:sin:DVR0)'
      length = (x_max - x_min)!*real(nDVR,rk)/(nDVR+1)
      allocate(xx(nDVR), source=0.0_rk)
      do i = 1, nDVR
         xx(I) = x_a(i, length, nDVR, x_min)
      end do
      print*,'Get Trans matrix'
      xe = xx
      !$omp parallel do private(j,i) schedule(dynamic) 
      do j = 1, nDVR
         do i = 1, nDVR
            Uja0(i, j) = U_ja(nDVR, j, i)
         end do
      end do
      !$omp end parallel do 

      print*,'Get Mon Matrix '
      T = 0.0_rk
      !$omp parallel do private(i) schedule(dynamic) 
      do i = 1, nFBR
         T(i, i) = real(i*pi, rk)**2/Length**2/2.0_rk/mu
      end do
      !$omp end parallel do 

      allocate(TFBR(nDVR, nDVR), source=0.0_rk)
      TFBR = matmul(matmul(transpose(Uja0), T), (Uja0))
      print*,'Get DVR Matrix' 
      Hh = 0
      !$OMP parallel do private(j,i) schedule(dynamic) 
      do j = 1, nDVR
         do i = 1, nDVR
            Hh(i, j) = TFBR(i, j) + pot(i)*delta(i, j) + l*(l + 1.0_rk)/mu/2.0_rk/xx(i)**2*delta(i, j)
         end do
      end do
      !$omp end parallel do 

      deallocate (TFBR)
      ee = 0
      ierr = 0
      print*,'Get DVR Energy'
      allocate (fv1(nDVR), source=0.0_rk)
      allocate (fv2(nDVR), source=0.0_rk)
      call rs(nFBR, nDVR, Hh, ee, 1, WFs, fv1, fv2, ierr)
      deallocate (fv1, fv2)

   end subroutine

   function x_a(a, L, n, xmin) result(x)
      integer, intent(in)::a, n
      real(rk), intent(in):: xmin, L
      real(rk)::x
      x = xmin + a*L/(n + 1)
   end function
   function delta(i, j) result(d)
      integer, intent(in):: i, j
      integer:: d
      d = 0
      if (i == j) d = 1
   end function
   function U_ja(n, j, a) result(U)
      !j the index of basis
      !a the index of grids
      integer, intent(in)::n, j, a
      real(rk)::U
      u = sqrt(2.0_rk/(n + 1))*sin(j*a*pi/real(n + 1, rk))
   end function U_ja
   function basisSin(L, j, x, x_0)
      real(rk), intent(in)::L, x, x_0
      integer, intent(in)::j
      real(rk)::basisSin
      if (x .ge. x_0 .and. x .le. (L + x_0)) then
         basisSin = sqrt(2/L)*sin(j*pi*(x - x_0)/l)
      else
         basisSin = 0
      end if
   end function basisSin


   subroutine Solve_DVR(nDVR,nfbr,x_min,x_max,mu,pot,v,l,xe,WFs,ee,Uja0,T)
   !use mod_dsyev,only:DSYEV

   integer,intent(in)::nDVR,nfbr
   real(rk),intent(in)::x_min,x_max
   real(rk),intent(in)::mu,pot(nDVR)
   integer,intent(in)::v,l
   real(rk),intent(out)::T(nfbr,nfbr)
   real(rk),intent(out)::WFs(nDVR,nFBR),ee(nDVR),xe(nDVR),Uja0(nfbr,nDVR)

   integer,parameter::nn=1000
   real(rk),dimension(nDVR,nDVR)::Q_jk,TDVR,VDVR,Hh
   real(rk),dimension(nDVR)::xx,fv1,fv2
   real(rk)::xing,length
   integer::i,j,k,ierr

   length = x_max-x_min
   !$omp parallel do private(i) schedule(dynamic)
   do i=1,nDVR
      xx(i)=x_min+i*length/(ndvr+1)
   end do
   !$OMP end parallel do 

   Q_jk=0
   do j=1,nFBR
      do i=1,nFBR
         do k=1,nn
            xing = x_a(k,length,nn,x_min)
            Q_jk(i,j)=Q_jk(i,j)+basis(length,i,xing,x_min)*xing*basis(length,j,xing,x_min)*(length/(nn+1))
         end do
     end do
   end do
   print*,'v=',v,'l=',l


   xe=0

   call rs(nFBR,nDVR,Q_jk,xe,1,Uja0,fv1,fv2,ierr)


   if (ierr /= 0 ) stop 'stop at ierr!'
   do j=1,nFBR
      if (Uja0(1,j)<0) then
         do i=1,nDVR
            Uja0(i,j)=-Uja0(i,j)
         end do
      end if
   end do
  
   T=0
   
   do i=1,nFBR
      T(i,i)= real(i*pi,rk)**2/Length**2/2.0_rk/mu+l*(l+1)/mu/2
   end do
   close(66)

  
   TDVR=matmul(matmul(transpose(Uja0),T),(Uja0))

      VDVR =0.0_rk
   do i=1,nDVR
     VDVR(i,i)=pot(i)
   end do

   Hh = 0
   do j=1,nDVR
      do i=1,nDVR
         Hh(i,j)=TDVR(i,j)+VDVR(i,j)
      end do
   end do
   ee=0
   ierr=0
   fv1=0;fv2=0

   call rs(nFBR,nDVR,Hh,ee,1,WFs,fv1,fv2,ierr)

   end subroutine Solve_DVR 



 function basis(L,j,x,x_0)
   real(rk),intent(in)::L,x,x_0
   integer,intent(in)::j
   real(rk)::basis
      if (x .ge. x_0 .and. x .le. (L+x_0) ) then
      basis = sqrt(2/L) *sin(j*pi*(x-x_0)/l)
      else
      basis =0
      endif
   end function basis

   function Gamma(R)
   real(rk)::Gamma
   real(rk),intent(in)::R 
   Gamma =1.0e-10_rk * exp(-R/20.0_rk)
   end function Gamma


subroutine rbesjy (ell,x,z,zp,zu,zup)
      implicit double precision (a-h,o-z)
        !ell: l angular momenton
        !x: sqrt(2 mass E) * x
        !z real
        !zp dz/dx
        !zu imag
        !d zu /dx
!c          -----------------------------------------------------------------
!c    Riccati-Bessel functions of fractional order
!c          and their first derivatives with respect to x:
!c
!c          j(ell,x) = cj * exp(ej)
!c          y(ell,x) = cy * exp(ey)
!c          d/dx j(ell,x) = dj * exp(ej)
!c          d/dx y(ell,x) = dy * exp(ey)
!c          -----------------------------------------------------------------
!c
      if (x.le.0.0d0 .or. ell.lt.-0.5d0) stop 'rbesjy 0'
      v = ell+0.5d0
      call bessjy (v,x,cj,dj,ej,cy,dy,ey)
      ! pi = acos(-1.d0)
      ex = 0.5d0*dlog(pi*x/2.d0)
      dj = dj+cj/(2.d0*x)
      sj = sqrt(cj*cj+dj*dj)
      cj = cj/sj
      dj = dj/sj
      ej = ej+dlog(sj)+ex
      dy = dy+cy/(2.d0*x)
      sy = sqrt(cy*cy+dy*dy)
      cy = cy/sy
      dy = dy/sy
      ey = ey+dlog(sy)+ex

!c add by ZDH
      z = cj * exp(ej)
      zu = -cy * exp(ey)
      zp = dj * exp(ej)
      zup = dy * exp(ey)

      return
      end

      subroutine bessjy (v,x,cj,dj,ej,cy,dy,ey)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
!c
!c          -----------------------------------------------------------------
!c          This subroutine uses a combination of methods (mostly due
!c          to Temme) to calculate the Ordinary Bessel functions

!c          J(v,x) = cj * exp(ej)
!c          Y(v,x) = cy * exp(ey)

!c          and their first derivatives with respect to x

!c          d/dx J(v,x) = dj * exp(ej)
!c          d/dx Y(v,x) = dy * exp(ey)

!c          for a given real order v >= 0 and real argument x > 0.
!c          Note the exponential scaling, which is used to avoid
!c          overflow of Y(v,x) and underflow of J(v,x) for v >> x.
!c          -----------------------------------------------------------------

      double precision, parameter :: eps = 1.d-15! consistent with rgamma
      integer,parameter :: maxit = 1000

      if (v.lt.0.d0 .or. x.le.0.d0) stop 'bessjy 0'
      ! pi = acos(-1.d0)
      xmin = 3.d0
      xmax = 5.d0-dlog10(eps)

!c          begin by calculating Y(a,x) and Y(a+1,x) for |a| <= 1/2

      na = int(v+0.5d0)
      a = v-na
      if (x .lt. xmin) then

!c             using Temme's series (bessya) for small x
!c             [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]

         b = x/2.d0
         d = -dlog(b)
         e = a*d
         if (abs(a) .lt. eps) then
            c = 1.d0/pi
         else
            c = a/sin(a*pi)
         endif
         if (abs(e) .lt. eps) then
            s = 1.d0
         else
            s = sinh(e)/e
         endif
         e = exp(e)
         g = e*rgamma(a,p,q)
         e = (e+1.d0/e)/2.d0
         f = 2*c*(p*e+q*s*d)
         e = a*a
         p = g*c
         q = 1.d0/g/pi
         c = a*pi/2.d0
         if (abs(c) .lt. eps) then
            r = 1.d0
         else
            r = sin(c)/c
         endif
         r = pi*c*r*r
         c = 1.d0
         d = -b*b
         ya = f+r*q
         ya1 = p
         do n = 1,maxit
            f = (f*n+p+q)/(n*n-e)
            c = c*d/n
            p = p/(n-a)
            q = q/(n+a)
            g = c*(f+r*q)
            h = c*p-n*g
            ya = ya+g
            ya1 = ya1+h
            del = abs(g)/(1.d0+abs(ya))
            del1 = abs(h)/(1.d0+abs(ya1))
            if (del+del1 .lt. eps) go to 1
         enddo
         stop 'bessjy 1'
   1     f = -ya
         g = -ya1/b
      else if (x.ge.xmin .and. x.lt.xmax) then

!c             Temme's PQ method (besspqa) for intermediate x
!c             [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]

         c = 0.25d0-a*a
         b = x+x
         p = pi
         e = (x*cos(a*pi)/pi/eps)**2
         p = 1.d0
         q = -x
         r = 1.d0+x*x
         s = r
         do n = 2,maxit
            d = (n-1+c/n)/s
            p = (2*n-p*d)/(n+1)
            q = (-b+q*d)/(n+1)
            s = p*p+q*q
            r = r*s
            if (r*n*n .gt. e) go to 2
         enddo
         stop 'bessjy 2'
   2     p = p/s
         f = p
         q = -q/s
         g = q
         do m = n,1,-1
            r = (m+1)*(2.d0-p)-2.d0
            s = b+(m+1)*q
            d = (m-1+c/m)/(r*r+s*s)
            p = d*r
            q = d*s
            e = f+1.d0
            f = p*e-g*q
            g = q*e+p*g
         enddo
         f = 1.d0+f
         d = f*f+g*g
         pa = f/d
         qa = -g/d
         d = a+0.5d0-p
         q = q+x
         pa1 = (pa*q-qa*d)/x
         qa1 = (qa*q+pa*d)/x
         b = x-pi*(a+0.5d0)/2.d0
         c = cos(b)
         s = sin(b)
         d = sqrt(2.d0/x/pi)
         f = d*(pa*s+qa*c)
         g = d*(qa1*s-pa1*c)
      else if (x .ge. xmax) then

!c             and Hankel's asymptotic expansions for large x
!c             [ Abramowitz and Stegun, Section 9.2 ]

         p = 0.d0
         q = 0.d0
         do ia = 0,1
            pa = p
            qa = q
            y = 4.d0*(a+ia)**2
            z = 8.d0*x
            d = 0.d0
            w = -1.d0
            p = 1.d0
            q = 0.d0
            tp = 1.d0
            do k = 1,maxit
               d = d+z
               w = w+2.d0
               tq = +tp*(y-w*w)/d
               q = q+tq
               d = d+z
               w = w+2.d0
               tp = -tq*(y-w*w)/d
               p = p+tp
               if (abs(tp)+abs(tq) .lt. eps) go to 3
            enddo
            stop 'bessjy 3'
   3        p = p-0.5d0*tp
            q = q-0.5d0*tq
         enddo
         pa1 = p
         qa1 = q
         b = x-pi*(a+0.5d0)/2.d0
         c = cos(b)
         s = sin(b)
         d = sqrt(2.d0/x/pi)
         f = d*(pa*s+qa*c)
         g = d*(qa1*s-pa1*c)
      endif

!c          now recur upwards from Y(a,x) to Y(v,x),
!c          scaling to avoid overflow along the way

      p = 0.d0
      if (na .gt. 0) then
         y = 2.d0/x
         do n = 1,na
            h = y*(a+n)*g-f
            f = g
            g = h
   4        if (abs(f) .gt. 4.d0) then
               p = p+1.d0
               f = 0.0625d0*f
               g = 0.0625d0*g
               go to 4
            endif
         enddo
      endif
      cy = f
      dy = (v/x)*f-g
      sy = sqrt(cy*cy+dy*dy)
      cy = cy/sy
      dy = dy/sy
      ey = dlog(sy)+p*dlog(16.d0)

!c          finally, calculate J(v,x) and dJ(v,x)/dx

      vv = max(xmin,v)
      if (x .ge. vv) then

!c             using upward recursion in the classically allowed region

         f = d*(pa*c-qa*s)
         g = d*(qa1*c+pa1*s)
         if (na .gt. 0) then
            y = 2.d0/x
            do n = 1,na
               h = y*(a+n)*g-f
               f = g
               g = h
            enddo
         endif
         cj = f
         dj = (v/x)*f-g
         sj = sqrt(cj*cj+dj*dj)
         cj = cj/sj
         dj = dj/sj
         ej = dlog(sj)
      else

!c             and CF1 in the classically forbidden region
!c             [ Numerical Recipes, 2nd Edition, Section 6.7 ]

         ap = 1.d0
         a = v/x
         bp = 0.d0
         b = 1.d0
         f = 0.d0
         g = 0.d0
         y = 2.d0/x
         w = y/pi
         do n = 1,maxit
            an = y*(v+n)*a-ap
            ap = a
            a = an
            bn = y*(v+n)*b-bp
            bp = b
            b = bn
            if (abs(b) .gt. abs(a)) then
               ap = ap/b
               a = a/b
               bp = bp/b
               b = 1.d0
               if (abs(a-f) .lt. eps*abs(f)) then
                  cj = w/(dy-cy*a)
                  dj = a*cj
                  go to 5
               endif
               f = a
            else
               bp = bp/a
               b = b/a
               ap = ap/a
               a = 1.d0
               if (abs(b-g) .lt. eps*abs(g)) then
                  dj = w/(dy*b-cy)
                  cj = b*dj
                  go to 5
               endif
               g = b
            endif
         enddo
         stop 'bessjy 4'
   5     sj = sqrt(cj*cj+dj*dj)
         cj = cj/sj
         dj = dj/sj
         ej = dlog(sj)-ey
      endif
      return
      end

      function rgamma(x,odd,even)
      implicit double precision (a-h,o-z)
      implicit integer  (i-n)

!c          -----------------------------------------------------------------
!c          Direct fortran translation of Temme's algol routine for computing
!c          rgamma = 1/Gamma(1-x), along with its odd and even parts, for
!c          abs(x) .le. 0.5. [ N.M.Temme, J Comput Phys 19 (1975) 324-337 ]
!c          -----------------------------------------------------------------

      dimension b(12)
      data b / -0.283876542276024d0, -0.076852840844786d0, &
              +0.001706305071096d0, +0.001271927136655d0,&
             +0.000076309597586d0, -0.000004971736704d0,&
              -0.000000865920800d0, -0.000000033126120d0,&
              +0.000000001745136d0, +0.000000000242310d0,&
              +0.000000000009161d0, -0.000000000000170d0 /
      save b

      x2 = x*x*8.d0
      alfa = -0.000000000000001d0
      beta = 0.d0
      do i = 12,2,-2
         beta = -(2*alfa+beta)
         alfa = -beta*x2-alfa+b(i)
      enddo
      even = (beta/2.d0+alfa)*x2-alfa+0.921870293650453d0
      alfa = -0.000000000000034d0
      beta = 0.d0
      do i = 11,1,-2
         beta = -(2*alfa+beta)
         alfa = -beta*x2-alfa+b(i)
      enddo
      odd = 2*(alfa+beta)
      rgamma = odd*x+even
      return
      end


subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
  integer n,nm,ierr,matz
  real(rk) a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
  !c
  !c     this subroutine calls the recommended sequence of
  !c     subroutines from the eigensystem subroutine package (eispack)
  !c     to find the eigenvalues and eigenvectors (if desired)
  !c     of a real symmetric matrix.
  !c
  !c     on input
  !c
  !c        nm  must be set to the row dimension of the two-dimensional
  !c        array parameters as declared in the calling program
  !c        dimension statement.
  !c
  !c        n  is the order of the matrix  a.
  !c
  !c        a  contains the real symmetric matrix.
  !c
  !c        matz  is an integer variable set equal to zero if
  !c        only eigenvalues are desired.  otherwise it is set to
  !c        any non-zero integer for both eigenvalues and eigenvectors.
  !c
  !c     on output
  !c
  !c        w  contains the eigenvalues in ascending order.
  !c
  !c        z  contains the eigenvectors if matz is not zero.
  !c
  !c        ierr  is an integer output variable set equal to an error
  !c           completion code described in the documentation for tqlrat
  !c           and tql2.  the normal completion code is zero.
  !c
  !c        fv1  and  fv2  are temporary storage arrays.
  !c
  !c     questions and comments should be directed to burton s. garbow,
  !c     mathematics and computer science div, argonne national laboratory
  !c
  !c     this version dated august 1983.
  !c
  !c     ------------------------------------------------------------------
  !c
  if (n .le. nm) go to 10
  ierr = 10 * n
  go to 50
  !c
  10 if (matz .ne. 0) go to 20
  !c     .......... find eigenvalues only ..........
  call  tred1(nm,n,a,w,fv1,fv2)
  !*  tqlrat encounters catastrophic underflow on the Vax
  !*     call  tqlrat(n,w,fv2,ierr)
  call  tql1(n,w,fv1,ierr)
  go to 50
  !c     .......... find both eigenvalues and eigenvectors ..........
  20 call  tred2(nm,n,a,w,fv1,z)
  call  tql2(nm,n,w,fv1,z,ierr)

  50 return
  end
  !!!***********************************
  subroutine tred1(nm,n,a,d,e,e2)

  integer i,j,k,l,n,ii,nm,jp1
  real(rk) a(nm,n),d(n),e(n),e2(n)
  real(rk) f,g,h,scale
  !c
  !c     this subroutine is a translation of the algol procedure tred1,
  !c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
  !c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
  !c
  !c     this subroutine reduces a real symmetric matrix
  !c     to a symmetric tridiagonal matrix using
  !c     orthogonal similarity transformations.
  !c
  !c     on input
  !c
  !c        nm must be set to the row dimension of two-dimensional
  !c          array parameters as declared in the calling program
  !c          dimension statement.
  !c
  !c        n is the order of the matrix.
  !c
  !c        a contains the real symmetric input matrix.  only the
  !c          lower triangle of the matrix need be supplied.
  !c
  !c     on output
  !c
  !c        a contains information about the orthogonal trans-
  !c          formations used in the reduction in its strict lower
  !c          triangle.  the full upper triangle of a is unaltered.
  !c
  !c        d contains the diagonal elements of the tridiagonal matrix.
  !c
  !c        e contains the subdiagonal elements of the tridiagonal
  !c          matrix in its last n-1 positions.  e(1) is set to zero.
  !c
  !c        e2 contains the squares of the corresponding elements of e.
  !c          e2 may coincide with e if the squares are not needed.
  !c
  !c     questions and comments should be directed to burton s. garbow,
  !c     mathematics and computer science div, argonne national laboratory
  !c
  !c     this version dated august 1983.
  !c
  !c     ------------------------------------------------------------------
  !c
  do 100 i = 1, n
     d(i) = a(n,i)
     a(n,i) = a(i,i)
  100 continue
  !c     .......... for i=n step -1 until 1 do -- ..........
  do 300 ii = 1, n
     i = n + 1 - ii
     l = i - 1
     h = 0.0d0
     scale = 0.0d0
     if (l .lt. 1) go to 130
  !c     .......... scale row (algol tol then not needed) ..........
     do 120 k = 1, l
  120    scale = scale + dabs(d(k))
  !c
     if (scale .ne. 0.0d0) go to 140
  !c
     do 125 j = 1, l
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = 0.0d0
  125    continue
  !c
  130    e(i) = 0.0d0
     e2(i) = 0.0d0
     go to 300
  !c
  140    do 150 k = 1, l
        d(k) = d(k) / scale
        h = h + d(k) * d(k)
  150    continue
  !c
     e2(i) = scale * scale * h
     f = d(l)
     g = -dsign(dsqrt(h),f)
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g
     if (l .eq. 1) go to 285
  !c     .......... form a*u ..........
     do 170 j = 1, l
  170    e(j) = 0.0d0
  !c
     do 240 j = 1, l
        f = d(j)
        g = e(j) + a(j,j) * f
        jp1 = j + 1
        if (l .lt. jp1) go to 220
  !c
        do 200 k = jp1, l
           g = g + a(k,j) * d(k)
           e(k) = e(k) + a(k,j) * f
  200       continue
  !c
  220       e(j) = g
  240    continue
  !c     .......... form p ..........
     f = 0.0d0
  !c
     do 245 j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
  245    continue
  !c
     h = f / (h + h)
  !c     .......... form q ..........
     do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
  !c     .......... form reduced a ..........
     do 280 j = 1, l
        f = d(j)
        g = e(j)
  !c
        do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
  !c
  280    continue
  !c
  285    do 290 j = 1, l
        f = d(j)
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = f * scale
  290    continue
  !c
  300 continue
  !c
  return
  end

  subroutine tred2(nm,n,a,d,e,z)

  integer i,j,k,l,n,ii,nm,jp1
  real(rk) a(nm,n),d(n),e(n),z(nm,n)
  real(rk) f,g,h,hh,scale
  !c
  !c     this subroutine is a translation of the algol procedure tred2,
  !c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
  !c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
  !c
  !c     this subroutine reduces a real symmetric matrix to a
  !c     symmetric tridiagonal matrix using and accumulating
  !c     orthogonal similarity transformations.
  !c
  !c     on input
  !c
  !c        nm must be set to the row dimension of two-dimensional
  !c          array parameters as declared in the calling program
  !c          dimension statement.
  !c
  !c        n is the order of the matrix.
  !c
  !c        a contains the real symmetric input matrix.  only the
  !c          lower triangle of the matrix need be supplied.
  !c
  !c     on output
  !c
  !c        d contains the diagonal elements of the tridiagonal matrix.
  !c
  !c        e contains the subdiagonal elements of the tridiagonal
  !c          matrix in its last n-1 positions.  e(1) is set to zero.
  !c
  !c        z contains the orthogonal transformation matrix
  !c          produced in the reduction.
  !c
  !c        a and z may coincide.  if distinct, a is unaltered.
  !c
  !c     questions and comments should be directed to burton s. garbow,
  !c     mathematics and computer science div, argonne national laboratory
  !c
  !c     this version dated august 1983.
  !c
  !c     ------------------------------------------------------------------
  !c
  do 100 i = 1, n
  !c
     do 80 j = i, n
  80    z(j,i) = a(j,i)
  !c
     d(i) = a(n,i)
  100 continue

  if (n .eq. 1) go to 510
  !c     .......... for i=n step -1 until 2 do -- ..........
  do 300 ii = 2, n
     i = n + 2 - ii
     l = i - 1
     h = 0.0d0
     scale = 0.0d0
     if (l .lt. 2) go to 130
  !c     .......... scale row (algol tol then not needed) ..........
     do 120 k = 1, l
  120    scale = scale + dabs(d(k))

     if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)

     do 135 j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0d0
        z(j,i) = 0.0d0
  135    continue

     go to 290

  140    do 150 k = 1, l
        d(k) = d(k) / scale
        h = h + d(k) * d(k)
  150    continue

     f = d(l)
     g = -dsign(dsqrt(h),f)
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g
  !     .......... form a*u ..........
     do 170 j = 1, l
  170    e(j) = 0.0d0

     do 240 j = 1, l
        f = d(j)
        z(j,i) = f
        g = e(j) + z(j,j) * f
        jp1 = j + 1
        if (l .lt. jp1) go to 220

        do 200 k = jp1, l
           g = g + z(k,j) * d(k)
           e(k) = e(k) + z(k,j) * f
  200       continue

  220       e(j) = g
  240    continue
  !c     .......... form p ..........
     f = 0.0d0

     do 245 j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
  245    continue

     hh = f / (h + h)
  !c     .......... form q ..........
     do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
  !c     .......... form reduced a ..........
     do 280 j = 1, l
        f = d(j)
        g = e(j)
  !c
        do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)

        d(j) = z(l,j)
        z(i,j) = 0.0d0
  280    continue

  290    d(i) = h
  300 continue
  !c     .......... accumulation of transformation matrices ..........
  do 500 i = 2, n
     l = i - 1
     z(n,l) = z(l,l)
     z(l,l) = 1.0d0
     h = d(i)
     if (h .eq. 0.0d0) go to 380

     do 330 k = 1, l
  330    d(k) = z(k,i) / h

     do 360 j = 1, l
        g = 0.0d0

        do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)

        do 360 k = 1, l
           z(k,j) = z(k,j) - g * d(k)
  360    continue

  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0

  500 continue

  510 do 520 i = 1, n
     d(i) = z(n,i)
     z(n,i) = 0.0d0
  520 continue

  z(n,n) = 1.0d0
  e(1) = 0.0d0
  return
  end



  subroutine tql1(n,d,e,ierr)
  integer i,j,l,m,n,ii,l1,l2,mml,ierr
  real(rk) d(n),e(n)
  real(rk) c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
  !c
  !c     this subroutine is a translation of the algol procedure tql1,
  !c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
  !c     wilkinson.
  !c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).

  !c     this subroutine finds the eigenvalues of a symmetric
  !c     tridiagonal matrix by the ql method.

  !c     on input

  !c        n is the order of the matrix.

  !c        d contains the diagonal elements of the input matrix.

  !c        e contains the subdiagonal elements of the input matrix
  !c          in its last n-1 positions.  e(1) is arbitrary.

  !c      on output

  !c        d contains the eigenvalues in ascending order.  if an
  !c          error exit is made, the eigenvalues are correct and
  !c          ordered for indices 1,2,...ierr-1, but may not be
  !c          the smallest eigenvalues.

  !c        e has been destroyed.

  !c        ierr is set to
  !c          zero       for normal return,
  !c          j          if the j-th eigenvalue has not been
  !c                     determined after 30 iterations.

  !c     calls pythag for  dsqrt(a*a + b*b) .

  !c     questions and comments should be directed to burton s. garbow,
  !c     mathematics and computer science div, argonne national laboratory

  !c     this version dated august 1983.

  !c     ------------------------------------------------------------------

  ierr = 0
  if (n .eq. 1) go to 1001

  do 100 i = 2, n
  100 e(i-1) = e(i)

  f = 0.0d0
  tst1 = 0.0d0
  e(n) = 0.0d0

  do 290 l = 1, n
     j = 0
     h = dabs(d(l)) + dabs(e(l))
     if (tst1 .lt. h) tst1 = h
  !c     .......... look for small sub-diagonal element ..........
     do 110 m = l, n
        tst2 = tst1 + dabs(e(m))
        if (tst2 .eq. tst1) go to 120
  !c     .......... e(n) is always zero, so there is no exit
  !c                through the bottom of the loop ..........
  110    continue

  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
     j = j + 1
  !c     .......... form shift ..........
     l1 = l + 1
     l2 = l1 + 1
     g = d(l)
     p = (d(l1) - g) / (2.0d0 * e(l))
     r = pythag(p,1.0d0)
     d(l) = e(l) / (p + dsign(r,p))
     d(l1) = e(l) * (p + dsign(r,p))
     dl1 = d(l1)
     h = g - d(l)
     if (l2 .gt. n) go to 145

     do 140 i = l2, n
  140    d(i) = d(i) - h

  145    f = f + h
  !c     .......... ql transformation ..........
     p = d(m)
     c = 1.0d0
     c2 = c
     el1 = e(l1)
     s = 0.0d0
     mml = m - l
  !c     .......... for i=m-1 step -1 until l do -- ..........
     do 200 ii = 1, mml
        c3 = c2
        c2 = c
        s2 = s
        i = m - ii
        g = c * e(i)
        h = c * p
        r = pythag(p,e(i))
        e(i+1) = s * r
        s = e(i) / r
        c = p / r
        p = c * d(i) - s * g
        d(i+1) = h + s * (c * g + s * d(i))
  200    continue

     p = -s * s2 * c3 * el1 * e(l) / dl1
     e(l) = s * p
     d(l) = c * p
     tst2 = tst1 + dabs(e(l))
     if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
  !c     .......... order eigenvalues ..........
     if (l .eq. 1) go to 250
  !c     .......... for i=l step -1 until 2 do -- ..........
     do 230 ii = 2, l
        i = l + 2 - ii
        if (p .ge. d(i-1)) go to 270
        d(i) = d(i-1)
  230    continue
  !c
  250    i = 1
  270    d(i) = p
  290 continue
  !c
  go to 1001
  !c     .......... set error -- no convergence to an
  !c                eigenvalue after 30 iterations ..........
  1000 ierr = l
  1001 return
  end


  subroutine tql2(nm,n,d,e,z,ierr)
  integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
  real(rk) d(n),e(n),z(nm,n)
  real(rk) c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
  !c
  !c     this subroutine is a translation of the algol procedure tql2,
  !c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
  !c     wilkinson.
  !c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
  !c
  !c     this subroutine finds the eigenvalues and eigenvectors
  !c     of a symmetric tridiagonal matrix by the ql method.
  !c     the eigenvectors of a full symmetric matrix can also
  !c     be found if  tred2  has been used to reduce this
  !c     full matrix to tridiagonal form.
  !c
  !c     on input
  !c
  !c        nm must be set to the row dimension of two-dimensional
  !c          array parameters as declared in the calling program
  !c          dimension statement.
  !c
  !c        n is the order of the matrix.
  !c
  !c        d contains the diagonal elements of the input matrix.
  !c
  !c        e contains the subdiagonal elements of the input matrix
  !c          in its last n-1 positions.  e(1) is arbitrary.
  !c
  !c        z contains the transformation matrix produced in the
  !c          reduction by  tred2, if performed.  if the eigenvectors
  !c          of the tridiagonal matrix are desired, z must contain
  !c          the identity matrix.
  !c
  !c      on output
  !c
  !c        d contains the eigenvalues in ascending order.  if an
  !c          error exit is made, the eigenvalues are correct but
  !c          unordered for indices 1,2,...,ierr-1.
  !c
  !c        e has been destroyed.
  !c
  !c        z contains orthonormal eigenvectors of the symmetric
  !c          tridiagonal (or full) matrix.  if an error exit is made,
  !c          z contains the eigenvectors associated with the stored
  !c          eigenvalues.
  !c
  !c        ierr is set to
  !c          zero       for normal return,
  !c          j          if the j-th eigenvalue has not been
  !c                     determined after 30 iterations.
  !c
  !c     calls pythag for  dsqrt(a*a + b*b) .
  !c
  !c     questions and comments should be directed to burton s. garbow,
  !c     mathematics and computer science div, argonne national laboratory
  !c
  !c     this version dated august 1983.
  !c
  !c     ------------------------------------------------------------------

  ierr = 0
  if (n .eq. 1) go to 1001
  !c
  do 100 i = 2, n
  100 e(i-1) = e(i)
  !c
  f = 0.0d0
  tst1 = 0.0d0
  e(n) = 0.0d0
  !c
  do 240 l = 1, n
     j = 0
     h = dabs(d(l)) + dabs(e(l))
     if (tst1 .lt. h) tst1 = h
  !c     .......... look for small sub-diagonal element ..........
     do 110 m = l, n
        tst2 = tst1 + dabs(e(m))
        if (tst2 .eq. tst1) go to 120
  !c     .......... e(n) is always zero, so there is no exit
  !c                through the bottom of the loop ..........
  110    continue
  !c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
     j = j + 1
  !c     .......... form shift ..........
     l1 = l + 1
     l2 = l1 + 1
     g = d(l)
     p = (d(l1) - g) / (2.0d0 * e(l))
     r = pythag(p,1.0d0)
     d(l) = e(l) / (p + dsign(r,p))
     d(l1) = e(l) * (p + dsign(r,p))
     dl1 = d(l1)
     h = g - d(l)
     if (l2 .gt. n) go to 145
  !c
     do 140 i = l2, n
  140    d(i) = d(i) - h
  !c
  145    f = f + h
  !c     .......... ql transformation ..........
     p = d(m)
     c = 1.0d0
     c2 = c
     el1 = e(l1)
     s = 0.0d0
     mml = m - l
  !c     .......... for i=m-1 step -1 until l do -- ..........
     do 200 ii = 1, mml
        c3 = c2
        c2 = c
        s2 = s
        i = m - ii
        g = c * e(i)
        h = c * p
        r = pythag(p,e(i))
        e(i+1) = s * r
        s = e(i) / r
        c = p / r
        p = c * d(i) - s * g
        d(i+1) = h + s * (c * g + s * d(i))
  !c     .......... form vector ..........
        do 180 k = 1, n
           h = z(k,i+1)
           z(k,i+1) = s * z(k,i) + c * h
           z(k,i) = c * z(k,i) - s * h
  180       continue
  !c
  200    continue
  !c
     p = -s * s2 * c3 * el1 * e(l) / dl1
     e(l) = s * p
     d(l) = c * p
     tst2 = tst1 + dabs(e(l))
     if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
  !c     .......... order eigenvalues and eigenvectors ..........
  do 300 ii = 2, n
     i = ii - 1
     k = i
     p = d(i)
  !c
     do 260 j = ii, n
        if (d(j) .ge. p) go to 260
        k = j
        p = d(j)
  260    continue
  !c
     if (k .eq. i) go to 300
     d(k) = d(i)
     d(i) = p
  !c
     do 280 j = 1, n
        p = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = p
  280    continue
  !c
  300 continue
  !c
  go to 1001
  !c     .......... set error -- no convergence to an
  !c                eigenvalue after 30 iterations ..........
  1000 ierr = l
  1001 return
  end

   function pythag(a,b)
   real(rk)::pythag
   real(rk):: a,b
   !c
   !c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
   !c
   real(rk) p,r,s,t,u
   p = dmax1(dabs(a),dabs(b))
   if (p .eq. 0.0d0) go to 20
   r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
   t = 4.0d0 + r
   if (t .eq. 4.0d0) go to 20
   s = r/t
   u = 1.0d0 + 2.0d0*s
   p = u*p
   r = (s/u)**2 * r
   go to 10
   20 pythag = p
   return
   end function pythag


end module
