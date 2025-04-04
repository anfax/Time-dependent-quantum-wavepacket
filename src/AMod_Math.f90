
module math

implicit none
integer,parameter::rk=8 
integer,parameter::dp = rk
private
public::FGHLD,RS,tred1,tred2,tql1,tql2,pythag,normalize,lengdrePolynomial,QsortC,matmat_real8
contains
! 子程序：diagonalize_matrix
! 描述：计算一般复数矩阵的特征值和特征向量，并打印结果
! 输入参数：
!     a: 输入复数矩阵
! 输出参数：
!     w: 特征值
!     vr: 右特征向量
! 内部变量：
!     n: 矩阵大小
!     lwork: 工作空间大小
!     work: 工作数组
!     rwork: 实数工作数组
!     info: 返回信息


! 辅助子程序：print_vector
! 描述：打印复数向量
! 输入参数：
!     v: 复数向量

subroutine print_vector(v)
   complex(rk), dimension(:) :: v
   integer :: i
   do i = 1, size(v)
       write(*, '("(",F6.2,",",F6.2,")")') v(i)
   end do
end subroutine print_vector

! 辅助子程序：print_matrix
! 描述：打印复数矩阵
! 输入参数：
!     m: 复数矩阵

subroutine print_matrix(m)
   complex(rk), dimension(:, :) :: m
   integer :: i, j
   do i = 1, size(m, 1)
       write(*, '(4("(",F6.2,",",F6.2,")"))') m(i, :)
   end do
end subroutine print_matrix

function lengdrePolynomial(l_max,x) result(f)
   integer,intent(in)::l_max
   integer::l
   real(rk),intent(in):: x
   real(rk):: f
   real(rk):: p(l_max)
   p=0
   if (l_max ==1 )then
      f = sqrt(0.5_rk)
   else
      p(1) = 1
      p(2) = X
      do l=2,l_max-1
         p(l+1)=real((2*l-1)*x*p(l)-(l-1)*p(l-1),rk)/real(l,rk)
      enddo
      f = sqrt(l_max-0.5_rk)*p(l_max)
   end if
end function lengdrePolynomial

function factorial(n) result(f)
integer,Intent(in)::n
integer::i,f
f = 1
if (n == 0) then
 f = 1
else
   do i=1,n,1
      f= f*i
   end do
end if
end function factorial

subroutine normalize(f,ff)
   complex(rk):: f(:),ff(:)
   real(rk)::N
   N = sum(abs(f)**2)
   ff = f / sqrt(N)
end subroutine

subroutine diag(V,vv)
   real(rk)::v(:),vv(:,:)
   integer::n,i
   v= size(V)
   vv = 0
   do i=1,N
      vv(i,i) = v(i)
   enddo
end subroutine diag

SUBROUTINE FGHLD(NDIM,NX,ZMU,ZL,VV,ZR,WCH,FV1,FV2,AR,IWRIT,NJJ)
  !C.....PROGRAM TO PERFORM FGH CALCULATION TO SET UP A 1 DIMENSIONAL
  !C.....WAVEFUNCTION ON A GRID.
  !C.....BASED ON C.C. MARSTON AND G.G. BALINT-KURTI, J. CHEM. PHYS.,
  !C.....91, 3571 (1989).
  !C.....THIS PROGRAM USES AN EVEN NUMBER OF GRID POINTS, AND THE
  !C.....THEORY OF THE ABOVE PAPER HAS BEEN EXTENDED TO ACCOMMODATE
  !C.....THIS.THE PROGRAM IS TO BE USED EITHER TO
  !C.....SET UP AN INITIAL VIBRATIONAL WAVEFUNCTION FOR A TIME DEPENEDENT
  !C.....QUANTUM DYNAMICAL CALCULATION, OR IT MAY BE USED TO SET UP THE
  !C.....COEFFICIENT ARRAY FOR THE FINAL STATE ANALYSIS.
  !C.....  THE FOLLOWING ARRAYS REPRESENT :-
  !C.....NDIM -- DIMENSION OF ARRAYS.
  !C.....AR  - THE HAMILTONIAN MATRIX
  !C.....NX -- NUMBER OF POINTS IN GRID, THIS IS 1 GREATER THAN THAT
  !C.....        USED IN THE REST OF THE PROGRAM TO MAKE IT ODD.
  !C.....ZMU -- REDUCED MASS.
  !C.....ZL -- LENGTH OF INTERVAL (ONE INCREMENT LONGER THAN IN MAIN
  !C.....      PART OF CODE.
  !C.....VV -- POTENTIAL VECTOR.
  !C.....ZR  - THE EIGENVECTOR (X,Y) : X = WAVEFUNCTION
  !C.....                              Y = THE ENERGg LEVEL
  !C.....WCH - THE EIGENVALUES FOR RESPECTIVE ENERGY LEVELS
  !C.....     FV1 - (2/NX)*T(L)  WHERE L = AN INTEGER
  !C.....                              NX = NUMBER OF GRID POINTS
  !C.....     FV2 - COS((2*PI*L)/NX)
  !use molecule, only : mrx

  IMPLICIT real(dp)(A-H,O-Z)
  implicit integer(i-n)
  DIMENSION AR(NDIM,NX),WCH(NX),ZR(NDIM,NX),VV(NX)
  DIMENSION FV1(NX),FV2(NX)
  INTEGER MATZ
  !C.....  THE FOLLOWING DATA CONSTANTS REPRESENT :-
  DATA PI/3.141592653589793D0/
  !C.....     NPRIN    - 0 OR 1 FOR EIGENVALUE OR EIGENVECTOR RESPECTIVELY
  DATA NPRIN/1/
  DATA MATZ/1/
  !C.....  TITLE FOR DATA TO BE PRINTED
  IF(IWRIT.GE.10) THEN
  !C      WRITE(6,*)' LOWEST 10 ENERGY LEVELS FOR THE DIATOMIC FRAGMENT'
  END IF
  !C.....  TEST THAT NX IS EVEN.
  ITEST=MOD(NX,2)
  IF(ITEST.NE.0) THEN
  !C         WRITE(6,*)' **** NX MUST BE EVEN __ FATAL ERROR  ****'
     STOP
  END IF

  NHALF=NX/2
  NHAM1=NHALF-1
  !C.....  NOW SET UP HAMILTONIAN MATRIX
  DARG=2.0D0*PI/DFLOAT(NX)
  TARG=4.0D0*((PI/ZL)**2.0D0)/(ZMU*DFLOAT(NX))
  !C.....  PRECALCULATE KINETIC ENERGY FACTOR.
  DO 001 L=1,NX
     FV1(L)=TARG*DFLOAT(L**2)
     FV2(L)=DCOS(DARG*DFLOAT(L))
  001   CONTINUE
  !C.....  INITIALISE VARIABLES
  !C.....    CONST - ((-1)**(I-J))*T(NX/2)
  !C.....    EQUIJ - WHEN I=J THEN SUMMATION ONLY THAT OF T(L)
  !C.....    INIJ  - (I-J)
  !C.....    SUM   - TOTAL SUMMATION OF EQUATION WITHIN SUM
  CONST=0.0D0
  EQUIJ=0.0D0

  !C.....  PRECALCULATE THE SUM OF THE T(L) WHEN I=J
  DO 002 L=1,NHAM1
     EQUIJ=EQUIJ+FV1(L)
  002   CONTINUE
  DO 003 I=1,NX
     DO 004 J=1,I
        IJ=(I-J)
        INIJ=IJ
        SUM=0.0D0
        IF(IJ.EQ.0) THEN
           SUM=EQUIJ
        ELSE
           DO 005 L=1,NHAM1
              SUM=SUM+FV1(L)*FV2(IJ)
              IJ=IJ+INIJ
  !C.....  COSINE IS PERIODIC, SO ONLY NEED VALUE WITHIN ONE CYCLE/PERIOD
              IF(IJ.GT.NX)IJ=MOD(IJ,NX)
  005            CONTINUE
        END IF
  !C.....  ADD IN T(NX/2)
        CONST=((-1)**INIJ)*FV1(NHALF)/2
        AR(I,J)=SUM+CONST
  004      CONTINUE
  !c      write(*,*)i

  !C.....  ADD THE POTENTIAL VALUE WHEN KNONICLEAR DELTA FUNCTION
  !C.....     EQUALS 1, I.E. WHEN I AND J ARE EQUAL
     AR(I,I)=AR(I,I)+VV(I)!+mrx(i)*(NJJ*(NJJ+1))
  003   CONTINUE

  !C.....  NOW FILL OUT HAMILTONIAN MATRIX.
  DO 006 I=1,NX
     DO 007 J=1,I
        AR(J,I)=AR(I,J)
  007      CONTINUE
  006   CONTINUE

  !C.....  NOW CALL EIGENVALUE SOLVER

  CALL RS(NDIM,NX,AR,WCH,MATZ,ZR,FV1,FV2,IERR)

  I=30
  !C.....  NOW WRITE OUT LOWEST 10 EIGENVALUES AND EIGENVECTORS.
  DO 009 I=1,10
     IF(IWRIT.GE.10)THEN
  !C         WRITE(6,*)' ENERGY LEVEL NO ',I,' EIGENVALUE= ',WCH(I)
     IF (NPRIN.EQ.1) THEN
        DO 010 J=1,NX
  !            WRITE(6,*)J,' WAVEFUNCTION=',ZR(J,I)
  010        CONTINUE
     ENDIF
     END IF
  009   CONTINUE
  DO I=1,NX
   AR(I,I)=AR(I,I)-VV(I)
  END DO
  do j=1,nx
   if(zr(1,j).lt.0.d0)then
      do k=1,nx
         zr(k,j)=-zr(k,j)
      end do
   end if
  end do
  RETURN
  END

  subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
  integer n,nm,ierr,matz
  real(dp) a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
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
  real(dp) a(nm,n),d(n),e(n),e2(n)
  real(dp) f,g,h,scale
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
  real(dp) a(nm,n),d(n),e(n),z(nm,n)
  real(dp) f,g,h,hh,scale
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
  real(dp) d(n),e(n)
  real(dp) c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
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
  real(dp) d(n),e(n),z(nm,n)
  real(dp) c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
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
   real(dp)::pythag
   real(dp):: a,b
   !c
   !c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
   !c
   real(dp) p,r,s,t,u
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

   function  matrix_T(A) result(B)
   real(rk),intent(in)::A(:,:)
   real(rk)::B(size(A(:,1)),size(A(1,:)))
   integer:: i,j
   do j = 1, size(A(1,:))
   do i = 1, size(A(:,1))
   B(i,j) = A(j,i)
   enddo
   enddo
   end function matrix_T

   function matmat_real8(a,b) result(c)
      real(8),intent(in)::a(:,:),b(:,:)
      real(8)::c(size(a(:,1)),size(b(1,:)))
      integer::na1,na2,nb1,nb2,i,j,k,l
         na1 = size(a(:,1))
         na2 = size(a(1,:))
         nb1 = size(b(:,1))
         nb2 = size(b(1,:))

         do j=1,nb2
            do i=1,nb1
               c(i,j)=sum(a(i,:)*b(:,j))
            end do
         end do

   end function


   recursive subroutine QsortC(A)
    real(rk), intent(in out), dimension(:) :: A
    integer :: iq
    if(size(A) > 1) then
      call Partition(A, iq)
      call QsortC(A(:iq-1))
      call QsortC(A(iq:))
    endif
  end subroutine QsortC

  subroutine Partition(A, marker)
   !for Qsort'C'
    real(rk), intent(in out), dimension(:) :: A
    integer, intent(out) :: marker
    integer :: i, j
    real(rk) :: temp
    real(rk) :: x      ! pivot point
    x = A(1)
    i= 0
    j= size(A) + 1
    do
      j = j-1
      do
        if (A(j) <= x) exit
        j = j-1
      end do
      i = i+1
      do
        if (A(i) >= x) exit
        i = i+1
      end do
      if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
      elseif (i == j) then
        marker = i+1
        return
      else
        marker = i
        return
      endif
    end do
  end subroutine Partition

function MTM(A,B) result(C)
    real(rk),intent(in)::A(:,:),B(:,:)
    real(Rk)::C(size(A(:,1)),size(B(1,:)))
    integer::m,n,p,g,h,i,j,k,l
    m=size(A(:,1))
    n=size(A(1,:))
    p=size(B(1,:))
    !print*,m,n,p
    if (n /= size(B(:,1))) stop "The size must matching!"
    c=0
    do h=1,m
      do g=1,p
        do l=1,p
          do k=1,n
            do j=1,n
              do i=1,m
        c(g,h)=c(g,h)+A(i,j)*B(l,k)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    end function

end module math
