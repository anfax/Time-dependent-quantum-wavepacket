!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  gausdriv                                                *
!  *   Function   :  the drive routine for gauss quadrature                  *
!  *                                                                         *
!  ***************************************************************************
        module GaussianDrive
        private
        public::GAUSDRIV
        contains
       subroutine GAUSDRIV(kind,nodes,xnode,wnode)
       implicit none
c
c......declare variables

       integer :: memstatus,kind,nodes,kpts
       real*8  :: xnode(nodes),wnode(nodes),alpha,beta,endpts(2)

       real*8,dimension(:),allocatable :: dworks
       allocate(dworks(nodes),stat=memstatus)
       if (memstatus.ne.0) stop 'failed to allocate memory in gausdriv'
c
c......for gauss-legendre and gauss-chebyshev quadrature

       if (kind.ne.1.and.kind.ne.2)
     &    stop 'this value of kind is not availble,use GAUSSQ directly'

       kpts=0; alpha=0.0d0; beta=0.0d0
       call GAUSSQ(kind,nodes,alpha,beta,kpts,endpts,dworks,xnode,wnode)
c
c......release memory

       deallocate(dworks,stat=memstatus)
       if (memstatus.ne.0) stop 'failed to release memory in gausdriv'
       return
       end
c-------------------------------------------------------------------------
c
c     Gauss quadrature
c
c-------------------------------------------------------------------------
      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)
      IMPLICIT REAL*8 (A-H,O-Z)
      implicit integer (i-n) 
C
C           THIS SET OF ROUTINES COMPUTES THE NODES X(I) AND WEIGHTS
C        C(I) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE
C
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX
C
C                              N
C        BY                   SUM C  F(X )
C                             I=1  I    I
C
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.
C
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF
C        ORTHOGONAL POLYNOMIALS.  THE NODES X(I) ARE JUST THE ZEROES
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.
C
C     INPUT PARAMETERS
C
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF
C                 QUADRATURE RULE
C
C        KIND = 1=  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)
C        KIND = 2=  CHEBYSHEV QUADRATURE OF THE FIRST KIND
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)
C        KIND = 3=  CHEBYSHEV QUADRATURE OF THE SECOND KIND
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)
C        KIND = 4=  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON
C                   (-INFINITY, +INFINITY)
C        KIND = 5=  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.
C                   NOTE= KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.
C        KIND = 6=  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1
C
C        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE
C        ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-
C                 LAGUERRE QUADRATURE (OTHERWISE USE 0.).
C        BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE--
C                 (OTHERWISE USE 0.).
C        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-
C                 POINT (OR BOTH) OF THE INTERVAL IS REQUIRED TO BE A
C                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO
C                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED
C                 ENDPOINTS (1 OR 2).
C        ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF
C                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.
C        B        REAL SCRATCH ARRAY OF LENGTH N
C
C     OUTPUT PARAMETERS (BOTH ARRAYS OF LENGTH N)
C
C        T        WILL CONTAIN THE DESIRED NODES X(1),,,X(N)
C        W        WILL CONTAIN THE DESIRED WEIGHTS C(1),,,C(N)
C
C     SUBROUTINES REQUIRED
C
C        GBSLVE, CLASS, AND GBTQL2 ARE PROVIDED. UNDERFLOW MAY SOMETIMES
C        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE
C        TURNED OFF AS THEY ARE ON THIS MACHINE.
C
C     ACCURACY
C
C        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,
C        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP
C        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,
C        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT
C        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR
C        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME
C        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,
C        COMPLETELY HARMLESS.
C
C     METHOD
C
C           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION
C        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE
C        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE
C        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH
C        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF
C        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,
C        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A
C        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.
C        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF
C        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.
C
C     REFERENCES
C
C        1.  GOLUB, G. H., AND WELSCH, J. H.,  CALCULATION OF GAUSSIAN
C            QUADRATURE RULES,  MATHEMATICS OF COMPUTATION 23 (APRIL,
C            1969), PP. 221-230.
C        2.  GOLUB, G. H.,  SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,
C            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).
C        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE-
C            HALL, ENGLEWOOD CLIFFS, N.J., 1966.
C
C     ..................................................................
C
      REAL*8  MUZERO
      DIMENSION  B(N),T(N),W(N),ENDPTS(2)
C
      CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)
C
C           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.
C           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY
C           B THE OFF-DIAGONAL ELEMENTS.
C           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2
C           SUBMATRIX.
C
      IF (KPTS.EQ.0)  GO TO 100
      IF (KPTS.EQ.2)  GO TO  50
C
C           IF KPTS=1, ONLY T(N) MUST BE CHANGED
C
      T(N) =GBSLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)
      GO TO 100
C
C           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED
C
   50 GAM =GBSLVE(ENDPTS(1), N, T, B)
      T1 = ((ENDPTS(1) - ENDPTS(2))/(GBSLVE(ENDPTS(2), N, T, B) - GAM))
      B(N-1) =  SQRT(T1)
      T(N) = ENDPTS(1) + GAM*T1
C
C           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
C           AND THUS THE VALUE OF B(N) IS ARBITRARY.
C           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL
C           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.
C           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING
C
  100 W(1) = 1.0D0
      DO 105 I = 2, N
  105    W(I) = 0.0D0
C
      CALL GBTQL2 (N, T, B, W, IERR)
      DO 110 I = 1, N
  110    W(I) = MUZERO * W(I) * W(I)
C
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION GBSLVE(SHIFT, N, A, B)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION
C
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,
C
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN
C       THE N-TH POSITION.
C
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.
C
C
      DIMENSION  A(N),B(N)
C
      ALPHA = A(1) - SHIFT
      NM1 = N - 1
      DO 10 I = 2, NM1
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA
      GBSLVE = 1.0D0  /ALPHA
      RETURN
      END
C
C
C
      SUBROUTINE CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)
      IMPLICIT REAL*8 (A-H,O-Z)
      implicit integer (i-n) 
C
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE
C        RECURRENCE RELATION
C
C             B P (X) = (X - A ) P   (X) - B   P   (X)
C              J J            J   J-1       J-1 J-2
C
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS,
C        AND THE ZERO-TH MOMENT
C
C             MUZERO = INTEGRAL W(X) DX
C
C        OF THE GIVEN POLYNOMIAL   WEIGHT FUNCTION W(X).  SINCE THE
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS
C        GUARANTEED TO BE SYMMETRIC.
C
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS
C        REQUIRE THE GAMMA FUNCTION.
C
C     ..................................................................
C
      DIMENSION  A(N),B(N)
      REAL*8  MUZERO
      DATA PI / 3.141592653589793D0  /
C
      NM1 = N - 1
      GO TO (10, 20, 30, 40, 50, 60), KIND
C
C              KIND = 1=  LEGENDRE POLYNOMIALS P(X)
C              ON (-1, +1), W(X) = 1.
C
   10 MUZERO = 2.0D0
      DO 11 I = 1, NM1
         A(I) = 0.0D0
         ABI = I
   11    B(I) = ABI/ SQRT(4*ABI*ABI - 1.0D0  )
      A(N) = 0.0D0
      RETURN
C
C              KIND = 2=  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)
C
   20 MUZERO = PI
      DO 21 I = 1, NM1
         A(I) = 0.0D0
   21    B(I) = 0.5D0
      B(1) =  SQRT(0.5D0  )
      A(N) = 0.0D0
      RETURN
C
C              KIND = 3=  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)
C              ON (-1, +1), W(X) = SQRT(1 - X*X)
C
   30 MUZERO = PI/2.0D0
      DO 31 I = 1, NM1
         A(I) = 0.0D0
   31    B(I) = 0.5D0
      A(N) = 0.0D0
      RETURN
C
C              KIND = 4=  HERMITE POLYNOMIALS H(X) ON (-INFINITY,
C              +INFINITY), W(X) = EXP(-X**2)
C
   40 MUZERO =  SQRT(PI)
      DO 41 I = 1, NM1
         A(I) = 0.0D0
   41    B(I) =  SQRT(I/2.0D0  )
      A(N) = 0.0D0
      RETURN
C
C              KIND = 5=  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND
C              BETA GREATER THAN -1
C
   50 AB = ALPHA + BETA
      ABI = 2.0D0   + AB
      MUZERO = 2.0D0   ** (AB + 1.0D0  ) * DGAMMA(ALPHA + 1.0D0  ) * DGA
     VMMA(
     X BETA + 1.0D0  ) / DGAMMA(ABI)
      A(1) = (BETA - ALPHA)/ABI
      B(1) =  SQRT(4.0D0  *(1.0D0   + ALPHA)*(1.0D0   + BETA)/((ABI + 1.
     V0D0  )*
     1  ABI*ABI))
      A2B2 = BETA*BETA - ALPHA*ALPHA
      DO 51 I = 2, NM1
         ABI = 2.0D0  *I + AB
         A(I) = A2B2/((ABI - 2.0D0  )*ABI)
   51    B(I) =  SQRT (4.0D0  *I*(I + ALPHA)*(I + BETA)*(I + AB)/
     1   ((ABI*ABI - 1)*ABI*ABI))
      ABI = 2.0D0  *N + AB
      A(N) = A2B2/((ABI - 2.0D0  )*ABI)
      RETURN
C
C              KIND = 6=  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER
C              THAN -1.
C
   60 MUZERO = DGAMMA(ALPHA + 1.0D0  )
      DO 61 I = 1, NM1
         A(I) = 2.0D0  *I - 1.0D0   + ALPHA
   61    B(I) =  SQRT(I*(I + ALPHA))
      A(N) = 2.0D0  *N - 1 + ALPHA
      RETURN
      END
C     ------------------------------------------------------------------
C
      SUBROUTINE GBTQL2(N, D, E, Z, IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL
C     METHOD, AND IS ADAPTED FROM THE EISPAK ROUTINE IMTQL2
C
C     ON INPUT=
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;
C
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.
C
C      ON OUTPUT=
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;
C
C        E HAS BEEN DESTROYED;
C
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES;
C
C        IERR IS SET TO
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     ------------------------------------------------------------------
C
      INTEGER I, J, K, L, M, N, II, MML, IERR
      DIMENSION  D(N),E(N),Z(N)
      REAL*8  MACHEP
C
C     ========== MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ==========
       MACHEP=1.0D-14
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      E(N) = 0.0D0
      DO 240 L = 1, N
         J = 0
C     ========== LOOK FOR SMALL SUB-DIAGONAL ELEMENT ==========
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            IF ( ABS(E(M)) .LE. MACHEP * ( ABS(D(M)) +  ABS(D(M+1))))
     X         GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ========== FORM SHIFT ==========
         G = (D(L+1) - P) / (2.0D0   * E(L))
         R =  SQRT(G*G+1.0D0  )
         G = D(M) - P + E(L) / (G +  SIGN(R, G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     ========== FOR I=M-1 STEP -1 UNTIL L DO -- ==========
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF ( ABS(F) .LT.  ABS(G)) GO TO 150
            C = G / F
            R =  SQRT(C*C+1.0D0  )
            E(I+1) = F * R
            S = 1.0D0   / R
            C = C * S
            GO TO 160
  150       S = F / G
            R =  SQRT(S*S+1.0D0  )
            E(I+1) = G * R
            C = 1.0D0   / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0   * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
C     ========== FORM FIRST COMPONENT OF VECTOR ==========
            F = Z(I+1)
            Z(I+1) = S * Z(I) + C * F
            Z(I) = C * Z(I) - S * F
C
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0
         GO TO 105
  240 CONTINUE
C     ========== ORDER EIGENVALUES AND EIGENVECTORS ==========
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         P = Z(I)
         Z(I) = Z(K)
         Z(K) = P
C
  300 CONTINUE
C
      GO TO 1001
C     ========== SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ==========
 1000 IERR = L
 1001 RETURN
C     ========== LAST CARD OF GBTQL2 ==========
      END

      DOUBLE PRECISION FUNCTION  DGAMMA(Z)
      IMPLICIT REAL*8 (A-H,O-Z)
C  THIS IS A PROCEDURE THAT EVALUATES GAMMA(Z) FOR
C     0 LT Z LE 3 TO 16 SIGNIFICANT FIGURES
C    IT IS BASED ON A CHEBYSHEV-TYPE POLYNOMIAL
C   APPROXIMATION GIVEN IN H. WERNER AND R. COLLINGE, MATH. COMPUT.
C    15 (1961), PP. 195-97.
C   APPROXIMATIONS TO THE GAMMA FUNCTION, ACCURATE UP TO 18 SIGNIFICANT
C   DIGITS, MAY BE FOUND IN THE PAPER QUOTED ABOVE
C
C
C
      DIMENSION  A(18)
C
       A(1)=1.0D0
       A(2)=.4227843350984678D0
       A(3)=.4118403304263672D0
      A(4)=.0815769192502609D0
      A(5)=.0742490106800904D0
      A(6)=-.0002669810333484D0
      A(7)=.0111540360240344D0
      A(8)=-.0028525821446197D0
      A(9)=.0021036287024598D0
      A(10)=-.0009184843690991D0
      A(11)=.0004874227944768D0
      A(12)=-.0002347204018919D0
      A(13)=.0001115339519666D0
      A(14)=-.0000478747983834D0
      A(15)=.0000175102727179D0
      A(16)=-.0000049203750904D0
      A(17)=.0000009199156407D0
      A(18)=-.0000000839940496D0
C
C
C
      IF(Z.LE.1.0D0  ) GO TO 10
      IF(Z.LE.2.0D0  ) GO TO 20
      T=Z-2.0D0
      GO TO 30
10    T=Z
      GO TO 30
20    T=Z-1.0D0
30    P=A(18)
      DO 40 K1=1,17
      K=18-K1
      P=T*P+A(K)
40    CONTINUE
C
      IF(Z.GT.2.0D0  ) GO TO 50
      IF(Z.GT.1.0D0  ) GO TO 60
      DGAMMA=P/(Z*(Z+1.0D0  ))
      RETURN
60    DGAMMA=P/Z
      RETURN
50    DGAMMA=P
      RETURN
      END
      end module GaussianDrive


