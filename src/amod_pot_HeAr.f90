module amod_pot_HeAr 
    implicit none 
    integer,parameter::rk = selected_real_kind(15,307)
    private
    public:: pot_HBePlus,pot_Hestar_Ar,pot_He_ArPlus,Gamma_R_approx,Gamma_R

    contains 
function f_E(R,R0,beta,B)
  !PhysRevA.18.1435
  real(rk),intent(in)::R,R0,beta,B
  real(rk)::f_E
  f_E = (1+exp(-(R-R0)/B))**(-1) &
  *(1+beta*(1+exp((R-R0)/B))**(-1))
end function f_E 

function V_lr(C6,C4,R)
  !PhysRevA.18.1435
  !long range potential 
    real(rk), intent(in) :: R, C6,C4
    real(rk):: V_lr
    ! Check for non-positive distance
    if (R <= 0.0_rk) then
        print *, "Distance must be positive."
        V_lr = 0.0_rk
        return
    endif
    ! Calculate the potential
    V_lr = -C6 / R**6-C4/R**4
end function 

function pot_Hestar_Ar(R)
  !PhysRevA.18.1435
  !He^* + Ar
  real(rk),intent(in)::R
  real(rk):: pot_Hestar_Ar,R0 
  real(rk),parameter::  A_au = 4.39678_rk,&
                        b_a0 = 0.9675_rk,&
                        epsilon_m_au = 0.1583e-3_rk,&
                        R_m_a0 = 9.4488_rk,&
                        C_6_aua0s6 = 226.0_rk,&
                        C_4_aua0s4 = 0.0_rk,& 
                        beta=2.81505_rk 

  pot_Hestar_Ar = 0.0_rk
  R0 = R_m_a0 - log((beta+1)/(beta-1))*b_a0
  pot_Hestar_Ar = A_au*exp(-R/b_a0) &
   + V_lr(C_6_aua0s6,C_4_aua0s4,R)*f_E(R,R0,beta,b_a0)
end function 

function pot_He_ArPlus(R)
  !PhysRevA.18.1435
  !He + Ar^+
  real(rk),intent(in)::R
  real(rk):: pot_He_ArPlus,R0
  real(rk),parameter::  A_au = 44.6527_rk,&
                        b_a0 = 0.532_rk,&
                        epsilon_m_au = 0.945e-3_rk,&
                        R_m_a0 = 5.936_rk,&
                        C_6_aua0s6 = 7.94_rk,&
                        C_4_aua0s4 = 0.6904_rk,& 
                        beta=6.42317_rk 
  pot_He_ArPlus = 0.0_rk
  R0 = R_m_a0 - log((beta+1)/(beta-1))*b_a0
  pot_He_ArPlus = A_au*exp(-R/b_a0) &
  + V_lr(C_6_aua0s6,C_4_aua0s4,R)*f_E(R,R0,beta,b_a0)

end function 

function Gamma_R(R)
  !PhysRevA.18.1435
  !He^* + Ar 
 
  real(rk),intent(in):: R 
  real(rk):: Gamma_R
  real(rk),parameter:: A_G = 1e-10_rk!13400.0_rk!1.0_rk! 73000_rk !4000.0_rk!  7400_rk !
  real(rk),parameter:: aG = 10.0_rk! 0.357_rk!0.667_rk!0.34    !0.36_rk! 0.357_rk  !
  Gamma_R = A_G*exp(-R/aG)!*EV2AU
end function Gamma_R

function Gamma_R_approx(R)
  !PhysRevA.18.1435
  !He^* + Ar 
  real(rk),intent(in):: R 
  real(rk):: Gamma_R_approx
  real(rk),parameter:: A_m =100.0_rk!23.41_rk! 50.0_rk! 31.85_rk!  1.0_rk
  real(rk),parameter:: am = 0.68_rk!0.72_rk ! 0.7_rk! 0.714_rk!  1.334_rk!
  real(rk),parameter:: epsilon_0 = 0.1505_rk 

  Gamma_R_approx = (4/sqrt(2*epsilon_0))*A_m**2*exp(-2*R/am)
end function Gamma_R_approx


subroutine pot_HBePlus(nx,rpts,pot1,pot2,pot3,u12,u23,u13,dipx,dipa,dipb)
implicit none
integer:: ir,nx,spline,irr
real(8):: rpts(nx),pot1(nx),pot2(nx),pot3(nx),u12(nx),u23(nx),u13(nx)
real(8):: x_pot1(25), y_pot1(25), x_pot2(25),y_pot2(25), x_pot3(25),y_pot3(25),&
& x_dip12(23),y_dip12(123), x_dip23(23),y_dip23(23),x_dip13(20),y_dip13(20)
real(8):: dip_x(25),dip_a(25),dip_b(23),x_dip_x(25),x_dip_b(23)
real(8),dimension(nx):: dipa,dipb,dipx

open(77,file='dip_x.txt')
do irr=1,25
    read(77,*)x_dip_x(irr),dip_x(irr)
enddo
close(77)

do irr=1,nx
    spline=-1
    call enspl(x_dip_x,dip_x,25,spline,rpts(irr),dipx(irr))
enddo

open(66,file='dip_a.txt')
do irr=1,25
    read(66,*)dip_a(irr)
enddo
close(66)
do irr=1,nx
    spline=-1
    call enspl(x_dip_x,dip_a,25,spline,rpts(irr),dipa(irr))
enddo

open(66,file='dip_b.txt')
do irr=1,23
    read(66,*)x_dip_b(irr),dip_b(irr)
enddo
  close(66)
do irr=1,nx
    spline=-1
    call enspl(x_dip_b,dip_b,23,spline,rpts(irr),dipb(irr))
enddo

open(11,file='G.txt')
do irr=1,25
    read(11,*)x_pot1(irr),y_pot1(irr)
enddo
close(11)
do irr=1,nx
    spline=-1
    call enspl(x_pot1,y_pot1,25,spline,rpts(irr),pot1(irr))
!write(11,*) rpts(irr),pot1(irr)
enddo

open(22,file='A.txt')
do irr=1,25
    read(22,*)x_pot2(irr),y_pot2(irr)
enddo
close(22)
do irr=1,nx
    spline=-1
    call enspl(x_pot2,y_pot2,25,spline,rpts(irr),pot2(irr))
    write(11,*) rpts(irr),pot1(irr)
enddo

open(33,file='B.txt')
do irr=1,25
    read(33,*)x_pot3(irr),y_pot3(irr)
enddo
close(33)
do irr=1,nx
    spline=-1
    call enspl(x_pot3,y_pot3,25,spline,rpts(irr),pot3(irr))
write(22,*) rpts(irr),pot2(irr)
enddo

open(44,file='GA.txt')
do irr=1,23
    read(44,*)x_dip12(irr),y_dip12(irr)
enddo
close(44)
do irr=1,nx
    spline=-1
    call enspl(x_dip12,y_dip12,23,spline,rpts(irr),u12(irr))
enddo


open(66,file='AB.txt')
do irr=1,23
    read(66,*)x_dip23(irr),y_dip23(irr)
enddo
close(66)
do irr=1,nx
    spline=-1
    call enspl(x_dip23,y_dip23,23,spline,rpts(irr),u23(irr))
enddo

open(66,file='GB.txt')
do irr=1,20
    read(66,*)x_dip13(irr),y_dip13(irr)
enddo
close(66)
do irr=1,nx
    spline=-1
    call enspl(x_dip13,y_dip13,20,spline,rpts(irr),u13(irr))
enddo

return
endsubroutine
        
SUBROUTINE ENSPL(X,Y,N,K,T,Z)
implicit real(rk)(A-H,O-Z)
implicit integer(I-N)
DIMENSION X(N),Y(N)
real(8) X,Y,T,Z,G1,G2,U1,U2,U3,U4,U5,A,B,C,D,S
Z=0.0
IF (N.LE.0) RETURN
IF (N.EQ.1) THEN
  Z=Y(1)
  K=0
  A=Y(1)
  B=0.0
  C=0.0
  D=0.0
  RETURN
END IF
IF (N.EQ.2) THEN
  Z=(Y(1)*(T-X(2))-Y(2)*(T-X(1)))/(X(1)-X(2))
  K=1
  A=Y(1)
  B=(Y(2)-Y(1))/(X(2)-X(1))
  C=0.0
  D=0.0
  RETURN
END IF
S=1.0
IF (K.LE.0) THEN
  S=-1.0
  IF (T.LE.X(2)) THEN
    K=1
  ELSE IF (T.GE.X(N)) THEN
    K=N-1
  ELSE
    K=1
    M=N
10  IF (IABS(K-M).NE.1) THEN
      L=(K+M)/2
      IF (T.LT.X(L)) THEN
        M=L
      ELSE
        K=L
      END IF
      GOTO 10
    END IF
  END IF
END IF
IF (K.GE.N) K=N-1
U3=(Y(K+1)-Y(K))/(X(K+1)-X(K))
IF (N.EQ.3) THEN
  IF (K.EQ.1) THEN
    U4=(Y(3)-Y(2))/(X(3)-X(2))
    U5=2.0*U4-U3
    U2=2.0*U3-U4
    U1=2.0*U2-U3
  ELSE IF (K.EQ.2) THEN
    U2=(Y(2)-Y(1))/(X(2)-X(1))
    U1=2.0*U2-U3
    U4=2.0*U3-U2
    U5=2.0*U4-U3
  END IF
ELSE
  IF (K.LE.2) THEN
    U4=(Y(K+2)-Y(K+1))/(X(K+2)-X(K+1))
    IF (K.EQ.2) THEN
      U2=(Y(2)-Y(1))/(X(2)-X(1))
      U1=2*U2-U3
      IF (N.EQ.4) THEN
        U5=2.0*U4-U3
      ELSE
        U5=(Y(5)-Y(4))/(X(5)-X(4))
      END IF
    ELSE
      U2=2*U3-U4
      U1=2*U2-U3
      U5=(Y(4)-Y(3))/(X(4)-X(3))
    END IF
  ELSE IF (K.GE.(N-2)) THEN
    U2=(Y(K)-Y(K-1))/(X(K)-X(K-1))
    IF (K.EQ.(N-2)) THEN
      U4=(Y(N)-Y(N-1))/(X(N)-X(N-1))
      U5=2*U4-U3
      IF (N.EQ.4) THEN
        U1=2.0*U2-U3
      ELSE
        U1=(Y(K-1)-Y(K-2))/(X(K-1)-X(K-2))
      END IF
    ELSE
      U4=2*U3-U2
      U5=2*U4-U3
      U1=(Y(K-1)-Y(K-2))/(X(K-1)-X(K-2))
    END IF
  ELSE
    U2=(Y(K)-Y(K-1))/(X(K)-X(K-1))
    U1=(Y(K-1)-Y(K-2))/(X(K-1)-X(K-2))
    U4=(Y(K+2)-Y(K+1))/(X(K+2)-X(K+1))
    U5=(Y(K+3)-Y(K+2))/(X(K+3)-X(K+2))
  END IF
END IF
A=ABS(U4-U3)
B=ABS(U1-U2)
IF ((A+1.0.EQ.1.0).AND.(B+1.0.EQ.1.0)) THEN
  G1=(U2+U3)/2.0
ELSE
  G1=(A*U2+B*U3)/(A+B)
END IF
A=ABS(U4-U5)
B=ABS(U3-U2)
IF ((A+1.0.EQ.1.0).AND.(B+1.0.EQ.1.0)) THEN
  G2=(U3+U4)/2.0
ELSE
  G2=(A*U3+B*U4)/(A+B)
END IF
A=Y(K)
B=G1
D=X(K+1)-X(K)
C=(3*U3-2*G1-G2)/D
D=(G2+G1-2*U3)/(D*D)
IF (S.LT.0.0) THEN
  S=T-X(K)
  Z=A+B*S+C*S*S+D*S*S*S
END IF
RETURN
END
         
         
end module amod_pot_HeAr