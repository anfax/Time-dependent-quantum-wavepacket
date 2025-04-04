module scatter1_sudden 
    !======================================================================
    !This is the module for the 2D Time-Dependent Quantum Wavapacket Scattering 
    !======================================================================
    use omp_lib 
    implicit none 
    
    integer,parameter:: rk = 8 
    !nature constants 
    real(rk),parameter:: pi=3.141592653589_rk
    complex(rk),parameter:: eye=cmplx(0.0_rk,1.0_rk,rk)
    !atomic masses
    character(5)::mame_A, name_B, name_C ! the name of the atoms
    real(rk):: mass_A,mass_B, mass_C ! the mass of the atoms 
    real(rk):: MassABC,MassBC ! the reduced mass of the atoms about three and two 
    !Grid set
    integer:: nx,na !number of grid points in radial and angle
    integer:: NE ! the number of energy points 
    real(rk):: xmin,xmax,r0,len_x,flux_x0 ! the range of the grid in x
    real(rk):: emin,emax,e0,gau_delta ! the range of the energy grid 
    real(rk):: dx,da,de ! the step size in x, angle and energy
    real(rk),allocatable:: x(:),ang_cos(:),e(:),ang_weight(:),ang(:) ! the grid points in x, angle and energy
    complex(rk),allocatable:: WFing(:,:,:),WFtemp(:,:,:),ae(:) ! the initial and temporary wave function
    complex(rk),allocatable:: WFRE(:,:,:,:) ! the final wave function in x, angle and energy
    real(rk),allocatable:: U_ang(:,:),U_rad(:,:) ! the transport matrix of  angle and radial direction 
    real(rk),allocatable:: T_rad(:) ! the kinetic energy matrix of angle and radial direction
    real(rk),allocatable:: Vs(:,:),T_BC(:),T_ABC(:,:) ! the potential energy matrix of angle and radial direction
    real(rk):: tmin,dt,tlim !the begin time and delta time and the limit of time step in the propagation
    real(rk):: absorp_a,absorp_n,absorp_len  ! the absorption parameters cofficient, power number and length
    !Quantum number
    integer:: nstate,v0,j0,k0 ! the total rotational quantum number, the vibrational quantum number, the rotational quantum number and magnetic quantum number 
    !others
    character(10),parameter:: dir='./output/' ! the directory for the output files
    
    private
    public:: TDQW2D 
    contains
    subroutine TDQW2D()!(out_WFRE,out_x,out_a,out_e)
    !======================================================================
    !This is the program for the 2D TDQW
    !======================================================================
    integer:: i,j,ie 
    ! complex(rk),intent(out)::out_WFRE(nx,na,ne) ! the final wave function in x, angle and energy for output 
    ! real(rk),intent(out)::out_x(:),out_a(:),out_e(:) !the grid points in x, angle and energy for output 
    complex(rk),allocatable::out_WFRE(:,:,:,:) 
    
    
    call readparamters()
    
    call grid_set()
    
    call WF_initial() 
    ! out_x = x 
    ! out_a = ang_cos 
    ! out_e = e 
    call getpot() 
    allocate(out_WFRE(nx,na,nstate,ne))
    call propagation(out_WFRE)
    
    call LAR(out_WFRE)
    open(66,file=trim(adjustl(dir))//'WFRE1.dat')
    do j=1,na 
    do i=1,nx 
    write(66,'(*(f15.7,x))') x(i),ang_cos(j),(abs(WFRE(i,j,2,ie))**2,ie=1,ne)
    enddo 
    enddo 
    close(66)
    
    open(66,file=trim(adjustl(dir))//'WFRE2.dat')
    do j=1,na 
    do i=1,nx 
    write(66,'(*(f15.7,x))') x(i),ang_cos(j),(abs(WFRE(i,j,2,ie))**2,ie=1,ne)
    enddo 
    enddo 
    close(66)
    end subroutine  
    
    subroutine LAR(out_WFRE)
    !======================================================================
    !This is the program for the LAR
    !======================================================================
    use Math_Functions,only:Gamma
    complex(rk),intent(inout)::out_WFRE(:,:,:,:)
    integer:: i,j,k,ie
    complex(rk):: LARs(ne) 
    do ie =1,ne 
      do k=1,nstate 
        do j=1,na
          do i=1,nx
            LARs(ie) = LARs(ie) + out_WFRE(i,j,k,ie)*sqrt(Gamma(x(i)))*dconjg(out_WFRE(i,j,k,ie))
          enddo 
        enddo
      enddo
    enddo
    open(66,file=trim(adjustl(dir))//'LAR.dat')
    do ie =1,ne
    write(66,'(3(e15.7,x))') e(ie),real(LARs(ie),rk),aimag(LARs(ie))
    enddo
    close(66)
    end subroutine
    
    
    subroutine readparamters 
    use AtomMassInau,only:ATOMASS
    !======================================================================
    !This subroutine reads the parameters from the input file
    !======================================================================
    print*,'Reading parameters from the input file' 
    open(66,file='./input/1.in',status='old')
    read(66,*) nx,na,NE
    read(66,*) name_B, name_C
    read(66,*) xmin,xmax
    read(66,*) r0,e0,gau_delta 
    read(66,*) emin,emax,flux_x0
    read(66,*) tmin,dt,tlim
    read(66,*) absorp_n,absorp_a,absorp_len
    read(66,*) nstate,v0,j0,k0 
    close(66) 
    na=1
    write(109,*) 'The parameters are:'
    write(109,*) 'nx,na,NE:',nx,na,NE
    write(109,*) 'name_B,name_C:',name_B,name_C
    write(109,*) 'xmin,xmax:',xmin,xmax
    write(109,*) 'r0,e0,gau_delta:',r0,e0,gau_delta
    write(109,*) 'emin,emax:',emin,emax
    write(109,*) 'tmin,dt,tlim:',tmin,dt,tlim
    write(109,*) 'absorp_n,absorp_a,absorp_len:',absorp_n,absorp_a,absorp_len
    write(109,*) 'nstate,v0,j0,k0:',nstate,v0,j0,k0
    ! mass_A= 7295.54494_rk !He atom
    ! mass_B= 1836.9626949398998_rk 
    ! mass_C= 1836.9626949398998_rk 
    call ATOMASS(name_B,mass_B)
    call ATOMASS(name_C,mass_C)
    
    write(109,*)  'The masses are:',mass_B,mass_C
    MassBC = mass_B*mass_C/(mass_B+mass_C)
    ! MassABC = mass_A*(mass_B+Mass_C)/(mass_A+mass_B+mass_C)
    len_x = xmax-xmin 
    write(109,*) 'The length of the grid is:',len_x
    write(109,*) 'The reduced mass of BC is:',MassBC
    write(109,*) 'The parameters are:'
    write(109,*) 'nx,na,NE:',nx,na,NE
    write(109,*) 'name_B,name_C:',name_B,name_C
    write(109,*) 'xmin,xmax:',xmin,xmax
    write(109,*) 'r0,e0,gau_delta:',r0,e0,gau_delta
    write(109,*) 'emin,emax:',emin,emax
    write(109,*) 'tmin,dt,tlim:',tmin,dt,tlim
    write(109,*) 'absorp_n,absorp_a,absorp_len:',absorp_n,absorp_a,absorp_len
    write(109,*) 'nstate,v0,j0,k0:',nstate,v0,j0,k0
    write(109,*)  'The masses are:',mass_A,mass_B,mass_C
    write(109,*) 'The length of the grid is:',len_x
    write(109,*) 'The reduced mass of BC is:',MassBC
    write(109,*) 'The reduced mass of ABC is:',MassABC
    
    end subroutine
    
    
    subroutine grid_set()
    !======================================================================
    !This subroutine sets the grid
    !======================================================================
    integer:: i,j
    
    print*,'Setting the grid' 
    dx=(len_x)/(nx+1)
    de=(emax-emin)/(ne-1) 
    allocate(x(nx),ang_cos(na),e(ne),ae(ne))
    allocate(ang(na),source=0.0_rk) 

    !$omp parallel do private(i) 
    do i=1,nx
    x(i)=xmin+i*dx 
    end do
    !$omp end parallel do 
    
    !$omp parallel do private(i) 
    do i=1,ne 
    e(i)=emin+(i-1) *de
    end do
    !$omp end parallel do
    
    allocate(WFing(nx,na,nstate),source=(0.0_rk,0.0_rk))
    allocate(WFtemp(nx,na,nstate),source=(0.0_rk,0.0_rk)) 
    allocate(WFRE(nx,na,nstate,ne),source=(0.0_rk,0.0_rk))
    
    
    allocate(T_ABC(nx,na),source=0.0_rk)
    allocate(T_BC(na),source=0.0_rk)
    !T_total =0 
    !l(l+1)/(2 mBC R^2) 
    !$omp parallel do private(i)
    do j=1,na 
    do i=1,nx 
    T_ABC(i,j)=j0*(j0-1)/2.0_rk/massBC/x(i)**2
    end do
    end do
    !$omp end parallel do 
    
    end subroutine 
    
    
    subroutine WF_initial()
    !======================================================================
    !This subroutine initializes the wave function
    !======================================================================
    use GaussianDrive,only:GAUSDRIV
    use P_lm,only:yjm0
    complex(rk),allocatable::WF_X(:),WF_A(:)
    integer:: i,j 
    print*,'Initializes the wave function!'
    
    allocate(WF_X(nx)) 
    call wavepacket_gaussian(nx,x,r0,len_x,gau_delta,ne,e,e0,MassBC,-1.0_rk,WF_X,ae)
    print*,'The population of the initial radial wave function is:',sum(abs(WF_X)**2)
    !d^2/dR^2 
    allocate(T_rad(nx))
    !$omp parallel do private(i)
    do i=1,nx
    T_rad(i)=(i*pi)**2/(len_x)**2/2.0_rk/massBC
    end do
    !$omp end parallel do 
    !R:DVR -> FBR 
    allocate(U_rad(nx,nx))
    !$omp parallel do private(i,j) schedule(dynamic)
    do i=1,nx
    do j=1,nx
    U_rad(j,i)=U_ja(nx,i,j)!sqrt(dx)*sqrt(2/len_x)*sin(i*pi*x(j)/len_x)
    end do
    end do
    !$omp end parallel do 
    
    allocate(ang_weight(na))
    allocate(U_ang(na,na))
    
    call GAUSDRIV(1,na,ang_cos,ang_weight)
    do j=1,na 
        ang(j)=acos(ang_cos(j))
    end do

    !Angle: DVR -> FBR 
    !$omp parallel do private(i,j)  schedule(dynamic)
    do i=1,na
    do j=1,na 
    U_ang(j,i)=sqrt(ang_weight(i))*yjm0(j-1,k0,ang_cos(i))
    end do
    end do
    !$omp end parallel do 
    
    !the initial wavefuntion in angle dimension 
    allocate(WF_A(na))
    !$omp parallel do private(i) schedule(dynamic)
    do i=1,na 
    WF_A(i) = 1.0_rk!sqrt(ang_weight(i))*yjm0(j0,k0,ang_cos(i))
    enddo 
    !$omp end parallel do 
    print*,'The population of the initial angular wave function is:',sum(abs(WF_A)**2)
    
    !the initial wave function in 2D 
    
    !$omp parallel do private(i,j) schedule(dynamic)
    do i=1,nx
    do j=1,na
    WFing(i,j,1)=WF_X(i)*WF_A(j)
    end do
    end do
    !$omp end parallel do 
    
    print*,'The population of the initial wave function is:',sum(abs(WFing)**2)
    
    open(66,file=trim(dir)//'wf_t0.dat')
    
    do j=1,na 
    do i=1,nx 
    write(66,'(*(f15.7,x))') x(i)*ang_cos(j),x(i)*sqrt(1-ang_cos(j)**2),real(WFing(i,j,1),rk)/sqrt(ang_weight(j)),aimag(WFing(i,j,1))/sqrt(ang_weight(j))
    end do
    end do
    close(66) 
    
    end subroutine 
    
    subroutine Flux_Opertaor_sin(n,r_min,Len_r,r_0,m_r,F)
    ! the flux operator in the radial direction ( use FBR basis wavefuntion in the radial direction)
    integer,intent(in):: n
    real(rk),intent(in):: Len_r,r_0,m_r,r_min
    
    complex(rk),intent(out)::F(n,n)
    real(rk),parameter::PI=3.14159265358979323846_rk 
    
    integer::i,j
    
    F = 0.0_rk
    do j=1,n
    do i=1,n
    F(i,j) = -PI/(m_r*Len_r**2)*(j*sin(i*pi*(r_0-r_min)/Len_r)*cos(j*pi*(r_0-r_min)/Len_r) &
                          -i*sin(j*pi*(r_0-r_min)/Len_r)*cos(i*pi*(r_0-r_min)/Len_r))
    end do
    end do
    F=F*cmplx(0.0_rk,1.0_rk,rk) 
    end subroutine
    
    subroutine wavepacket_gaussian(n,rs,x_0,length,width,n_e,Es,E_0,MASS,arrow,wf,As)
    !The gaussian function as the initial wave function in radial direction 
    integer,intent(in)::n,n_e  !the number of grids and energy points 
    real(rk),intent(in)::x_0,rs(n),length,width,E_0,MASS,arrow,Es(n_e) !the parameters of the gaussian function
    complex(rk),intent(out)::wf(n),As(n_e) !the wave function and the normalization factor of energy 
    complex(rk):: wf_temp(n) 
    real(rk)::Sum_wf,MV
    integer::i
    MV = sqrt(2*E_0*MASS)
    do i=1,n
    wf_temp(i) = sqrt(length/(n+1)) /(2.0_rk*pi*width**2)**0.25_rk &
    *exp(-(rs(i)-x_0)**2/(4*width**2))*exp(arrow*eye*mv*rs(i))
    end do
    
    sum_wf=0
    Sum_wf = sum(abs(wf_temp)**2)
    wf = wf_temp!/sqrt(Sum_wf)
    ! wf = wf_temp 
    open(66,file = trim(dir)//'/gaussian.dat')
    do i=1,n
    write(66,'(*(f15.7,x))') rs(i),real(wf(i),rk),imag(wf(i))
    end do
    close(66)
    
    write(*,*) 'The population of Gaussian wavepacket is:',sum(abs(wf)**2)
    
    call A_E1(j0*1.0_rk,n_e,Es,MASS,n,rs,wf,As)
    
    end subroutine
    
    subroutine A_E1(ll,n_e,Es,Mass,nr,R,Psi,Aoe)
    ! the hankel function for get the a(E) for psi 
    use Math_Functions,only:rbesjy
    integer,intent(in)::n_e,nr
    real(rk),intent(in)::ll, Es(n_e),Mass,R(nr)
    complex(rk),intent(in)::Psi(:)
    complex(rk),intent(out)::Aoe(N_E)
    real(rk),allocatable,dimension(:)::P
    real(rk)::dR
    real(rk):: nl,jl,nlp,jlp,c_mr
    integer:: i,j,k,l
    
    allocate(P(N_E),source=0.0_rk)
    Aoe=0.0_rk
    dR=R(2)-R(1)
    c_mr = sqrt(dr)*sqrt(mass/2.0_rk/pi)
    
    do j=1,N_E
    P(j)=sqrt(2*Mass*Es(j))
    do i=1,Nr 
    call rbesjy((ll*1.0_rk),r(i)*p(j),nl,nlp,jl,jlp)
    if (jl >=5.0_rk ) jl =0.0_rk
    aoe(j)=aoe(j)+c_mr/sqrt(p(j))*psi(i)*dconjg(cmplx(jl,-nl,rk))
    end do
    end do
    
    
    open(66,file=trim(adjustl(dir))//'Aoe.dat')
    do i=1,N_E
    write(66,'(3(f15.7,x))')Es(i),real(aoe(i),rk),aimag(aoe(i))
    end do
    close(66)
    
    deallocate(P)
    end subroutine A_E1
    
    function U_ja(n,j,a) result(U)
    !j the index of basis
    !a the index of grids
    integer,intent(in)::n,j,a
    real(rk)::U
    u=sqrt(2.0_rk/(n+1))*sin(j*a*pi/real(n+1,rk))
    end function U_ja
    
    function basisSin(L,j,x,x_0)
    ! the basis function in the radial direction (sine function)
    real(rk),intent(in)::L,x,x_0
    integer,intent(in)::j
    real(rk)::basisSin
    if (x .ge. x_0 .and. x .le. (L+x_0) ) then
    basisSin = sqrt(2/L) *sin(j*pi*(x-x_0)/l)
    else
    basisSin =0
    endif
    end function basisSin
    
    function Vabs(n,a,x,x0,x1) result(vv)   
    ! the absorption potential in the radial direction 
    real(rk),intent(in):: n,a,x,x0,x1
    complex(rk):: vv
    if (x<x0) then
    vv =0
    else
    vv = - eye*a *(abs(x-x0)/(x1-x0))**n
    end if
    end function
    
    subroutine getpot()
    ! use amod_pot_KrRb,only:pot_ini=>KrStar_Rb,pot_fin=>Kr_RbPlus
      use amod_pot_HeAr,only:pot_ini=>pot_Hestar_Ar,pot_fin=>pot_He_ArPlus
    ! the potential energy matrix of angle and radial direction
    integer:: i,j 
    allocate(Vs(nx,nstate))
    
    do i=1,nx 
    ! call pot_ini(x(i),Vs(i,1))
    ! call pot_fin(x(i),Vs(i,2))
      vs(i,1) = pot_ini(x(i))
      vs(i,2) = pot_fin(x(i))
    end do 
    
    open(66,file=trim(adjustl(dir))//'Vs.dat')
    
    do i=1,nx
    write(66,'(3(f15.7,x))') x(i),(Vs(i,j),j=1,nstate)
    end do
    
    close(66)
    
    end subroutine 
    
    subroutine propagation(out_WFRE)
    use amod_propatation,only:split_s,T2E,U_state 
    use Math_Functions,only:Gamma
    integer::i,j,k,it,ie
    real(rk)::ting,poping
    complex(rk),allocatable:: Tx(:),TBC(:),TABC(:,:),VsC(:,:)
    complex(rk),allocatable:: Cin(:,:),Cout(:,:)
    real(rk),allocatable:: U_s(:,:,:)
    complex(rk),intent(out):: out_WFRE(nx,na,nstate,ne)
    character(20)::tlabel
    allocate(Cin(nstate,ne))
    allocate(Cout(nstate,ne))
    
    !$omp parallel do private(iE) schedule(dynamic)
    do ie=1,nE
    !Cin(ie) = eye*(E(iE)-minval(Vs(NX,:))*0.0_rk)!(Es(ie))
    Cin(1,ie) = eye*(E(iE)+(Vs(NX/2,1))*1.0_rk)
    Cout(1,ie)= dt/(sqrt(2.0_rk*pi)*ae(ie))
    Cin(2,ie) = eye*(E(iE)+(Vs(NX/2,2))*1.0_rk)
    Cout(2,ie)= dt/(sqrt(2.0_rk*pi)*ae(ie))
    enddo
    !$omp end parallel do
    ! exp(-i V dt /2 )
    allocate(VsC(nx,nstate))
    !$omp parallel do private(i,j) schedule(dynamic)
    do j=1,nstate 
    do i=1,nx 
    VsC(i,j)=exp(-eye*dt*( & 
    +Vabs(absorp_n,absorp_a,x(i),xmax-absorp_len,xmax) &
    +Vabs(absorp_n,absorp_a*0.1,x(i),xmax-absorp_len,xmax) &
    +eye*Gamma(x(i))*0.0_rk)*0.5_rk)
    end do
    end do
    !$omp end parallel do
    ! exp(-i T dt )
    allocate(Tx(nx)) 
    !$omp parallel do private(i) schedule(dynamic)
    do i=1,nx 
    Tx(i)=exp(-eye*dt*T_rad(i))
    end do
    !$omp end parallel do
    
    ! exp(-i l(l+1)/(2 mu R^2) dt )
    allocate(TABC(nx,na))
    !$omp parallel do private(i,j) schedule(dynamic)
    do j=1,na
    do i=1,nx
    TABC(i,j)=exp(-eye*dt*T_ABC(i,j))
    end do
    end do
    !$omp end parallel do 
    ! U_state(nstate,n,x,v,u_s)
    allocate(U_s(nstate,nstate,nx),source=0.0_rk)
    call U_state(nstate,nx,x,Vs,U_s)
    
    ting =0.0_rk 
    it=1 
    poping = sum(abs(wfing)**2)
    time: do while(poping>1e-9)
      it=it+1
      ting=ting+dt 
      output: if (mod(it,10000)==0 .or. (it<10000 .and.  mod(it,100)==0)) then
        poping = sum(abs(wfing)**2) 
        print*,'time:',ting,'pop:',sum(abs(wfing(:,:,1))**2),sum(abs(wfing(:,:,2))**2)
        write(tlabel,'(i20)')it 
        open(66,file='wf/wfra_'//trim(adjustl(tlabel))//'.dat')
        do j=1,na
        do i=1,nx 
        write(66,'(*(f15.7,x))') x(i)*ang_cos(j),x(i)* sqrt(1-ang_cos(j)**2),(abs(wfing(i,j,k))**2/ang_weight(j),k=1,nstate)
        enddo 
        enddo 
        call getflux()
        call LAR(out_WFRE)
    
      end if output 
    
      !subroutine split_s(nstate,n,l,dt,T_rad,T_rot,Vs,VsC,U_rad,U_ang,wfd,wff,U_s)
      call split_s(nstate,nx,na,dt,Tx,TABC,Vs,VsC,U_rad,U_ang,wfing,WFtemp,u_s)
      !T2E(nstate,n,l,m,ting,Cin,Cout,wf,wfre)
      call T2E(nstate,nx,na,ne,ting,Cin,Cout,wfing,out_WFRE)
      call T2E(nstate,nx,na,ne,ting,Cin,Cout,WFtemp,WFRE)
    
    enddo time
    
    
    end subroutine 
    
    
    subroutine getflux()
      complex(rk),allocatable:: FF(:,:)
      real(rk),allocatable:: flux(:,:),flux_ang(:,:,:)
      integer::i,j,k,ie,ii 
      allocate(ff(nx,nx),source=(0.0_rk,0.0_rk))
      call Flux_Opertaor_sin(nx,xmin,len_x,flux_x0,massBC,FF)
      !    do j=1,nx
      !    do i=1,nx 
      !    FF(i,j) = (2.0_rk/massABC/len_x**2 &
      !    *(pi*j*sin(i*pi*((xmax-absorp_len)-xmin)/len_x)*cos(j*pi*((xmax-absorp_len)-xmin)/len_x)))
      !    end do
      !    end do
    
      allocate(flux(nstate,ne),source=0.0_rk) 
      
      do ii=1,nstate
      !$omp parallel do private(ie,j,i,k) schedule(dynamic)
      do ie=1,ne
      do k=1,na 
      do j=1,nx 
      do i=1,nx 
      Flux(ii,ie) = Flux(ii,ie) + (dconjg(WFRE(i,k,ii,ie)) * ff(i,j) * WFRE(j,k,ii,ie))
      enddo 
      enddo 
      enddo 
      enddo 
      !$omp end parallel do 
      enddo 
    
      open(66,file=trim(adjustl(dir))//'Flux.dat')
      do ie=1,ne 
      write(66,'(*(e15.7,x))') E(ie),(flux(ii,ie),ii=1,nstate)
      enddo 
      close(66)
    
      deallocate(flux)
    
      allocate(flux_ang(na,nstate,ne),source=0.0_rk)
      do ie=1,ne
      do ii=1,nstate 
        do k=1,na 
          do j=1,nx  
            do i=1,nx 
        flux_ang(k,ii,ie) = flux_ang(k,ii,ie) + (dconjg(WFRE(i,k,ii,ie)) * ff(i,j) * WFRE(j,k,ii,ie))
            enddo 
          enddo 
        enddo 
      enddo 
      enddo
      deallocate(ff)
      
        open(66,file=trim(adjustl(dir))//'Flux_ang.dat')
        ! write(66,'(*(e15.7,x))')0.0_rk, (E(ie),ie=1,ne)
        ! do k=1,na
        ! write(66,'(*(e15.7,x))') ang_cos(k),(flux_ang(k,1,ie),ie=1,ne)
        ! enddo
        do j=1,na 
            do ie =1,ne 
            write(66,'(*(e15.7,x))') e(ie)*ang_cos(j),E(ie)*sqrt(1-ang_cos(j)**2),(flux_ang(j,ii,ie)/ang_weight(j),ii=1,nstate)
            enddo 
        enddo
        close(66) 
        
    
      deallocate(flux_ang) 
    
      end subroutine 
    
    end module scatter1_sudden
  