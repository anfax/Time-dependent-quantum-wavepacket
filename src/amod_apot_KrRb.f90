!======================================================================
module amod_pot_KrRb 
    implicit none
    private 
    public :: KrStar_Rb,Kr_RbPlus  
    contains 
  
  !!!======================================================================== 
  !!!  Kr* + Rb PES
  subroutine KrStar_Rb(r,v)
      real*8, parameter :: au2cm=219474.d0
      real*8, parameter :: au2ev=27.2114d0
      real*8, parameter :: pi=dacos(-1.d0)
      real*8, parameter :: eps=1.d-8
      real(8),intent(in) :: r
      real(8),intent(out) :: v 
      real*8 :: dv(1,1),v0(1),r0(1,1)
      r0(1,1)=r 
      call potSigma(1,1,r0,V0,dv)
      v =  v0(1) 
  end subroutine KrStar_Rb 
  
  subroutine potSigma(ndim,ntot,r0,vx,dv)
      implicit none
      integer,intent(in) :: ndim,ntot
      real*8,intent(in) :: r0(ndim,ntot)
      real*8,intent(out) :: vx(ntot),dv(ndim,ntot)
    
      integer,parameter :: s0=1,s1=5,s2=5
      character*99,save :: wfile='./input/sigma.txt'
    
      real*8 :: dvx(ndim),tmp
      integer :: i
      integer,save :: init=0
      integer,save :: fid
      real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3
      real*8,save :: rg(2,s0),vg(2)
      real*8, allocatable :: rx(:,:)
    
      allocate(rx(ndim,ntot))
      rx=exp(-0.2d0*r0)
      if(ndim.ne.s0) stop "ndim .ne. s0"
      if (init.eq.0) then
        fid=123456
        open(fid,file=trim(wfile),status='old',action='read')
        read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
        close(fid)
        init=1
      endif
      if(ntot.ge.24) then
        call nsimx(ntot,rx,vx,1,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
      else
        do i=1,ntot
        call nsim(rx(1,i),vx(i),1,dvx,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
        dv(:,i)=dvx
        enddo
      endif
    
      vx=(vx-0.3396043773d-4)/27.2114d0
    
      !dv
      do i=1,ntot
        dv(:,i)=-dv(:,i)*rx(1,i)*0.2d0
      end do
      dv=dv/27.2114d0
    
      return
    end subroutine
    !======================================================================
    subroutine nsim(r0,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
      implicit none
      integer,intent(in) :: n0,n1,n2,idv
      real*8,intent(in)  :: r0(n0),rg(2,n0),vg(2)
      real*8,intent(in)  :: w1(n0,n1),b1(n1)
      real*8,intent(in)  :: w2(n1,n2),b2(n2)
      real*8,intent(in)  :: w3(n2),b3
      real*8,intent(out) :: v,dv(n0)
      integer :: i,j,k
      real*8  :: r(n0),rgg(n0),vgg,ax(n1),bx(n2)
      real*8  :: dvtm,rtmp,rt1(n1),rt2(n2)
      real*8,external :: ddot
      v=0.d0
      r=r0
      vgg=vg(2)-vg(1)
      ! mapminmax [-1,1]
      do i=1,n0
        rgg(i)=rg(2,i)-rg(1,i)
        r(i)=2.d0*(r(i)-rg(1,i))/rgg(i)-1.d0
      end do
      ! 1st layer
      rt1=b1
      call dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1)
      ax=dtanh(rt1)
      ! 2nd layer
      rt2=b2
      call dgemv('t',n1,n2,1.d0,w2,n1,ax,1,1.d0,rt2,1)
      bx=dtanh(rt2)
      ! output layer
      v=b3+ddot(n2,w3,1,bx,1)
      !reverse map
      v=vgg*(v+1.d0)/2+vg(1)
      if(idv.ne.1) return
      ! calculate first derivatives, dv(i)=dv/dr(i)
      dv=0.d0
      do i=1,n0
        do k=1,n2
          dvtm=0.d0
          do j=1,n1
            dvtm=dvtm+w2(j,k)*w1(i,j)*(1-ax(j)**2)
          enddo
          dv(i)=dv(i)+w3(k)*dvtm*(1-bx(k)**2)
        enddo
        dv(i)=dv(i)*vgg/rgg(i)
      enddo
      return
    end subroutine
    ! ====================================================================
    subroutine nsimx(nt,r,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
      implicit none
      integer,intent(in) :: nt,idv,n0,n1,n2
      real*8,intent(in) :: r(n0,nt)
      real*8,intent(out) :: v(nt),dv(n0,nt)
      real*8,intent(in) :: w1(n0,n1),b1(n1),w2(n1,n2),b2(n2),w3(n2),b3
      real*8,intent(in) :: rg(2,n0),vg(2)
      real*8 :: x0(n0,nt),x1(nt,n1),x2(nt,n2),tmp,rgg(n0),vgg
      integer :: i,j,k,n
      v=0.d0
      do i=1,n0
        rgg(i)=rg(2,i)-rg(1,i)
      enddo
      vgg=vg(2)-vg(1)
      x0=r
      do i=1,n0
        x0(i,1:nt)=(x0(i,1:nt)-rg(1,i))/rgg(i)*2.d0-1.d0
      enddo
      do i=1,n1
       x1(1:nt,i)=b1(i)
      enddo
      call dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
      x1=dtanh(x1)
      do i=1,n2
       x2(1:nt,i)=b2(i)
      enddo
      call dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
      x2=dtanh(x2)
      v(1:nt)=b3
      call dgemv('n',nt,n2,1.d0,x2,nt,w3,1,1.d0,v,1)
      v=(v+1.d0)/2.d0*vgg+vg(1)
      if (idv.ne.1) return
      dv=0.d0
      do n=1,nt
        do i=1,n0
          do k=1,n2
            tmp=0.d0
            do j=1,n1
              tmp=tmp+w2(j,k)*w1(i,j)*(1-x1(n,j)**2)
            enddo
            dv(i,n)=dv(i,n)+w3(k)*tmp*(1-x2(n,k)**2)
          enddo
          dv(i,n)=dv(i,n)*vgg/(rgg(i))
        enddo
      enddo
      return
    end subroutine
    ! ====================================================================
    !!!======================================================================== 
  !!!  Kr + Rb+ PES
  subroutine Kr_RbPlus(r,v)
      real*8, parameter :: bohr2angs=0.5291772109d0
      real*8, parameter :: eshift=487.4129897322d0
      real*8, parameter :: au2cm=219474.6313632d0
      real*8 :: r,v,v0
      call Vtot(r*bohr2angs,v0)
      v = v0/au2cm 
  end subroutine Kr_RbPlus 
  !==============================================================================
  subroutine Vtot(r,v)
      implicit none
      real*8, intent(in) :: r
      real*8, intent(out) :: v
      real*8, parameter :: bohr2angs=0.5291772109d0
      real*8, parameter :: au2cm=219474.6313632d0
      real*8, parameter :: Z=1.055d0
      real*8, parameter :: ad=2.48d0 
      real*8, parameter :: aq=3.97d0
      real*8, parameter :: ao=16.35d0
      real*8, parameter :: B=-6.53d0
      real*8, parameter :: gam=30.2d0
      real*8, parameter :: C6=83.31d0
      real*8, parameter :: C8=1936d0
      real*8, parameter :: a1=1.078d8
      real*8, parameter :: b1=3.412d0
      real*8 :: Vrpl, Vindd, Vindq, V6, V7, Vindo, Vgam, V8
      real*8 :: f4, f6, f7, f8
    
      call  TTdamp(4,b1,r,f4)
      call  TTdamp(6,b1,r,f6)
      call  TTdamp(7,b1,r,f7)
      call  TTdamp(8,b1,r,f8)
    
      Vrpl=a1*exp(-b1*r)
    
      Vindd=-ad*Z**2/(2.d0*r**4)*bohr2angs*au2cm
      Vindq=-aq*Z**2/(2.d0*r**6)*bohr2angs*au2cm
      V6=-C6/(r/bohr2angs)**6*au2cm
      V7=B*Z**3/(2.d0*r**7)*bohr2angs*au2cm
      Vindo=-ao*Z**2/(2.d0*r**8)*bohr2angs*au2cm
      Vgam=-gam*Z**4/(24.d0*r**8)*bohr2angs*au2cm
      V8=-C8/(r/bohr2angs)**8*au2cm
    
      V=Vrpl + f4*Vindd + f6*Vindq + f6*V6 + f7*V7 + f8*Vindo + f8*Vgam + f8*V8
    
      return
    end
    !==============================================================================
    subroutine TTdamp(n,b,r,f)
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: b, r
      real*8, intent(out) :: f
      real*8 :: fac, factorial, t
      integer :: i,j
    
      fac=exp(-b*r)
      f=1.d0
    
      do i=0,n
        factorial=1.d0
        do j=1,i
          factorial=factorial*dble(j)
        end do
        f=f-fac*(b*r)**i/factorial
      end do
    
      return
    end subroutine 
  
    end module amod_pot_KrRb 