module amod_propatation 
  implicit none
  integer,parameter:: rk = selected_real_kind(15) 
  real(rk),parameter:: pi = 4.0_rk*atan(1.0_rk) 
  complex(rk),parameter:: iunit = (0.0_rk,1.0_rk) 
  private
  public:: split, T2E ,split_s, U_state 
  contains 
  subroutine split(nstate,n,l,dt,T_rad,T_rot,Vs,VsC,U_rad,U_ang,wfd,wff)
    integer,intent(in):: n,l,nstate
    complex(rk),intent(in):: T_rad(n),T_rot(n,l),VsC(n,nstate)
    real(rk),intent(in):: Vs(n,nstate),U_rad(n,n),U_ang(l,l),dt
    complex(rk),intent(inout):: wfd(n,l,nstate),wff(n,l,nstate)

    complex(rk),allocatable:: WFtemp(:,:,:)
    integer:: i,j,k

    complex(kind=8), parameter :: one = (1.0d0, 0.0d0), zero = (0.0d0, 0.0d0)
    complex(rk),allocatable:: wfni(:),wfnf(:),wfli(:),wflf(:) 
    allocate(WFtemp(n,l,nstate)) 
    allocate(wfni(n),wfnf(n),wfli(l),wflf(l))

    

    !$omp parallel do private(i,j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l 
    do i=1,n
      WFtemp(i,j,k)=wfd(i,j,k)*exp(-iunit*0.5_rk*dt*Vs(i,k))*Vsc(i,k)
    end do
    end do
    end do 
    !$omp end parallel do 



    !print*,'ang:DVR->FBR',sum(abs(wfd)**2)
    !$omp parallel do private(i,k) schedule(dynamic)
    do k=1,nstate 
    do i=1,n
    wfd(i,:,k) =  matmul(U_ang,WFtemp(i,:,k))
    enddo
    enddo
    !$omp end parallel do 

    !print*,'mul T_rot',sum(abs(wfd)**2)
    !$omp parallel do private(i,j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l
    do i=1,n
    WFtemp(i,j,k) = T_rot(i,j) * wfd(i,j,k)
    enddo 
    enddo 
    enddo 
    !$omp end parallel do 

    !print*,'rad:DVR->FBR',sum(abs(wfd)**2)
    !$omp parallel do private(j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l
    wfd(:,j,k) = matmul(U_rad,WFtemp(:,j,k))
    enddo
    enddo 
    !$omp end parallel do

    !print*,'mul T_rad',sum(abs(wfd)**2)
    !$omp parallel do private(i,j,k) schedule(dynamic)
    do k=1,nstate
    do j=1,l 
    do i=1,n
    WFtemp(i,j,k) = T_rad(i) * wfd(i,j,k)
    enddo
    enddo 
    enddo 
    !$omp end parallel do
    
   
    !print*,'ang:FBR->DVR',sum(abs(wfd)**2)
    !$omp parallel do private(i,k,wfli,wflf) schedule(dynamic)
    do k=1,nstate 
    do i=1,n
    wfd(i,:,k) = matmul(transpose(U_ang),WFtemp(i,:,k))
    enddo 
    enddo 
    !$omp end parallel do 
  
    wff = wfd
    !print*,'rad:FBR->DVR',sum(abs(wfd)**2)
    !$omp parallel do private(j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l
    WFtemp(:,j,k) = matmul(U_rad,wfd(:,j,k))
      ! call zgemv('N', n, n, one, U_rad, n, wfd(:, j, k), 1, zero, WFtemp(:, j, k), 1)
    end do
    end do 
    !$omp end parallel do 



    ! print*,'mul V/2',sum(abs(wfd)**2)
    !$omp parallel do private(i,j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l
    do i=1,n
    wfd(i,j,k) =WFtemp(i,j,k) *exp(-iunit*0.5_rk*dt*Vs(i,k))*Vsc(i,k)
    enddo 
    enddo 
    enddo 
    !$omp end parallel do


   deallocate(WFtemp) 
  end subroutine split
 
  subroutine U_state(nstate,n,x,v,u_s)
    use Math_Functions,only:Gamma
    use math,only:rs 
    integer,intent(in):: nstate,n 
    real(rk),intent(in):: x(n)
    real(rk),intent(out):: u_s(nstate,nstate,n)
    real(rk),intent(inout):: V(n,nstate)
    real(rk),allocatable:: H(:,:),fv1(:),fv2(:),eigenv(:),eigenf(:,:)
    integer:: i, ierr 
    allocate(H(nstate,nstate)) 
    allocate(fv1(nstate),fv2(nstate))
    allocate(eigenv(nstate),eigenf(nstate,nstate))
    do i=1,n 
      H =0.0_rk 
      H(1,1) = V(i,1)
      H(2,2) = V(i,2)
      H(1,2) = sqrt(Gamma(x(i)))
      H(2,1) = H(1,2)
      fv1 = 0.0_rk 
      fv2 = 0.0_rk 
      eigenf=0.0_rk
      eigenv=0.0_rk 
      call rs(nstate,nstate,H,eigenv,1,eigenf,fv1,fv2,ierr)
      u_s(:,:,i) = eigenf
      V(i,:) = eigenv 
    end do

    deallocate(H,fv1,fv2)
  end subroutine U_state 
    


  subroutine split_s(nstate,n,l,dt,T_rad,T_rot,Vs,VsC,U_rad,U_ang,wfd,wff,U_s)
    integer,intent(in):: n,l,nstate
    complex(rk),intent(in):: T_rad(n),T_rot(n,l),VsC(n,nstate)
    real(rk),intent(in):: Vs(n,nstate),U_rad(n,n),U_ang(l,l),dt,u_s(nstate,nstate,n)
    
    complex(rk),intent(inout):: wfd(n,l,nstate),wff(n,l,nstate)
    complex(rk),allocatable:: WFtemp(:,:,:)
    integer:: i,j,k

    complex(kind=8), parameter :: one = (1.0d0, 0.0d0), zero = (0.0d0, 0.0d0)
    complex(rk),allocatable:: wfni(:),wfnf(:),wfli(:),wflf(:) 
    allocate(WFtemp(n,l,nstate)) 
    allocate(wfni(n),wfnf(n),wfli(l),wflf(l))
    
    !print*,'mul V/2',sum(abs(wfd)**2)

    
    !$omp parallel do private(i,j) schedule(dynamic)
    do j=1,l 
      do i=1,n 
        WFtemp(i,j,:) =  matmul(u_s(:,:,i),wfd(i,j,:))
      end do
    end do
    !$omp end parallel do


    !$omp parallel do private(i,j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l 
    do i=1,n
    wfd(i,j,k)=WFtemp(i,j,k)*exp(-iunit*0.5_rk*dt*Vs(i,k))*Vsc(i,k)
    end do
    end do
    end do 
    !$omp end parallel do 

    !$omp parallel do private(i,j) schedule(dynamic)
    do j=1,l 
      do i=1,n 
        WFtemp(i,j,:) =  matmul(transpose(u_s(:,:,i)),wfd(i,j,:))
      end do
    end do
    !$omp end parallel do

    !print*,'ang:DVR->FBR',sum(abs(wfd)**2)
    !$omp parallel do private(i,k) schedule(dynamic)
    do k=1,nstate 
    do i=1,n
    wfd(i,:,k) =  matmul(U_ang,WFtemp(i,:,k))
    enddo
    enddo
    !$omp end parallel do 

    !print*,'mul T_rot',sum(abs(wfd)**2)
    !$omp parallel do private(i,j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l
    do i=1,n
    WFtemp(i,j,k) = T_rot(i,j) * wfd(i,j,k)
    enddo 
    enddo 
    enddo 
    !$omp end parallel do 

    !print*,'rad:DVR->FBR',sum(abs(wfd)**2)
    !$omp parallel do private(j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l
    wfd(:,j,k) = matmul(U_rad,WFtemp(:,j,k))
    enddo
    enddo 
    !$omp end parallel do

    !print*,'mul T_rad',sum(abs(wfd)**2)
    !$omp parallel do private(i,j,k) schedule(dynamic)
    do k=1,nstate
    do j=1,l 
    do i=1,n
    WFtemp(i,j,k) = T_rad(i) * wfd(i,j,k)
    enddo
    enddo 
    enddo 
    !$omp end parallel do
    
   
    !print*,'ang:FBR->DVR',sum(abs(wfd)**2)
    !$omp parallel do private(i,k,wfli,wflf) schedule(dynamic)
    do k=1,nstate 
    do i=1,n
    wfd(i,:,k) = matmul(transpose(U_ang),WFtemp(i,:,k))
    enddo 
    enddo 
    !$omp end parallel do 
  
    wff = wfd 
    !print*,'rad:FBR->DVR',sum(abs(wfd)**2)
    !$omp parallel do private(j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l
    WFtemp(:,j,k) = matmul(U_rad,wfd(:,j,k))
    end do
    end do 
    !$omp end parallel do 

    ! print*,'mul V/2',sum(abs(wfd)**2)
    !$omp parallel do private(i,j) schedule(dynamic)
    do j=1,l 
      do i=1,n 
        wfd(i,j,:) =  matmul(u_s(:,:,i),WFtemp(i,j,:))
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(i,j,k) schedule(dynamic)
    do k=1,nstate 
    do j=1,l 
    do i=1,n
      WFtemp(i,j,k)=wfd(i,j,k)*exp(-iunit*0.5_rk*dt*Vs(i,k))*Vsc(i,k)
    end do
    end do
    end do 
    !$omp end parallel do 
    !$omp parallel do private(i,j) schedule(dynamic)
    do j=1,l 
      do i=1,n 
        wfd(i,j,:) =  matmul(transpose(u_s(:,:,i)),WFtemp(i,j,:))
      end do
    end do
    !$omp end parallel do

   deallocate(WFtemp) 

  end subroutine split_s 

  subroutine T2E(nstate,n,l,m,ting,Cin,Cout,wf,wfre)
    integer,intent(in):: nstate,n,l,m
    real(rk),intent(in):: ting
    complex(rk),intent(in):: wf(n,l,nstate),Cin(nstate,m),Cout(nstate,m)
    complex(rk),intent(out):: wfre(n,l,nstate,m)
    integer:: i,j,k,iE 
    complex(rk):: CC(nstate,m)
    
    do iE = 1,m 
        do k=1,nstate 
            CC(k,ie) = exp((Cin(k,ie))*ting)*Cout(k,ie)
            !! $omp parallel do private(i,j) schedule(dynamic)
            do j=1,l 
                do i=1,n
                    WFRE(i,j,k,iE) = WFRE(i,j,k,iE) + WF(i,j,k)*CC(k,ie)
                end do 
            end do 
            !!$omp end parallel do
        end do 
    end do 
  
  end subroutine T2E

end module amod_propatation