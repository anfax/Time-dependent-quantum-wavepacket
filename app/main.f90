program main
  ! use scatter2, only: TDQW2D
  use scatter1_sudden,only: TDQW2D 
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  
  open(109,file='ALog.md')
  call TDQW2D()
  close(109)
end program main
