module amod_initial_wf
    implicit none
    integer, parameter :: rk = selected_real_kind(15)
    real(rk), parameter :: pi = 4.0_rk * atan(1.0_rk)
    complex(rk), parameter :: eye = (0.0_rk, 1.0_rk)

contains
    subroutine initial_wf(mass,Es, E0, sigmaE, k0, a_E, xs, psi)
        real(rk), intent(in) ::mass, Es(:), E0, sigmaE, k0, xs(:)
        complex(rk), intent(out) :: a_E(:), psi(:)
        integer :: i, j, m, n
        real(rk) :: dx, dE, norm_factor

        real(rk),allocatable:: Ks(:)

        
        m = size(Es)
        n = size(xs)

        if (size(a_E) /= m .or. size(psi) /= n) stop "Array size mismatch"
        
        allocate(Ks(size(Es))) 

        Ks = sqrt(2.0_rk * mass * Es)

       
        do j = 1, m
            a_E(j) = (1.0_rk / (2.0_rk * pi * sigmaE**2))**0.25_rk * &
                     exp(-(Es(j) - E0)**2 / (4.0_rk * sigmaE**2)) * &
                     exp(eye * k0 * (Es(j) - E0))
        end do

        
        dE = sqrt(( Es(2) - Es(1))*2*mass) 

        do i = 1, n
            psi(i) = (0.0_rk, 0.0_rk)
            do j = 1, m
                psi(i) = psi(i) + a_E(j) * exp(-eye * xs(i) * Ks(j)) * dE
            end do
        end do

       
        dx = xs(2) - xs(1)
        norm_factor = sqrt(sum(abs(psi)**2))
        if (norm_factor <= tiny(1.0_rk)) stop "Norm factor too small"
        psi = psi / norm_factor

       
        do j = 1, m
            a_E(j) = (0.0_rk, 0.0_rk)
            do i = 1, n
                a_E(j) = a_E(j) + psi(i) * exp(eye * xs(i) * Ks(j)) 
            end do
        end do

        
    end subroutine initial_wf
end module amod_initial_wf
