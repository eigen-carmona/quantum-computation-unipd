module FFTW3
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
end module

module split_op
! Compile as ```gfortran -c split_op.f90```
use FFTW3
contains
    subroutine split_op_1D(Psi,int_V_t,dt,N_grid,p_grid,timesteps,dx)
    implicit none
    integer, intent(in) :: timesteps, N_grid
    integer :: ii
    integer*8 :: plan_psi_0, plan_psi_1, plan_psi_2, plan_psi_3
    real*8 , intent(in):: dx, dt
    real*8 :: norm_sq
    complex*16, dimension(:,0:), intent(out) :: Psi! Evolution of target state. Assumed to have t=0 state at Psi(:,0)
    real*8, dimension(:), intent(in) :: p_grid! momentum grid values
    real*8, dimension(:,:), intent(in) :: int_V_t! integrated potential evaluated at each interval timestep.
    complex*16, allocatable, dimension(:) :: T_exop, V_exop, psi_x_0, psi_k_0, psi_x_1, psi_k_1, psi_x_2, psi_k_2, psi_k_3, psi_x_3

        ! TODO: validation for memory in allocation
        allocate(T_exop(1:N_grid))
        ! Isn't the potential time-dependent? how are we accounting for this at each step?
        allocate(V_exop(1:N_grid))

        ! Populate the operators
        T_exop = exp(dcmplx(0.d0,-(dt/2)*(p_grid**2/2)))
        ! start fft plans
        allocate(psi_x_0(1:N_grid))
        allocate(psi_k_0(1:N_grid))
        allocate(psi_x_1(1:N_grid))
        allocate(psi_k_1(1:N_grid))
        allocate(psi_x_2(1:N_grid))
        allocate(psi_k_2(1:N_grid))
        allocate(psi_x_3(1:N_grid))
        allocate(psi_k_3(1:N_grid))
        call dfftw_plan_dft_1d(plan_psi_0,N_grid,psi_x_0,psi_k_0,FFTW_FORWARD,FFTW_MEASURE)
        call dfftw_plan_dft_1d(plan_psi_1,N_grid,psi_k_1,psi_x_1,FFTW_BACKWARD,FFTW_MEASURE)
        call dfftw_plan_dft_1d(plan_psi_2,N_grid,psi_x_2,psi_k_2,FFTW_FORWARD,FFTW_MEASURE)
        call dfftw_plan_dft_1d(plan_psi_3,N_grid,psi_k_3,psi_x_3,FFTW_BACKWARD,FFTW_MEASURE)

        do ii = 1, timesteps
            psi_x_0 = Psi(:,ii-1)
            ! First take the state to momentum space
            call dfftw_execute_dft(plan_psi_0,psi_x_0,psi_k_0)
            ! Then apply the kinetic operator half step
            psi_k_1 = T_exop*psi_k_0
            ! Take back to position space
            call dfftw_execute_dft(plan_psi_1,psi_k_1,psi_x_1)
            ! Then apply the integrated potential step
            V_exop = exp(dcmplx(0.d0,-int_V_t(:,ii)))
            psi_x_2 = V_exop*psi_x_1
            ! Take to momentum space
            call dfftw_execute_dft(plan_psi_2,psi_x_2,psi_k_2)
            ! Apply once more the kinetic operator half step
            psi_k_3 = T_exop*psi_k_2
            ! Take back to position space
            call dfftw_execute_dft(plan_psi_3,psi_k_3,psi_x_3)
            ! Normalize for stability
            norm_sq = sum(abs(psi_x_3)**2)*dx
            ! Apply the V operator once more
            Psi(:,ii) = psi_x_3/sqrt(norm_sq)
        end do
    end subroutine
end module