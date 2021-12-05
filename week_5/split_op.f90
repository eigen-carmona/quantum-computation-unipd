module FFTW3
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
end module

module split_op
! Compile as ```gfortran -c split_op.f90```
use FFTW3
contains
    subroutine split_op_1D(Psi,V_t,dt,N_grid,p_grid,timesteps,dx)
    implicit none
    integer, intent(in) :: timesteps, N_grid
    integer :: ii
    integer*8 :: plan_psi_x, plan_psi_k
    real*8 , intent(in):: dx, dt
    real*8 :: norm_sq
    complex*16, dimension(:,0:), intent(out) :: Psi! Evolution of target state. Assumed to have t=0 state at Psi(:,0)
    real*8, dimension(:), intent(in) :: p_grid! momentum grid values
    real*8, dimension(:,0:), intent(in) :: V_t! potential evaluated in at least t=timesteps times.
    complex*16, allocatable, dimension(:) :: T_exop, V_exop, psi_x_0, psi_k_0, psi_x_1, psi_k_1, psi_x_2

        ! TODO: validation for memory in allocation
        allocate(T_exop(1:N_grid))
        ! Isn't the potential time-dependent? how are we accounting for this at each step?
        allocate(V_exop(1:N_grid))

        ! Populate the operators
        T_exop = exp(dcmplx(0.d0,-dt)*dcmplx(p_grid**2/2,0))
        ! start fft plans
        allocate(psi_x_0(1:N_grid))
        allocate(psi_k_0(1:N_grid))
        allocate(psi_x_1(1:N_grid))
        allocate(psi_k_1(1:N_grid))
        allocate(psi_x_2(1:N_grid))
        call dfftw_plan_dft_1d(plan_psi_x,N_grid,psi_x_0,psi_k_0,FFTW_FORWARD,FFTW_MEASURE)
        call dfftw_plan_dft_1d(plan_psi_k,N_grid,psi_k_1,psi_x_1,FFTW_BACKWARD,FFTW_MEASURE)
        do ii = 1, timesteps
            V_exop = exp(dcmplx(0.d0,-dt/2)*dcmplx((V_t(:,ii)),0))
            ! Apply the V operator to the last step configuration
            psi_x_0 = V_exop*Psi(:,ii-1)
            ! Take to momentum space
            call dfftw_execute_dft(plan_psi_x,psi_x_0,psi_k_0)
            ! Apply the T operator to the momentum space half step
            psi_k_1 = T_exop*psi_k_0
            ! Take back to position space
            call dfftw_execute_dft(plan_psi_k,psi_k_1,psi_x_1)
            psi_x_2 = V_exop*psi_x_1
            ! Normalize for stability
            norm_sq = sum(abs(psi_x_2)**2)*dx
            ! Apply the V operator once more
            Psi(:,ii) = psi_x_2/sqrt(norm_sq)
        end do
    end subroutine
end module