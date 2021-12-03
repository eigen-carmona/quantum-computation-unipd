module FFTW3
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
end module

program phasing_oscillator
    ! Compile as ```gfortran time_dependent.f90 -o time_dependent finite_difference.o -llapack```
    ! After compiling ```gfortran -c finite_difference.f90 -llapack```
    use finite_diff
    
    use FFTW3
    implicit none
    
    interface
        function potential_V(N_grid,x_pos,election)
        implicit none
        integer     :: N_grid
        real*8      :: x_pos(1:N_grid)
        character*2 :: election
        real*8      :: potential_V(1:N_grid)
        end function
    end interface
    
    integer :: NN, ii, n_eigen, timesteps, INFO, N_grid
    integer*8 :: plan_psi_x, plan_psi_k
    real*8, parameter :: L = 10, pi = acos(-1.0)
    real*8 :: dx, dp, dt = 0.05, norm_sq
    complex*16, allocatable, dimension(:,:) :: Psi
    real*8, allocatable, dimension(:,:) :: psi_states
    real*8, allocatable, dimension(:) :: x_grid, p_grid, energies, V_pot
    complex*16, allocatable, dimension(:) :: T_exop, V_exop, psi_x_0, psi_k_0, psi_x_1, psi_k_1, psi_x_2
    character*2 :: election
    
    write(*,*) "Currently enabled potential choices:"
    write(*,'(A,/,A,/,A,/,A)') " - 'ho': harmonic oscillator"," - 'qs': V=5 step at x = 0",&
    " - 'qb': V = 5 barrier when x in (-5,5)", " - 'qw': V = -5 well when x in (-5,5)"
    write(*,*) "Any option not listed leads to the (-10,10) quantum box."
    write(*,*) "Enter the potential type, grid size and the number of desired energy levels"
    read(*,*) election, NN, n_eigen

    write(*,*) "Evolve for how many timesteps?"
    read(*,*) timesteps
    
    if (NN.le.n_eigen) stop "The grid size must be higher than the number of energy levels"
    
    ! The grid size determines the increment value
    dx = 2*L/NN
    dp = pi/L
    
    if (dx.ge.1) stop "Grid size too small. Minumum should be 20 (advised 400 and above)."
    
    ! Since we assume x_0=x_n=0, we reduce the problem to a (NN-1)x(NN-1) matrix
    N_grid = NN - 1
    
    ! Compute the potential energy array
    allocate(V_pot(1:N_grid))
    allocate(x_grid(1:N_grid))
    allocate(p_grid(1:N_grid))
    x_grid = (/(ii*dx - L, ii = 1,N_grid)/)
    p_grid(1:int(N_grid/2)) = (/(ii*dp, ii = 1,int(N_grid/2))/)
    p_grid(int(N_grid/2)+1:N_grid) = (/(-int(N_grid/2)*dp+ii*dp, ii = 0,int(N_grid/2))/)

    V_pot = potential_V(N_grid,x_grid,election)
    
    ! Calling single particle routine for diagonal potentials and vanishing boundaries
    call finite_diff_1D(-L,L,NN,n_eigen,energies,psi_states,V_pot,INFO)
    
    ! Validate there are no routine errors
    if (INFO.ne.0) then
        write(*,'(A,I4)') "Unsuccessful calculation. Status = ",INFO
    stop "Not storing results. Exiting program..."
    end if
    
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

    allocate(Psi(1:N_grid,0:timesteps))
    ! Use the eigenstate of interest
    norm_sq = sum(abs(psi_states(:,n_eigen))**2)*dx
    Psi(:,0) = psi_states(:,n_eigen)/sqrt(norm_sq)
    ! evolve the state
    open(12, file = 'states_ev.dat')
    write(12,'(*(F0.8,SP,F0.8,"j",","))') Psi(:,0)
    do ii = 1, timesteps
        V_exop = exp(dcmplx(0.d0,-dt/2)*dcmplx((x_grid-ii*dt/5)**2/2,0))
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
        write(12,'(*(F0.8,SP,F0.8,"j",","))') Psi(:,ii)
    end do
    close(12)
    end program phasing_oscillator



! Catallogue of potentials
function potential_V(N_grid,x_pos,election)
implicit none
integer     :: N_grid, ii
real*8      :: x_pos(1:N_grid)
character*2 :: election
real*8      :: potential_V(1:N_grid)
    select case (election)
    case ('ho')
        ! Harmonic oscillator potential
        potential_V = x_pos**2/2
    case ('qs')
        ! Step potential
        do ii = 1, N_grid
            if (x_pos(ii).lt.0) then
                potential_V(ii) = 0
            else
                potential_V(ii) = 1
            end if
        end do
    case ('qb','qw')
        ! Barrier potential
        do ii = 1, N_grid
            if ((x_pos(ii).lt.-5).or.(x_pos(ii).gt.5)) then
                potential_V(ii) = 0
            else
                potential_V(ii) = 1
            end if
        end do
        ! Potential well
        if (election.eq.'we') then
            potential_V = -potential_V
        end if
    case default
        potential_V = (/(0,ii = 1,N_grid)/)
        election = 'bx'
    end select
end function