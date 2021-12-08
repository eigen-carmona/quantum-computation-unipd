program time_dependent
    ! Compile as ```gfortran time_dependent.f90 -o time_dependent finite_difference.o split_op.o -llapack -Wall -lfftw3 -frecursive```
    ! After compiling ```gfortran -c finite_difference.f90 -llapack``` and ```gfortran -c split_op.f90```
    use finite_diff
    use split_op
    
    integer :: NN, ii, n_eigen, timesteps, INFO, N_grid
    real*8, parameter :: L = 10, pi = acos(-1.0)
    real*8 :: dx, dp, dt, norm_sq, T
    complex*16, allocatable, dimension(:,:) :: Psi
    real*8, allocatable, dimension(:,:) :: psi_states, int_V_t
    real*8, allocatable, dimension(:) :: x_grid, p_grid, energies, V_0
    character*42 :: states_file
    

    write(*,*) "Enter the target state to evolve, the timesteps and the constant T"
    read(*,*) n_eigen, timesteps, T
    dt = T/timesteps

    n_eigen = n_eigen + 1! In our setting, the ground state is shifted to n = 1

    NN = 2**12 ! Just because...

    if (NN.le.n_eigen) stop "The grid size must be higher than the number of energy levels"

    ! The grid size determines the increment value
    dx = 2*L/NN
    dp = pi/L
    
    ! Since we assume x_0=x_n=0, we reduce the problem to a (NN-1)x(NN-1) matrix
    N_grid = NN - 1
    
    ! Compute the potential energy array
    allocate(x_grid(1:N_grid))
    allocate(p_grid(1:N_grid))
    x_grid = (/(ii*dx - L, ii = 1,N_grid)/)
    p_grid(1:int(N_grid/2)) = (/(ii*dp, ii = 1,int(N_grid/2))/)
    p_grid(int(N_grid/2)+1:N_grid) = (/(-int(N_grid/2)*dp+ii*dp, ii = 0,int(N_grid/2))/)

    ! Prepare t=0 potential
    allocate(V_0(1:N_grid))
    V_0 = x_grid**2/2
    ! Prepare time-evolving integrated potential
    allocate(int_V_t(1:N_grid,1:timesteps))
    do ii = 1, timesteps
        int_V_t(:,ii) = ((x_grid*T-(ii-1)*dt)**3-(x_grid*T-ii*dt)**3)/(6*T**2)
    end do

    ! Calling single particle routine for diagonal potentials and vanishing boundaries
    call finite_diff_1D(-L,L,NN,n_eigen,energies,psi_states,V_0,INFO)
    ! Validate there are no routine errors
    if (INFO.ne.0) then
        write(*,'(A,I4)') "Unsuccessful calculation. Status = ",INFO
    stop "Not storing results. Exiting program..."
    end if

    allocate(Psi(1:N_grid,0:timesteps))
    ! Use the eigenstate of interest
    norm_sq = sum(abs(psi_states(:,n_eigen))**2)*dx
    Psi(:,0) = psi_states(:,n_eigen)/sqrt(norm_sq)
    call split_op_1D(Psi,int_V_t,dt,N_grid,p_grid,timesteps,dx)
    ! evolve the state
    write(states_file,'(A,I0,A,I0,A,I0,A)') "data/state_",n_eigen-1,"_evo_T",int(T),"_",timesteps,".dat"
    open(12, file = states_file)
    do ii = 0, timesteps
        write(12,'(*(F0.8,SP,F0.8,"j",","))') Psi(:,ii)
    end do
    close(12)
end program time_dependent
