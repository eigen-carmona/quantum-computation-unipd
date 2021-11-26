program phasing_oscillator
    implicit none
    integer :: NN, jj, ii, n_eigen
    real*8, parameter :: x_0 = 0, x_n = 1
    real*8 :: dx, dt = 0.0005, x_i = x_0
    real*8, allocatable :: T_exop(:,:), V_exop(:,:), Psi(:,:), x_grid(:), p_grid(:), psi_0(:)

    ! DSTEIN-Specific variables
    integer :: M, LDZ, INFO
    integer, allocatable :: IBLOCK(:), ISPLIT(:), IWORK(:), IFAIL(:)
    real*8, allocatable :: D(:), E(:), W(:), Z(:,:), WORK(:)
    
    ! TODO: shouldn't this be the interval a,b instead?
    write(*,*) "How many points should be calculated?"
    read(*,*) NN
    
    !write(*,*) "How many energy values should be obtained?"
    !read(*,*) n_eigen
    !
    !! NN-1 should be at least n_eigen
    !if (NN-1.lt.n_eigen) stop "The grid size should be bigger than the desired number of energy levels."
    
    ! The number of points fixes the dx
    dx = (x_n-x_0)/NN

    ! As soon as dx is defined, we may build x_grid
    allocate(x_grid(1:NN-1))
    do ii = 1, NN-1
        x_grid(ii) = ii*dx
    end do

    ! We compute the momentum, that is, the fourier transform
    allocate(p_grid(1:NN-1))
    p_grid = fft(x_grid)

    ! TODO: validation for memory in allocation
    allocate(T_exop(1:NN-1,1:NN-1))
    ! Isn't the potential time-dependent? how are we accounting for this at each step?
    allocate(V_exop(1:NN-1,1:NN-1))

    ! Populate the operators
    do ii = 1, NN-1
        T_exop(ii,ii) = exp(dcmplx(0.d0,-dt)*p_grid(ii)**2/(2*m))
        V_exop(ii,ii) = exp(dcmplx(0.d0,-dt/2)*x_grid(ii)**2/2)
    end do

    ! evolve the state
    allocate(Psi(1:NN-1,0:NN-1))
    Psi(:,0) = psi_0
    do ii = 1, NN-1
        ! Apply the V operator to the last step configuration
        ! Take to momentum space
        Psi(:,ii)= fft(matmul(V_exop,Psi(:,ii-1)))
        ! Apply the T operator to the transformed state
        ! Go back to coordinate space
        Psi(:,ii) = invfft(matmul(T_exop,Psi(:,ii)))
        ! Apply the V operator once more
        Psi(:,ii) = matmul(V_exop,Psi(:,ii))
    end do

    end program phasing_oscillator
    