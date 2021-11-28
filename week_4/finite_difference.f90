module finite_diff
! Compile as ```gfortran -c finite_difference.f90 -llapack```
contains
    subroutine finite_diff_1D(x_0,x_n,dx,M,W,Z,V,INFO)
    ! Computes the first M energies and eigenstates of a 1-D system under potential V
    ! Arguments:
    ! x_0 - Lower boundary of x grid. Psi is assumed to vanish here
    ! x_n - Upper boundary of x grid. Psi is assumed to vanish here, too
    ! dx - increments, defines the coarseness of the grid
    ! M - Number of states to be obtained
    ! W - Energies of the states
    ! Z - Array of column vectors describing the states
    ! V - Array with the potential energy evaluated from x_0+dx to x_n-dx
    implicit none
    real*8, intent(in) :: x_0, x_n, dx, V(:)
    integer :: M, LDZ, LWORK, LIWORK,INFO, NN
    integer, allocatable :: ISUPPZ(:), ISPLIT(:), IWORK(:), IFAIL(:)
    real*8, allocatable :: D(:), E(:), WORK(:)
    real*8, intent(out), allocatable :: W(:), Z(:,:)
    logical :: TRYRAC = .True.

        ! We'll tackle the problem by means of finite difference, with a second derivative:
        ! \frac{d^2 f_{j}}{dx^2} = \frac{f_{j+1} - 2f_{j} + f_{j-1}}{h^2}.
        ! A key assumption in the boundary conditions is that \psi_{0} = \psi_{N+1} = 0,
        ! this immediately leads to a convenient tridiagonal matrix

        NN = int((x_n-x_0)/dx)-1 ! we start at x_0 + dx and end at x_N - dx

        ! Since we know eigenvectors will also be computed,
        ! Optimal WORK and IWORK dimensions can be provided a priori
        LWORK = 18*NN
        LIWORK = 10*NN

        ! DSTEMR-specific allocations
        allocate(D(1:NN))
        allocate(E(1:NN-1))
        allocate(W(1:NN))
        allocate(ISUPPZ(1:2*NN))
        allocate(ISPLIT(1:NN))
        LDZ = NN
        allocate(Z(1:LDZ,1:M))
        allocate(WORK(1:LWORK))
        allocate(IWORK(1:LIWORK))
        allocate(IFAIL(1:M))

        ! Build the hamiltonian matrix
        ! Diagonal elements
        D = 2/dx**2/2 + V
        ! Subdiagonal elements
        E = E  - 1/dx**2/2        

        ! USE LAPACK DSTEMR
        call dstemr(&
            'V',&! We wish for both eigenvalues and eigenvectors
            'I',&! We aim to obtain the first M eigenvalues
            NN,&! This is the order of the matrix
            D,&! The diagonal elements of the matrix. Will be overwritten
            E,&! Sub-diagonal elements of the matrix. Also to be overwritten
            0.d0,&! Irrelevant, since we're obtaining eigenvalues by index
            0.d0,&! Idem
            1,&! We want to start at the very first eigenvalue...
            M,&! and go up to the Mth (M<=NN) eigenvalue
            M,&! The total number of eigenvalues to be found
            W,&! The first M elements will be the first M eigenvalues
            Z,&! Its columns will store the first M eigenvectors
            LDZ,&! The leading dimension of Z. In general, NN
            M,&! M columns in Z are required to store the eigenvectors
            ISUPPZ,&! Support of the eigenvectors
            TRYRAC,&! Validate if the matrix defines its e-values w/high relative accuracy
            WORK,&! Work vector for dstmr
            LWORK,&! Work array dimension
            IWORK,&
            LIWORK,&! Dimension of IWORK array
            INFO)! Status info
        
    end subroutine
end module