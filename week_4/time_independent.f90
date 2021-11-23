program harmonic_oscillator
implicit none
integer*4 :: NN, jj, ii
real*8, parameter :: x_0 = 0, x_n = 1
real*8 :: dx, x_i = x_0
real*8, allocatable :: H_op(:,:)

! We'll tackle the problem by means of finite difference, with a second derivative:
! \frac{d^2 f_{j}}{dx^2} = \frac{f_{j+1} - 2f_{j} + f_{j-1}}{h^2}.
! A key assumption in the boundary conditions is that \psi_{0} = \psi_{N+1} = 0,
! this immediately leads to a convenient tridiagonal matrix

! TODO: shouldn't this be the interval a,b instead?
write(*,*) "How many points should be calculated?"
read(*,*) NN

! The number of points fixes the dx
dx = (x_n-x_0)/NN

! TODO: validation for memory in allocation
allocate(H_op(1:NN-1,1:NN-1))

! Build the hamiltonian matrix
do ii = 1, NN-1
    ! Increase the value of x_i
    x_i = x_i + dx
    ! Now evaluate the diagonal element of the matrix
    ! TODO: generalize potential
    H_op(ii,ii) = 2/dx**2 + x_i**2/2!+ V(x_i)

    ! Avoid trying to fill an NN,NN-1 entry in the matrix
    if (ii.lt.NN-1) then
        H_op(ii+1,ii) = -1/dx**2
    end if

    ! Avoid trying to fill an 0,1 entry in the matrix
    if (ii.gt.1) then
        H_op(ii-1,ii) = -1/dx**2
    end if
end do

! TODO: debug mode
do ii = 1, NN-1
    write(*,*) H_op(ii,:)
end do

end program harmonic_oscillator