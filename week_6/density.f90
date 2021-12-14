program density_matrix
! executed as ```gfortran density.f90 -o density_matrix definitions.o -llapack```

use definitions
implicit none

interface
    function rand_wf(NN)
    implicit none
    integer :: NN
    double complex, dimension(1:NN) :: rand_wf
    end function
end interface

double complex, allocatable :: wvfn(:), density(:,:), red_dens(:,:), matsie(:,:), wvfn_2(:), two_qubit(:,:)
logical :: separable
integer :: ii, jj, hilbert_dim, DD, NN, ok
complex*16, allocatable :: org_hermitian(:,:), work(:), eigenvectors(:,:)
real*8,allocatable :: eigenvalues(:), rwork(:)


write(*,*) "What's the dimension D of the Hilbert space?"
read(*,*) DD
write(*,*) "How many subsystems compose the system?"
read(*,*) NN
write(*,*) "Is the system separable?"
read(*,*) separable

! Describe a separable and a non-separable state
if (separable) then
    write(*,*) "Building representation of separable state"
    hilbert_dim = DD*NN
    allocate(wvfn(1:hilbert_dim))
    do ii = 1, NN
        wvfn((ii-1)*DD+1:ii*DD) = rand_wf(DD)
    end do
else
    write(*,*) "Building representation of non-separable state"
    hilbert_dim = DD**NN
    allocate(wvfn(1:hilbert_dim))
    wvfn = rand_wf(hilbert_dim)
end if

! Build the density matrix for an N = 2 pure system
hilbert_dim = DD**2
allocate(wvfn_2(1:hilbert_dim))
allocate(density(1:hilbert_dim,1:hilbert_dim))
! Generate the pure system
wvfn_2 = rand_wf(hilbert_dim)
! Build the density matrix (outer product of wavefunction)
do jj = 1, hilbert_dim
    do ii = 1, hilbert_dim
        density(ii,jj) = wvfn_2(ii)*conjg(wvfn_2(jj))
    end do
end do

! get reduced matrix for a generic N = 2 density matrix
! Theorem: Any positive trace-1 matrix is a density matrix 
density = rand_positive(hilbert_dim)
density = density/mat_trace(density)
allocate(red_dens(1:DD,1:DD))
allocate(matsie(1:DD,1:DD))
do jj = 1, DD
    do ii = 1, DD
        matsie = density((ii-1)*DD+1:ii*DD,(jj-1)*DD+1:jj*DD)
        red_dens(ii,jj) = mat_trace(matsie)
    end do
end do

! apply to a two-qubit system
allocate(two_qubit(1:2**2,1:2**2))
two_qubit = rand_positive(2**2)
two_qubit = two_qubit/mat_trace(two_qubit)

! Obtain eigenvalues
allocate(org_hermitian(1:2**2,1:2**2))
allocate(eigenvalues(1:2**2))
allocate(work(1:max(5000,2*2**2-1)))
allocate(rwork(1:3*2**2-2))

! copy the hermitian matrix to avoid modifying it
eigenvectors = two_qubit

! Obtain an optimal lwork
call zheev('N','U',2**2,eigenvectors,2**2,eigenvalues,work,-1,rwork,ok)
! the first element of work is the optimal lwork
call zheev('N','U',2**2,eigenvectors,2**2,eigenvalues,work,min(5000,int(work(1))),rwork,ok)

if (ok.ne.0) stop "Unsuccessful eigenvalue calculation for two-qubit system"

write(*,*) "Eigenvalues of the two-qubit density matrix"
do ii = 1, 2**2
    write(*,*) eigenvalues(ii)
end do


end program density_matrix


function rand_wf(NN)
implicit none
integer :: NN, ii
double complex, dimension(1:NN) :: rand_wf
real*8, dimension(1:NN) :: r_wf, i_wf

    call random_number(r_wf)
    call random_number(i_wf)

    r_wf = 2*r_wf - 1
    i_wf = 2*i_wf - 1

    do ii = 1, NN
        rand_wf(ii) = complex(r_wf(ii),i_wf(ii))
    end do

    ! Normalized state
    rand_wf = rand_wf/dot_product(rand_wf,rand_wf)

end function rand_wf

