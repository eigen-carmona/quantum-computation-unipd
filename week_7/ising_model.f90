program ising_model
implicit none
interface
    function tens_prod(A,B,MM,NN,OO,PP)
        implicit none
        integer :: MM,NN,OO,PP
        double complex, dimension(:,:) :: A,B
        double complex, dimension(1:NN*OO,1:MM*PP) :: tens_prod
    end function
end interface

interface
    function identity(NN)
        implicit none
        integer :: NN
        double complex, dimension(1:NN,1:NN) :: identity
    end function
end interface

double complex, dimension(1:2,1:2) :: sigma_x, sigma_y, sigma_z
double complex, dimension(1:2**2,1:2**2) :: sigma_int
integer ii, jj, NN, allocate_status
double complex, dimension(:,:), allocatable :: hamiltonian,mean_int,neigh_int
real*8 :: lambda


sigma_x = dcmplx(reshape((/0,1,&
                           1,0/),&
                    (/2,2/)),0)

sigma_y = dcmplx(0,reshape((/0,-1,&
                             1,0/),&
                    (/2,2/)))

sigma_z = dcmplx(reshape((/1,0,&
                           0,-1/),&
                    (/2,2/)),0)

sigma_int = tens_prod(sigma_x,sigma_x,2,2,2,2)

write(*,*) "How many spins compose the system? Whats the field intensity?"
read(*,*) NN, lambda

if (NN.lt.3) stop "The system must have at least 3 spins"

! BRUTE FORCE TENSOR PRODUCT
! Can be simplified noting that most matrices in the tensor product are the identity
allocate(hamiltonian(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"
allocate(mean_int(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"
allocate(neigh_int(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"


! Transverse field interaction
! Setting everything for ii = 1
mean_int = tens_prod(&
    sigma_z,&
    identity(2**(NN-1)),&
    2,2,2**(NN-1),2**(NN-1))

hamiltonian = hamiltonian + mean_int

! Setting for states inbetween
do ii = 2, NN-1
    mean_int(1:2**ii,1:2**ii) = tens_prod(&
            identity(2**(ii-1)),&
            sigma_z,&
            2**(ii-1),2**(ii-1),2,2)
    mean_int = tens_prod(&
    mean_int(1:2**ii,1:2**ii),&
    identity(2**(NN-ii)),&
    2**ii,2**ii,2**(NN-ii),2**(NN-ii))

    hamiltonian = hamiltonian + mean_int

end do

! Setting everything for ii = NN
mean_int = tens_prod(&
    identity(2**(NN-1)),&
    sigma_z,&
    2**(NN-1),2**(NN-1),2,2)

hamiltonian = hamiltonian + mean_int
hamiltonian = lambda*hamiltonian

! Nearest neighbours interaction
! Prepare interaction of the first two
neigh_int = tens_prod(&
                sigma_int,&
                identity(2**(NN-2)),&
                2**2,2**2,2**(NN-2),2**(NN-2))

hamiltonian = hamiltonian + neigh_int

do ii = 2, NN-2
    neigh_int(1:2**(ii+1),1:2**(ii+1)) = tens_prod(&
                    identity(2**(ii-1)),&
                    sigma_int,&
                    2**(ii-1),2**(ii-1),2**2,2**2)
    neigh_int = tens_prod(&
                    neigh_int(1:2**(ii+1),1:2**(ii+1)),&
                    identity(2**(NN-ii-1)),&
                    2**(ii+1),2**(ii+1),2**(NN-ii-1),2**(NN-ii-1))
    hamiltonian = hamiltonian + neigh_int
end do

! Prepare interaction of the last two
neigh_int = tens_prod(&
                identity(2**(NN-2)),&
                sigma_int,&
                2**(NN-2),2**(NN-2),2**2,2**2)

hamiltonian = hamiltonian + neigh_int


end program ising_model


function tens_prod(A,B,MM,NN,OO,PP)
implicit none
double complex, dimension(:,:) :: A,B
double complex, dimension(1:MM*OO,1:NN*PP) :: tens_prod
integer :: ii, jj, MM, NN, OO, PP

    do ii = 1, MM
        do jj = 1, NN
            tens_prod((ii-1)*OO+1:ii*OO,&
            (jj-1)*PP+1:jj*PP) = A(ii,jj)*B
        end do
    end do
end function

function identity(NN)
implicit none
integer :: NN, ii
double complex, dimension(1:NN,1:NN) :: identity

identity = 0*identity
do ii = 1, NN
    identity(ii,ii) = dcmplx(1,0)
end do

end function
