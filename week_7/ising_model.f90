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

double complex, dimension(1:2,1:2) :: sigma_x, sigma_y, sigma_z, mat
double complex, dimension(1:2**2,1:2**2) :: sigma_int
integer ii, jj, NN, allocate_status
double complex, dimension(:,:), allocatable :: hamiltonian,mean_int,neigh_int,hamiltonian_1
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
allocate(hamiltonian_1(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"



! Transverse field interaction
do ii = 1, NN
    if (ii.eq.1) then
        mean_int(1:2,1:2) = sigma_z
    else
        mean_int(1:2,1:2) = identity(2)
    end if

    do jj = 2, NN
        if (jj.eq.ii) then
            mat = sigma_z
        else
            mat = identity(2)
        end if
        mean_int(1:2**jj,1:2**jj) = tens_prod(&
            mean_int(1:2**(jj-1),1:2**(jj-1)),&
            mat,&
            2**(jj-1),2**(jj-1),2,2)
    end do
    hamiltonian = hamiltonian + mean_int

end do

! Setting everything for ii = 1
mean_int = tens_prod(&
    sigma_z,&
    identity(2**(NN-1)),&
    2,2,2**(NN-1),2**(NN-1))

hamiltonian_1 = hamiltonian_1 + mean_int

! Setting for states inbetween
do ii = 2, NN-1
    ! Initialize whole matrix as identity
    mean_int = identity(2**NN)
    mean_int(1:2**ii,1:2**ii) = tens_prod(&
            identity(2**(ii-1)),&
            sigma_z,&
            2**(ii-1),2**(ii-1),2,2)
    mean_int = tens_prod(&
    mean_int(1:2**ii,1:2**ii),&
    identity(2**(NN-ii)),&
    2**ii,2**ii,2**(NN-ii),2**(NN-ii))

    hamiltonian_1 = hamiltonian_1 + mean_int

end do

! Setting everything for ii = NN
mean_int = tens_prod(&
    identity(2**(NN-1)),&
    sigma_z,&
    2**(NN-1),2**(NN-1),2,2)

hamiltonian_1 = hamiltonian_1 + mean_int

hamiltonian = lambda*hamiltonian
hamiltonian_1 = lambda*hamiltonian_1

! Nearest neighbours interaction
do ii = 1, NN-1
    if (ii.eq.1) then
        neigh_int(1:2,1:2) = sigma_x
    else
        neigh_int(1:2,1:2) = identity(2)
    end if

    do jj = 2, NN
        if ((jj.eq.ii+1).or.(jj.eq.ii)) then
            mat = sigma_x
        else
            mat = identity(2)
        end if
        neigh_int(1:2**jj,1:2**jj) = tens_prod(&
            neigh_int(1:2**(jj-1),1:2**(jj-1)),&
            mat,&
            2**(jj-1),2**(jj-1),2,2)
    end do
    hamiltonian = hamiltonian + neigh_int
end do

! Nearest neighbours interaction
! Prepare interaction of the first two
neigh_int = tens_prod(&
                sigma_int,&
                identity(2**(NN-2)),&
                2**2,2**2,2**(NN-2),2**(NN-2))

hamiltonian_1 = hamiltonian_1 + neigh_int

do ii = 2, NN-2
    neigh_int(1:2**(ii+1),1:2**(ii+1)) = tens_prod(&
                    identity(2**(ii-1)),&
                    sigma_int,&
                    2**(ii-1),2**(ii-1),2**2,2**2)
    neigh_int = tens_prod(&
                    neigh_int(1:2**(ii+1),1:2**(ii+1)),&
                    identity(2**(NN-ii-1)),&
                    2**(ii+1),2**(ii+1),2**(NN-ii-1),2**(NN-ii-1))
    hamiltonian_1 = hamiltonian_1 + neigh_int
end do

! Prepare interaction of the last two
neigh_int = tens_prod(&
                identity(2**(NN-2)),&
                sigma_int,&
                2**(NN-2),2**(NN-2),2**2,2**2)

hamiltonian_1 = hamiltonian_1 + neigh_int

!do ii = 1, 2**NN
!    do jj = 1, 2**NN
!        if ((real(real(hamiltonian_1(ii,jj))).ne.real(real(hamiltonian(ii,jj))))&
!        .or.((real(aimag(hamiltonian_1(ii,jj))).ne.real(aimag(hamiltonian(ii,jj))))))&
!        stop "They're different :("
!
!    end do
!end do
!
!
!do ii = 1, 2**NN
!    write(*,'(*(F0.0,SP,F0.0,"i",", "))') hamiltonian(ii,:)
!    write(*,'(*(F0.0,SP,F0.0,"i",", "))') hamiltonian_1(ii,:)
!    write(*,*) ""
!end do

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
