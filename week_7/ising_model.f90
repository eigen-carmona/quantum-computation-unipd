module ising_functs
contains

subroutine print_mat(mat_a)
implicit none
double complex, dimension(:,:) :: mat_a
integer :: ii, shape_a(1:2)
shape_a = shape(mat_a)
do ii = 1, shape_a(1)
    write(*,'(*(F0.0,SP,F0.0,"i",","))') mat_a(ii,:)
end do
write(*,*) "---------------------------------------"
end subroutine

subroutine compare_mat(mat_a,mat_b)
implicit none
double complex, dimension(:,:) :: mat_a,mat_b
integer :: ii, jj, shape_a(1:2), shape_b(1:2)
    shape_a = shape(mat_a)
    shape_b = shape(mat_b)
    do ii = 1,2
        if (shape_a(ii).ne.shape_b(ii)) stop "Matrices differ in shape"
    end do
    do jj = 1, shape_a(2)
        do ii = 1, shape_a(1)
            if(&
                (real(real(mat_a(ii,jj))).ne.real(real(mat_b(ii,jj))))&
                .or.&
                (real(aimag(mat_a(ii,jj))).ne.real(aimag(mat_b(ii,jj))))&
                ) then
                write(*,'(A,I0)') "matrix A, row ", ii
                call print_mat(mat_a(ii:ii,:))
                write(*,'(A,I0)') "matrix B, row ", ii
                call print_mat(mat_b(ii:ii,:))
                write(*,'(A,I0)') "different in element ",jj
                stop
            end if
        end do
    end do    
end subroutine

function tens_prod(A,B,MM,NN,OO,PP)
implicit none
double complex, dimension(:,:) :: A,B
double complex, dimension(1:MM*OO,1:NN*PP) :: tens_prod
integer :: ii, jj, MM, NN, OO, PP

    !tens_prod = 0
    do jj = 1, NN
        do ii = 1, MM
            !if(A(ii,jj).eq.0) continue
            tens_prod((ii-1)*OO+1:ii*OO,&
            (jj-1)*PP+1:jj*PP) = A(ii,jj)*B
        end do
    end do
end function

function identity(NN)
implicit none
integer :: NN, ii
double complex, dimension(1:NN,1:NN) :: identity

    identity = 0
    do ii = 1, NN
        identity(ii,ii) = dcmplx(1,0)
    end do

end function

function diag_val(OO,val)
implicit none
integer :: OO, ii
double complex :: val, diag_val(1:OO,1:OO)

    diag_val = 0
    !if (val.eq.dcmplx(0,0)) return
    do ii = 1, OO
        diag_val(ii,ii) = val
    end do

end function

function tens_id(mat_a,MM,NN,id_nn)
    ! Obtains the tensor product of a matrix A and the id_nn identity matrix
    implicit none
    double complex, dimension(:,:) :: mat_a
    double complex, dimension(1:MM*id_nn,1:NN*id_nn) :: tens_id
    integer :: ii, jj, MM, NN, id_nn
    
        !tens_id = 0
        do jj = 1, NN
            do ii = 1, MM
                !if(mat_a(ii,jj).eq.0) continue
                tens_id((ii-1)*id_nn+1:ii*id_nn,&
                (jj-1)*id_nn+1:jj*id_nn) = diag_val(id_nn,mat_a(ii,jj))
            end do
        end do
end function    

end module

program ising_model
use ising_functs
implicit none

double complex, dimension(1:2,1:2) :: sigma_x, sigma_y, sigma_z
double complex, dimension(1:2**2,1:2**2) :: sigma_int
integer ii, jj, NN, allocate_status
double complex, dimension(:,:), allocatable :: hamiltonian,mean_int,neigh_int, hamiltonian_1
real*8 :: lambda, start, end

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



call cpu_time(start)
!
do jj = 1, 100
! Transverse field interaction
! Setting everything for ii = 1
mean_int = tens_prod(&
    sigma_z,&
    identity(2**(NN-1)),&
    2,2,2**(NN-1),2**(NN-1))

hamiltonian = mean_int


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

end do
call cpu_time(end)
write(*,*) "Timing first method:", end-start

! **** ALTERNATIVE ****
call cpu_time(start)
do jj = 1, 100
! Transverse field interaction
! Setting everything for ii = 1
mean_int = tens_id(&
    sigma_z,&
    2,2,2**(NN-1))

hamiltonian_1 = mean_int


! Setting for states inbetween
do ii = 2, NN-1
    mean_int(1:2**ii,1:2**ii) = tens_prod(&
            identity(2**(ii-1)),&
            sigma_z,&
            2**(ii-1),2**(ii-1),2,2)
    mean_int = tens_id(&
    mean_int(1:2**ii,1:2**ii),&
    2**ii,2**ii,2**(NN-ii))

    hamiltonian_1 = hamiltonian_1 + mean_int

end do

! Setting everything for ii = NN
mean_int = tens_prod(&
    identity(2**(NN-1)),&
    sigma_z,&
    2**(NN-1),2**(NN-1),2,2)

hamiltonian_1 = hamiltonian_1 + mean_int

! Nearest neighbours interaction
! Prepare interaction of the first two
neigh_int = tens_id(&
                sigma_int,&
                2**2,2**2,2**(NN-2))

hamiltonian_1 = hamiltonian_1 + neigh_int

do ii = 2, NN-2
    neigh_int(1:2**(ii+1),1:2**(ii+1)) = tens_prod(&
                    identity(2**(ii-1)),&
                    sigma_int,&
                    2**(ii-1),2**(ii-1),2**2,2**2)
    neigh_int = tens_id(&
                    neigh_int(1:2**(ii+1),1:2**(ii+1)),&
                    2**(ii+1),2**(ii+1),2**(NN-ii-1))
    hamiltonian_1 = hamiltonian_1 + neigh_int
end do

! Prepare interaction of the last two
neigh_int = tens_prod(&
                identity(2**(NN-2)),&
                sigma_int,&
                2**(NN-2),2**(NN-2),2**2,2**2)

hamiltonian_1 = hamiltonian_1 + neigh_int
end do
call cpu_time(end)
write(*,*) "Timing second method:", end-start


call compare_mat(hamiltonian,hamiltonian_1)

end program ising_model
