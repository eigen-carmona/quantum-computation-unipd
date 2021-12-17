module ising_functs
contains

subroutine diagonalization(org_hermitian, eigenvalues, NN)
! Obtains the e-values of an hermitian matrix
! Modifies org_hermitian
implicit none
integer :: NN, ok, allocate_status, lwork
complex*16 :: org_hermitian(1:NN,1:NN)
real*8 :: eigenvalues(1:NN)
complex*16, allocatable :: work(:)
real*8, allocatable :: rwork(:)

    ! we don't want to saturate the memory
    allocate(work(1:2*NN-1), stat = allocate_status)
    if (allocate_status .ne. 0) stop "***Not enough memory to allocate diagonalization arrays***"

    ! Obtain an optimal lwork
    call zheev('N','U',NN,org_hermitian,NN,eigenvalues,work,-1,rwork,ok)

    ! the first element of work is the optimal lwork
    lwork = int(work(1))
    if (lwork.le.0) stop "overflow in lwork"
    ! Resize the work vector and validate memory
    deallocate(work)
    allocate(work(1:lwork), stat = allocate_status)
    if (allocate_status .ne. 0) stop "***Not enough memory to allocate diagonalization arrays***"

    allocate(rwork(1:3*NN-2), stat = allocate_status)
    if (allocate_status .ne. 0) stop "***Not enough memory to allocate diagonalization arrays***"

    call zheev('N','U',NN,org_hermitian,NN,eigenvalues,work,lwork,rwork,ok)

    ! Releasing memory
    deallocate(work,rwork)

    if (ok.ne.0) stop "Unsuccessful eigenvalue calculation"

end subroutine

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

function tens_prod_2(A,B,MM,NN,OO,PP)
    implicit none
    double complex, dimension(:,:) :: A,B
    double complex, dimension(1:MM*OO,1:NN*PP) :: tens_prod_2
    integer :: ii, jj, MM, NN, OO, PP, base(1:MM), out(1:MM)
    logical :: mask(1:MM)

        base = (/(ii,ii=1,MM)/)
        tens_prod_2 = 0
        do jj = 1, NN
            mask = A(:,jj)/=0
            out = pack(base,mask)
            do ii = 1, count(mask)
                !if(A(ii,jj).eq.0) continue
                tens_prod_2((out(ii)-1)*OO+1:out(ii)*OO,&
                (jj-1)*PP+1:jj*PP) = A(out(ii),jj)*B
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


function tens_id_2(mat_a,MM,NN,id_nn)
    ! Obtains the tensor product of a matrix A and the id_nn identity matrix
    implicit none
    double complex, dimension(:,:) :: mat_a
    double complex, dimension(1:MM*id_nn,1:NN*id_nn) :: tens_id_2
    integer :: ii, jj, MM, NN, id_nn, base(1:MM), out(1:MM)
    logical :: mask(1:MM)
    
        base = (/(ii,ii=1,MM)/)
        tens_id_2 = 0
        do jj = 1, NN
            mask = mat_a(:,jj) /= 0
            out = pack(base,mask)
            do ii = 1, count(mask)
                !if(mat_a(ii,jj).eq.0) continue
                tens_id_2((out(ii)-1)*id_nn+1:out(ii)*id_nn,&
                (jj-1)*id_nn+1:jj*id_nn) = diag_val(id_nn,mat_a(out(ii),jj))
            end do
        end do
end function


end module

program ising_model
use ising_functs
implicit none

double complex, dimension(1:2,1:2) :: sigma_x, sigma_y, sigma_z
double complex, dimension(1:2**2,1:2**2) :: sigma_int
integer ii, NN, allocate_status
double complex, dimension(:,:), allocatable :: hamiltonian,holder_int
real*8, dimension(:), allocatable :: energies
real*8 :: lambda, start,end
character*42 :: energy_file

sigma_x = dcmplx(reshape((/0,1,&
                           1,0/),&
                    (/2,2/)),0)

sigma_y = dcmplx(0,reshape((/0,-1,&
                             1,0/),&
                    (/2,2/)))

sigma_z = dcmplx(reshape((/1,0,&
                           0,-1/),&
                    (/2,2/)),0)

sigma_int = tens_prod_2(sigma_x,sigma_x,2,2,2,2)

write(*,*) "How many spins compose the system? What's the field intensity?"
read(*,*) NN, lambda

! TODO:
! diagonalization
! Lambda variation
! Energy level plotting
! NN timing and fitting
! memory efficiency study (non-zero entries)
call cpu_time(start)
if ((NN.lt.3).or.(NN.gt.30)) stop "The system must have at least 3 spins and no more than 15 spins"


allocate(hamiltonian(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"
allocate(holder_int(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"

! Transverse field interaction
! Setting everything for ii = 1
holder_int = tens_id_2(&
    sigma_z,&
    2,2,2**(NN-1))

hamiltonian = -holder_int


! Setting for states inbetween
do ii = 2, NN-1
    holder_int(1:2**ii,1:2**ii) = tens_prod_2(&
            identity(2**(ii-1)),&
            sigma_z,&
            2**(ii-1),2**(ii-1),2,2)
    holder_int = tens_id_2(&
    holder_int(1:2**ii,1:2**ii),&
    2**ii,2**ii,2**(NN-ii))

    hamiltonian = hamiltonian - holder_int

end do

! Setting everything for ii = NN
holder_int = tens_prod_2(&
    identity(2**(NN-1)),&
    sigma_z,&
    2**(NN-1),2**(NN-1),2,2)

hamiltonian = hamiltonian - holder_int

hamiltonian = lambda*hamiltonian

! Nearest neighbours interaction
! Prepare interaction of the first two
holder_int = tens_id_2(&
                sigma_int,&
                2**2,2**2,2**(NN-2))

hamiltonian = hamiltonian - holder_int

do ii = 2, NN-2
    holder_int(1:2**(ii+1),1:2**(ii+1)) = tens_prod_2(&
                    identity(2**(ii-1)),&
                    sigma_int,&
                    2**(ii-1),2**(ii-1),2**2,2**2)
    holder_int = tens_id_2(&
                    holder_int(1:2**(ii+1),1:2**(ii+1)),&
                    2**(ii+1),2**(ii+1),2**(NN-ii-1))
    hamiltonian = hamiltonian - holder_int
end do

! Prepare interaction of the last two
holder_int = tens_prod_2(&
                identity(2**(NN-2)),&
                sigma_int,&
                2**(NN-2),2**(NN-2),2**2,2**2)

hamiltonian = hamiltonian - holder_int

! Representation efficiency check
write(*,*) count(hamiltonian/=0), "non-zero elements out of ", 2**NN*2**NN

call cpu_time(end)
write(*,*) "Hamiltonian building time in seconds: ", end-start

! Once the hamiltonian has been computed, we can release the memory from holder_int
deallocate(holder_int)
! Now we may allocate the energy array
allocate(energies(1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"

call diagonalization(hamiltonian,energies,2**NN)

call cpu_time(end)
write(*,*) "Total time in seconds: ", end-start

write(*,*) "****ENERGIES****"
write(*,*) energies

if (lambda.lt.1) then
    write(energy_file,'(A,I0,A,f0.2,A)') "data/energies_",NN,"_spins_lambda_0",lambda,".dat"
else
    write(energy_file,'(A,I0,A,f0.2,A)') "data/energies_",NN,"_spins_lambda_",lambda,".dat"
end if
open(12,file = energy_file)
do ii = 1,2**NN
    write(12,*) energies(ii)
end do
close(12)
write(*,*) "Successfully wrote to ", energy_file

end program ising_model
