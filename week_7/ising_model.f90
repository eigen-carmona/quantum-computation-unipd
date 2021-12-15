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

double complex, dimension(1:2,1:2) :: sigma_x, sigma_y, sigma_z, sigma_t, mat
integer ii, jj, NN, allocate_status
double complex, dimension(:,:), allocatable :: hamiltonian,mean_int,neigh_int

sigma_t = dcmplx(reshape((/1,0,&
                           0,1/),&
                    (/2,2/)),0)

sigma_x = dcmplx(reshape((/0,1,&
                           1,0/),&
                    (/2,2/)),0)

sigma_y = dcmplx(0,reshape((/0,-1,&
                             1,0/),&
                    (/2,2/)))

sigma_z = dcmplx(reshape((/1,0,&
                           0,-1/),&
                    (/2,2/)),0)

write(*,*) "How many spins compose the system?"
read(*,*) NN

! BRUTE FORCE TENSOR PRODUCT
! Can be simplified noting that most matrices in the tensor product are the identity
allocate(hamiltonian(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"
allocate(mean_int(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"
allocate(neigh_int(1:2**NN,1:2**NN),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate hamiltonian***"

! Transverse field interaction
do ii = 1, NN
    if (ii.eq.1) then
        mean_int(1:2,1:2) = sigma_z
    else
        mean_int(1:2,1:2) = sigma_t
    end if

    do jj = 2, NN
        if (jj.eq.ii) then
            mat = sigma_z
        else
            mat = sigma_t
        end if
        mean_int(1:2**jj,1:2**jj) = tens_prod(&
            mean_int(1:2**(jj-1),1:2**(jj-1)),&
            mat,&
            2**(jj-1),2**(jj-1),2,2)
    end do
    hamiltonian = hamiltonian + mean_int
end do

! Nearest neighbours interaction
do ii = 1, NN-1
    if (ii.eq.1) then
        neigh_int(1:2,1:2) = sigma_x
    else
        neigh_int(1:2,1:2) = sigma_t
    end if
    
    do jj = 2, NN
        if ((jj.eq.ii+1).or.(jj.eq.ii)) then
            mat = sigma_x
        else
            mat = sigma_t
        end if
        neigh_int(1:2**jj,1:2**jj) = tens_prod(&
            neigh_int(1:2**(jj-1),1:2**(jj-1)),&
            mat,&
            2**(jj-1),2**(jj-1),2,2)
    end do
    hamiltonian = hamiltonian + neigh_int
end do

do ii = 1, 2**NN
    write(*,'(*(F0.2,SP,F0.2,"i",", "))') hamiltonian(ii,:)
end do

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

