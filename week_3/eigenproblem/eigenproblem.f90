program eigenproblem
implicit none

interface
    function random_hermitian(N)
    integer :: N
    complex*16, dimension(1:N,1:N) :: random_hermitian
    end function
end interface

integer :: N, ii, ok
complex*16, allocatable :: hermitian(:,:), work(:)
real*16,allocatable :: eigenvalues(:), rwork(:), spacings(:)

write(*,*) "How many rows should the hermitian matrix have?"
read(*,*) N

allocate(hermitian(1:N,1:N))
allocate(eigenvalues(1:N))
allocate(spacings(1:N-1))
allocate(work(1:2*N-1))
allocate(rwork(1:3*N-2))

hermitian = random_hermitian(N)

do ii = 1,N
    write(*,*) hermitian(ii,:)
end do

call zheev('N','U',N,hermitian,N,eigenvalues,work,2*N,work,ok)

do ii = 1, N-1
    spacings(ii) = eigenvalues(ii+1)-eigenvalues(ii)
end do

spacings = spacings*(N-1)/sum(spacings)

do ii = 1, N-1
    write(*,*) spacings(ii)
end do


end program eigenproblem

function random_hermitian(N)
implicit none
integer :: N, ii, jj
complex*16, dimension(1:N,1:N) :: random_hermitian
real*16, dimension(1:N,1:N) :: A, B

    call random_number(A)
    A = 2*A - 1
    call random_number(B)
    B = 2*B - 1
    do ii = 1, N-1
        do jj = ii+1, N
            A(jj,ii) = A(ii,jj)
            B(jj,ii) = -B(ii,jj)
        end do
    end do

    do ii = 1, N
        do jj = 1, N
            random_hermitian(ii,jj) = complex(A(ii,jj),B(ii,jj))
        end do
    end do

end function random_hermitian