program eigenproblem
implicit none
interface
    function random_hermitian(N)
    integer*4 :: N
    complex*8, dimension(1:N,1:N) :: random_hermitian
    end function
end interface

integer*4 :: N, ii
complex*8, allocatable :: hermitian(:,:)

write(*,*) "How many rows should the hermitian matrix have?"
read(*,*) N

allocate(hermitian(1:N,1:N))

hermitian = random_hermitian(N)

do ii = 1,N
    write(*,*) hermitian(ii,:)
end do

end program eigenproblem

function random_hermitian(N)
implicit none
integer*4 :: N, ii, jj
complex*8, dimension(1:N,1:N) :: random_hermitian
real*8, dimension(1:N,1:N) :: A, B

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