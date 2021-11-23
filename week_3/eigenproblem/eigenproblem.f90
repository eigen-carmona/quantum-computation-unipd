program eigenproblem
implicit none

interface
    function random_hermitian(N)
    integer :: N
    complex*16, dimension(1:N,1:N) :: random_hermitian
    end function
end interface

integer :: N, ii, ok
complex*16, allocatable :: org_hermitian(:,:), work(:), eigenvectors(:,:)
real*8,allocatable :: eigenvalues(:), rwork(:), spacings(:)
real*8 :: average_spacing
character (len = 15) :: filename

write(*,*) "How many rows should the hermitian matrix have?"
write(*,*) "Type in the name of the file to store the normalized spacings. Type '' if no storing is desired."
read(*,*) N, filename

allocate(org_hermitian(1:N,1:N))
allocate(eigenvalues(1:N))
allocate(spacings(1:N-1))
allocate(work(1:max(5000,2*N-1)))
allocate(rwork(1:3*N-2))

! generate random hermitian matrix
org_hermitian = random_hermitian(N)

! copy the hermitian matrix to avoid modifying it
eigenvectors = org_hermitian

! Obtain an optimal lwork
call zheev('N','U',N,eigenvectors,N,eigenvalues,work,-1,rwork,ok)
! the first element of work is the optimal lwork
call zheev('N','U',N,eigenvectors,N,eigenvalues,work,min(5000,int(work(1))),rwork,ok)

if (ok.ne.0) stop "Unsuccessful eigenvalue calculation"

! obtaing the normalized spacings
average_spacing = (eigenvalues(N)-eigenvalues(1))/(N-1)
do ii = 1, N-1
    spacings(ii) = (eigenvalues(ii+1)-eigenvalues(ii))/average_spacing
end do

! log or store the normalized spacings
if (filename .ne. '') then
    open(15, file =filename)
    do ii = 1, N-1
        write(15,*) spacings(ii)
    end do
    close(15)
else
    write(*,*) "******SPACINGS******"
    do ii = 1, N-1
        write(*,*) spacings(ii)
    end do
end if

end program eigenproblem


function random_hermitian(N)
! Generates an order N hermitian matrix with random entries x+iy
! where each x and y are uniformly distributed between -1 and 1
implicit none
integer :: N, ii, jj
complex*16, dimension(N,N) :: random_hermitian
real*8, dimension(N,N) :: A, B

    call random_number(A)
    ! rescaling into an element in (-1,1)
    A = 2*A - 1
    call random_number(B)
    ! rescaling into an element in (-1,1)
    B = 2*B - 1

    ! populating the matrix
    do ii = 1,N
        do jj = 1,N
            ! replacing entries in the lower triangle
            if ((jj.gt.ii).and.(ii.lt.N)) then
                A(jj,ii) = A(ii,jj)
                ! The imaginary part is negative, to achieve hermicity
                B(jj,ii) = -B(ii,jj)
            end if
            ! cast into a complex entry
            random_hermitian(ii,jj) = complex(A(ii,jj),B(ii,jj))
        end do
    end do

end function random_hermitian
