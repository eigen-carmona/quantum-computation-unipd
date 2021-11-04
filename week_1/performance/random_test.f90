program test_mat
implicit none

real*4, allocatable:: A(:,:)
integer:: N,ii,jj

write(*,*) "give me the maximum number of rows for the square matrix"
read(*,*) N

do ii = 2,N
	allocate(A(1:ii,1:ii))
	call random_number(A)
	do jj = 1,ii
		write(*,*) A(jj,1:ii)
	end do
	write(*,*) ""
	deallocate(A)
end do

end program test_mat
