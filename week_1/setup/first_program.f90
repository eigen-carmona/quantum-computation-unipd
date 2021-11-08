
program era_sieve

implicit none
integer*4 :: N ! The natural number up to which prime numbers will be calculated
integer*4 :: ii, jj ! indexing variables
integer*4 :: allocate_status
logical, allocatable :: is_prime(:)
logical :: true = 1.gt.0, false = 0.gt.1

write(*,*) "Up to which number do you wish to get primes? Enter an integer N>=2."
read(*,*) N

! N should be at least 2
if (N.lt.2) stop "***The upper limit should be at least 2***"

allocate(is_prime(2:N),stat = allocate_status)
if (allocate_status .ne. 0) stop "***Not enough memory to allocate sieve vector***"
is_prime = (/(true, ii =2, N)/)

ii = 2

! Eratosthenes Sieve
do while (ii .le. sqrt(real(N,4))) ! Casting N as real to calculate sqrt(N)
	if (is_prime(ii)) then
		jj = ii**2
		do while (jj .le. N) ! Spotting multiples of ii starting from ii**2
			is_prime(jj) = false
			jj = jj + ii
		end do
	end if
	ii = ii + 1
end do


do ii = 2,N
	if (is_prime(ii)) then
		write(*,*) ii
	end if
end do

write(*,*) "Enter any number to normally exit the program."
read(*,*) N

end program era_sieve
