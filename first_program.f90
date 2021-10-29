program era_sieve

implicit none
integer*2 :: N ! The natural number up to which prime numbers will be calculated
integer*2 :: ii, jj ! indexing variables
logical, allocatable :: is_prime(:)
logical :: index_value, true = 1.gt.0, false = 0.gt.1

write(*,*) "Up to which number do you wish to get primes?"
read(*,*) N

allocate(is_prime(2:N))
is_prime = (/(true, ii =2, N)/)

ii = 2

! Erathostenes Sieve
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
