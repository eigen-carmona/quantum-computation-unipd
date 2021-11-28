program harmonic_oscillator
! Compile as ```gfortran time_independent.f90 -o harmonic_oscillator finite_difference.o -llapack```
! After compiling ```gfortran -c finite_difference.f90 -llapack```
use finite_diff
implicit none

integer :: NN, ii, n_eigen, INFO
real*8, parameter :: x_0 = -10, x_n = 10
real*8 :: dx=0.0005, x_i = x_0
real*8, allocatable :: energies(:), psi_states(:,:), V_pot(:), x_pos(:)

! Compute the potential energy array
NN = int((x_n-x_0)/dx)-1
allocate(V_pot(1:NN))
allocate(x_pos(1:NN))

do ii = 1, NN
    x_i = x_i + dx
    x_pos(ii) = x_i
end do

V_pot = x_pos**2/2

write(*,*) "How many energy values should be calculated?"
read(*,*) n_eigen

call finite_diff_1D(x_0,x_n,dx,n_eigen,energies,psi_states,V_pot,INFO)

write(*,*) "status:",INFO, "(0 means successful)"
open(15, file = 'energies3.dat')
open(16, file = 'psi_states3.dat')
do ii = 1, n_eigen
    write(15,*) energies(ii)
    write(16,*) psi_states(:,ii)! Storing as row vectors
end do
close(16)
close(15)

end program harmonic_oscillator
