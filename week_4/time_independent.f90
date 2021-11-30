program harmonic_oscillator
! Compile as ```gfortran time_independent.f90 -o harmonic_oscillator finite_difference.o -llapack```
! After compiling ```gfortran -c finite_difference.f90 -llapack```
use finite_diff
implicit none

integer :: NN, ii, n_eigen, INFO, N_grid
real*8, parameter :: x_0 = -10, x_n = 10
real*8 :: dx
real*8, allocatable :: energies(:), psi_states(:,:), V_pot(:), x_pos(:)
character*25 :: energy_file, states_file

write(*,*) "Enter the grid size and the number of desired energy levels"
read(*,*) NN, n_eigen

if (NN.le.n_eigen) stop "The grid size must be higher than the number of energy levels"

! The grid size determines the increment value
dx = (x_n-x_0)/NN

if (dx.ge.1) stop "Grid size too small. Minumum should be 20 (advised 400 and above)."

! Since we assume x_0=x_n=0, we reduce the problem to a (NN-1)x(NN-1) matrix
N_grid = NN - 1

! Compute the potential energy array
allocate(V_pot(1:N_grid))
allocate(x_pos(1:N_grid))
x_pos = (/(ii*dx + x_0, ii = 1,N_grid)/)
V_pot = x_pos**2/2

! Calling single particle routine for diagonal potentials and vanishing boundaries
call finite_diff_1D(x_0,x_n,NN,n_eigen,energies,psi_states,V_pot,INFO)

! Validate there are no routine errors
if (INFO.ne.0) then
    write(*,'(A,I4)') "Unsuccessful calculation. Status = ",INFO
stop "Not storing results. Exiting program..."
end if

write(energy_file,'(A,I4,A)') "energies",NN,".dat"
write(states_file,'(A,I4,A)') "psi_states",NN,".dat"
open(15, file = energy_file)
open(16, file = states_file)
do ii = 1, n_eigen
    write(15,*) energies(ii)
    write(16,*) psi_states(:,ii)! Storing as row vectors
end do
close(16)
close(15)

end program harmonic_oscillator
