program time_independent
! Compile as ```gfortran time_independent.f90 -o time_independent finite_difference.o -llapack```
! After compiling ```gfortran -c finite_difference.f90 -llapack```
use finite_diff
implicit none

interface
    function potential_V(N_grid,x_pos,election)
    implicit none
    integer     :: N_grid
    real*8      :: x_pos(1:N_grid)
    character*2 :: election
    real*8      :: potential_V(1:N_grid)
    end function
end interface

integer :: NN, ii, n_eigen, INFO, N_grid
real*8, parameter :: x_0 = -10, x_n = 10
real*8 :: dx
real*8, allocatable :: energies(:), psi_states(:,:), V_pot(:), x_pos(:)
character*42 :: energy_file, states_file
character*2 :: election

write(*,*) "Currently enabled potential choices:"
write(*,'(A,/,A,/,A,/,A)') " - 'ho': harmonic oscillator"," - 'qs': V=5 step at x = 0",&
" - 'qb': V = 5 barrier when x in (-5,5)", " - 'qw': V = -5 well when x in (-5,5)"
write(*,*) "Any option not listed leads to the (-10,10) quantum box."
write(*,*) "Enter the potential type, grid size and the number of desired energy levels"
read(*,*) election, NN, n_eigen

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
V_pot = potential_V(N_grid,x_pos,election)

! Calling single particle routine for diagonal potentials and vanishing boundaries
call finite_diff_1D(x_0,x_n,NN,n_eigen,energies,psi_states,V_pot,INFO)

! Validate there are no routine errors
if (INFO.ne.0) then
    write(*,'(A,I4)') "Unsuccessful calculation. Status = ",INFO
stop "Not storing results. Exiting program..."
end if

write(energy_file,'(A,A,A,I0,A)') "energies_",election,"_",NN,".dat"
write(states_file,'(A,A,A,I0,A)') "psi_states_",election,"_",NN,".dat"
open(15, file = energy_file)
open(16, file = states_file)
do ii = 1, n_eigen
    write(15,*) energies(ii)
    write(16,*) psi_states(:,ii)! Storing as row vectors
end do
close(16)
close(15)

end program time_independent



! Catallogue of potentials
function potential_V(N_grid,x_pos,election)
implicit none
integer     :: N_grid, ii
real*8      :: x_pos(1:N_grid)
character*2 :: election
real*8      :: potential_V(1:N_grid)
    select case (election)
    case ('ho')
        ! Harmonic oscillator potential
        potential_V = x_pos**2/2
    case ('qs')
        ! Step potential
        do ii = 1, N_grid
            if (x_pos(ii).lt.0) then
                potential_V(ii) = 0
            else
                potential_V(ii) = 1
            end if
        end do
    case ('qb','qw')
        ! Barrier potential
        do ii = 1, N_grid
            if ((x_pos(ii).lt.-5).or.(x_pos(ii).gt.5)) then
                potential_V(ii) = 0
            else
                potential_V(ii) = 1
            end if
        end do
        ! Potential well
        if (election.eq.'we') then
            potential_V = -potential_V
        end if
    case default
        potential_V = (/(0,ii = 1,N_grid)/)
        election = 'bx'
    end select
end function
