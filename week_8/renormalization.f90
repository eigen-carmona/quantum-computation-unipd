module renormalization
    contains


    function mat_adj(entries) result(adj_entries)
        implicit none
        complex*16, intent(in) :: entries(:,:)
        complex*16, allocatable :: adj_entries(:,:)
        integer*4 :: ii, jj, NN, MM
        integer*4 :: mat_shape(1:2)

        mat_shape = shape(entries)
        NN = mat_shape(1)
        MM = mat_shape(2)
        allocate(adj_entries(1:MM,1:NN))
        do ii = 1, NN
            do jj = 1, MM
                adj_entries(jj, ii) = conjg(entries(ii,jj))
            end do
        end do
    end function mat_adj


    subroutine zallocate_1(array,MM)
    implicit none
    integer :: allocate_status, MM
    double complex, allocatable :: array(:)

        allocate(array(1:MM), stat = allocate_status)
        if (allocate_status .ne. 0) stop "***Not enough memory to allocate arrays***"

    end subroutine

    
    subroutine zallocate_2(array,MM,NN)
    implicit none
    integer :: allocate_status, MM,NN
    double complex, allocatable :: array(:,:)

        allocate(array(1:MM,1:NN), stat = allocate_status)
        if (allocate_status .ne. 0) stop "***Not enough memory to allocate arrays***"

    end subroutine


    subroutine rallocate_1(array,MM)
    implicit none
    integer :: allocate_status, MM
    real*8, allocatable :: array(:)

        allocate(array(1:MM), stat = allocate_status)
        if (allocate_status .ne. 0) stop "***Not enough memory to allocate arrays***"

    end subroutine


    subroutine rallocate_2(array,MM,NN)
    implicit none
    integer :: allocate_status, MM,NN
    real*8, allocatable :: array(:,:)

        allocate(array(1:MM,1:NN), stat = allocate_status)
        if (allocate_status .ne. 0) stop "***Not enough memory to allocate arrays***"

    end subroutine


    subroutine diagonalization(org_hermitian, eigenvalues, NN)
    ! Obtains the e-values of an hermitian matrix
    ! Modifies org_hermitian
    implicit none
    integer :: NN, ok, lwork
    complex*16 :: org_hermitian(1:NN,1:NN)
    real*8 :: eigenvalues(1:NN)
    complex*16, allocatable :: work(:)
    real*8, allocatable :: rwork(:)
    
        ! we don't want to saturate the memory
        call zallocate_1(work,2*NN-1)

        ! Obtain an optimal lwork
        call zheev('V','U',NN,org_hermitian,NN,eigenvalues,work,-1,rwork,ok)
    
        ! the first element of work is the optimal lwork
        lwork = int(work(1))
        if (lwork.le.0) stop "overflow in lwork"
        ! Resize the work vector and validate memory
        deallocate(work)
        call zallocate_1(work,lwork)

        call rallocate_1(rwork,3*NN-2)
   
        call zheev('V','U',NN,org_hermitian,NN,eigenvalues,work,lwork,rwork,ok)
    
        ! Releasing memory
        deallocate(work,rwork)
    
        if (ok.ne.0) stop "Unsuccessful eigenvalue calculation"
    
    end subroutine


    subroutine print_mat(mat_a)
    implicit none
    double complex, dimension(:,:) :: mat_a
    integer :: ii, shape_a(1:2)
        shape_a = shape(mat_a)
        do ii = 1, shape_a(1)
            write(*,'(*(F0.8,SP,F0.8,"i",","))') mat_a(ii,:)
        end do
        write(*,*) "---------------------------------------"
    end subroutine


    subroutine compare_mat(mat_a,mat_b)
    implicit none
    double complex, dimension(:,:) :: mat_a,mat_b
    integer :: ii, jj, shape_a(1:2), shape_b(1:2)
        shape_a = shape(mat_a)
        shape_b = shape(mat_b)
        do ii = 1,2
            if (shape_a(ii).ne.shape_b(ii)) stop "Matrices differ in shape"
        end do
        do jj = 1, shape_a(2)
            do ii = 1, shape_a(1)
                if(&
                    (real(real(mat_a(ii,jj))).ne.real(real(mat_b(ii,jj))))&
                    .or.&
                    (real(aimag(mat_a(ii,jj))).ne.real(aimag(mat_b(ii,jj))))&
                    ) then
                    write(*,'(A,I0)') "matrix A, row ", ii
                    call print_mat(mat_a(ii:ii,:))
                    write(*,'(A,I0)') "matrix B, row ", ii
                    call print_mat(mat_b(ii:ii,:))
                    write(*,'(A,I0)') "different in element ",jj
                    stop
                end if
            end do
        end do    
    end subroutine


    function tens_prod_2(A,B,MM,NN,OO,PP)
        implicit none
        double complex, dimension(:,:) :: A,B
        double complex, dimension(1:MM*OO,1:NN*PP) :: tens_prod_2
        integer :: ii, jj, MM, NN, OO, PP, base(1:MM), out(1:MM)
        logical :: mask(1:MM)
    
            base = (/(ii,ii=1,MM)/)
            tens_prod_2 = 0
            do jj = 1, NN
                mask = A(:,jj)/=0
                out = pack(base,mask)
                do ii = 1, count(mask)
                    !if(A(ii,jj).eq.0) continue
                    tens_prod_2((out(ii)-1)*OO+1:out(ii)*OO,&
                    (jj-1)*PP+1:jj*PP) = A(out(ii),jj)*B
                end do
            end do
    end function

    
    function identity(NN)
    implicit none
    integer :: NN, ii
    double complex, dimension(1:NN,1:NN) :: identity
    
        identity = 0
        do ii = 1, NN
            identity(ii,ii) = dcmplx(1,0)
        end do
    
    end function


    function diag_val(OO,val)
    implicit none
    integer :: OO, ii
    double complex :: val, diag_val(1:OO,1:OO)
    
        diag_val = 0
        !if (val.eq.dcmplx(0,0)) return
        do ii = 1, OO
            diag_val(ii,ii) = val
        end do

    end function
    
    
    function tens_id_2(mat_a,MM,NN,id_nn)
    ! Obtains the tensor product of a matrix A and the id_nn identity matrix
    implicit none
    double complex, dimension(:,:) :: mat_a
    double complex, dimension(1:MM*id_nn,1:NN*id_nn) :: tens_id_2
    integer :: ii, jj, MM, NN, id_nn, base(1:MM), out(1:MM)
    logical :: mask(1:MM)
        
        base = (/(ii,ii=1,MM)/)
        tens_id_2 = 0
        do jj = 1, NN
            mask = mat_a(:,jj) /= 0
            out = pack(base,mask)
            do ii = 1, count(mask)
                !if(mat_a(ii,jj).eq.0) continue
                tens_id_2((out(ii)-1)*id_nn+1:out(ii)*id_nn,&
                (jj-1)*id_nn+1:jj*id_nn) = diag_val(id_nn,mat_a(out(ii),jj))
            end do
        end do

    end function


    function outer(vec_1,vec_2,MM,NN)
    implicit none
    integer :: MM, NN, ii, jj
    double complex :: vec_1(1:MM), vec_2(1:NN), outer(1:MM,1:NN)

        do jj = 1, NN
            do ii = 1, MM
                outer(ii,jj) = vec_1(ii)*conjg(vec_2(jj))
            end do
        end do

    end function


    function ising_model_1D(NN,lambda) result(hamiltonian)
    ! Builds hamiltonian for the 1/2 spin chain with nearest-neighbours interaction
    implicit none
    double complex, dimension(1:2,1:2) :: sigma_x, sigma_y, sigma_z
    double complex, dimension(1:2**2,1:2**2) :: sigma_int
    integer ii, NN
    double complex, dimension(:,:), allocatable :: hamiltonian,holder_int
    real*8 :: lambda

        sigma_x = dcmplx(reshape((/0,1,&
                                   1,0/),&
                            (/2,2/)),0)
        
        sigma_y = dcmplx(0,reshape((/0,-1,&
                                     1,0/),&
                            (/2,2/)))
        
        sigma_z = dcmplx(reshape((/1,0,&
                                   0,-1/),&
                            (/2,2/)),0)
        
        sigma_int = tens_prod_2(sigma_x,sigma_x,2,2,2,2)
        call zallocate_2(hamiltonian,2**NN,2**NN)
        call zallocate_2(holder_int,2**NN,2**NN)
        
        ! Transverse field interaction
        ! Setting everything for ii = 1
        holder_int = tens_id_2(&
            sigma_z,&
            2,2,2**(NN-1))
        
        hamiltonian = holder_int
        
        
        ! Setting for states inbetween
        do ii = 2, NN-1
            holder_int(1:2**ii,1:2**ii) = tens_prod_2(&
                    identity(2**(ii-1)),&
                    sigma_z,&
                    2**(ii-1),2**(ii-1),2,2)
            holder_int = tens_id_2(&
            holder_int(1:2**ii,1:2**ii),&
            2**ii,2**ii,2**(NN-ii))
            hamiltonian = hamiltonian + holder_int
        end do
        
        ! Setting everything for ii = NN
        holder_int = tens_prod_2(&
            identity(2**(NN-1)),&
            sigma_z,&
            2**(NN-1),2**(NN-1),2,2)
        
        hamiltonian = hamiltonian + holder_int
        
        hamiltonian = lambda*hamiltonian
        
        ! Nearest neighbours interaction
        ! Prepare interaction of the first two
        holder_int = tens_id_2(&
                        sigma_int,&
                        2**2,2**2,2**(NN-2))
        
        hamiltonian = hamiltonian - holder_int
        
        do ii = 2, NN-2
            holder_int(1:2**(ii+1),1:2**(ii+1)) = tens_prod_2(&
                            identity(2**(ii-1)),&
                            sigma_int,&
                            2**(ii-1),2**(ii-1),2**2,2**2)
            holder_int = tens_id_2(&
                            holder_int(1:2**(ii+1),1:2**(ii+1)),&
                            2**(ii+1),2**(ii+1),2**(NN-ii-1))
            hamiltonian = hamiltonian - holder_int
        end do
        
        ! Prepare interaction of the last two
        holder_int = tens_prod_2(&
                        identity(2**(NN-2)),&
                        sigma_int,&
                        2**(NN-2),2**(NN-2),2**2,2**2)
        
        hamiltonian = hamiltonian - holder_int
        
    end function


    subroutine allocate_group(holder_ham,H_int_A,H_int_B,eigenvalues,projector_0,red_proj,curr_size,M_size)
    implicit none
    integer :: curr_size, M_size
    double complex, allocatable :: holder_ham(:,:), projector_0(:,:), red_proj(:,:), H_int_A(:,:), H_int_B(:,:)
    real*8, allocatable :: eigenvalues(:)

        call zallocate_2(holder_ham,curr_size,curr_size)
        call zallocate_2(H_int_A,curr_size,curr_size)
        call zallocate_2(H_int_B,curr_size,curr_size)
        call rallocate_1(eigenvalues,curr_size)
        call zallocate_2(projector_0,curr_size,curr_size)
        call zallocate_2(red_proj,curr_size,M_size)

    end subroutine


    subroutine real_space_rg(&
        N_size,&! Number of particles described by the initial hamiltonian
        d_states,&! Number of states for the individual systems
        hamiltonian_N,&! The analytical hamiltonian for the N particle system
        boundary_A,&! The leftmost boundary interaction term
        boundary_B,&! The rightmost boundary interaction term
        M_size,&! dimension of the truncated hamiltonian
        iterations,& ! using instead of target size
        !target_size,&! number of particles of the final system. Must be N_size*2**w_it for some non negative integer w_it
        output_hamiltonian&
        )
    implicit none
    integer :: N_size, M_size, d_states, iterations, curr_size, ii, jj
    double complex :: output_hamiltonian(1:M_size,1:M_size), H_int(1:M_size**2,1:M_size**2)
    double complex :: hamiltonian_N(1:d_states**N_size,1:d_states**N_size), hamiltonian_M(1:M_size**2,1:M_size**2)
    double complex :: boundary_A(1:d_states,1:d_states), boundary_B(1:d_states,1:d_states)
    double complex :: H_int_A(1:M_size,1:M_size), H_int_B(1:M_size,1:M_size)
    double complex, allocatable :: holder_ham(:,:), projector_0(:,:), red_proj(:,:)
    double complex, allocatable :: H_A(:,:), H_B(:,:)
    real*8, allocatable :: eigenvalues(:)

        curr_size = d_states**N_size
        !output_hamiltonian = 0
        ! The renormalization group dimension must be lower or equal to the initial configuration dimension.
        if (curr_size.lt.M_size) stop "Dimension of hilbert space for initial configuration too small."

        call allocate_group(holder_ham,H_A,H_B,eigenvalues,projector_0,red_proj,curr_size,M_size)
        holder_ham = hamiltonian_N

        ! Diagonalize the initial hamiltonian
        call diagonalization(holder_ham, eigenvalues, curr_size)

        ! Project on the basis for the first M eigenstates
        projector_0 = holder_ham
        !projector_0 = 0
        !do jj = 1, M_size
        !    projector_0 = projector_0 + outer(holder_ham(:,jj),holder_ham(:,jj),curr_size,curr_size)
        !end do
        red_proj = projector_0(1:curr_size,1:M_size)
        write(*,*) eigenvalues(1:M_size)
        output_hamiltonian = matmul(mat_adj(red_proj),matmul(hamiltonian_N,red_proj))
        H_A = tens_prod_2(&
                identity(d_states**(N_size-1)),&
                boundary_A,&
                d_states**(N_size-1),d_states**(N_size-1),d_states,d_states)
        H_B = tens_id_2(&
                boundary_B,&
                d_states,d_states,d_states**(N_size-1))
        H_int_A = matmul(mat_adj(red_proj),matmul(H_A,red_proj))
        H_int_B = matmul(mat_adj(red_proj),matmul(H_B,red_proj))

        ! Build interaction hamiltonian
        H_int = tens_prod_2(&
            H_int_A,&
            H_int_B,&
            M_size,M_size,M_size,M_size)

        ! Resize the dummy variables
        curr_size = M_size*M_size
        deallocate(holder_ham,H_A,H_B,eigenvalues,projector_0,red_proj)
        call allocate_group(holder_ham,H_A,H_B,eigenvalues,projector_0,red_proj,curr_size,M_size)

        ! repeat until the desired number of particles is reached
        do ii = 1, iterations

            ! build the 2N-particle hamiltonian as tens(H,id)+tens(id,H)+H_int
            hamiltonian_M = tens_id_2(output_hamiltonian,M_size,M_size,M_size)&
            + tens_prod_2(identity(M_size),output_hamiltonian,M_size,M_size,M_size,M_size)&
            + H_int! Projected interaction term

            ! diagonalize N-particle system hamiltonian
            holder_ham = hamiltonian_M
            call diagonalization(holder_ham, eigenvalues, curr_size)

            ! project only the first m e-values
            projector_0 = holder_ham
            !projector_0 = 0
            !do jj = 1, M_size
            !    projector_0 = projector_0 + outer(holder_ham(:,jj),holder_ham(:,jj),curr_size,curr_size)
            !end do
            red_proj = projector_0(1:curr_size,1:M_size)
            output_hamiltonian = matmul(mat_adj(red_proj),matmul(hamiltonian_M,red_proj))
            ! Build interaction hamiltonian
            ! Fixed interaction terms
            H_A = tens_prod_2(&
                    identity(M_size),&
                    H_int_A,&
                    M_size,M_size,M_size,M_size)
            H_B = tens_id_2(&
                    H_int_B,&
                    M_size,M_size,M_size)
            H_int_A = matmul(mat_adj(red_proj),matmul(H_A,red_proj))
            H_int_B = matmul(mat_adj(red_proj),matmul(H_B,red_proj))
            H_int = tens_prod_2(&
                H_int_A,&
                H_int_B,&
                M_size,M_size,M_size,M_size)

        end do

        ! Releasing memory
        deallocate(holder_ham,eigenvalues,projector_0,red_proj)

    end subroutine


end module

program bullshit
use renormalization
implicit none
integer :: NN, MM, iterations
double complex, allocatable :: hamiltonian_N(:,:), hamiltonian_M(:,:)
double complex :: sigma_x(1:2,1:2)
real*8 :: lambda
real*8, allocatable :: eigenvalues(:)

sigma_x = dcmplx(reshape((/0,1,&
                           1,0/),&
                    (/2,2/)),0)


write(*,*) "Enter NN, MM, lambda, iterations"
read(*,*) NN, MM, lambda, iterations

call zallocate_2(hamiltonian_N,2**NN,2**NN)
call zallocate_2(hamiltonian_M,MM,MM)
call rallocate_1(eigenvalues,MM)

hamiltonian_N = ising_model_1D(NN,lambda)

call real_space_rg(NN,2,hamiltonian_N,sqrt(lambda)*sigma_x,-sqrt(lambda)*sigma_x,MM,iterations,hamiltonian_M)

call diagonalization(hamiltonian_M,eigenvalues,MM)

write(*,*) eigenvalues

end program