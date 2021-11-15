program checkpoints_test
implicit none

interface
    subroutine checkpoint(debug,variable,message)
        logical :: debug
        real, optional :: variable
        character (len = *), optional :: message
    end subroutine
end interface

logical :: debug
real :: variable = 0.0
character (len = *), parameter :: message = "Message placeholder"

debug = .True.
write(*,*) "---Debugging without optional arguments---"
call checkpoint(debug)
write(*,*) "---Debugging with message---"
call checkpoint(debug, message = message)
write(*,*) "---Debugging with variable---"
call checkpoint(debug, variable = variable)
write(*,*) "---Debugging with both---"
call checkpoint(debug, variable = variable, message = message)
write(*,*) "---Not debugging---"
call checkpoint(.False.)



end program checkpoints_test


subroutine checkpoint(debug,variable,message)
implicit none
logical :: debug
real, optional :: variable
character (len = *), optional :: message

if (debug .eqv. .True.) then
    write(*,*) "***DEBUG MODE***"
    if (present(message)) then
        write(*,*) "Message:",message
    end if
    if (present(variable)) then
        write(*,*) "Variable:",variable
    end if
end if

end subroutine