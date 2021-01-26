subroutine error_handling()

    external f_abort
    integer :: f_abort 
    integer :: iret1, iret2, procnum
    integer, parameter :: SIGABRT = 6
    integer, parameter :: SIGSEGV =174

    iret1 = SIGNAL(SIGSEGV, f_abort)
    write(*,*) 'Set signal handler #1. Return = ', iret1
    iret2 = kill(getpid(), SIGSEGV)
    write(*,*) 'Raised signal. Return = ', iret2

end subroutine error_handling

function f_abort(sig_num) result(abort)

    use MessageHandling
    use iso_c_utils

    integer, intent(in) :: sig_num
    integer :: abort

    call mess(LEVEL_ERROR, 'In signal handler function h_abort for SIG$ABORT')
    WRITE(*,*) 'signum = ', sig_num
    abort = 1

end function f_abort
