subroutine print_unix(string, var1, var2)

    character(LEN=1024) :: string
    integer, optional :: var1
    integer, optional :: var2

    if (present(var2)) then
        write(*,*) trim(string), var1, var2
    elseif(present(var1)) then
        write(*,*) trim(string), var1, var2
    else
        write(*,*) trim(string)
    endif
end subroutine print_unix
