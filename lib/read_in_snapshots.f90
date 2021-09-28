subroutine read_in_trajectory(filein, N, snaps, r)
    implicit none 
    integer :: i, j
    character(128), intent(in) :: filein
    integer, intent(in) :: N, snaps  
    integer, dimension(snaps,N) :: type 
    real(8),intent(out),dimension(snaps,N,3) :: r 
    character(64) :: dumb1,dumb2 

    open(10, file = filein)
    do j = 1, snaps 
        read(10,*) dumb1 
        read(10,*) dumb2 
        do i = 1, N 
            read(10,*) type(j,i), r(j,i,1), r(j,i,2), r(j,i,3)
        enddo 
    enddo 

    close(10)

end subroutine 