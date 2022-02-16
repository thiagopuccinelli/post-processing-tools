subroutine print_dimer_center_of_mass2d(N,Nmolecule,x,y,snaps,fileout)
    implicit none 
    integer, intent(in) :: N 
    integer, intent(in) :: Nmolecule
    integer, intent(in) :: snaps
    real(8), intent(in) :: x(snaps,N), y(snaps,N)
    integer i,j,k 
    real(8) :: xcm, ycm
    character(128), intent(in) :: fileout 

    open(unit = 2, file = fileout, action="write")
    do j = 1, snaps 
        write(2,*) Nmolecule
        write(2,*) "Atoms. Timestep: ", j
        do i = 1, N - 1, 2 
            xcm = (x(j,i) + x(j,i+1)) / 2.
            ycm = (y(j,i) + y(j,i+1)) / 2.
            write(2,*) 1, xcm, ycm, 0.0 
        enddo 
    enddo 

    close(2)

end subroutine
    