subroutine rdf_3d(N, x, y, z, Lx, Ly, Lz, snaps, avg, r)
    implicit none 
 
    integer, intent(in) :: N
    integer, intent(in) :: snaps  
    real(8), intent(in) :: Lx,Ly, Lz 
    real(8), intent(in) :: x(snaps,N), y(snaps,N), z(snaps,N)
    integer, parameter :: Nhis = 2.**8. 
    real(8), intent(out) :: avg(Nhis),r(Nhis)
    integer i, j, k, ig
    DOUBLE PRECISION::rr,delg,pi,xr,yr,zr,r2,vb,nid,rho
    DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr

    delg = Lx / (2.* Nhis)
    rho = N / (Lx * Ly * Lz)
    pi = 4. * ATAN(1.)
    gr = 0.0d0
    avg(:) = 0.d0

    do k = 1, snaps 
        do i = 1, N - 1 
            do j = i+1, N
                xr = x(k,i) - x(k,j)
                yr = y(k,i) - y(k,j)
                zr = z(k,i) - z(k,j)

                xr = xr - Lx * (NINT(xr / Lx))
                yr = yr - Ly * (NINT(yr / Ly))
                zr = zr - Lz * (NINT(zr / Lz))

                r2 = xr * xr + yr * yr + zr * zr 
                rr = sqrt(r2)
                if (rr .lt. Lx / 2.0) then 
                    ig = ceiling(rr / delg)
                    gr(k,ig) = gr(k,ig) + 2.
                endif 
            enddo 
        enddo 
    enddo 

    DO j=1,Nhis
        DO i=1,snaps
           r(j)=delg*(j+0.5)
           vb=((j+1)**3.-j**3.)*delg**3.
           nid = (4./3.)*pi*vb*rho
           gr(i,j)=gr(i,j)/(N*nid)
           avg(j)=avg(j)+gr(i,j)/snaps/N
        END DO
    END DO



end subroutine