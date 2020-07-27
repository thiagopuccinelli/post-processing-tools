subroutine radial_distribution_2d(N,filename,snaps,LL,avg,r)
    implicit none 
    integer, intent(in) :: N 
    character(len=64), intent(in) :: filename
    integer, intent(in)  :: snaps 
    real(8), intent(in) :: LL 
    integer, parameter :: Nhis = 2.**8.
    real(8), intent(out) :: avg(Nhis), r(Nhis)
    real(8) :: x(snaps, N),y(snaps, N),z(snaps, N)
    integer :: type(snaps,N)
    integer :: i, j, k, ig 
    integer :: dumb1 
    character(len=128) :: dumb2 
    DOUBLE PRECISION::rr,delg,pi,xr,yr,zr,r2,vb,nid,rho
    DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr   

	delg = LL / (2.* Nhis)
	rho = N / (LL ** 2.)
	pi=4*ATAN(1.)
	gr = 0.0d0 
	avg(:) = 0.d0
	
	open(10, file=filename)

    do j = 1, snaps
        read(10,*) dumb1 
        read(10,*) dumb2 
        do i = 1, N 
            read(10,*) type(j,i), x(j,i), y(j,i), z(j,i)
        enddo 
	enddo 
	
	close(10)

	do k = 1, snaps 
		do i = 1, N-1 
			do j = i + 1, N 
				xr = x(k,i) - x(k,j)
				yr = y(k,i) - y(k,j)
				zr = z(k,i) - z(k,j)

				xr = xr - LL * (NINT(xr / LL))
				yr = yr - LL * (NINT(yr / LL))
				zr = zr - LL * (NINT(zr / LL))
				r2 = xr * xr + yr * yr 
				rr = sqrt(r2)

				if (rr .lt. LL/2.d0) then 
					ig = ceiling(rr/delg)
					gr(k,ig) = gr(k,ig) + 2. 
				endif 
			enddo 
		enddo 
	enddo 


    DO j=1,Nhis
        DO i=1,snaps
           r(j)=delg*(j+0.5)
           vb=((j+1)**2.-j**2.)*delg**2.
           nid = pi*vb*rho
           gr(i,j)=gr(i,j)/(N*nid)
        END DO
    END DO

	DO i=1,Nhis
        DO j=1,snaps
           avg(i)=avg(i)+gr(j,i)
        END DO
    END DO

	OPEN(unit=2,file='rdf_'//trim(filename)//'.dat',action="write")
	
    DO i=1,Nhis
        WRITE(2,'(2(f17.10,1X))')r(i),avg(i)/snaps
    END DO

	close(2)

	DO i=1,Nhis
        avg(i) = avg(i)/snaps
    END DO


end subroutine

subroutine msd_evaluation(N,filename,snaps,LL,time,msdvec,timevec)
    implicit none 
    integer, intent(in) :: N 
    character(len=64), intent(in) :: filename
    integer, intent(in)  :: snaps 
    real(8), intent(in) :: LL 
    integer, parameter :: Nhis = 2.**8.
    real(8) :: x(snaps, N),y(snaps, N),z(snaps, N), ix(snaps,N), iy(snaps,N), iz(snaps,N)
    integer :: type(snaps,N)
    integer :: i, j, k
    integer :: dumb1 
    character(len=128) :: dumb2 
	real(8) :: dr(snaps, N), x0,y0, z0, msd
    real(8), intent(out) :: msdvec(snaps), timevec(snaps)
    real(8), intent(in) :: time  !number of steps between each measure x time step

	open(10, file=filename)

    do j = 1, snaps
        read(10,*) dumb1 
        read(10,*) dumb2
        do i = 1, N 
            read(10,*) type(j,i),x(j,i), y(j,i), z(j,i)
        enddo 
	enddo 

	close(10)

	dr = 0.0d0
     do i=2,snaps
        do j=1,N
           x0 = x(i-1,j)
           y0 = y(i-1,j) 
           z0 = z(i-1,j) 
           x(i,j) = (x(i,j)-x0)
           x(i,j) = x(i,j) - anint((x(i,j))/LL)*LL
           x(i,j) = (x(i,j)+x0)
           y(i,j) = (y(i,j)-y0)
           y(i,j) = y(i,j) - anint((y(i,j))/LL)*LL
           y(i,j) = (y(i,j)+y0)
           z(i,j) = (z(i,j)-z0)
           z(i,j) = z(i,j) - anint((z(i,j))/LL)*LL
           z(i,j) = (z(i,j)+z0)
           dr(i,j) = (x(i,j)-x0)**2.0d0+(y(i,j)-y0)**2.0d0+(z(i,j)-z0)**2.0d0
        enddo
     enddo

	 OPEN(unit=2,file='msd'//trim(filename)//'.dat',action="write")
     msdvec = 0.0d0
     timevec = 0.0d0
     msd = 0.0d0 
     do j=2,snaps
        do i=1, N
           msd = msd + dr(j,i)
        enddo
        write(2,*) (j-2)*time , msd/N
	 enddo
	 close(2)
	 msdvec = 0.0d0
     timevec = 0.0d0
     msd = 0.0d0 
     do j=2,snaps
        do i=1, N
           msd = msd + dr(j,i)
        enddo
        msdvec(j) = msd/N
        timevec(j) = (j-2) * time 
     enddo
end subroutine 
subroutine radial_distribution_3d(N,filename,snaps,LL,avg,r)
    implicit none 
    integer, intent(in) :: N 
    character(len=64), intent(in) :: filename
    integer, intent(in)  :: snaps 
    real(8), intent(in) :: LL 
    integer, parameter :: Nhis = 2.**8.
    real(8), intent(out) :: avg(Nhis), r(Nhis)
    real(8) :: x(snaps, N),y(snaps, N),z(snaps, N)
    integer :: type(snaps,N)
    integer :: i, j, k, ig 
    integer :: dumb1 
    character(len=128) :: dumb2 
    DOUBLE PRECISION::rr,delg,pi,xr,yr,zr,r2,vb,nid,rho
    DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr   

	delg = LL / (2.* Nhis)
	rho = N / (LL ** 3.0d0)
	pi=4*ATAN(1.)
	gr = 0.0d0 
	avg(:) = 0.d0
	
	open(10, file=filename)

    do j = 1, snaps
        read(10,*) dumb1 
        read(10,*) dumb2 
        do i = 1, N 
            read(10,*) type(j,i), x(j,i), y(j,i), z(j,i)
        enddo 
	enddo 
	
	close(10)

	do k = 1, snaps 
		do i = 1, N-1 
			do j = i + 1, N 
				xr = x(k,i) - x(k,j)
				yr = y(k,i) - y(k,j)
				zr = z(k,i) - z(k,j)

				xr = xr - LL * (NINT(xr / LL))
				yr = yr - LL * (NINT(yr / LL))
				zr = zr - LL * (NINT(zr / LL))
				r2 = xr * xr + yr * yr + zr * zr 
				rr = sqrt(r2)

				if (rr .lt. LL/2.d0) then 
					ig = ceiling(rr/delg)
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
        END DO
    END DO

	DO i=1,Nhis
        DO j=1,snaps
           avg(i)=avg(i)+gr(j,i)
        END DO
    END DO

	OPEN(unit=2,file='rdf_'//trim(filename)//'.dat',action="write")
	
    DO i=1,Nhis
        WRITE(2,'(2(f17.10,1X))')r(i),avg(i)/snaps
    END DO

	close(2)

	DO i=1,Nhis
        avg(i) = avg(i)/snaps
    END DO


end subroutine
