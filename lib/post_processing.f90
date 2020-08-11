! MODULE FOR SOME USEFULL COMPUTATION IN MOLECULAR DYNAMICS 
! POST PROCESSING RESULTS PRODUCTION FROM SNAPSHOTS 
!
! This module treats LAMMPS .xyz file. 

! 2d standard rdf production for one style particle: 
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
! 3d standard rdf production for one style particle:
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

! 2d dimer simulation:
subroutine radial_distribution_center_of_mass(N,filename,snaps,Lx,Ly,avg,r) 
    implicit none 
 
    integer, intent(in) :: N
    character(64), intent(in) :: filename 
    integer, intent(in) :: snaps  
    real(8), intent(in) :: Lx,Ly
    integer, parameter :: Nhis = 2.**8. 
    real(8), parameter :: time = 1.0
    real(8), intent(out) :: avg(Nhis),r(Nhis)
    real(8) :: x(snaps,N), y(snaps,N), z(snaps,N), xcm(snaps,N), ycm(snaps,N), zcm(snaps,N)
    integer :: type(snaps,N), id(snaps,N)
    integer i, j, k,ig,N1,nd
    character(64) :: dumb1
    character(64) :: dumb2
    DOUBLE PRECISION::rr,delg,pi,xr,yr,zr,r2,vb,nid,rho
    DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr

    delg=10.0/(Nhis)
    pi=4*ATAN(1.)
    N1 = N/2
    rho = N1/((Lx*Ly)**2.0d0)
    gr = 0.0d0
    avg(:)=0.d0
    
    open(10,file=trim(filename)//'.snap')
 
    do j = 1, snaps 
        read(10,*) dumb1
        read(10,*) dumb2
        do i = 1, N 
           read(10,*) type(j,i), x(j,i), y(j,i), z(j,i)
        enddo
        nd = 0 
        do i=1, N , 2 
           nd = nd + 1 
           xcm(j,nd) = (x(j,i+1) + x(j,i))/2.
           ycm(j,nd) = (y(j,i+1) + y(j,i))/2.
           zcm(j,nd) = (z(j,i+1) + z(j,i))/2.
        enddo 
    enddo 
 
    close(10)

    DO k=1,snaps
        DO i=1,N1 - 1 
           DO j=i+1,N1 
           xr=xcm(k,i)-xcm(k,j)
           yr=ycm(k,i)-ycm(k,j)
           
           
           xr=xr-Lx*(NINT(xr/Lx))
           yr=yr-Ly*(NINT(yr/Ly))
           
           r2=xr*xr+yr*yr
           rr=SQRT(r2)
  
           IF(rr.LT.10.d0)THEN
                 ig=ceiling(rr/delg)
                 gr(k,ig)=gr(k,ig)+2.
           END IF
           END DO
        END DO
    END DO
  
     
  
     
    DO j=1,Nhis
        DO i=1,snaps
           r(j)=delg*(j+0.5)
           vb=((j+1)**2.-j**2.)*delg**2.
           nid = pi*vb/(Lx*Ly)
           gr(i,j)=gr(i,j)/(N1*nid)
        END DO
    END DO
  
     
     
  
     
    DO i=1,Nhis
        DO j=1,snaps
           avg(i)=avg(i)+gr(j,i)
        END DO
    END DO
     
    OPEN(unit=2,file='rdf_'//trim(filename)//'_cm.dat',action="write")
     
    DO i=1,Nhis
        WRITE(2,'(2(f17.10,1X))')r(i),avg(i)/snaps/N1
    END DO

    close(2)
 
end subroutine 

subroutine radial_distribution_monomer_to_monomer(N,filename,snaps,Lx,Ly,avg,r)
    implicit none 
    integer, intent(in) :: N
    character(64), intent(in) :: filename 
    integer, intent(in) :: snaps  
    real(8), intent(in) :: Lx,Ly
    integer, parameter :: Nhis = 2.**8. 
    real(8), parameter :: time = 1.0
    real(8), intent(out) :: avg(Nhis),r(Nhis)
    real(8) :: x(snaps,N), y(snaps,N), z(snaps,N), xcm(snaps,N), ycm(snaps,N), zcm(snaps,N)
    integer :: type(snaps,N), id(snaps,N)
    integer i, j, k,ig,N1,nd
    character(64) :: dumb1
    character(64) :: dumb2
    DOUBLE PRECISION::rr,delg,pi,xr,yr,zr,r2,vb,nid,rho
    DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr


    delg=10.0/(Nhis)
    pi=4*ATAN(1.)
    N1 = N
    rho = N1/((Lx*Ly)**2.0d0)
    gr = 0.0d0
    avg(:)=0.d0
    
    open(10,file=trim(filename)//'.snap')
 
    do j = 1, snaps 
        read(10,*) dumb1
        read(10,*) dumb2
        do i = 1, N 
           read(10,*) type(j,i), x(j,i), y(j,i), z(j,i)
        enddo
    enddo     

    close(10)

    DO k=1,snaps
        DO i=1,N1 - 1 
           DO j=i+1,N1 
           xr=x(k,i)-x(k,j)
           yr=y(k,i)-y(k,j)
           
           
           xr=xr-Lx*(NINT(xr/Lx))
           yr=yr-Ly*(NINT(yr/Ly))
           
           r2=xr*xr+yr*yr
           rr=SQRT(r2)
  
           IF(rr<=10.0)THEN
                if ((-1.0)**i == 1.0) then !se i é par cacula sempre
                    ig=ceiling(rr/delg)
                    gr(k,ig)=gr(k,ig)+2.
                else if ((-1.0)**i == -1.0 .and. j /= i+1) then !se i é impar e j diferente de i+1
                    ig=ceiling(rr/delg)
                    gr(k,ig)=gr(k,ig)+2.
                endif
            END IF
           END DO
        END DO
    END DO
  
     
  
     
    DO j=1,Nhis
        DO i=1,snaps
           r(j)=delg*(j+0.5)
           vb=((j+1)**2.-j**2.)*delg**2.
           nid = pi*vb/(Lx*Ly)
           gr(i,j)=gr(i,j)/(N1*nid)
        END DO
    END DO
  
     
     
  
     
    DO i=1,Nhis
        DO j=1,snaps
           avg(i)=avg(i)+gr(j,i)
        END DO
    END DO
     
    OPEN(unit=2,file='rdf_'//trim(filename)//'_mon_mon.dat',action="write")
     
    DO i=1,Nhis
        WRITE(2,'(2(f17.10,1X))')r(i),avg(i)/snaps/N1
    END DO

    close(2)


end subroutine

subroutine print_center_of_mass_XYZ(N,filename,snaps,Lx,Ly,avg,r) 
    implicit none 
 
    integer, intent(in) :: N
    character(64), intent(in) :: filename 
    integer, intent(in) :: snaps  
    real(8), intent(in) :: Lx,Ly
    integer, parameter :: Nhis = 2.**8. 
    real(8), parameter :: time = 1.0
    real(8), intent(out) :: avg(Nhis),r(Nhis)
    real(8) :: x(snaps,N), y(snaps,N), z(snaps,N), xcm(snaps,N), ycm(snaps,N), zcm(snaps,N)
    integer :: type(snaps,N), id(snaps,N)
    integer i, j, k,ig,N1,nd
    character(64) :: dumb1
    character(64) :: dumb2
    DOUBLE PRECISION::rr,delg,pi,xr,yr,zr,r2,vb,nid,rho
    DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr

    N1 = N/2
    
    open(10,file=trim(filename)//'.snap')
 
    do j = 1, snaps 
        read(10,*) dumb1
        read(10,*) dumb2
        do i = 1, N 
           read(10,*) type(j,i), x(j,i), y(j,i), z(j,i)
        enddo
        nd = 0 
        do i=1, N , 2 
           nd = nd + 1 
           xcm(j,nd) = (x(j,i+1) + x(j,i))/2.
           ycm(j,nd) = (y(j,i+1) + y(j,i))/2.
           zcm(j,nd) = (z(j,i+1) + z(j,i))/2.
        enddo 
    enddo 

    close(10)
 
    OPEN(unit=2,file=trim(filename)//'_cm.snap',action="write")


    DO k=1,snaps
        write(2,*) 2000
        write(2,*) dumb2
        DO i=1,N1  
           write(2,*) 1, xcm(k,i),ycm(k,i), 0.0 
        END DO
    END DO
  
    close(2)
 
end subroutine 

! post processing radial distribution functions: 
! compute excess entropy 
subroutine excess_entropy_from_rdf(filein, fileout, sex)
    implicit none 
    integer :: i 
    character(64), intent(in) :: filein 
    character(64), intent(in) :: fileout  
    integer, parameter :: Nhis = 2.**8.
    real(8) :: gor(Nhis), r(Nhis), csum 
    real(8), intent(out) :: sex 

    sex = 0.d0
    csum = 0.d0 
    OPEN(unit=110,file=trim(fileout),action="write")
    OPEN(unit=120,file=trim(filein),action="read")

    do i = 1, Nhis 
        read(120,*) r(i), gor(i)
        if (gor(i) > 0.d0) then 
            csum = csum - (gor(i)*log(gor(i)) - gor(i) +1)*r(i) !for excess entropy
            if (gor(i) /= 0.0d0) sex = sex + (gor(i)*log(gor(i)) - gor(i) +1) !for excess entropy
            write(110,'(f12.6,x,f12.6)') r(i), -3.14*csum
        endif 
    enddo 

    close(110)
    close(120)

end subroutine
subroutine radial_distribution_3max(filein, max1, max1r, max2, max2r, max3, max3r)
    implicit none 
    integer :: i, j, k 
    character(64), intent(in) :: filein 
    integer, parameter :: Nhis = 2.**8.
    real(8) :: gor(Nhis), r(Nhis)
    real(8), intent(out) ::  max1, max1r, max2, max2r, max3, max3r 

    OPEN(unit=120,file=trim(filein),action="read")

    do i = 1, Nhis 
        read(120,*) r(i), gor(i)
    enddo 

    close(120)
    i = 3 
    do while (i <= Nhis-2)
        if (gor(i) > 1.0) then
            if ( ((gor(i+2)-gor(i+1))/(r(i+2)-r(i+1))) < 0.0d0 .and. ((gor(i-2)-gor(i-1))/(r(i-2)-r(i-1))) > 0.0d0 ) then 
                max1 = gor(i)
                max1r = r(i)
                j = i + 3 
                i = Nhis+1 
            endif 
        endif 
        i = i + 1 
    enddo 

    do while (j<= Nhis-2)
        if  (gor(j) > 1.0) then
            if ( ((gor(j+2)-gor(j+1))/(r(j+2)-r(j+1))) < 0.0d0 .and. ((gor(j-2)-gor(j-1))/(r(j-2)-r(j-1))) > 0.0d0 ) then
                max2 = gor(j)
                max2r = r(j)
                k = j + 3 
                j = Nhis + 1 
            endif 
        endif
        j = j + 1  
    enddo 

    do while (k<= Nhis-2)
        if  (gor(k) > 1.0) then
            if ( ((gor(k+2)-gor(k+1))/(r(k+2)-r(k+1))) < 0.0d0 .and. ((gor(k-2)-gor(k-1))/(r(k-2)-r(k-1))) > 0.0d0 ) then
                max3 = gor(k)
                max3r = r(k)
                k = Nhis + 1  
            endif 
        endif
        k = k + 1  
    enddo 


    

end subroutine