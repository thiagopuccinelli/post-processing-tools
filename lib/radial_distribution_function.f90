subroutine radial_distribution_one_particle_species(N,filename,snaps,LL,avg,r) 
   implicit none 

   integer, intent(in) :: N
   character(64), intent(in) :: filename 
   integer, intent(in) :: snaps  
   real(8), intent(in) :: LL
   integer, parameter :: Nhis = 2.**8. 
   real(8), parameter :: time = 1.0
   real(8), intent(out) :: avg(Nhis),r(Nhis)
   real(8) :: x(snaps,N), y(snaps,N), z(snaps,N), vx(snaps,N), vy(snaps,N), vz(snaps,N)
   integer :: type(snaps,N), id(snaps,N)
   integer i, j, k,ig,N1
   integer dumb1,dumb2,dumb3
   DOUBLE PRECISION::rr,delg,pi,xr,yr,zr,r2,vb,nid,rho
   DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr
        
   delg=LL/(2.*Nhis)
   pi=4*ATAN(1.)
   N1 = N
   rho = N1/(LL**3.0d0)
   gr = 0.0d0
   avg(:)=0.d0

   open(10, file=filename)

   do j = 1, snaps 
      read(10,*) dumb1 
      read(10,*) dumb2 
      read(10,*) dumb3
      do i = 1, N 
         read(10,*) id(j,i), type(j,i), x(j,i), y(j,i), z(j,i),vx(j,i), vy(j,i), vz(j,i)
      enddo 
   enddo 

   DO k=1,snaps
      DO i=1,N-1
         DO j=i+1,N
            xr=x(k,i)-x(k,j)
            yr=y(k,i)-y(k,j)
            zr=z(k,i)-z(k,j)
               
            xr=xr-LL*(NINT(xr/LL))
            yr=yr-LL*(NINT(yr/LL))
            zr=zr-LL*(NINT(zr/LL))
            r2=xr*xr+yr*yr+zr*zr
            rr=SQRT(r2)

            IF(rr.LT.LL/2.d0)THEN
               ig=ceiling(rr/delg)
               gr(k,ig)=gr(k,ig)+2.
            END IF
         END DO
      END DO
   END DO




   DO j=1,Nhis
      DO i=1,snaps
         r(j)=delg*(j+0.5)
         vb=((j+1)**3.-j**3.)*delg**3.
         nid=(4./3.)*pi*vb*rho
         gr(i,j)=gr(i,j)/(N1*nid)
      END DO
   END DO





   DO i=1,Nhis
      DO j=1,snaps
         avg(i)=avg(i)+gr(j,i)
      END DO
   END DO

   OPEN(unit=2,file='rdf_one_species.dat',action="write")

   DO i=1,Nhis
      WRITE(2,'(2(f17.10,1X))')r(i),avg(i)/snaps
      avg(i) = avg(i)/snaps
   END DO

end subroutine 

