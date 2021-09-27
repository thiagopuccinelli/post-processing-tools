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

subroutine radial_distribution_close_to_particle(N,filename,snaps,LL,density,avgc,rc,avgf,rf)
   implicit none 

   integer, intent(in) :: N
   character(64), intent(in) :: filename 
   integer, intent(in) :: snaps  
   real(8), intent(in) :: LL
   character(64), intent(in) :: density
   character(10) ::  flag
   integer, parameter :: Nhis = 2.**8. 
   real(8), parameter :: time = 1.0
   real(8), intent(out) :: avgc(Nhis),rc(Nhis),avgf(Nhis),rf(Nhis)

   real(8) :: xc(snaps,N), yc(snaps,N), zc(snaps,N), xcc(snaps,N), ycc(snaps,N), zcc(snaps,N)
   real(8) :: xcf(snaps,N), ycf(snaps,N), zcf(snaps,N)
   real(8) :: xp(snaps,N), yp(snaps,N), zp(snaps,N)
   real(8) :: x(snaps,N), y(snaps,N), z(snaps,N), vx(snaps,N), vy(snaps,N), vz(snaps,N)
   integer type(snaps,N), id(snaps, N), ncc(snaps), ncf(snaps), mark(snaps,N)
   integer dumb1, dumb2, dumb3
   integer i,j,k, nc, np, nncc, nncf, ig,N1
   DOUBLE PRECISION::rr,xr,yr,zr,r2,delg,pi,vb,nid,rho
   DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr

   open(10, file=filename)
   
   do j=1,snaps
      read (10,*) dumb1
      read (10,*) dumb2
      read (10,*) dumb3
      nc = 0
      np = 0
      do i= 1, N
         read(10,*) id(j,i), type(j,i), x(j,i), y(j,i), z(j,i), vx(j,i), vy(j,i), vz(j,i) 
         if (type(j,i) == 1) then
            np = np+1
            xp(j,np) = x(j,i)*LL
            yp(j,np) = y(j,i)*LL
            zp(j,np) = z(j,i)*LL
         else
            nc = nc+1
            xc(j,nc) = x(j,i)*LL
            yc(j,nc) = y(j,i)*LL
            zc(j,nc) = z(j,i)*LL
         endif
      enddo
   enddo


   ncc = 0
   ncf = 0
   nncf = 0
   nncc = 0
   mark = 0
   do k = 1, snaps
      do i = 1, np
         do j = 1, nc
            xr=xp(k,i)-xc(k,j)
            yr=yp(k,i)-yc(k,j)
            zr=zp(k,i)-zc(k,j)
            
            xr=xr-LL*(NINT(xr/LL))
            yr=yr-LL*(NINT(yr/LL))
            zr=zr-LL*(NINT(zr/LL))
            
            r2=xr*xr+yr*yr+zr*zr
            rr=SQRT(r2)
            
            if (rr < (2.0**(1./6.)) .and. mark(k,j) == 0) then
               nncc = nncc+1
               xcc(k,nncc) = xc(k,j)
               ycc(k,nncc) = yc(k,j)
               zcc(k,nncc) = zc(k,j)
               mark(k,j) = 1
            endif
         enddo
      enddo
      ncc(k) = nncc
      ncf(k) = nc - nncc
      j = 0
      do i=1,nc
         if (mark(k,i) == 0) then
            j = j+1
            xcf(k,j) = xc(k,i)
            ycf(k,j) = yc(k,i)
            zcf(k,j) = zc(k,i)
         endif
      enddo
      nncc = 0
      nncf = 0
   enddo
   flag = "close"
   if (flag == "close") then 
      delg = LL / (2. * Nhis)
      pi = 4*ATAN(1.)
      N1 = 0
      do k = 1, snaps 
         N1 = N1 + ncc(k)
      end do 
      N1 = (N1 / snaps)
      rho = N1 / (LL ** 3.d0)
      gr = 0.d0 
      avgc(:) = 0.d0 
      do k = 1, snaps 
         do i = 1, ncc(k) - 1 
            do j = i+1, ncc(k)
               xr = xcc(k,i) - xcc(k,j)
               yr = ycc(k,i) - ycc(k,j)
               zr = zcc(k,i) - zcc(k,j)

               xr=xr-LL*(NINT(xr/LL))
               yr=yr-LL*(NINT(yr/LL))
               zr=zr-LL*(NINT(zr/LL))
               r2=xr*xr+yr*yr+zr*zr
               rr=SQRT(r2)

               IF(rr.LT.LL/2.d0)THEN
                  ig=ceiling(rr/delg)
                  gr(k,ig)=gr(k,ig)+2.
               END IF
            enddo 
         enddo 
      enddo 

      DO j=1,Nhis
         DO i=1,snaps
            rc(j)=delg*(j+0.5)
            vb=((j+1)**3.-j**3.)*delg**3.
            nid=(4./3.)*pi*vb*rho
            gr(i,j)=gr(i,j)/(N1*nid)
         END DO
      END DO

      DO i=1,Nhis
         DO j=1,snaps
            avgc(i)=avgc(i)+gr(j,i)
         END DO
      END DO
      
      OPEN(unit=2,file='rdf_close_phif_'//trim(density)//'.dat',action="write")
      
      DO i=1,Nhis
         WRITE(2,'(2(f17.10,1X))')rc(i),avgc(i)/snaps
         avgc(i) = avgc(i)/snaps
      END DO
   end if 

   flag = "far"

   if (flag == "far") then 
      delg = LL / (2. * Nhis)
      pi = 4*ATAN(1.)
      N1 = 0
      do k = 1, snaps 
         N1 = N1 + ncf(k)
      end do 
      N1 = (N1 / snaps)
      rho = N1 / (LL ** 3.d0)
      gr = 0.d0 
      avgf(:) = 0.d0 
      
      do k = 1, snaps 
         do i = 1, ncf(k) - 1 
            do j = i+1, ncf(k)
               xr = xcf(k,i) - xcf(k,j)
               yr = ycf(k,i) - ycf(k,j)
               zr = zcf(k,i) - zcf(k,j)

               xr=xr-LL*(NINT(xr/LL))
               yr=yr-LL*(NINT(yr/LL))
               zr=zr-LL*(NINT(zr/LL))
               r2=xr*xr+yr*yr+zr*zr
               rr=SQRT(r2)

               IF(rr.LT.LL/2.d0)THEN
                  ig=ceiling(rr/delg)
                  gr(k,ig)=gr(k,ig)+2.
               END IF
            enddo 
         enddo 
      enddo 

      DO j=1,Nhis
         DO i=1,snaps
            rf(j)=delg*(j+0.5)
            vb=((j+1)**3.-j**3.)*delg**3.
            nid=(4./3.)*pi*vb*rho
            gr(i,j)=gr(i,j)/(N1*nid)
         END DO
      END DO

      DO i=1,Nhis
         DO j=1,snaps
            avgf(i)=avgf(i)+gr(j,i)
         END DO
      END DO
      
      OPEN(unit=2,file='rdf_far_phif_'//trim(density)//'.dat',action="write")
      
      DO i=1,Nhis
         WRITE(2,'(2(f17.10,1X))')rf(i),avgf(i)/snaps
         avgf(i) = avgf(i)/snaps
      END DO
   end if 

end subroutine