subroutine rdf_monomeric2d(N,x,y,Lx,Ly,snaps,avg,r) 
    implicit none 
 
    integer, intent(in) :: N
    integer, intent(in) :: snaps  
    real(8), intent(in) :: Lx,Ly
    real(8), intent(in) :: x(snaps,N), y(snaps,N)
    integer, parameter :: Nhis = 2.**8. 
    real(8), intent(out) :: avg(Nhis),r(Nhis)
    integer i, j, k, ig
    DOUBLE PRECISION::rr,delg,pi,xr,yr,r2,vb,nid,rho
    DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr

    delg=10.0/(Nhis)
    pi=4*ATAN(1.)
    rho = N/((Lx*Ly))
    gr = 0.0d0
    avg(:)=0.d0
    

    DO k=1,snaps
        DO i=1,N - 1 
           DO j=i+1,N 
           xr=x(k,i)-x(k,j)
           yr=y(k,i)-y(k,j)
           
           
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
           gr(i,j)=gr(i,j)/(N*nid)
           avg(j)=avg(j)+gr(i,j)/snaps/N
        END DO
    END DO
  
    
end subroutine 
