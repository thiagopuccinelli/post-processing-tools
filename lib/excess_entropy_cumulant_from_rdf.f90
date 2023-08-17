subroutine excess_entropy_cumulant_from_rdf(gor,r,nbins,csum,newR)
    implicit none 
    integer :: i, ndim
    integer, intent(in) :: nbins
    real(8) :: cum
    real(8),intent(in) :: gor(nbins), r(nbins)  
    real(8), intent(out) :: csum(nbins),newR(nbins)
    DOUBLE PRECISION :: pi

    pi=4*ATAN(1.)
    cum = 0.d0 
    csum = 0.d0 
    newR = 0.d0
    ndim = 0
    do i = 1, nbins 
        if (gor(i) > 0.d0) then
            ndim = ndim + 1  
            cum = cum - (gor(i)*log(gor(i)) - gor(i) +1)*r(i) 
            csum(ndim) = - pi * cum
            newR(ndim) = r(i) 
        endif 
    enddo 


end subroutine