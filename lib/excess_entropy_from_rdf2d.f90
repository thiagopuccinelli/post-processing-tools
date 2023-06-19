subroutine excess_entropy_from_rdf2d(gor,r,sex,csum,newR)
    implicit none 
    integer :: i, ndim
    integer, parameter :: Nhis = 2**9
    real(8) :: cum
    real(8),intent(in) :: gor(Nhis), r(Nhis)  
    real(8), intent(out) :: sex
    real(8), intent(out) :: csum(Nhis),newR(Nhis)

    sex = 0.d0
    cum = 0.d0 
    csum = 0.d0 
    newR = 0.d0
    ndim = 0
    do i = 1, Nhis 
        if (gor(i) > 0.d0) then
            ndim = ndim + 1  
            cum = cum - (gor(i)*log(gor(i)) - gor(i) +1)*r(i) !for excess entropy
            if (gor(i) /= 0.0d0) sex = sex + (gor(i)*log(gor(i)) - gor(i) +1) !for excess entropy
            csum(ndim) = -3.14 * cum
            newR(ndim) = r(i) 
        endif 
    enddo 


end subroutine