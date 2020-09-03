MODULE ice_fldread_std


use ice_kinds_mod

    contains
          subroutine Tf_nonlin(Tf,sss,n)
          implicit none
          real (kind=dbl_kind), intent(inout), dimension(n) :: Tf
          real (kind=dbl_kind), intent(in),    dimension(n) :: sss
          integer (kind=int_kind), intent(in) :: n
    
          integer (kind=int_kind) :: i
          real (kind=dbl_kind) :: S
    ! Millero, 1979
          do i=1,n
            S=sss(i)
            Tf(i) =-0.0575*S+1.710523e-3*S**1.5-2.154996e-4*S*S
          enddo
    
          return
          end subroutine Tf_nonlin
    
    END MODULE ice_fldread_std
    