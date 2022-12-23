      module ice_nemo_concepts

      use ice_kinds_mod
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: alloc_nemo_concepts, diagnostic_wind_components

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         u_zref , & ! wind at zref, x-direction (m/s)
         v_zref     ! wind at zref, y-direction (m/s)

!=======================================================================

      contains

!=======================================================================

      subroutine alloc_nemo_concepts

      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks

      integer (int_kind) :: ierr

      allocate( &
         u_zref      (nx_block,ny_block,max_blocks), &
         v_zref      (nx_block,ny_block,max_blocks), &
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_nemo_concepts): Out of memory')

      end subroutine alloc_nemo_concepts

!=======================================================================

!
! Compute components of diagnostic winds
!

      subroutine diagnostic_wind_components(iblk)

      use ice_blocks, only: nx_block, ny_block, get_block, block
      use ice_domain, only: blocks_ice
      use ice_fileunits, only: nu_diag
      use ice_flux, only: uatm, vatm, wind, Uref
      use ice_grid, only: tmask
      use icepack_intfc, only: icepack_query_parameters, icepack_warnings_flush

      integer (kind=int_kind), intent(in) :: &
         iblk            ! block index

      ! local variables

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i,j                ! horizontal indices

      type (block) :: &
         this_block         ! block information for current block

      real (kind=dbl_kind)    :: &
         umin               ! min wind speed for turbulent fluxes (m/s)

      character(len=*), parameter :: subname = '(diagnostic_wind_components)'

      call icepack_query_parameters(umin_out=umin)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      this_block = get_block(blocks_ice(iblk),iblk)
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      do j = jlo, jhi
      do i = ilo, ihi
         ! compute x and y components of diag wind
         u_zref(i,j,iblk) = uatm(i,j,iblk) * Uref(i,j,iblk) / max(umin, wind(i,j,iblk))
         v_zref(i,j,iblk) = vatm(i,j,iblk) * Uref(i,j,iblk) / max(umin, wind(i,j,iblk))
      enddo
      enddo

      end subroutine diagnostic_wind_components

!=======================================================================

      end module ice_nemo_concepts

!=======================================================================
