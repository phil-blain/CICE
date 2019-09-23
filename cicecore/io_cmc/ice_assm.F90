MODULE ice_assm
#ifdef CMCform
   !!======================================================================
   !!                       ***  MODULE  ice_assm  ***
   !! Ice Assimilatio:  performs ice direct insertion of data
   !!                 (std files)
   !!=====================================================================
   !! History :  9.0  !  Feb 2014  (F. Roy) 
   !!----------------------------------------------------------------------

   use ice_kinds_mod
   use ice_blocks, only: nx_block, ny_block
   use ice_domain_size
   use ice_communicate, only: my_task, master_task
   use ice_constants
   use ice_calendar, only: sec, idate
   use ice_domain
   use ice_gather_scatter
   use ice_fileunits
   use ice_broadcast, only: broadcast_scalar
   use ice_state
   use ice_flux, only: Tair, Tf
   use ice_grid, only: tmask

   use mod_gemdate

   use ice_exit, only: abort_ice

   IMPLICIT NONE
   PRIVATE   

! -------------------------------------------------
! Direct insertion file management
! -------------------------------------------------

   INTEGER, SAVE                 ::   nu_din
   INTEGER, PARAMETER            ::   nfilemax = 300 ! Max number of direct insertion files
   INTEGER (kind=int_kind), SAVE, &
            DIMENSION(nfilemax)  ::   idate1_in, &   ! yyyymmdd for which direct insertion is done
                                      idate2_in      ! hhmmss for which direct insertion is done
   CHARACTER(LEN=200), SAVE, &
            DIMENSION(nfilemax)  ::   cdndin         ! related file name for direct insertion
   INTEGER, SAVE                 ::   nfile_din

   REAL(kind=real_kind)  , &
     DIMENSION(nx_global,ny_global)                     :: glbtab
   REAL(kind=dbl_kind)  , &
     DIMENSION(nx_global,ny_global)                     :: glbtab8
   INTEGER(kind=int_kind), &
     DIMENSION(nx_global,ny_global)                     :: glbtabm
! aggregated conc thickness
   REAL(kind=dbl_kind),    &
     DIMENSION(nx_block,ny_block,     max_blocks), SAVE :: stdtab_gl, stdtab_ge
   INTEGER(kind=int_kind), &
     DIMENSION(nx_block,ny_block,     max_blocks), SAVE :: stdtabm
! partial conc
   REAL(kind=dbl_kind),    &
     DIMENSION(nx_block,ny_block,ncat,max_blocks), SAVE :: stdtab_glp

   logical (kind=log_kind), public :: dirassm_ice      ! if true, direct insertion

   PUBLIC   ice_Direct_Assm

CONTAINS

   SUBROUTINE ice_Direct_Assm
      implicit none

      !!-------------------------------------------------------------------
      !!                    ***  ROUTINE ice_Direct_Assm  ***
      !!
      !! ** Purpose :   Update the sea-ice state with std files
      !!                created with readdex and interpolated on this grid
      !!
      !! ** Method  :   read
      !!--------------------------------------------------------------------
      INTEGER  ::   ifile
      INTEGER  ::   i, j, k, icat, iblk       ! dummy loop indices
      REAL(kind=dbl_kind) ::   zidto          ! temporary scalar
      CHARACTER(LEN=12)   :: cdindate
      LOGICAL, SAVE       :: first=.true.
      integer        :: hh,mm,ss,idate1_now,idate2_now
      !--------------------------------------------------------------------

      IF (first) THEN
        ! Initialisation
        call get_fileunit(nu_din)
        open(nu_din, file='./cice_din.sequence', &
                     status='OLD', &
                     form='FORMATTED', &
                     access='SEQUENTIAL')
        read(nu_din,*) nfile_din
        IF (nfile_din == 0) THEN
          dirassm_ice = .FALSE.
          return
        ENDIF
        IF (nfile_din > nfilemax) call abort_ice('ice_Direct_Assm: not enough storage')
        cdndin(:) = 'NONE'
        idate1_in(:) = -1
        idate2_in(:) = -1
        if (my_task == master_task) then
          write(nu_diag,*) 'ice_Direct_Assm: nfile_din = ', nfile_din
        endif
        DO ifile = 1, nfile_din
          read(nu_din,'(i8.8,x,i6.6,x,a)') idate1_in(ifile), idate2_in(ifile), cdndin(ifile)
          if (my_task == master_task) write(nu_diag,*) &
               'ice_Direct_Assm:: ifile,idate1_in,idate2_in,cdndin = ',&
                           ifile,idate1_in(ifile),idate2_in(ifile),cdndin(ifile)
        ENDDO
        first = .false.
      ENDIF

      DO ifile = 1, nfile_din

        hh=sec/3600
        mm=(sec-hh*3600)/60
        ss=sec-hh*3600-mm*60
        idate1_now=idate
        idate2_now=hh*10000+mm*100+ss

        IF( idate1_now == idate1_in(ifile) .and. idate2_now == idate2_in(ifile) ) THEN
          IF ( my_task == master_task ) THEN
            write(nu_diag,*) 'ice_Direct_Assm: ==============================='
            write(nu_diag,*) 'ice_Direct_Assm: Ice Direct Insertion of file:'
            write(nu_diag,*) 'ice_Direct_Assm: ',TRIM(cdndin(ifile))
            write(nu_diag,*) 'ice_Direct_Assm: at idate1_now= ', idate1_now
            write(nu_diag,*) 'ice_Direct_Assm:    idate2_now= ', idate2_now
            write(nu_diag,*) 'ice_Direct_Assm: ==============================='
          ENDIF
          CALL ice_Direct_Assm_read(cdndin(ifile))
        ENDIF

      ENDDO

      RETURN

   END SUBROUTINE ice_Direct_Assm

   SUBROUTINE ice_Direct_Assm_read(cdname)

      USE ice_itd 

      IMPLICIT NONE
      CHARACTER(LEN=200), INTENT(IN) :: cdname

      INTEGER(kind=int_kind)                        :: stdunit
      INTEGER(kind=int_kind)                        :: icat,iblk,ierr
      INTEGER(kind=int_kind)                        :: i,j,it
      CHARACTER(LEN=4)               :: nomvar
      CHARACTER(LEN=2)               :: atmp2
      LOGICAL                        :: vreset
      REAL(kind=real_kind), PARAMETER   &
                                     :: epsila = 1.e-5

      integer  fnom,fclos,fstouv,fstfrm,fstlir
      external fnom,fclos,fstouv,fstfrm,fstlir

      if (my_task == master_task) then
        stdunit=0
        ierr=fnom(stdunit,cdname,'STD+RND+R/O',0)
        if (ierr.lt.0) then
          call abort_ice('ice_Direct_Assm_read: error unit assingment: '//TRIM(cdname))
        endif
        ierr = fstouv(stdunit,'RND')
        if (ierr.lt.0) then
          call abort_ice('ice_Direct_Assm_read: error opening file (1): '//TRIM(cdname))
        endif
      endif

      ! ... GE
      nomvar='GE'

      CALL cice_din_getstdtab(stdunit,nomvar)

      ! ... GL
      nomvar='GL'
      CALL cice_din_getstdtab(stdunit,nomvar)
 
      ! ... GL01, GL02, ... (partial concentrations)
      do icat = 1, 8
        write(atmp2,'(i2.2)') icat
        nomvar(1:2)='GL'
        nomvar(3:4)=atmp2
        CALL cice_din_getstdtab(stdunit,nomvar,icat)
      enddo
      do icat = 9, ncat
        stdtab_glp(:,:,icat,:) = c0
      enddo

      if (my_task == master_task) then
        ierr = fstfrm(stdunit)
        if (ierr.lt.0) then
          call abort_ice('ice_Direct_Assm_read: error closing std file unit '//TRIM(cdname))
        endif
        ierr = fclos(stdunit)
        if (ierr.lt.0) then
          call abort_ice('ice_Direct_Assm_read: error closing std file '//TRIM(cdname))
        endif
      endif

      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block
        IF (.not.tmask(i,j,iblk)) stdtabm(i,j,iblk) = 0
        IF (stdtabm(i,j,iblk) == 1) THEN
          IF (stdtab_gl(i,j,iblk) < epsila) stdtab_gl(i,j,iblk) = 0.
        ENDIF
      enddo
      enddo
      enddo

      do iblk = 1, nblocks

         call set_state_var_cice_din(nx_block,     ny_block,           &
                             stdtabm   (:,:,    iblk),                 &
                             stdtab_gl (:,:,    iblk),                 &
                             stdtab_ge (:,:,    iblk),                 &
                             stdtab_glp(:,:,  :,iblk),                 &
                             tmask(:,:,    iblk),                      &
                             Tair (:,:,    iblk),                      &
                             Tf   (:,:,    iblk),                      &
                             aicen(:,:,  :,iblk), trcrn(:,:,:,:,iblk), &
                             vicen(:,:,  :,iblk), vsnon(:,:,  :,iblk), &
                             eicen(:,:,  :,iblk), esnon(:,:,  :,iblk))

         aice(:,:,iblk) = c0
         vice(:,:,iblk) = c0
         vsno(:,:,iblk) = c0
         eice(:,:,iblk) = c0
         esno(:,:,iblk) = c0
         do it = 1, ntrcr
            trcr(:,:,it,iblk) = c0
         enddo

         call aggregate (nx_block, ny_block,  &
                         aicen(:,:,:,iblk),   &
                         trcrn(:,:,:,:,iblk), &
                         vicen(:,:,:,iblk),   &
                         vsnon(:,:,:,iblk),   &
                         eicen(:,:,:,iblk),   &
                         esnon(:,:,:,iblk),   &
                         aice (:,:,  iblk),   &
                         trcr (:,:,:,iblk),   &
                         vice (:,:,  iblk),   &
                         vsno (:,:,  iblk),   &
                         eice (:,:,  iblk),   &
                         esno (:,:,  iblk),   &
                         aice0(:,:,  iblk),   &
                         t4sf (:,:,  iblk),   &
                         tmask(:,:,  iblk),   &
                         trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo                     ! iblk

      !
   END SUBROUTINE ice_Direct_Assm_read

   SUBROUTINE cice_din_getstdtab(stdunit,nomvar,icat_pass)
   IMPLICIT NONE
   INTEGER(kind=int_kind), INTENT(IN)             :: stdunit
   CHARACTER(LEN=4), INTENT(IN)    :: nomvar
   INTEGER(kind=int_kind)                         :: infogd, NI,NJ,NK
   INTEGER(kind=int_kind), INTENT(IN), OPTIONAL   :: icat_pass
   INTEGER(kind=int_kind)                         :: icat_pass_in,icat
   REAL(kind=dbl_kind),    &
     DIMENSION(nx_block,ny_block    ,max_blocks) :: wrk_din

   integer(kind=int_kind)  fnom,fclos,fstouv,fstfrm,fstlir,ierr
   external fnom,fclos,fstouv,fstfrm,fstlir

   icat_pass_in=-1
   if (present(icat_pass)) icat_pass_in=icat_pass
 
   if (my_task == master_task) then
     ierr = fstlir(glbtabm, stdunit, NI, NJ, NK, &
                       -1, ' ', -1, -1, -1, &
                       '@@', nomvar)
     if (ierr.lt.0) then
       call abort_ice('cice_din_getstdtab: error reading (@@): '//TRIM(nomvar))
     endif
     IF (NI /= nx_global) call abort_ice('cice_din_getstdtab NI wrong dimension (@@): '//trim(nomvar))
     IF (NJ /= ny_global) call abort_ice('cice_din_getstdtab NJ wrong dimension (@@): '//trim(nomvar))
     IF (NK /= 1)         call abort_ice('cice_din_getstdtab NK wrong dimension (@@): '//trim(nomvar))
     glbtab8(:,:) = glbtabm(:,:)
   endif
   call scatter_global(wrk_din, glbtab8, &
                       master_task, distrb_info, &
                       field_loc_center, field_type_scalar)
   stdtabm(:,:,:) = wrk_din(:,:,:)

   if (my_task == master_task) then
     ierr = fstlir(glbtab, stdunit, NI, NJ, NK, &
                       -1, ' ', -1, -1, -1, &
                       'A@', nomvar)
     if (ierr.lt.0) then
       call abort_ice('cice_din_getstdtab: error reading: '//TRIM(nomvar))
     endif
     IF (NI /= nx_global) call abort_ice('cice_din_getstdtab NI wrong dimension: '//trim(nomvar))
     IF (NJ /= ny_global) call abort_ice('cice_din_getstdtab NJ wrong dimension: '//trim(nomvar))
     IF (NK /= 1)         call abort_ice('cice_din_getstdtab NK wrong dimension: '//trim(nomvar))
     glbtab8(:,:) = glbtab(:,:)
   endif
   if (TRIM(nomvar).eq.'GL') then
     call scatter_global(stdtab_gl, glbtab8, &
                         master_task, distrb_info, &
                         field_loc_center, field_type_scalar)
    elseif (TRIM(nomvar).eq.'GE') then
     call scatter_global(stdtab_ge, glbtab8, &
                         master_task, distrb_info, &
                         field_loc_center, field_type_scalar)
    elseif ( icat_pass_in /= -1 ) then
     call scatter_global(wrk_din, glbtab8, &
                         master_task, distrb_info, &
                         field_loc_center, field_type_scalar)
     stdtab_glp(:,:,icat_pass_in,:) = wrk_din(:,:,:)
    else
     call abort_ice('cice_din_getstdtab: nomvar not adequate')
   endif

   END SUBROUTINE cice_din_getstdtab


!=======================================================================
!BOP
!
! !IROUTINE: set_state_var_cice_din - initialize multi-category state variables
!                                     overwritten by ice direct insertion
!
! !INTERFACE:
!
      subroutine set_state_var_cice_din (nx_block, ny_block, &
                                glm,             &
                                gl,              &
                                ge,              &
                                glp,             &
                                tmask,           &
                                Tair,            &
                                Tf,              &
                                aicen,    trcrn, &
                                vicen,    vsnon, &
                                eicen,    esnon) 
!
! !DESCRIPTION:
!
! Initialize state in each ice thickness category
! for direct ice insertion purpose
!
! !REVISION HISTORY:
!
! authors: F. Roy, inspired by set_state_var original
!
! !USES:
!
      use ice_state, only: nt_Tsfc
      use ice_therm_vertical, only: heat_capacity, calc_Tsfc, Tmlt
      use ice_itd, only: ilyr1, slyr1, hin_max
      use ice_grid, only: grid_type
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions

      integer (kind=int_kind), dimension (nx_block,ny_block), intent(in) :: &
         glm         ! Assimilated ice fraction mask (0 or 1)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         gl      , & ! Assimilated ice fraction
         ge          ! Assimilated ice thickness in meters

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         glp         ! Assimilated partion concentration of ice

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         tmask      ! true for ice/ocean cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tair    , & ! air temperature  (K)
         Tf          ! freezing temperature (C) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
         intent(out) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(out) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(out) :: &
         esnon     ! energy of melting for each ice layer  (J/m^2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k           , & ! ice layer index
         n           , & ! thickness category index
         it          , & ! tracer index
         icells          ! number of cells initialized with ice

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind) :: &
         slope, Ti, hbar, asum, &
         hinit(ncat)

      indxi(:) = 0
      indxj(:) = 0

      ! Overwrites some state variables.

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
          if (glm(i,j) > 0) then
            trcrn(i,j,nt_Tsfc,n) = Tf(i,j)  ! surface temperature
            eicen(i,j,n) = c0
            esnon(i,j,n) = c0
            if (ntrcr >= 2) then
               do it = 2, ntrcr
                  trcrn(i,j,it,n) = c0
               enddo
            endif
          endif
         enddo
         enddo
      enddo

      ! initial category areas in cells with ice
      hbar = c3  ! initial ice thickness with greatest area
                    ! Note: the resulting average ice thickness 
                    ! tends to be less than hbar due to the
                    ! nonlinear distribution of ice thicknesses 
      do n = 1, ncat
         if (n < ncat) then
            hinit(n) = p5*(hin_max(n-1) + hin_max(n)) ! m
         else                ! n=ncat
            hinit(n) = (hin_max(n-1) + c1) ! m
         endif
      enddo

      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (glm(i,j) > 0) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
         endif               ! tmask
      enddo                  ! i
      enddo                  ! j

      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
        asum=0.
        do n = 1, ncat
              asum=asum+glp(i,j,n)
              if (glp(i,j,n).lt.0.)  &
     &          call abort_ice('set_state_var_cice_din: negative partial concentration')
        enddo
        if (asum.gt.1.01) then
              write(*,*) 'set_state_var_cice_din: total concentration too large i,j=',i,j
              write(*,*) 'set_state_var_cice_din: total concentration too large asum=',asum
              call abort_ice('set_state_var_cice_din: total concentration too large')
        endif
        if (asum.gt.1.) then
              do n=1,ncat
                glp(i,j,n)=glp(i,j,n)/asum
              enddo
        endif
      enddo
 
      do n = 1, ncat

         ! ice volume, snow volume
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
!froy
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            aicen(i,j,n) = glp(i,j,n)
            vicen(i,j,n) = hinit(n) * glp(i,j,n)
            vsnon(i,j,n) = c0
         enddo               ! ij

         ! surface temperature
         if (calc_Tsfc) then
     
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               trcrn(i,j,nt_Tsfc,n) = min(Tsmelt, Tair(i,j) - Tffresh) !deg C
            enddo

         else    ! Tsfc is not calculated by the ice model

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               trcrn(i,j,nt_Tsfc,n) = Tf(i,j)   ! not used
            enddo

         endif       ! calc_Tsfc

         ! other tracers (none at present)

         if (heat_capacity) then

            ! ice energy
            do k = 1, nilyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  ! assume linear temp profile and compute enthalpy
                  slope = Tf(i,j) - trcrn(i,j,nt_Tsfc,n)
                  Ti = trcrn(i,j,nt_Tsfc,n) &
                     + slope*(real(k,kind=dbl_kind)-p5) &
                             /real(nilyr,kind=dbl_kind)

!froy - prevent overflow if Ti=0 
!       (trcrn=0 and Tair=Tf ==> fresh water under ice + Tair > 0)
                  if (Ti /= c0) then
                    eicen(i,j,ilyr1(n)+k-1) = &
                       -(rhoi * (cp_ice*(Tmlt(k)-Ti) &
                       + Lfresh*(c1-Tmlt(k)/Ti) - cp_ocn*Tmlt(k))) &
                       * vicen(i,j,n)/real(nilyr,kind=dbl_kind)
                   else
                    eicen(i,j,ilyr1(n)+k-1) = &
                       -(rhoi * (-cp_ice*Ti &
                       + Lfresh)) &
                       * vicen(i,j,n)/real(nilyr,kind=dbl_kind)
                  endif
               enddo            ! ij
            enddo               ! nilyr

            ! snow energy
            do k = 1, nslyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  Ti = min(c0, trcrn(i,j,nt_Tsfc,n))
                  esnon(i,j,slyr1(n)+k-1) = -rhos*(Lfresh - cp_ice*Ti) &
                                            *vsnon(i,j,n) &
                                            /real(nslyr,kind=dbl_kind)
               enddo            ! ij
            enddo               ! nslyr

         else  ! one layer with zero heat capacity

            ! ice energy
            k = 1

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               eicen(i,j,ilyr1(n)+k-1) = &
                       - rhoi * Lfresh * vicen(i,j,n)
            enddo            ! ij

            ! snow energy
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               esnon(i,j,slyr1(n)+k-1) = & 
                       - rhos * Lfresh * vsnon(i,j,n)
            enddo            ! ij

         endif               ! heat_capacity
      enddo                  ! ncat

      end subroutine set_state_var_cice_din

#endif
END MODULE ice_assm
