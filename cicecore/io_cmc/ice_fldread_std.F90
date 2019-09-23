MODULE ice_fldread_std
#ifdef key_rpne
   !!======================================================================
   !!                       ***  MODULE  ice_fldread_std  ***
   !! Ice state:  read input field
   !!                 (std files)
   !!=====================================================================
   !! History :  9.0  !  Feb 2017  (F. Roy) 
   !!----------------------------------------------------------------------

   use ice_kinds_mod
   use ice_blocks, only: nx_block, ny_block
   use ice_domain_size
   use ice_communicate, only: my_task, master_task
   use ice_constants
   use ice_calendar, only: istep1, sec, dt, yday, &
                           dats_zero
   use ice_domain
   use ice_gather_scatter
   use ice_fileunits, only: nu_diag
   use ice_broadcast, only: broadcast_scalar

   use mod_gemdate

   use ice_exit, only: abort_ice

   IMPLICIT NONE
   PRIVATE   

! -------------------------------------------------
! General options for atm. and oceanic forcing
! -------------------------------------------------

   CHARACTER(LEN=16), PUBLIC ::  cmc_dateo  
             !: date corresponding to initial state

! -------------------------------------------------
! General options for reading atm. forcing
! -------------------------------------------------

   TYPE, PUBLIC :: cmc_ipt                         ! CMC input attributes for namelist
      character (len=10)       :: long_name        ! uniform long name defined in CICE 
      character (len=4)        :: nomvar           ! nomvar provided by user
      character (len=1)        :: structure        ! data structure : field (F) or constant (C)
                                                   ! (note constant is provided in namelist)
      real (kind=real_kind)    :: const_value      ! if data structure is constant, the value
      character (len=1)        :: field_format     ! field format if so: binary (B) or integer (I) 
                                                   ! (note, binary assumed real_kind)
      integer (kind=int_kind)  :: field_ip1_1      ! ip1 first  value to select field
      integer (kind=int_kind)  :: field_ip1_2      ! ip1 second value to select field
      real (kind=real_kind)    :: scaling          ! field scaling value if so 
      real (kind=real_kind)    :: offset           ! field offset  value if so 
                                                   !   field=scaling*field + offset)
   END TYPE cmc_ipt

   TYPE(cmc_ipt), &   !CMC field or constant forcing input specifications
      PUBLIC :: &
      cmc_ulat,cmc_ulon,cmc_kmt,cmc_hte,cmc_htn, &  !GRID
      cmc_tlat,cmc_tlon,                         &  !GRID
      cmc_gl,cmc_i8,cmc_snow, &
      cmc_uice0,cmc_vice0,cmc_tair

   LOGICAL , PUBLIC :: &
              cmc_ln_wmocat, &  ! Apply reading of ice categories according to wmo defs. (T)
              cmc_ln_avgcat, &  ! Apply reading of ice conc. (GL) and thickness (I8=vice/GL) (T)
              cmc_ln_itdspr     ! Apply ice re-distribution (T)

   integer(kind=int_kind), public :: &
      cmc_tt_dateo, cmc_tt_deet, cmc_tt_npas, &
      cmc_tt_ip1, cmc_tt_ip2, cmc_tt_ip3, &
      cmc_tt_ig1, cmc_tt_ig2, cmc_tt_ig3, cmc_tt_ig4, &
      cmc_tt_datyp
   character(LEN=1), public :: cmc_tt_typvar, cmc_tt_grtyp

   real(kind=real_kind), public, &
                         dimension(nx_global,ny_global) :: &
                         cmclat, cmclon

   real(kind=real_kind), public, &
                         dimension(nx_global,1) :: &
                         cmclx

   real(kind=real_kind), public, &
                         dimension(1,ny_global) :: &
                         cmcly

   logical             , public :: &
                         cmc_ln_gemgrid, cmc_ln_nemogrid, cmc_ln_bdyopn

   integer(kind=int_kind), public :: &
                         cmc_nbuff_zone    ! Buffer zone dimension

! -------------------------------------------------
! Working variables
! -------------------------------------------------

   TYPE, PUBLIC ::   FLD                 ! Input field related variables
      CHARACTER(len=4) ::   nomvar       ! Name of the variable in STD files
      REAL(kind=real_kind), DIMENSION(nx_block,ny_block,max_blocks,2) ::   fdta 
   END TYPE FLD

   INTEGER , PUBLIC, PARAMETER ::   cmcnfld = 2  ! maximum number of fields to read
   TYPE(FLD), PUBLIC, DIMENSION(cmcnfld) :: cmcdta
                    ! structure of input fields (variables information, fields read)

   INTEGER :: iun_std(1000), &  ! Unit number vector
  &           nfile_std
   LOGICAL :: frc_file_L=.false.

   INTEGER, DIMENSION(2) :: cmc_lev_code

   PUBLIC  ice_fldread_state, &
  &        ice_fldread_grid_g, ice_fldread_grid_l, &
  &        bcast_cmc_options, Tf_nonlin, &
  &        ice_fldread_grid_att

CONTAINS

   SUBROUTINE ice_fldread_state(nxblk,nyblk,n,nblk, &
                                aicen,vicen,vsnon,Tair, &
                                uvel,vvel, &
                                restart,hin_max)
      implicit none
      integer (kind=int_kind), intent(in) :: nxblk,nyblk,n,nblk
      real    (kind=dbl_kind), dimension(nxblk,nyblk,n,nblk), &
            intent(inout) :: aicen,vicen,vsnon
      real    (kind=dbl_kind), dimension(nxblk,nyblk,nblk), &
            intent(inout) :: Tair,uvel,vvel
      logical, intent(in) :: restart
      real (kind=dbl_kind), intent(in) :: &
         hin_max(0:n) ! category limits (m)


      character*16 datev
      character*2  cnn
      character*4  nvmcat 
      integer iblk,i,j,icat,jf,ierror
      real    (kind=dbl_kind), dimension(n) :: ainit,hinit,sinit
!jlei For MCH spreading
      real (kind=dbl_kind) :: &
         zfbar, zhbar, zhsno, zhinf, &
         za, zb, zhm1, zhp1, zjk, zfk, zhk, zhup, zfseuil=0.85
!jlei end

      real    (kind=dbl_kind) :: asum, h1, h2

      TYPE(cmc_ipt) :: cmc_wrk

      aicen=c0
      vicen=c0
      vsnon=c0
      uvel=c0
      vvel=c0
      Tair = Tffresh - c5

      if (restart) then
        write(nu_diag,*) 'ice_fldread_state: WILL SKIP READING OF VARIABLES restart=',restart
        call flush(nu_diag)
        return
       else
        write(nu_diag,*) 'ice_fldread_state: WILL READ VARIABLES restart=',restart
        call flush(nu_diag)
      endif

      if (nxblk.ne.nx_block.or.nyblk.ne.ny_block.or.  &
     &    n.ne.ncat.or.nblk.ne.nblocks) then
        write(nu_diag,*) 'ice_fldread_state inconsistent dimensions'
        write(nu_diag,*) 'nxblk,nyblk,n,nblk',nxblk,nyblk,n,nblk
        call abort_ice('ice_fldread_state inconsistent dimensions') 
      endif

      do icat = 1, n
         if (icat < n) then
            hinit(icat) = p5*(hin_max(icat-1) + hin_max(icat)) ! m
         else                ! icat=n
            hinit(icat) = (hin_max(icat-1) + c1) ! m
         endif
      enddo

      if (my_task == master_task .and. .not.frc_file_L) call frc_openf

      DO jf = 1, SIZE( cmcdta )         
        cmcdta(jf)%fdta(:,:,:,:) = 0.0 
      ENDDO

      datev=' '

! ------------------------------------------------------------------
! First read in partial concentrations and snow for all categories
! (use cmcdta temporarily, will then be reserved for forcing data)
! ------------------------------------------------------------------

      if (my_task == master_task ) &
     &  write(nu_diag,*) 'ICE STATE: '

      if (cmc_ln_wmocat) then

        if (n.ne.10) &
         call abort_ice('ice_fldread_state, ncat must be equal to 10')
        if (my_task == master_task ) &
     &    write(nu_diag,*) '===>WMO CATEGORIES '
        call flush(nu_diag)
        if (cmc_gl%structure.eq.'C') then
          call flush(nu_diag)
          call abort_ice('ice_fldread_state, option not managed cmc_gl%structure=C')
        endif
        if (trim(cmc_gl%long_name).ne.'null') then
          if (cmc_gl%nomvar(3:4).ne.'') then
            call flush(nu_diag)
            call abort_ice('ice_fldread_state, WMO name for cmc_gl exceed 2 characters')
          endif
        endif
        if (trim(cmc_i8%long_name).ne.'null') then
          if (cmc_i8%nomvar(3:4).ne.'') then
            call flush(nu_diag)
            call abort_ice('ice_fldread_state, WMO name for cmc_i8 exceed 2 characters')
          endif
          if (cmc_i8%nomvar(1:1).eq.''.or.cmc_i8%nomvar(2:2).eq.'') then
            call abort_ice('ice_fldread_state, WMO name for cmc_i8 should have 2 characters')
          endif
          write(nu_diag,*) 'WARNING, ice_fldread_state:cmc_ln_wmocat=T'
          write(nu_diag,*) 'WARNING, ice_fldread_state:cmc_i8 used here represents GE, vicen...'
          call flush(nu_diag)
        endif
        if (trim(cmc_snow%long_name).ne.'null') then
          if (cmc_snow%nomvar(3:4).ne.''.and.trim(cmc_snow%long_name).ne.'null') then
            call flush(nu_diag)
            call abort_ice('ice_fldread_state, WMO name for cmc_snow exceed 2 characters')
          endif
          if (cmc_snow%nomvar(1:1).eq.''.or.cmc_snow%nomvar(2:2).eq.'') then
            call flush(nu_diag)
            call abort_ice('ice_fldread_state, WMO name for cmc_snow should have 2 characters')
          endif
        endif
        if (cmc_gl%nomvar(1:1).eq.''.or.cmc_gl%nomvar(2:2).eq.'') then
          call flush(nu_diag)
          call abort_ice('ice_fldread_state, WMO name for cmc_gl should have 2 characters')
        endif
! ... first deals with concentration
        if (trim(cmc_gl%long_name).ne.'null') then
          do icat=1,n
            write(cnn,'(i2.2)') icat
            nvmcat(3:4)=cnn
            cmcdta(1)%nomvar='GL'  !Tmp
            nvmcat(1:2)=trim(cmc_gl%nomvar)
            cmc_lev_code(1)=cmc_gl%field_ip1_1
            cmc_lev_code(2)=cmc_gl%field_ip1_2
            call readfst (cmcdta(1)%fdta (:,:,:,1), &
     &                       nvmcat, &
     &                       datev,cmc_lev_code, &
     &                       cmc_gl%scaling, &
     &                       cmc_gl%offset, &
     &                       cmc_gl%field_format)
            do iblk=1,nblk
            do j=1,nyblk
            do i=1,nxblk
              aicen(i,j,icat,iblk)=cmcdta(1)%fdta (i,j,iblk,1)
            enddo
            enddo
            enddo
          enddo
          do iblk=1,nblk
          do j=1,nyblk
          do i=1,nxblk
            asum=0. 
            do icat=1,n
              asum=asum+aicen(i,j,icat,iblk)
              if (aicen(i,j,icat,iblk).lt.0.) then
                write(nu_diag,*) 'ice_fldread_state: negative concentration',i,j,icat,iblk
                write(nu_diag,*) 'ice_fldread_state: negative concentration asum=',asum
                call flush(nu_diag)
                call abort_ice('ice_fldread_state: negative partial concentration')
              endif
            enddo
            if (asum.gt.1.01) then
              write(nu_diag,*) 'ice_fldread_state: total concentration too large i,j,iblk=',i,j,iblk
              write(nu_diag,*) 'ice_fldread_state: total concentration too large asum=',asum
              write(nu_diag,*) 'ice_fldread_state: total concentration too large aicen=',aicen(i,j,:,iblk) ! FD
              call flush(nu_diag)
              call abort_ice('ice_fldread_state: total concentration too large')
            endif
            if (asum.gt.1.) then
              do icat=1,n
                aicen(i,j,icat,iblk)=aicen(i,j,icat,iblk)/asum
              enddo
            endif
          enddo
          enddo
          enddo
         else
          if (my_task == master_task ) &
     &      write(nu_diag,*) '===> ALL WMO CATEGORIES NOT DEFINED AND SET TO ZERO'
          call flush(nu_diag)
        endif
! ... now deals with volume and snow
! ... Careful notation i8 used here as GE
        do icat=1,n
          write(cnn,'(i2.2)') icat
          nvmcat(3:4)=cnn
          vicen(:,:,icat,:) = p5*(hin_max(icat)+hin_max(icat-1))*aicen(:,:,icat,:)
          if (trim(cmc_i8%long_name).ne.'null') then
            if (cmc_i8%structure.eq.'C') then
              call flush(nu_diag)
              call abort_ice('ice_fldread_state, option not managed cmc_i8%structure=C') 
            endif
            cmcdta(1)%nomvar='I8'  !Tmp
            nvmcat(1:2)=trim(cmc_i8%nomvar)
            cmc_lev_code(1)=cmc_i8%field_ip1_1
            cmc_lev_code(2)=cmc_i8%field_ip1_2
            call readfst (cmcdta(2)%fdta (:,:,:,1), &
     &                       nvmcat, &
     &                       datev,cmc_lev_code, &
     &                       cmc_i8%scaling, &
     &                       cmc_i8%offset, &
     &                       cmc_i8%field_format)
            do iblk=1,nblk
            do j=1,nyblk
            do i=1,nxblk
              vicen(i,j,icat,iblk)=cmcdta(2)%fdta (i,j,iblk,1)
            enddo
            enddo
            enddo
           else
            if (my_task == master_task ) &
     &        write(nu_diag,*) '===> ALL WMO PARTIAL VOLUME NOT DEFINED AND SET TO DEFAULT'
          endif
          vsnon(:,:,icat,:) = c0
          if (trim(cmc_snow%long_name).ne.'null') then
            if (cmc_snow%structure.eq.'C') then
              call flush(nu_diag)
              call abort_ice('ice_fldread_state, option not managed cmc_snow%structure=C') 
            endif
            cmcdta(1)%nomvar='SNOW'  !Tmp
            nvmcat(1:2)=trim(cmc_snow%nomvar)
            cmc_lev_code(1)=cmc_snow%field_ip1_1
            cmc_lev_code(2)=cmc_snow%field_ip1_2
            call readfst (cmcdta(3)%fdta (:,:,:,1), &
     &                       nvmcat, &
     &                       datev,cmc_lev_code, &
     &                       cmc_snow%scaling, &
     &                       cmc_snow%offset, &
     &                       cmc_snow%field_format)
            do iblk=1,nblk
            do j=1,nyblk
            do i=1,nxblk
              vsnon(i,j,icat,iblk)=cmcdta(3)%fdta (i,j,iblk,1)
            enddo
            enddo
            enddo
           else
            if (my_task == master_task ) &
     &        write(nu_diag,*) '===> ALL WMO PARTIAL SNOW VOLUME NOT DEFINED AND SET TO ZERO'
            call flush(nu_diag)
          endif
        enddo
        write(nu_diag,*) '===> FINISHED READING WMO VARIABLES'
        call flush(nu_diag)
      endif  !cmc_ln_wmocat

      if (cmc_ln_avgcat) then

        if (my_task == master_task ) &
     &    write(nu_diag,*) '===>AVERAGED CONDITIONS '

        if (cmc_gl%structure.eq.'C') then
          cmcdta(1)%fdta (:,:,:,1)=cmc_gl%const_value
          if (my_task == master_task ) &
     &    write(nu_diag,*) '===>GL CONSTANT: ',cmc_gl%const_value
         elseif (cmc_gl%structure.eq.'F') then
          cmcdta(1)%nomvar='GL'  !Tmp
          cmc_lev_code(1)=cmc_gl%field_ip1_1
          cmc_lev_code(2)=cmc_gl%field_ip1_2
          call readfst (cmcdta(1)%fdta (:,:,:,1), &
     &                     cmc_gl%nomvar, &
     &                     datev,cmc_lev_code, &
     &                     cmc_gl%scaling, &
     &                     cmc_gl%offset, &
     &                     cmc_gl%field_format)
        endif

        if (trim(cmc_i8%long_name).ne.'null') then

          if (my_task == master_task ) &
     &    write(nu_diag,*) '===>WARNING, THICKNESS BASED ON I8'

          if (cmc_i8%structure.eq.'C') then
            cmcdta(2)%fdta (:,:,:,1)=cmc_i8%const_value
            if (my_task == master_task ) &
     &      write(nu_diag,*) '===>I8 CONSTANT: ',cmc_i8%const_value
           elseif (cmc_i8%structure.eq.'F') then
            cmcdta(2)%nomvar='I8'  !Tmp
            cmc_lev_code(1)=cmc_i8%field_ip1_1
            cmc_lev_code(2)=cmc_i8%field_ip1_2
            call readfst (cmcdta(2)%fdta (:,:,:,1), &
     &                       cmc_i8%nomvar, &
     &                       datev,cmc_lev_code, &
     &                       cmc_i8%scaling, &
     &                       cmc_i8%offset, &
     &                       cmc_i8%field_format)

          endif

! Update state variables
! jlei MCH ice spreading

          do iblk=1,nblk
          do j=1,nyblk
          do i=1,nxblk
            do icat=1,n
              aicen(i,j,icat,iblk) = c0
              vicen(i,j,icat,iblk) = c0
              vsnon(i,j,icat,iblk) = c0
            enddo
            if  (cmc_ln_itdspr) then
                zfbar = cmcdta(1)%fdta (i,j,iblk,1)
                zhbar = cmcdta(1)%fdta (i,j,iblk,1) * cmcdta(2)%fdta (i,j,iblk,1)
!                zhsno = work3(i,j,iblk)
                if ( zfbar <= zfseuil .and. zfbar > 0.01_dbl_kind .and. zhbar > 0.1_dbl_kind) then
                   ! Assume a linear distribution in the MIZ
                   zhinf = 3._dbl_kind*zhbar / zfbar
                   za=2._dbl_kind*zfbar / zhinf 
                   zb=za / zhinf
                   do icat=1,n
                      zhm1=hin_max(icat-1)
                      zhp1=hin_max(icat)
                      if (zhm1 .ge. zhinf ) then
                       zfk=0._dbl_kind
                         if (icat.eq.n ) then
                          zhk=(zhm1+1.)*zfk
                         else
                          zhk=0.5_dbl_kind*(zhp1+zhm1)*zfk
                         end if
                      else
                         if (zhinf.le.zhp1) then 
                            zhup=zhinf
                         else
                            zhup=zhp1
                         end if
                         zfk = za*(zhup - zhm1)-zb*(zhup*zhup-zhm1*zhm1)/2._dbl_kind
                         zhk = za*(zhup*zhup - zhm1*zhm1)/2._dbl_kind-zb*(zhup*zhup*zhup - zhm1*zhm1*zhm1)/3._dbl_kind
                      end if
                      aicen(i,j,icat,iblk)=zfk
                      vicen(i,j,icat,iblk)=zhk
                   end do
                 elseif ( zfbar > zfseuil .and. zhbar > 0.1_dbl_kind ) then
 ! Assume a parabolic distribution in compact pack ice
                   zhinf = 2._dbl_kind*zhbar / zfbar
                   za = 6._dbl_kind * zfbar / (zhinf*zhinf*zhinf) 
                   do icat=1,n
                      zhm1=hin_max(icat-1)
                      zhp1=hin_max(icat)
                      if (zhm1 .ge. zhinf ) then
                       zfk=0._dbl_kind
                         if (icat.eq.n) then
                          zhk=(zhm1+1.)*zfk
                         else
                           zhk=0.5_dbl_kind*(zhp1+zhm1)*zfk
                         end if
                       else
                         if (zhinf.le.zhp1) then
                            zhup=zhinf
                         else
                            zhup=zhp1
                         end if
                         zfk = za*( zhinf*(zhup*zhup - zhm1*zhm1) / 2._dbl_kind - ( zhup*zhup*zhup - zhm1*zhm1*zhm1 ) / 3._dbl_kind )
                         zhk = za*( zhinf*(zhup*zhup*zhup - zhm1*zhm1*zhm1) / 3._dbl_kind - ( zhup*zhup*zhup*zhup - zhm1*zhm1*zhm1*zhm1 ) / 4._dbl_kind )
                      end if
                      aicen(i,j,icat,iblk)=zfk
                      vicen(i,j,icat,iblk)=zhk
                   end do
                 else
                   do icat=1,n
                     !-----------------------------------------------------------------
                     ! Place ice at given thickness category
                     if ( hin_max(icat-1) <= cmcdta(2)%fdta (i,j,iblk,1) .and. &
                 &                        cmcdta(2)%fdta (i,j,iblk,1) < hin_max(icat) ) then
                        aicen(i,j,icat,iblk) = zfbar
                        vicen(i,j,icat,iblk) = zhbar ! m
                     endif
                   enddo
                 end if
!                vsnon(i,j,:,iblk) = zhsno * aicen(i,j,:,iblk)
             else   ! not cmc_ln_itdspr
              h1=cmcdta(2)%fdta (i,j,iblk,1)
              do icat=1,n 
                if (h1 > hin_max(icat-1) .and. h1 <= hin_max(icat)) then
                  aicen(i,j,icat,iblk) =  &
     &              max( cmcdta(1)%fdta (i,j,iblk,1), 0.0_real_kind )
                  vicen(i,j,icat,iblk) = h1 * aicen(i,j,icat,iblk) ! m
                endif
              enddo
            endif ! not cmc_ln_itdspr
          enddo
          enddo
          enddo

         else

          if (my_task == master_task ) &
     &    write(nu_diag,*) '===>WARNING, HARD CODED THICKNESS icat=7'

          do iblk=1,nblk
          do j=1,nyblk
          do i=1,nxblk
            do icat=1,n
              aicen(i,j,icat,iblk) = c0
              vicen(i,j,icat,iblk) = c0
              vsnon(i,j,icat,iblk) = c0
            enddo
            aicen(i,j,7,iblk) = cmcdta(1)%fdta (i,j,iblk,1)
            vicen(i,j,7,iblk) = hinit(7) * aicen(i,j,7,iblk) ! m
            !put all the ice in category 7
          enddo
          enddo
          enddo

        endif

      endif  !cmc_ln_avgcat

!... Ice velocity
      cmcdta(1)%nomvar='UICE'  !Tmp
      cmcdta(1)%fdta (:,:,:,1)=c0
      if (trim(cmc_uice0%long_name).ne.'null') then
        if (cmc_uice0%structure.eq.'C') then
          cmcdta(1)%fdta (:,:,:,1)=cmc_uice0%const_value
          if (my_task == master_task ) &
     &    write(nu_diag,*) '===>UICE CONSTANT: ',cmc_uice0%const_value
         elseif (cmc_uice0%structure.eq.'F') then
          if (my_task == master_task ) &
     &    write(nu_diag,*) '===>UICE FIELD: ',cmc_uice0%nomvar
          cmc_lev_code(1)=cmc_uice0%field_ip1_1
          cmc_lev_code(2)=cmc_uice0%field_ip1_2
          call readfst (cmcdta(1)%fdta (:,:,:,1), &
                           cmc_uice0%nomvar, &
     &                     datev, &
     &                     cmc_lev_code, &
     &                     cmc_uice0%scaling, &
     &                     cmc_uice0%offset, &
     &                     cmc_uice0%field_format )
        endif
      endif
      cmcdta(2)%nomvar='VICE'  !Tmp
      cmcdta(2)%fdta (:,:,:,1)=c0
      if (trim(cmc_vice0%long_name).ne.'null') then
        if (cmc_vice0%structure.eq.'C') then
          cmcdta(2)%fdta (:,:,:,1)=cmc_vice0%const_value
          if (my_task == master_task ) &
     &    write(nu_diag,*) '===>VICE CONSTANT: ',cmc_vice0%const_value
         elseif (cmc_vice0%structure.eq.'F') then
          if (my_task == master_task ) &
     &    write(nu_diag,*) '===>VICE FIELD: ',cmc_vice0%nomvar
          cmc_lev_code(1)=cmc_vice0%field_ip1_1
          cmc_lev_code(2)=cmc_vice0%field_ip1_2
          call readfst (cmcdta(2)%fdta (:,:,:,1), &
                           cmc_vice0%nomvar, &
     &                     datev, &
     &                     cmc_lev_code, &
     &                     cmc_vice0%scaling, &
     &                     cmc_vice0%offset, &
     &                     cmc_vice0%field_format )
        endif
      endif

      uvel=cmcdta(1)%fdta (:,:,:,1)
      vvel=cmcdta(2)%fdta (:,:,:,1)

!... Air temperature 
!... For ice temperature initialization 
!... and may have an impact on halo restoring 
!... of initial conditions when restarting
      if (trim(cmc_tair%long_name).eq.'null') then
        call abort_ice('ice_fldread_state:MUST PROVIDE AIR TEMPERATURE IF NOT RESTARTING')
      endif
      if (cmc_tair%structure.eq.'C') then
        Tair=cmc_tair%const_value
        if (my_task == master_task ) &
     &  write(nu_diag,*) '===>TAIR CONSTANT: ',cmc_tair%const_value
       elseif (cmc_tair%structure.eq.'F') then
        if (my_task == master_task ) &
     &  write(nu_diag,*) '===>TAIR FIELD: ',cmc_tair%nomvar
        cmcdta(1)%nomvar='TAIR'  !Tmp
        cmc_lev_code(1)=cmc_tair%field_ip1_1
        cmc_lev_code(2)=cmc_tair%field_ip1_2
        call readfst (cmcdta(1)%fdta (:,:,:,1), &
                         cmc_tair%nomvar, &
     &                   cmc_dateo, &
     &                   cmc_lev_code, &
     &                   cmc_tair%scaling, &
     &                   cmc_tair%offset, &
     &                   cmc_tair%field_format )
        Tair=cmcdta(1)%fdta (:,:,:,1)
      endif

      return

   END SUBROUTINE ice_fldread_state

   SUBROUTINE ice_fldread_grid_g(nomvar,ftype,scaling,offset, &
                                 ip1_1,ip1_2,work_g8,lx,ly)
      implicit none

      character(LEN=4), intent(in) :: nomvar
      character(LEN=1), intent(in) :: ftype
      real (kind=real_kind), &
            intent(in) :: scaling, offset
      integer (kind=int_kind), &
            intent(in) :: ip1_1, ip1_2
      real (kind=dbl_kind), &
            dimension(nx_global,ny_global), &
            intent(out) :: work_g8    ! output array (real, 8-byte)
      LOGICAL, optional                   :: lx,ly
      real (kind=real_kind), &
            dimension(nx_global,ny_global)  &
                        :: work_g     ! output array (real, 4-byte)

      character(len=16) datev
      integer (kind=int_kind), dimension(2) :: cmc_lev
      LOGICAL                             :: lxw,lyw

      cmc_lev(1)=ip1_1
      cmc_lev(2)=ip1_2

      if (my_task == master_task) then
        write(nu_diag,*) 'Trying to read global grid:',nomvar,ip1_1,ip1_2
        call flush(nu_diag)
      endif
      if (my_task == master_task ) then
        if (.not.frc_file_L) call frc_openf
        datev=' '
        if (present(lx).or.present(ly)) then
          lxw=.false.
          lyw=.false.
          if (present(lx)) lxw=lx
          if (present(ly)) lyw=ly
          call readfst_global(work_g,nomvar,datev,cmc_lev, &
                                 scaling,offset,ftype,lx=lxw,ly=lyw)
         else
          call readfst_global(work_g,nomvar,datev,cmc_lev, &
                                 scaling,offset,ftype)
        endif
        work_g8=work_g
      endif

   END SUBROUTINE ice_fldread_grid_g

   SUBROUTINE ice_fldread_grid_l(nomvar,ftype,scaling,offset, &
                                 ip1_1,ip1_2,work_l8)
      implicit none

      real (kind=dbl_kind), intent(out), &
            dimension (nx_block,ny_block,max_blocks) :: work_l8
      character(LEN=4), intent(in) :: nomvar
      character(LEN=1), intent(in) :: ftype
      real (kind=real_kind), &
            intent(in) :: scaling, offset
      integer (kind=int_kind), &
            intent(in) :: ip1_1, ip1_2

      real(kind=real_kind), dimension(nx_block,ny_block,max_blocks) &
                         :: work_l
      character*16 datev
      integer (kind=int_kind), dimension(2) :: cmc_lev

      cmc_lev(1)=ip1_1
      cmc_lev(2)=ip1_2

      if (my_task == master_task .and. .not.frc_file_L) call frc_openf

      datev=' '

      if (my_task == master_task) then
        write(nu_diag,*) 'Trying to read local grid:',nomvar,ip1_1,ip1_2
      endif
      call readfst (work_l,nomvar,datev,cmc_lev, &
                       scaling,offset,ftype)
      work_l8=work_l

   END SUBROUTINE ice_fldread_grid_l

   SUBROUTINE ice_fldread_grid_att(nomvar, dateo, deet, npas, &
                          ip1, ip2, ip3, ig1, ig2, ig3, ig4, &
                          datyp, typvar, grtyp)
      implicit none

      character(LEN=4), intent(in)        :: nomvar
      integer(kind=int_kind), intent(out) :: dateo, deet, npas, &
                                             ip1, ip2, ip3, &
                                             ig1, ig2, ig3, ig4, &
                                             datyp
      character(LEN=1), intent(out)       :: typvar, grtyp

!... Dummy variables for our purpose
      integer(kind=int_kind)              ::  NI, NJ, NK, NBITS, &
                                              SWA, LNG, DLTF, UBC, &
                                              EXTRA1, EXTRA2, EXTRA3
!
      integer(kind=int_kind)              ::  ierr,infon,nmax=1
      integer(kind=int_kind)              ::  fstinl,fstprm
      integer(kind=int_kind),dimension(1) ::  liste

      character(LEN=4)                    :: nomvarw,nomvarwo
      character(LEN=8)                    :: etiket

      if (my_task == master_task) then
        write(nu_diag,*) 'Trying to read Z grid attributes:',nomvar
        call flush(nu_diag)

        if (.not.frc_file_L) call frc_openf

        nomvarw=cmc_kmt%nomvar
        if (nomvarw.ne.'null') then
          ierr = fstinl(iun_std, NI, NJ, NK, -1, ' ', -1, -1, -1, ' ', &
                       nomvarw, LISTE, INFON, nmax)
          ierr = fstprm(liste(1), DATEO, DEET, NPAS, NI, NJ, NK, &
                     NBITS, DATYP, IP1, &
                     IP2, IP3, TYPVAR, NOMVARWO, &
                     ETIKET, GRTYP, IG1, IG2, IG3, IG4, &
                     SWA, LNG, DLTF, UBC, EXTRA1, EXTRA2, EXTRA3)
          if (grtyp.ne.'Z') &
           call abort_ice('NOT A GEM GRID and cmc_ln_gemgrid is true?')
         else
          IF(my_task == master_task) THEN
            WRITE(nu_diag,*) &
             'ice_fldread_grid_att: WARNING NO MASK IN INPUT FILE, NO CONTROL ON INPUT grtyp'
            WRITE(nu_diag,*) &
             'ice_fldread_grid_att: Make sure grtyp=Z ... '
          ENDIF
        endif

        ierr = fstinl(iun_std, NI, NJ, NK, -1, ' ', -1, -1, -1, ' ', &
                     nomvar, LISTE, INFON, nmax)
  
        ierr = fstprm(liste(1), DATEO, DEET, NPAS, NI, NJ, NK, &
                   NBITS, DATYP, IP1, &
                   IP2, IP3, TYPVAR, NOMVARWO, &
                   ETIKET, GRTYP, IG1, IG2, IG3, IG4, &
                   SWA, LNG, DLTF, UBC, EXTRA1, EXTRA2, EXTRA3)

        IF(my_task == master_task) THEN 
          write(nu_diag,*) 'GEM GRID ATTRIBUTES:'
          write(nu_diag,*) 'ip1=',ip1
          write(nu_diag,*) 'ip2=',ip2
          write(nu_diag,*) 'ip3=',ip3
          write(nu_diag,*) 'ig1=',ig1
          write(nu_diag,*) 'ig2=',ig2
          write(nu_diag,*) 'ig3=',ig3
          write(nu_diag,*) 'ig4=',ig4
          write(nu_diag,*) 'datyp=',datyp
          write(nu_diag,*) 'typvar=',typvar
          write(nu_diag,*) 'grtyp=',grtyp
        ENDIF

      endif

      return

   END SUBROUTINE ice_fldread_grid_att

   SUBROUTINE bcast_cmc_options
      implicit none

      IF(my_task == master_task) THEN
        WRITE(nu_diag,*)
        WRITE(nu_diag,*) 'cmc input files: '
        WRITE(nu_diag,*) '~~~~~~~~'
        WRITE(nu_diag,*) '          Namelist namcmc_nml'
        WRITE(nu_diag,*) '          Ice Analysis date       cmc_dateo       = ', cmc_dateo
        WRITE(nu_diag,*) '                                  cmc_ln_wmocat   = ', cmc_ln_wmocat
        WRITE(nu_diag,*) '                                  cmc_ln_avgcat   = ', cmc_ln_avgcat
        WRITE(nu_diag,*) '                                  cmc_ln_itdspr   = ', cmc_ln_itdspr
      ENDIF

      call broadcast_scalar(cmc_dateo,          master_task)
      call broadcast_scalar(cmc_ln_wmocat,      master_task)
      call broadcast_scalar(cmc_ln_avgcat,      master_task)
      call broadcast_scalar(cmc_ln_itdspr,      master_task)

      call broadcast_scalar(cmc_ln_gemgrid,     master_task)
      call broadcast_scalar(cmc_ln_nemogrid,    master_task)
      call broadcast_scalar(cmc_nbuff_zone,     master_task)
      call broadcast_scalar(cmc_ln_bdyopn ,     master_task)

      IF(my_task == master_task) THEN
        WRITE(nu_diag,*)
        WRITE(nu_diag,*) 'cmc_IO s and GRID specs: '
        WRITE(nu_diag,*) '~~~~~~~~'
        WRITE(nu_diag,*) '  Namelist namcmc_nml'
        WRITE(nu_diag,*) '  GEM grid (Z) or not       cmc_ln_gemgrid  = ', cmc_ln_gemgrid
        WRITE(nu_diag,*) '  NEMO grid (X) or not      cmc_ln_nemogrid = ', cmc_ln_nemogrid
        WRITE(nu_diag,*) '  N pts in buffer zone      cmc_nbuff_zone  = ', cmc_nbuff_zone
        WRITE(nu_diag,*) '  Type of CMC boundary      cmc_ln_bdyopn   = ', cmc_ln_bdyopn
      ENDIF

      if (my_task.eq.master_task) then
        write(nu_diag,*) '-----------------------------------------  '
        write(nu_diag,*) ' CMC INPUT FIELDS AND CONSTANT OPTIONS  :  '
      endif

      call bcast_cmc_ipt(cmc_ulat,master_task)
      call bcast_cmc_ipt(cmc_ulon,master_task)
      call bcast_cmc_ipt(cmc_tlat,master_task)
      call bcast_cmc_ipt(cmc_tlon,master_task)
      call bcast_cmc_ipt(cmc_kmt, master_task)
      call bcast_cmc_ipt(cmc_hte, master_task)
      call bcast_cmc_ipt(cmc_htn, master_task)
      call bcast_cmc_ipt(cmc_gl,  master_task)
      call bcast_cmc_ipt(cmc_i8,  master_task)
      call bcast_cmc_ipt(cmc_snow,master_task)
      call bcast_cmc_ipt(cmc_uice0,master_task)
      call bcast_cmc_ipt(cmc_vice0,master_task)
      call bcast_cmc_ipt(cmc_tair,master_task)

      if (my_task.eq.master_task) then
        write(nu_diag,*) '-----------------------------------------  '
      endif

      return

      CONTAINS

      SUBROUTINE bcast_cmc_ipt(cmc_var,bc_task)

      TYPE(cmc_ipt), INTENT(inout) :: cmc_var
      INTEGER(kind=int_kind), INTENT(in) :: bc_task

      call broadcast_scalar(cmc_var%long_name,   bc_task)
      call broadcast_scalar(cmc_var%nomvar,      bc_task)
      call broadcast_scalar(cmc_var%structure,   bc_task)
      call broadcast_scalar(cmc_var%const_value, bc_task)
      call broadcast_scalar(cmc_var%field_format,bc_task)
      call broadcast_scalar(cmc_var%field_ip1_1, bc_task)
      call broadcast_scalar(cmc_var%field_ip1_2, bc_task)
      call broadcast_scalar(cmc_var%scaling,     bc_task)
      call broadcast_scalar(cmc_var%offset,      bc_task)
      if (my_task.eq.master_task) then
        if (cmc_var%long_name.ne.'null') then
          write(nu_diag,*) '|-----------------------------------|  '
          write(nu_diag,*) ' VARIABLE:  ', &
                            cmc_var%long_name
          write(nu_diag,*) 'cmc_var%nomvar= ', &
                            cmc_var%nomvar
          write(nu_diag,*) 'cmc_var%structure= ', &
                            cmc_var%structure
          write(nu_diag,*) 'cmc_var%const_value= ', &
                            cmc_var%const_value
          write(nu_diag,*) 'cmc_var%field_format= ', &
                            cmc_var%field_format
          write(nu_diag,*) 'cmc_var%field_ip1_1= ', &
                            cmc_var%field_ip1_1
          write(nu_diag,*) 'cmc_var%scaling= ', &
                            cmc_var%scaling
          write(nu_diag,*) 'cmc_var%offset= ', &
                            cmc_var%offset
          write(nu_diag,*) '|-----------------------------------|  '
        endif
      endif

      return
      END SUBROUTINE bcast_cmc_ipt
      
   END SUBROUTINE bcast_cmc_options

    !--------------------------------------------
    ! Specific routine for STD file management  
    !-------------------------- -----------------
    
      SUBROUTINE frc_openf
      implicit none
      integer maxnfile
      parameter ( maxnfile=1000 )
      character*512 filename(maxnfile),fn
      integer  fnom,fstouv,fstlnk
      external fnom,fstouv,fstlnk
      integer err,err1,err2,i,cnt,unf,wkoffit
      nfile_std = 0
      cnt       = 0
      unf       = 0

#ifdef CICE_IN_NEMO
      if (fnom(unf,'liste_inputfiles_cice','SEQ+OLD',0).lt.0) &
         call abort_ice('frc_openf: liste_inputfiles_cice opening problems')
#else
      if (fnom(unf,'liste_inputfiles','SEQ+OLD',0).lt.0) &
         call abort_ice('frc_openf: liste_inputfiles opening problems')
#endif
 77   cnt=cnt+1
      if (cnt.gt.maxnfile) then 
      write(nu_diag,*) 'maxnfile not large enough in ???? ---ABORT---'
      call flush(nu_diag)
      call abort_ice('frc_openf: maxnfile not large enough in ???? ---ABORT---')
      endif
      read (unf, '(a)', end = 9120) filename(cnt)
      goto 77
 9120 nfile_std = cnt - 1
      close(unf)
      do cnt = 1, nfile_std
         fn  = filename(cnt)
         err = wkoffit(fn)
         if ((err.ne.1).and.(err.ne.33)) then
            filename(cnt) = '@#$%^&'
         endif
      end do
      i=0
      do cnt = 1, nfile_std
         if (filename(cnt).ne.'@#$%^&') then
            i = i+1
            filename(i) = filename(cnt)
         endif
      end do
      nfile_std = i
      if (nfile_std.lt.1) then 
         write(nu_diag,*) 'NO STD FILE available for ????? ---ABORT---'
         call flush(nu_diag)
         call abort_ice('frc_openf: NO STD FILE available for ???? ---ABORT---')
      endif
      do cnt = 1, nfile_std
         iun_std(cnt) = 0
         err1 = FNOM  (iun_std(cnt),trim(filename(cnt)),'RND+OLD',0)
         if (my_task == master_task) write (nu_diag,1001) trim(filename(cnt)),iun_std(cnt)
         err2 = FSTOUV(iun_std(cnt),'RND')
         if ((err1.lt.0).or.(err2.lt.0)) call  abort_ice('frc_openf: err1 err2 ---ABORT---')
      end do
      err = fstlnk (iun_std,nfile_std)
      frc_file_L = .true.
 1001 format ('Opening STD file: ',a,2x,'UNIT= ',i4)
 1002 format ('Opening BINARY file: ',a,2x,'UNIT= ',i4)
      return

      END SUBROUTINE frc_openf

      SUBROUTINE readfst (f,nomvar,dat,niv,factm,facta,ftype)
      implicit none

      real (kind=real_kind), dimension(nx_block,ny_block,max_blocks), &
            intent(inout) :: f              ! input array (real)
      real (kind=real_kind), intent(in) :: factm,facta
      character*(*), intent(in) :: nomvar,dat,ftype
      INTEGER (kind=int_kind), intent(in) :: niv(2)

      INTEGER fstinl,fstluk,ip1_all

      integer ni1,nj1,nk1,err,nlis,lislon,i,j,k,datm,mi,mj
      parameter (nlis = 1024)
      integer liste (nlis)
      real (kind=real_kind),  allocatable,dimension(:,:) :: wrk_g
      integer (kind=int_kind),allocatable,dimension(:,:) :: wrk_gi
      character*2 typvar(4)
      integer itpv

      typvar(1)=' '
      typvar(2)='P@'
      typvar(3)='A@'
      typvar(4)='C@'
      if (my_task == master_task) then
        allocate(wrk_g (nx_global,ny_global))
        if (trim(dat).eq.'') then
          datm=-1
         else
          call datp2f (datm,dat)
        endif
        if (niv(1).eq.-1) then
          lislon=0
          do itpv=1,4
            if (lislon.lt.1) & 
     &      err = fstinl (iun_std,ni1,nj1,nk1,datm,' ', &
     &                    niv(1),-1,-1,typvar(itpv),nomvar,  &
     &                    liste,lislon,nlis)
          enddo
         else
          lislon=0
          do itpv=1,4
            if (lislon.lt.1) &
     &      err = fstinl (iun_std,ni1,nj1,nk1,datm,' ', &
     &                    niv(1),-1,-1,typvar(itpv),nomvar,  &
     &                    liste,lislon,nlis)
            if (lislon.lt.1) &
     &      err = fstinl (iun_std,ni1,nj1,nk1,datm,' ', &
     &                    niv(2),-1,-1,typvar(itpv),nomvar,  &
     &                    liste,lislon,nlis)
          enddo
        endif
        if (lislon.lt.1) then
           write(nu_diag,*) 'Variable ',nomvar,' valid: ',dat,  &
     &             ' NOT found ---ABORT---'
           call flush(nu_diag)
           call abort_ice('readfst NOT found ---ABORT---')
        endif
        if (ni1.ne.nx_global.or.nj1.ne.ny_global.or.nk1.ne.1) then
           write(nu_diag,*) 'STD FILE: WRONG DIMENSIONS'
           write(nu_diag,*) ni1,nx_global
           write(nu_diag,*) nj1,ny_global
           write(nu_diag,*) nk1
           call flush(nu_diag)
           call abort_ice('readfst STD FILE: WRONG DIMENSIONS ---ABORT---')
        endif
        if (ftype.eq.'I') then
          allocate(wrk_gi(nx_global,ny_global))
          err = fstluk (wrk_gi,liste(1),ni1,nj1,nk1)
          wrk_g=wrk_gi
          deallocate(wrk_gi)
         elseif (ftype.eq.'B') then
          err = fstluk (wrk_g ,liste(1),ni1,nj1,nk1)
         else
          write(nu_diag,*) 'readfst: unknown data format:',ftype
          call flush(nu_diag)
          call abort_ice('readfst: unknown data format:')
        endif
! Spread interior data to buffer zone
        do j=1,ny_global
          do i=nx_global-cmc_nbuff_zone+1,nx_global
              wrk_g(i,j)=wrk_g(nx_global-cmc_nbuff_zone,j)
          enddo
          do i=1,cmc_nbuff_zone
              wrk_g(i,j)=wrk_g(cmc_nbuff_zone+1,j)
          enddo
        enddo
        do i=1,nx_global
          do j=ny_global-cmc_nbuff_zone+1,ny_global
              wrk_g(i,j)=wrk_g(i,ny_global-cmc_nbuff_zone)
          enddo
          do j=1,cmc_nbuff_zone
              wrk_g(i,j)=wrk_g(i,cmc_nbuff_zone+1)
          enddo
        enddo
       else
        allocate(wrk_g (1,1))
      endif

!     This operation scatter all data (water+land+ghost)
      call scatter_global(f, wrk_g, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      f   = f*factm + facta

      deallocate(wrk_g)
      return
      END SUBROUTINE readfst

      SUBROUTINE readfst_global (g,nomvar,dat,niv,factm,facta,ftype,lx,ly)
      implicit none

      real (kind=real_kind), dimension(nx_global,ny_global), &
            intent(inout) :: g              ! input array (real)
      real (kind=real_kind), intent(in)    :: factm,facta
      character*(*), intent(in) :: nomvar,dat,ftype
      INTEGER (kind=int_kind), intent(in) :: niv(2)
      LOGICAL, optional                   :: lx,ly

      integer (kind=int_kind), dimension(nx_global,ny_global)  &
                          :: gi             ! input array (integer)

      INTEGER fstinl,fstluk,ip1_all

      integer ni1,nj1,nk1,err,nlis,lislon,i,j,datm,mi,mj
      parameter (nlis = 1024)
      integer liste (nlis)
      LOGICAL                             :: lxw,lyw
      integer nxt,nyt

      if (my_task == master_task) then
        lxw=.false.
        lyw=.false.
        if (present(lx)) lxw=lx
        if (present(ly)) lyw=ly
        if (trim(dat).eq.'') then
          datm=-1
         else
          call datp2f (datm,dat)
        endif
        if (niv(1).eq.-1) then
          err = fstinl (iun_std,ni1,nj1,nk1,datm,' ',niv(1),-1,-1,' ',nomvar,  &
     &                                                   liste,lislon,nlis)
         else
          err = fstinl (iun_std,ni1,nj1,nk1,datm,' ',niv(1),-1,-1,' ',nomvar,  &
     &                                                   liste,lislon,nlis)
          if (lislon.lt.1) &
     &    err = fstinl (iun_std,ni1,nj1,nk1,datm,' ',niv(2),-1,-1,' ',nomvar,  &
     &                                                   liste,lislon,nlis)
        endif
        if (lislon.lt.1) then
           write(nu_diag,*) 'Variable ',nomvar,' valid: ',dat,  &
     &             ' NOT found ---ABORT---'
           call flush(nu_diag)
           call abort_ice('readfst_global nomvar NOT found ---ABORT---')
        endif
        nxt=nx_global
        nyt=ny_global
        if (lxw) nyt=1
        if (lyw) nxt=1
        if (ni1.ne.nxt.or.nj1.ne.nyt.or.nk1.ne.1) then
           write(nu_diag,*) 'STD FILE: WRONG DIMENSIONS'
           write(nu_diag,*) ni1,nxt,nx_global
           write(nu_diag,*) nj1,nyt,ny_global
           write(nu_diag,*) nk1
           call abort_ice('readfst_global : STD FILE: WRONG DIMENSIONS ---ABORT---')
        endif
        if (ftype.eq.'I') then
          if (lxw) then
            gi=0
            err = fstluk (gi(:,1),liste(1),ni1,nj1,nk1)
            g=gi
           elseif (lyw) then
            gi=0
            err = fstluk (gi(1,:),liste(1),ni1,nj1,nk1)
            g=gi
           else
            err = fstluk (gi,     liste(1),ni1,nj1,nk1)
            g=gi
          endif
         elseif (ftype.eq.'B') then
          if (lxw) then
            g=0.0_real_kind
            err = fstluk (g(:,1),liste(1),ni1,nj1,nk1)
           elseif (lyw) then
            g=0.0_real_kind
            err = fstluk (g(1,:),liste(1),ni1,nj1,nk1)
           else
            err = fstluk (g,     liste(1),ni1,nj1,nk1)
          endif
         else
          write(nu_diag,*) 'readfst_global: unknown data format:',ftype
          call flush(nu_diag)
          call abort_ice('readfst_global: unknown data format:')
        endif
! Spread interior data to buffer zone
! (only for mask kmt, the rest is model grid geometry)
        if (nomvar == cmc_kmt%nomvar) then
          do j=1,ny_global
            do i=nx_global-cmc_nbuff_zone+1,nx_global
                g(i,j)=g(nx_global-cmc_nbuff_zone,j)
            enddo
            do i=1,cmc_nbuff_zone
                g(i,j)=g(cmc_nbuff_zone+1,j)
            enddo
          enddo
          do i=1,nx_global
            do j=ny_global-cmc_nbuff_zone+1,ny_global
                g(i,j)=g(i,ny_global-cmc_nbuff_zone)
            enddo
            do j=1,cmc_nbuff_zone
                g(i,j)=g(i,cmc_nbuff_zone+1)
            enddo
          enddo
          if (cmc_ln_nemogrid) then
            if (trim(ew_boundary_type)=='closed') then
              do j=1,ny_global
                do i=nx_global-2+1,nx_global
                    g(i,j)=c0
                enddo
                do i=1,2
                    g(i,j)=c0
                enddo
              enddo
            endif
            if (trim(ns_boundary_type)=='closed') then
              do i=1,nx_global
                do j=ny_global-2+1,ny_global
                    g(i,j)=c0
                enddo
                do j=1,2
                    g(i,j)=c0
                enddo
              enddo
            endif
          endif
        endif

        g   = g*factm + facta

      endif

      return
      END SUBROUTINE readfst_global

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

#endif
END MODULE ice_fldread_std
