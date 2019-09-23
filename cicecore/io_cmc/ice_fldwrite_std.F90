MODULE ice_fldwrite_std
#ifdef FUTURE_MAYBE
   !!======================================================================
   !!                       ***  MODULE  ice_fldwrite_std  ***
   !! Ice outputs:  write in std files
   !!=====================================================================
   !! History :  9.0  !  Nov 2010  (F. Roy) 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   icecmc :        write output global fields
   !!----------------------------------------------------------------------

   use mod_gemdate
   use ice_constants
   use ice_calendar, only: istep1,istep0,dats_curr,dats_zero,dt,cmc_ofmt
   use ice_kinds_mod
   use ice_domain_size, only: nx_global, ny_global
   use ice_domain, only: ns_boundary_type
   use ice_exit, only: abort_ice

   IMPLICIT NONE
   PRIVATE   

   integer(kind=int_kind), public :: &
                         cmc_npako, &      ! npak for outputs
                         cmc_npakri, &     ! npak for restarts (analysis in)
                         cmc_npakro        ! npak for restarts (trials out)

   character(LEN=12), public :: cmc_etiket  ! etiket for outputs and restarts

   integer(kind=int_kind) :: &
                         nbuf

! -------------------------------------------------
! Working variables
! -------------------------------------------------

   real(kind=real_kind), &
                         dimension(0:nx_global+1,0:ny_global+1) :: &
                         fldw

   integer(kind=int_kind), &
                         dimension(0:nx_global+1,0:ny_global+1) :: &
                         mskw

   PUBLIC   putfld_cmc, wrap_tripole, wrap_tripole8

CONTAINS

   SUBROUTINE putfld_cmc(cmfile,nomvar,fld,msk,lat,lon,nx,ny,w_etiket,vec)
      implicit none

!I/O
      character(len=*), intent(in) :: cmfile
      character(len=4), intent(in) :: nomvar
      character(len=12) :: w_etiket
      logical, intent(in)          :: vec
      integer (kind=int_kind), intent(in) :: nx,ny
      real   (kind=real_kind), intent(inout), &
              dimension(nx,ny) :: fld
      real   (kind=real_kind), intent(in), &
              dimension(nx,ny) :: lat,lon
      integer(kind=int_kind),  intent(in), &
              dimension(nx,ny) :: msk

!File management
      character(len=300), SAVE :: cmfopn='null'
      logical :: tobeopn, tobecls

!STD file stuff
      character(len=1),       SAVE :: grtyp_o
      integer(kind=int_kind), SAVE :: iunit
      integer(kind=int_kind), SAVE :: ip1,ip2,ip3,ig1,ig2,ig3,ig4
      integer(kind=int_kind), SAVE :: datyp,deet,dateo,datev
      character(len=2),       SAVE :: typvar,typvarm
 
      integer(kind=int_kind) :: fnom,fstouv,fstecr,fclos,fstfrm
      external                  fnom,fstouv,fstecr,fclos,fstfrm
      integer(kind=int_kind) :: err,npas,npak
      real(kind=real_kind)   :: dumwrk
      logical                :: rewrit
      character(len=4)       :: nompp
      real(kind=real_kind)   :: zlev=c0
      character(len=8)       :: dumc
      character(len=1)       :: grtyp
      character(len=12)       :: etiket

      integer(kind=int_kind) :: i,j,nxb,nyb,nxo,nyo, &
                                    ixo1,iyo1,ixo2,iyo2
! masking land value
      real  (kind=real_kind) :: pavg
      integer(kind=int_kind) :: navg

      if (nx.ne.nx_global) call abort_ice('Wrong nx value in putfld_cmc')
      if (ny.ne.ny_global) call abort_ice('Wrong ny value in putfld_cmc')

      if (trim(cmfopn)=='null') then
        tobeopn=.true.
        tobecls=.false.
       elseif (trim(cmfopn).eq.trim(cmfile)) then
        tobeopn=.false.
        tobecls=.false.
       else
        tobeopn=.true.
        tobecls=.true.
      endif

      if (tobecls) then
        err=fstfrm(iunit)
        err=fclos (iunit)
        if (trim(cmfile).eq.'close') then
          cmfopn='null'
          return
        endif
      endif

      nbuf=cmc_nbuff_zone
      nxb=nx-2*nbuf
      nyb=ny-2*nbuf

      if (tobeopn) then

!... opening new file

        cmfopn=trim(cmfile)

        if (cmc_ln_gemgrid) then
          grtyp_o='Z'
         elseif (cmc_ln_nemogrid) then
          grtyp_o='X'
         else
          grtyp_o='X'
        endif

        dateo=dats_zero
        npas=0
        deet=0
        npak=0
        rewrit=.true.
        iunit=0
        err=fnom(iunit, cmfile, 'RND', 0)
        if (err < 0) call abort_ice('error in fnom')
        err=fstouv(iunit, 'RND')
        if (err < 0) call abort_ice('error in fstouv')

        if (grtyp_o.eq.'X') then
          ip1=1001
          ip2=1002
          ip3=1003
          typvar='PP'
          datyp=5
          grtyp='L'
          call cxgaig (grtyp, ig1, ig2, ig3, ig4, 0.0_real_kind, 0.0_real_kind, &
      &                                           1.0_real_kind, 1.0_real_kind)
          nompp='>>'
          etiket=w_etiket
          if (cmc_ln_nemogrid.and. &
              (trim(ns_boundary_type)=='tripole'.or. &
               trim(ns_boundary_type)=='tripoleT')) then
            if (nbuf /= 0) &
             call abort_ice('putfld_cmc:TRIPOLE GRID nbuf should be 0')
            ixo1=0
            ixo2=nxb+1
            iyo1=1
            iyo2=nyb+1
            nxo=nxb+2
            nyo=nyb+1
            call wrap_tripole(cmclon(1:nxb,1:nyb),msk(1:nxb,1:nyb), &
                          ns_boundary_type,nxb,nyb,vec, &
                          fldw(0:nxb+1,1:nyb+1),mskw(0:nxb+1,1:nyb+1))
           else
            ixo1=nbuf+1
            ixo2=nxb+nbuf
            iyo1=nbuf+1
            iyo2=nyb+nbuf
            nxo=nxb
            nyo=nyb
            fldw(    ixo1:ixo2,iyo1:iyo2) = &
              cmclon(ixo1:ixo2,iyo1:iyo2)
          endif
          err = fstecr (fldw(ixo1:ixo2,iyo1:iyo2),  &
              dumwrk, npak, iunit, &
              dateo, deet, npas, &
              nxo, nyo, 1, ip1, ip2, ip3, &
              typvar, nompp, etiket, &
              grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
          if (err < 0) call abort_ice('error in fstecr : Longitude')
          nompp='^^'
          etiket=w_etiket
          if (cmc_ln_nemogrid.and. &
              (trim(ns_boundary_type)=='tripole'.or. &
               trim(ns_boundary_type)=='tripoleT')) then
            ixo1=0
            ixo2=nxb+1
            iyo1=1
            iyo2=nyb+1
            nxo=nxb+2
            nyo=nyb+1
            call wrap_tripole(cmclat(1:nxb,1:nyb),msk(1:nxb,1:nyb), &
                          ns_boundary_type,nxb,nyb,vec, &
                          fldw(0:nxb+1,1:nyb+1),mskw(0:nxb+1,1:nyb+1))
           else
            ixo1=nbuf+1
            ixo2=nxb+nbuf
            iyo1=nbuf+1
            iyo2=nyb+nbuf
            nxo=nxb
            nyo=nyb
            fldw(    ixo1:ixo2,iyo1:iyo2) = &
              cmclat(ixo1:ixo2,iyo1:iyo2)
          endif
          err = fstecr (fldw(ixo1:ixo2,iyo1:iyo2),  &
              dumwrk, npak, iunit, &
              dateo, deet, npas, &
              nxo, nyo, 1, ip1, ip2, ip3, &
              typvar, nompp, etiket, &
              grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
          if (err < 0) call abort_ice('error in fstecr : Latitude')
        endif

        if (grtyp_o.eq.'Z') then
!          dateo= cmc_tt_dateo
!          deet=  cmc_tt_deet
!          npas=  cmc_tt_npas
          ip1=   cmc_tt_ip1
          ip2=   cmc_tt_ip2
          ip3=   cmc_tt_ip3
          ig1=   cmc_tt_ig1
          ig2=   cmc_tt_ig2
          ig3=   cmc_tt_ig3
          ig4=   cmc_tt_ig4
          datyp= cmc_tt_datyp
          typvar=cmc_tt_typvar
          grtyp= cmc_tt_grtyp
          
          nompp='>>'
          err = fstecr (cmclx(nbuf+1:nxb+nbuf,1),  &
              dumwrk, npak, iunit, &
              dateo, deet, npas, &
              nxb, 1, 1, ip1, ip2, ip3, &
              typvar, nompp, w_etiket, &
              grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
          if (err < 0) call abort_ice('error in fstecr : Longitude (Z)')
          nompp='^^'
          err = fstecr (cmcly(1,nbuf+1:nyb+nbuf),  &
              dumwrk, npak, iunit, &
              dateo, deet, npas, &
              1, nyb, 1, ip1, ip2, ip3, &
              typvar, nompp, w_etiket, &
              grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
          if (err < 0) call abort_ice('error in fstecr : Latitude (Z)')
        endif

!... preparing parameters for new fields
        call convip (ip1,zlev,0,+2,dumc,.false.)
        deet   =dt
!        datyp  =134
        datyp  =133
        if (npak.gt.-32) datyp=134 ! FD add extension for f16
        typvar ='P@'
        typvarm='@@'
        rewrit =.true.
        ip3=0
        if ( grtyp_o .eq. 'X') then
         ig1=1001
         ig2=1002
         ig3=1003
         ig4=0
        endif
        if ( grtyp_o .eq. 'Z') then
         ig1=cmc_tt_ip1
         ig2=cmc_tt_ip2
         ig3=cmc_tt_ip3
         ig4=0
        endif

      endif

      npas =istep1-istep0
      ip2  =(npas*deet)/c3600
      if (cmc_ofmt.eq.'H') then  ! Hindcast style
        dateo=dats_curr
        npas=0
        ip2=0
      endif

! set value over land
!froy
!Write field average on land to optimize compression
      pavg=0.0_real_kind
      navg=0
      do j=1,ny
      do i=1,nx
        if (msk(i,j).ne.0) then
          pavg = pavg + fld(i,j)
          navg = navg + 1
        endif
      enddo
      enddo
      if (navg.gt.0) pavg = pavg / navg
      do j=1,ny
      do i=1,nx
        if (msk(i,j).eq.0) fld(i,j) = pavg
      enddo
      enddo

      npak=cmc_npako
      grtyp=grtyp_o
      if (cmc_ln_nemogrid.and. &
          (trim(ns_boundary_type)=='tripole'.or. &
           trim(ns_boundary_type)=='tripoleT')) then
        ixo1=0
        ixo2=nxb+1
        iyo1=1
        iyo2=nyb+1
        nxo=nxb+2
        nyo=nyb+1
        call wrap_tripole(fld(1:nxb,1:nyb),msk(1:nxb,1:nyb), &
                          ns_boundary_type,nxb,nyb,vec, &
                          fldw(0:nxb+1,1:nyb+1),mskw(0:nxb+1,1:nyb+1))
       else
        ixo1=nbuf+1
        ixo2=nxb+nbuf
        iyo1=nbuf+1
        iyo2=nyb+nbuf
        nxo=nxb
        nyo=nyb
        fldw(    ixo1:ixo2,iyo1:iyo2) = &
          fld(   ixo1:ixo2,iyo1:iyo2)
        mskw(    ixo1:ixo2,iyo1:iyo2) = &
          msk(   ixo1:ixo2,iyo1:iyo2)
      endif
      err=fstecr(fldw(ixo1:ixo2,iyo1:iyo2), &
              dumwrk,npak,iunit,dateo,deet,npas, &
              nxo,nyo,1,ip1,ip2,ip3,typvar,nomvar,w_etiket, &
              grtyp,ig1,ig2,ig3,ig4,datyp,rewrit)
      if (err < 0) call abort_ice('error in fstecr : fld')

      err=fstecr(mskw(ixo1:ixo2,iyo1:iyo2), &
              dumwrk,-1,iunit,dateo,deet,npas, &
              nxo,nyo,1,ip1,ip2,ip3,typvarm,nomvar,w_etiket, &
              grtyp,ig1,ig2,ig3,ig4,2,rewrit)
      if (err < 0) call abort_ice('error in fstecr : msk')

   END SUBROUTINE putfld_cmc

   SUBROUTINE wrap_tripole(fld,msk,ns_boundary_type,nx,ny,v,fldw,mskw)

   character (char_len),    intent(in) :: ns_boundary_type
   integer (kind=int_kind), intent(in) :: nx,ny
   logical ,                intent(in) :: v
   real   (kind=real_kind), intent(in),  &
           dimension(nx,ny)            :: fld
   integer(kind=int_kind),  intent(in),  &
           dimension(nx,ny)            :: msk
   real   (kind=real_kind), intent(out), &
           dimension(0:nx+1,ny+1)      :: fldw
   integer(kind=int_kind),  intent(out), &
           dimension(0:nx+1,ny+1)      :: mskw

   integer i,j,isign,isrc,jsrc,xoffset,yoffset

   isign=1
   if (v) isign=-1  !vector case

   if (trim(ns_boundary_type)=='tripoleT') then
!froy - taken from ice_gather_scatter.F90
!!!!     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 3 !different logic than ice_gather_scatter.F90
!!!!     end select
    elseif (trim(ns_boundary_type)=='tripole') then
!froy - taken from ice_gather_scatter.F90
!!!!     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 2 !different logic than ice_gather_scatter.F90
!!!!     end select
    else
     call abort_ice('wrap_tripole: wrong type of boundary')
   endif

   fldw(1:nx,1:ny)=fld(1:nx,1:ny)
   mskw(1:nx,1:ny)=msk(1:nx,1:ny)

   fldw(0   ,1:ny)=fld(nx  ,1:ny)
   mskw(0   ,1:ny)=msk(nx  ,1:ny)
   fldw(nx+1,1:ny)=fld(1   ,1:ny)
   mskw(nx+1,1:ny)=msk(1   ,1:ny)

   jsrc = ny - yoffset
   do i=0,nx+1
     isrc = nx + xoffset - i
     if (isrc < 1 ) isrc = isrc + nx
     if (isrc > nx) isrc = isrc - nx
     fldw(i,ny+1) = isign * fld(isrc,jsrc)
     mskw(i,ny+1) =         msk(isrc,jsrc)
   end do

   END SUBROUTINE wrap_tripole

   SUBROUTINE wrap_tripole8(fld,msk,ns_boundary_type,nx,ny,v,fldw,mskw)

   character (char_len),    intent(in) :: ns_boundary_type
   integer (kind=int_kind), intent(in) :: nx,ny
   logical ,                intent(in) :: v
   real   (kind=dbl_kind), intent(in),  &
           dimension(nx,ny)            :: fld
   integer(kind=int_kind),  intent(in),  &
           dimension(nx,ny)            :: msk
   real   (kind=dbl_kind), intent(out), &
           dimension(0:nx+1,ny+1)      :: fldw
   integer(kind=int_kind),  intent(out), &
           dimension(0:nx+1,ny+1)      :: mskw

   integer i,j,isign,isrc,jsrc,xoffset,yoffset

   isign=1
   if (v) isign=-1  !vector case

   if (trim(ns_boundary_type)=='tripoleT') then
!froy - taken from ice_gather_scatter.F90
!!!!     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 3 !different logic than ice_gather_scatter.F90
!!!!     end select
    elseif (trim(ns_boundary_type)=='tripole') then
!froy - taken from ice_gather_scatter.F90
!!!!     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 2 !different logic than ice_gather_scatter.F90
!!!!     end select
    else
     call abort_ice('wrap_tripole: wrong type of boundary')
   endif

   fldw(1:nx,1:ny)=fld(1:nx,1:ny)
   mskw(1:nx,1:ny)=msk(1:nx,1:ny)

   fldw(0   ,1:ny)=fld(nx  ,1:ny)
   mskw(0   ,1:ny)=msk(nx  ,1:ny)
   fldw(nx+1,1:ny)=fld(1   ,1:ny)
   mskw(nx+1,1:ny)=msk(1   ,1:ny)

   jsrc = ny - yoffset
   do i=0,nx+1
     isrc = nx + xoffset - i
     if (isrc < 1 ) isrc = isrc + nx
     if (isrc > nx) isrc = isrc - nx
     fldw(i,ny+1) = isign * fld(isrc,jsrc)
     mskw(i,ny+1) =         msk(isrc,jsrc)
   end do

   END SUBROUTINE wrap_tripole8

#endif
END MODULE ice_fldwrite_std
