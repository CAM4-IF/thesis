
module convect_deep
!---------------------------------------------------------------------------------
! Purpose:
!
!  CAM interface to several deep convection interfaces. Currently includes:
!    IF (default)
!  removed ZM
!
! Author: D.B. Coleman, Sep 2004
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use cam_logfile,  only: iulog

   implicit none

   save
   private                         ! Make default type private to the module
!------------------------------------------------------------------------
!
!  Public methods
!
!  removed: convect_deep_tend_2,deep_scheme_does_scav_trans
!
!------------------------------------------------------------------------
   public ::&
      convect_deep_register,           &! register fields in physics buffer
      convect_deep_init,               &! initialize donner_deep module
      convect_deep_tend                 ! return tendencies
   
!------------------------------------------------------------------------
!
! Private module data
!
!------------------------------------------------------------------------
   character(len=16) :: deep_scheme    ! default set in phys_control.F90, use namelist to change

!=========================================================================================
  contains 
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  Removed:
!   - "function deep_scheme_does_scav_trans()" since if loop where
!     used in tphysbc was removed.
!   - any place where it is needed?
!
!  Function was called by tphysbc to determine if it needs to do scavenging and
!    convective transport or if those have been done by the deep convection
!    scheme.
!    Each scheme could have its own identical query function for a
!    less-knowledgable
!    interface but for now, we know that KE does scavenging & transport, and ZM
!    doesn't
!
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  Calls subroutine phys_getopts "get physics options"
!  - added 'IF' case
!  - probably only called at the very beginning of a model run
!
! Purpose: register fields with the physics buffer
!
!  This program is called from "subroutine initindx.F90".
!  However, "subroutine initindx.F90" does not appear to be called from
!    anywhere inside the physics package, so presumably earlier, in the
!    so called "build stage".
!
!------------------------------------------------------------------------
subroutine convect_deep_register

!------------------------------------------------------------------------
!
!   Use statements
!   changed zm -> if here
!
!------------------------------------------------------------------------
  use phys_buffer, only: pbuf_times, pbuf_add
  use if_conv_intr, only: if_conv_register
!------------------------------------------------------------------------
!
!  phys_getopts is where "deep_scheme_out" is defined:
!   - character(len=16), intent(out), optional :: deep_scheme_out
!   - need to define as 'IF'
!
!------------------------------------------------------------------------
  use phys_control, only: phys_getopts

  implicit none

  integer idx
!------------------------------------------------------------------------
!
!  Get deep_scheme: deep convection option
!
!  subroutine phys_getopts doesn't seem to do much, beyond making an
!   equivalence between deep_scheme_out and deep_scheme
!   however "deep_scheme" is initialized as:
!   character(len=16) :: deep_scheme     = unset_str  ! deep convection package
!   in top of "module phys_control", where unset_str = 'UNSET'.
!   "deep_scheme" is a namelist variable; see:
!   http://www.ccsm.ucar.edu/cgi-bin/eaton/namelist/nldef2html-pub
!   where it says: ZM is the only option supported in CAM4.
!      Default: 'ZM'
!   I have to introduce new 'IF' for 'deep_scheme'
!   May be possible to just hardwire 'IF' in rather than via namelist
!   Should run in "aqua planet" mode
!
! get "deep_scheme" setting from phys_control
!
!------------------------------------------------------------------------
  print*,'calling phys_getopts from convect_deep_register'
  call phys_getopts(deep_scheme_out = deep_scheme)

!------------------------------------------------------------------------
!
!  Replaced 'ZM' case with 'IF'
!  - if_conv_register is where I register anvil liquid, and anvil age.
!  - presumably done for the pbuf_add commands below
!
!------------------------------------------------------------------------
  deep_scheme = 'IF'
  select case ( deep_scheme )
  case('IF') !    IF convection
     print*,'calling if_conv_register from convect_deep_register'
     call if_conv_register
  end select

!------------------------------------------------------------------------
!
!  buffer calls added by Ian Folkins (IF)
!  advice on how to add is in:
!  http://www.cesm.ucar.edu/models/atm-cam/docs/phys-interface/phys-interface.pdf
!
!  These are all 2D fields
!------------------------------------------------------------------------
  call pbuf_add('PRECIPORG' , 'global', 1,1,     1, idx)
  call pbuf_add('MASSPREV' , 'global', 1,1,     1, idx)
  call pbuf_add('CM5PREV' , 'global', 1,1,     1, idx)

!  pbuf_add(name, scope, fdim, mdim, ldim, index)
!------------------------------------------------------------------------
!
!  name 'ICWMRDP' later associated with pointer ql in zm_conv_intr.F90
!  name 'RPRDDP' later associated with pointer rprd in zm_conv_intr.F90
!  I don't need either of these I don't think, but leave in for now.
!
!------------------------------------------------------------------------
  call pbuf_add('ICWMRDP' , 'physpkg', 1,pver,      1, idx)
  call pbuf_add('RPRDDP' , 'physpkg', 1,pver,      1, idx)
!------------------------------------------------------------------------
!  this seems to be a cam5 addition
!------------------------------------------------------------------------
  call pbuf_add('NEVAPR_DPCU' , 'physpkg', 1,pver,      1, idx)

end subroutine convect_deep_register

!=========================================================================================


!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
! Looks like a pretty simple subroutine
! Purpose:  declare output fields, initialize variables needed by convection
! - needs to know the value of 'deep_scheme', which should be 'IF'
!
!------------------------------------------------------------------------
subroutine convect_deep_init(hypi)

  use pmgrid,       only: plevp
  use spmd_utils,   only: masterproc
!  use zm_conv_intr, only: zm_conv_init
  use if_conv_intr, only: if_conv_init
  use abortutils,   only: endrun

  implicit none

  real(r8),intent(in) :: hypi(plevp)        ! reference pressures at interfaces

  integer k

!------------------------------------------------------------------------
!
!  Need to add 'IF' case : how to do this?
!  Why isn't Emanuel case here?
!  - removed ZM case, in case couldn't find subroutines
!
!------------------------------------------------------------------------
  select case ( deep_scheme )
  case('off') !     ==> no deep convection
     if (masterproc) write(iulog,*)'convect_deep: no deep convection selected'
  case('ZM') !
     print*,'case is ZM: may fail'
  case('IF') !
     if (masterproc) write(iulog,*)'convect_deep initializing IF convection'
     call if_conv_init(hypi)
  case default
     if (masterproc) write(iulog,*)'WARNING: convect_deep: no deep convection scheme. May fail.'
  end select

end subroutine convect_deep_init
!=========================================================================================
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  Original call was:
!
!  subroutine convect_deep_tend(prec    , &
!    pblh    ,mcon    ,cme     ,          &
!    tpert   ,dlf     ,pflx    ,zdu      , &
!    rliq    , &
!    ztodt   ,snow    ,&
!    state   ,ptend   ,landfrac ,pbuf  )
!
!
!------------------------------------------------------------------------
subroutine convect_deep_tend(rain_if  ,snow_if    ,&
            ztodt   ,state  ,ptend    ,landfrac   ,&
            pbuf)
!------------------------------------------------------------------------
!
!  Use statements:
!
!------------------------------------------------------------------------
   use physics_types, only: physics_state, physics_ptend, physics_tend, physics_ptend_init
   use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_get_fld_idx
   use constituents, only: pcnst
!------------------------------------------------------------------------
!  zm -> if
!------------------------------------------------------------------------
   use if_conv_intr, only: if_conv_tend
#if ( defined WACCM_PHYS )
   use gw_drag,         only: idx_zmdt
   use cam_history,     only: outfld
   use physconst,       only: cpair
#endif
!------------------------------------------------------------------------
!
!  IN/OUT Arguments
!  removed: pblh,tpert,mcon,dlf,pflx,cme,zdu,rliq
!
!------------------------------------------------------------------------
! Arguments
   real(r8), intent(out) :: rain_if(pcols)   ! total IF rain at surface
   real(r8), intent(out) :: snow_if(pcols)   ! total IF snow at surface
   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
   real(r8), intent(in) :: landfrac(pcols)                ! Land fraction
   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer
!------------------------------------------------------------------------
!
!  Pointers
!  IF - nothing removed here; evapcdp is a cam5 addition; can probably remove
!
!------------------------------------------------------------------------
   real(r8), pointer, dimension(:) :: jctop
   real(r8), pointer, dimension(:) :: jcbot
   real(r8), pointer, dimension(:,:,:) :: cld        
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation

!------------------------------------------------------------------------
!
!  IF added pointers 
!
!------------------------------------------------------------------------
   real(r8), pointer, dimension(:) :: preciporg
   real(r8), pointer, dimension(:) :: massprev
   real(r8), pointer, dimension(:) :: cm5prev

!------------------------------------------------------------------------
!
!  Other variables
!
!------------------------------------------------------------------------
  real(r8) zero(pcols, pver)

  integer i, k
  integer ifld

!---------------------------------------------------------------------------------
!
!  ??
!
!---------------------------------------------------------------------------------
#if ( defined WACCM_PHYS )
   real(r8), pointer, dimension(:,:) :: zmdt
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
#endif

!---------------------------------------------------------------------------------
!
!  Appears to make an association between jctop and 'CLDTOP'
!
!---------------------------------------------------------------------------------
   ifld = pbuf_get_fld_idx('CLDTOP')
   jctop => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)
   ifld = pbuf_get_fld_idx('CLDBOT')
   jcbot => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)

!---------------------------------------------------------------------------------
!
!   IF added 
!
!---------------------------------------------------------------------------------
   ifld = pbuf_get_fld_idx('PRECIPORG')
   preciporg => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)
   ifld = pbuf_get_fld_idx('MASSPREV')
   massprev => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)
   ifld = pbuf_get_fld_idx('CM5PREV')
   cm5prev => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)

!---------------------------------------------------------------------------------
!
!  Start cases:
!   - add 'IF' case; deep_scheme specified in namelist file
!
!---------------------------------------------------------------------------------
  select case ( deep_scheme )
!---------------------------------------------------------------------------------
!
!  when namelist variable deep_scheme = 'off'
!  - zeroing a bunch of arrays
!  - I removed: (mcon, dlf, pflx, cme, zdu, rliq) since I no longer
!    define these arrays at this level (only at one level higher).
!  - wanted to remove any reference to ZM at this level and lower.
!
!---------------------------------------------------------------------------------
  case('off') !    0 ==> no deep convection
    zero = 0     
    rain_if = 0
    snow_if = 0
!---------------------------------------------------------------------------------
!
! Associate pointers with physics buffer fields
! - think I should remove rprd and fracis; not defining
! - all for no convection case
! - evapcdp is a cam5 addition
!
!---------------------------------------------------------------------------------
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,state%lchnk,:) 
   ifld = pbuf_get_fld_idx('ICWMRDP')
   ql => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1)
   ifld = pbuf_get_fld_idx('RPRDDP')
   rprd => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1)
   ifld = pbuf_get_fld_idx('FRACIS')
   fracis  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1:pcnst)
   ifld = pbuf_get_fld_idx('NEVAPR_DPCU')
   evapcdp => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1)
!---------------------------------------------------------------------------------
!
!  remove some of these?
!  - Notice that the are all pointers; must be essentially saying that just
!  don't refer to anything.
!  - however odd that just associated above ...
!  (Remember in 'OFF' convection case still)
!
!---------------------------------------------------------------------------------
    jctop = 0
    jcbot = 0
    cld = 0
    ql = 0
    rprd = 0
    fracis = 0
    evapcdp = 0
!---------------------------------------------------------------------------------
!
!  zeroes tendency array ptend (tendencies from one physics param only)
!  (Remember in 'OFF' convection case still)
!
!---------------------------------------------------------------------------------
!  print*,'calling physics_ptend_init'
   call physics_ptend_init(ptend)
   ptend%name = "convect_deep"
!  print*,'finished calling physics_ptend_init'
!---------------------------------------------------------------------------------
!
!  ZM case: removed
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!
!  IF case:
!
!  original call to ZM was:
!
!    call zm_conv_tend(prec     , &
!         pblh    ,mcon    ,cme     ,          &
!         tpert   ,dlf     ,pflx    ,zdu      , &
!         rliq    , &
!         ztodt   ,snow    ,&
!         jctop, jcbot , &
!         state   ,ptend   ,landfrac ,pbuf  )
!
!  Note that jctop,jcbot are pointers here (not arrays as lower down).
!
!---------------------------------------------------------------------------------
  case('IF') !    IF convection
!   print*,'calling if_conv_tend'
    call if_conv_tend(rain_if     ,snow_if            ,&
         ztodt   ,jctop     ,jcbot     ,state         ,&
         ptend   ,landfrac  ,pbuf                     ,&
         preciporg, massprev, cm5prev) 
!   print*,'finished calling if_conv_tend'
!---------------------------------------------------------------------------------
!
!  end cases
!
!---------------------------------------------------------------------------------
  end select

!---------------------------------------------------------------------------------
!
!  Why outputting T/Q tendencies here?
!  commented this out since might cause problems, with 'ZMDT undefined
!
!---------------------------------------------------------------------------------
!#if ( defined WACCM_PHYS )
!   zmdt  => pbuf(idx_zmdt) %fld_ptr(1,:,:,state%lchnk,1)
!   ftem(:state%ncol,:pver) = ptend%s(:state%ncol,:pver)/cpair
!   zmdt(:state%ncol,:pver) = ftem(:state%ncol,:pver)
!   call outfld('ZMDT    ',ftem           ,pcols   ,state%lchnk   )
!   call outfld('ZMDQ    ',ptend%q(:,:,1) ,pcols   ,state%lchnk   )
!#endif
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!  end
!
!------------------------------------------------------------------------
!print*,'exiting convect_deep_tend'
end subroutine convect_deep_tend
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  removed:  subroutine convect_deep_tend_2( state,  ptend,  ztodt, pbuf  )
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
end module
