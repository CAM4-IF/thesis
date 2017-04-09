module if_conv
!---------------------------------------------------------------------------------
!
! Purpose:
!
! Interface with IF convection scheme.
!
! The arrays rvtend,rltend, and ritend 
! are sufficent to determine the change in surface pressure if
! conservation of dry air is invoked as an additional constraint. But it seems
! the changes in surface pressure from precipitation are not considered.
!
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
  use abortutils,      only: endrun
  use wv_saturation,   only: estblf, hlatv, tmin, hlatf, rgasv, pcf, &
                            cp, epsqs, ttrice
!------------------------------------------------------------------------
!
! IF added: to feed into if_conv_tend.f90 and output stuff on last timestep
!
!------------------------------------------------------------------------
  use time_manager,       only: get_nstep
!------------------------------------------------------------------------
!
! use statements to my modules 
! - should be made more specific at some point
!
!------------------------------------------------------------------------
  use if_conv_params          
  use if_conv_tend
  use if_conv_solvers
!------------------------------------------------------------------------
!
!  Added this so would know what cnst_get_ind  was. 
!  Took from check_energy.F90; may not need all.
!
!------------------------------------------------------------------------
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind
!------------------------------------------------------------------------
!
!  Should change over to using cam4 physical constants
!  - module needs access to physconst to avoid errors.
!
!------------------------------------------------------------------------
  use physconst,       only: cpair,gravit
  use cam_logfile,     only: iulog

  implicit none
!------------------------------------------------------------------------
!
!  not sure what this does
!
!------------------------------------------------------------------------
  save
!------------------------------------------------------------------------
!
!  I should retain this (but if make only subroutine public, probably 
!    doesn't do anything).
!
!------------------------------------------------------------------------
!  private                         ! Make default type private to the module
!------------------------------------------------------------------------
!
!  public interfaces
!
! The ZM scheme had:
!  public zm_convi                 ! ZM schemea
!  public zm_convr                 ! ZM schemea
!
!  This would allow other subroutines which had a "use zm_conv" command to
!    access variables in subroutines zm_convi, and zm_convr, but no others
!    (e.g. subroutine zm_conv_evap etc).
!
!  I think the default would be for me to do the same, though this is perhaps
!    silly because these are the only subroutines in this module anyway.
!
!  cam5 has a new interface: zmconv_readnl (Note the "l" is not a one)
!
!------------------------------------------------------------------------
!  public if_convi                 ! ZM schemea
!  public zmconv_readnl            ! read zmconv_nl namelist
                                   ! assume these three tuning parameters not
                                   ! needed elsewhere.
  public if_convr
! public cldwat               
!------------------------------------------------------------------------
!
! Private data
! "By default all module variables are available by all program units 
! which USE the module"
!
! moved from moistconvection.F90
!
!------------------------------------------------------------------------
   real(r8) :: rgrav       ! reciprocal of grav
   real(r8) :: grav        ! = gravit
!------------------------------------------------------------------------
!
!  convective switch
!  - this is used in cldwat.F90 to turn off stratiform snow production
!
!------------------------------------------------------------------------
integer :: ifconvection_activated(pcols), i_cf_all(pcols)
real(r8) :: cape_meann(pcols)
real(r8) :: rrain_if(pcols),conv_energy_all(pcols)
real(r8) :: cf_tar_all(pcols)
integer :: did_convect
!------------------------------------------------------------------------
!
!  limcnv: removed (along with msg)
!
!------------------------------------------------------------------------
contains
!------------------------------------------------------------------------
!************************************************************************
!
!   subroutine zmconv_readnl(nlfile)
!
!   - This was new to cam5
!   - This appears to define three tuning parameters 
!      c0_lnd = zmconv_c0_lnd
!      c0_ocn = zmconv_c0_ocn
!      ke = zmconv_ke
!   - I will assume I don't need and that no harm done in leaving out (would
!     not compile for me anyway, maybe since previously compiles in zm_conv?)
!
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
! removed: subroutine zm_convi
!
!  - This set resolution dependent ZM parameters I don't need.
!  - but may want something similar at some point.
!  - was called from subroutine zm_conv_init(hypi) in module convect_deep
!
!  - also set limcnv = limcnv_in : top interface level limit for convection
!    probably should do this somewhere else, since otherwise limcnv would
!    be undefined.
!
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  Subroutine if_convr:
!
!  Purpose: main driver for IF convection scheme 
! 
!   original call zm_convr from if_conv_intr.F90:
!   In cam5 landfrac is added as the last argument; Should I add? Its used to
!     calculate c0mask ....
!
! call zm_convr( lchnk   ,ncol    , &
!                state%t       ,state%q     ,prec    ,jctop   ,jcbot   , &
!                state%zm      ,state%phis    ,state%zi      ,ptend_loc%q(:,:,1)    , &
!                ptend_loc%s    ,state%pmid     ,state%pint    ,state%pdel     , &
!                .5_r8*ztodt    ,mcon    ,cme     , cape,      &
!                dlf     ,pflx    ,zdu     ,rprd    , &
!                lengath(lchnk) ,ql      ,rliq    )
!
!  Note that all the gathered quantities have the lchnk identifier
!  Gathered quantities are needed outside this subroutine
!  mainly because momtran and convtran use gathered mu/du, etc.
!  But since I don't call these, may not have to define these.
!
!  Here is call to convtran: the list of needed gathered quantities exactly
!    mirrors call zm_convr above, starting with mu and going to lengath, and the same
!    thing with momtran.
!
!   call convtran (lchnk,                                        &
!                  ptend%lq,state%q, pcnst,  mu(:,:,lchnk), md(:,:,lchnk),   &
!                  du(:,:,lchnk), eu(:,:,lchnk), ed(:,:,lchnk), dp(:,:,lchnk), dsubcld(:,lchnk),  &
!                  jt(:,lchnk),maxg(:,lchnk),ideep(:,lchnk), 1, lengath(lchnk),  &
!                  nstep,   fracis,  ptend%q, dpdry)
!
!  *************** HOW TO CALL if_convr:  *************************
!
!   Note that accessing pver both through use statement in module, and making
!     pver <-> nd association in argument of if_convr. Hope is OK.
!
! call if_convr( lchnk        ,ncol          , &
!                state%t      ,state%q       ,l_windt            ,winds    , &
!                2            ,state%pmid    ,state%pint         ,state%zm , &
!                state%phis   ,state%zi      ,ptend_loc%q(:,:,1) ,ptend%q  , &
!                ptend_loc%s  ,wind_tends    ,jctop              ,jcbot    , &
!                prec         ,.5_r8*ztodt   ,ptend%lq,          ,pcnst    , &
!                pver)  
!               
!    cam5 call has added landfrac as the last argument               
!
!------------------------------------------------------------------------
subroutine if_convr(lchnk       ,ncol        ,oomega        ,              &
                    tt          ,qh          ,domomtran    ,wind         , &
                    mcnst       ,pap         ,paph         ,zm           , &
                    geos        ,zi          ,qtnd         ,dqdt         , &
                    heat        ,dwdt        ,jctop        ,jcbot        , &
                    prec_if     ,delt        ,doconvtran   ,ncnst        , &
                    snow_if     ,rain_if     ,preciporg    ,massprev     , &
                    cm5prev     ,landfrac    )      
!----------------------------------------------------------------------- 
!
!----------------------------------------------------------------------- 
   use cam_history,  only: outfld
!----------------------------------------------------------------------- 
!
!
!
!----------------------------------------------------------------------- 
  use phys_grid, only: get_rlat_all_p
!----------------------------------------------------------------------- 
! 
!  Now using ncnst, so don't need pcnst
! 
!-----------------------------------------------------------------------
!   use constituents, only: pcnst
!-----------------------------------------------------------------------
! ************************ index of variables **********************
!
!  w  * cape     convective available potential energy.
!  i  * dpp      
!  ic  * delt     length of model time-step in seconds.
!  ic  * pver     number of model levels.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  w  * q        grid slice of mixing ratio.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  i/o * t       
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  O  * jcbot    row of base of cloud indices passed out.
!
!  Look at how momtran, or elsewhere determines u/v
!
!  i/o * u        grid slice of u-wind (real).
!  i/o * v        grid slice of v-wind (real).
!
! Need?
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!
!-----------------------------------------------------------------------
!
! multi-level i/o fields:
!
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!   INPUT FIELDS:
!   - do I still need nd as an input variable? association between nd and
!     pver should be being made in if_conv_params ...
!
!-----------------------------------------------------------------------
  integer, intent(in) :: lchnk         ! chunk identifier
  integer, intent(in) :: ncol          ! number of atmospheric columns
!-----------------------------------------------------------------------
!
!   - do I still need nd as an input variable? association between nd and
!     pver should be being made in if_conv_params ...
!   - shouldn't be any references to nd in this module ...
!
!-----------------------------------------------------------------------
!  integer, intent(in) :: nd            ! associated with pver
!-----------------------------------------------------------------------
!
!  FULL level t of model  
!  used to calculate km/hm
!  NEED (changed to tt)
!  state%t <-> tt
!
!-----------------------------------------------------------------------
  real(r8), intent(in) :: oomega(pcols,pver)     
  real(r8), intent(in) :: tt(pcols,pver)        ! grid slice of temperature at mid-layer.
!-----------------------------------------------------------------------
!
!   CONSTITUENT transport:
!
!  FULL level variable
!  calculate rv from qh(i,k,1)
!  calculate rl from qh(i,k,2)
!  - use subroutine initindx to define rl as second constituent.
!  
!  pcnst <-> ncnst
!  state%q <-> qh
!  ptend%lq <-> doconvtran (ncnst or pcnst as array size ???)
!
!  In zm_convir.F90, pcnst was accessed through the following:
!    use constituents, only: pcnst
!    Then dimension of qh was then pcnst
!  In subroutine convtran, pcnst was identified with ncnst in the argument of
!    convtran, and nnst used to define the size of q.
!  Maybe it doesn't matter, but using the convtran approach here (pver <-> nd also)
!
!-----------------------------------------------------------------------
  integer, intent(in) :: ncnst                   ! number of tracers to transport
  real(r8), intent(in) :: qh(pcols,pver,ncnst)   ! grid slice of constituents
  logical, intent(in) :: doconvtran(ncnst)       ! flag for doing convective transport
!-----------------------------------------------------------------------
!
!  MOMENTUM TRANSPORT (in):
!
!  Need three inputs from call momtran:
!  (i) l_windt <-> domomtran
!  (ii) winds <-> wind(pcols,pver,mcnst)
!  (iii) 2 <-> mcnst
!
!-----------------------------------------------------------------------
  logical, intent(in) :: domomtran(mcnst)         ! flag for doing convective transport
  real(r8), intent(in) :: wind(pcols,pver,mcnst)  ! Wind array
  integer, intent(in) :: mcnst                    ! number of mom tracers to transport
!-----------------------------------------------------------------------
!
! FULL pressure levels (hPa)
! 
! state%pmid <-> pap
!
!-----------------------------------------------------------------------
  real(r8), intent(in) :: pap(pcols,pver)     
!-----------------------------------------------------------------------
!
! INTERFACE pressure levels (Pa presumably)
! 
! state%pint  <-> paph
!
!-----------------------------------------------------------------------
  real(r8), intent(in) :: paph(pcols,pver+1)
!-----------------------------------------------------------------------
!
!  zm : FULL height above the local surface
!  z = zm + zs
!  in buoyan hmn = cp*t + grav*z + rl*q, so this is also the height I
!    should use to define my hm.
!  
!  state%zm  <-> zm
!
!-----------------------------------------------------------------------
  real(r8), intent(in) :: zm(pcols,pver)
!-----------------------------------------------------------------------
!
!  surface geopotential.
!  - used to define zs
!
!  state%phis <-> geos
!
!-----------------------------------------------------------------------
  real(r8), intent(in) :: geos(pcols)
!-----------------------------------------------------------------------
!
!  INTERFACE heights
!  presumably zi(i,pver+1) = 0.
!  since zf = zi + zs
!  - used to define zh
!
!  state%zi
!
!-----------------------------------------------------------------------
  real(r8), intent(in) :: zi(pcols,pver+1)
!-----------------------------------------------------------------------
!
! OUTPUT arguments:
!
!  Q tendency
!  in call zm_convr: qtnd(pcols,pver) <-> ptend_loc%q(:,:,1)
!  calculate from my rvtend
!
!  ptend_loc%q(:,:,1) <-> qtnd
!
!-----------------------------------------------------------------------
  real(r8), intent(out) :: qtnd(pcols,pver)   ! specific humidity tendency (kg/kg/s)
!-----------------------------------------------------------------------
!
!  Constituent tendency:
!  Note : in original program dqdt was gathered version of qtnd
!  - this terminology taken from convtran
!  - Need for CLDLIQ
!
!  ptend%q <-> dqdt
!
!-----------------------------------------------------------------------
  real(r8), intent(out) :: dqdt(pcols,pver,ncnst)
!-----------------------------------------------------------------------
!
!  Heating rate
!
!  ptend_loc%s <-> heat
!
!-----------------------------------------------------------------------
  real(r8), intent(out) :: heat(pcols,pver)   ! heating rate (dry static energy tendency, W/kg)
!-----------------------------------------------------------------------
!
!  MOMENTUM TRANSPORT (out):
!
!  one output: 
!  (i) wind_tends <-> dwdt(pcols,pver,ncnst)
!
!-----------------------------------------------------------------------
!  real(r8), intent(out) :: dwdt(pcols,pver,ncnst) ! Wind tendency
  real(r8), intent(out) :: dwdt(pcols,pver,mcnst) ! Wind tendency: changed to mcnst (Nov 2010)
!-----------------------------------------------------------------------
!
!  Should remove:
!
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  zm defines a parameter limcnv, jctop must be less than limcnv+1
!   integer, intent(in) :: limcnv                 ! convection limiting level
!   lel here is the index of highest theoretical convective plume, I guess the undilute
!      detrainment level.
!      jt(i) = max(lel(i),limcnv+1)
!      jt(i) = min(jt(i),pver)
!  scattered from gathered:
!  jctop(ideep(i)) = jt(i)
!  jctop allocated:
!  jctop => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)
!  needed on output? yes need to define gathered version jt for defining pcont:
!  pcont(ideep(i,lchnk)) = state%pmid(ideep(i,lchnk),jt(i,lchnk))  ! gathered array (or jctop ungathered)
!  Could easilly calculate.
!  Used as follows in if_convect_deep.F90:
!    ifld = pbuf_get_fld_idx('CLDTOP')
!    jctop => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)
!   ifld = pbuf_get_fld_idx('CLDBOT')
!   jcbot => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)
!  DEFINE
!  jctop
!
!-----------------------------------------------------------------------
  real(r8), intent(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
!-----------------------------------------------------------------------
!
!  should remove
!
! SCATTERED version of maxg index (gathered maximum MSE I think)
! jcbot => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)
! jcbot(ideep(i)) = maxg(i)
!  DEFINE? This would not be relevant to me; have convection from multiple
!    levels.
! searched if_conv_intr.F90: no evidence where used
! in if_convect_deep.F90:
!    ifld = pbuf_get_fld_idx('CLDBOT')
!    jcbot => pbuf(ifld)%fld_ptr(1,1:pcols,1,state%lchnk,1)
! Looks like jcbot identified with CLDBOT, but not clear why care if defined.
! NEED?
! jcbot
!
!-----------------------------------------------------------------------
  real(r8), intent(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
!-----------------------------------------------------------------------
!
!  - prec is the total precipitation
!  - in the zm scheme, snow is produced by subroutine zm_conv_evap, so
!    precip would have referred to rain only, at this stage.
!  - however, in the zm scheme, prec is modified in "zm_conv_evap" where
!    it is an inout variable, where it seems to end up referring to the
!    total precipitation (rain + snow)
!  - from "subroutine zm_conv_tend" snow is passed up in call "if_conv_tend"
!    as snow.
!  - from "subroutine convect_deep_tend" it is passed up
!    from snow to "prec_zmc" in "call convect_deep_tend".
!  - "prec_zmc" is then used as a water flux condition before the energy check,
!    with the command:
!    flx_cnd(:ncol) = prec_zmc(:ncol) + rliq(:ncol)
!  - on entry into "subroutine check_energy_chng", flx_cnd is described:
!    flx_cnd(pcols)      ! boundary flux of liquid+ice    (m/s)
!  - i.e. as rain+snow (+rliq in the past)
!  - since I will be passing "prec" up to "prec_zmc", the simplest thing to do
!    is is to define prec as the total.
!  - the current check on energy conservation ignores issue of temperature of
!    the exiting precipitation.
!
! prec <-> prec
!
!  Only prec_if and snow_if are needed higher up.
!  Rest are being passed up higher for diagnostic purposes
!
!-----------------------------------------------------------------------
  real(r8), intent(out) :: prec_if(pcols)
  real(r8), intent(out) :: rain_if(pcols)
  real(r8), intent(out) :: snow_if(pcols)
!-----------------------------------------------------------------------
!
!  nonlinear prev variables
!
!-----------------------------------------------------------------------
  real(r8), intent(inout) :: preciporg(pcols) 
  real(r8), intent(inout) :: massprev(pcols) 
  real(r8), intent(inout) :: cm5prev(pcols) 
  real(r8), intent(in) :: landfrac(pcols) 

!-----------------------------------------------------------------------
!
!   END of IN/OUT variables
!
!
! Local variables
!
!-----------------------------------------------------------------------
  real(r8) zs(pcols)

!-----------------------------------------------------------------------
!
!  1. 3D Rain diagnostics 
!
!-----------------------------------------------------------------------
  real(r8) uprain_1_diag(pcols)
  real(r8) uprain_2_diag(pcols)
  real(r8) uprain_surf_rlp_diag(pcols)
  real(r8) uprain_surf_rv_diag(pcols)
  real(r8) uprain_start_rlp_diag(pcols)
  real(r8) uprain_rlp_rvv_diag(pcols)
  real(r8) uprain_evap_diag(pcols)
  real(r8) anrain_down_diag(pcols)
  real(r8) anrain_evap_diag(pcols)
  real(r8) ansnow_melt_diag(pcols)
  real(r8) ansnow_subl_diag(pcols)
  real(r8) uprain_down_diag(pcols)

  real(r8) ansnow_conv_diag(pcols)
  real(r8) ansnow_strat_diag(pcols)
  real(r8) ansnow_strat_rv_diag(pcols)
  real(r8) ansnow_strat_ri_diag(pcols)
  real(r8) ansnow_surf_diag(pcols)
  real(r8) anrain_surf_diag(pcols)

  real(r8) rel_err_diag(pcols)
  real(r8) abs_err_diag(pcols)

  real(r8) mass_inc_diag(pcols)
  real(r8) mass_rat_diag(pcols)

!-----------------------------------------------------------------------
!
!  2. 3D Column Water
!
!-----------------------------------------------------------------------
  real(r8) colwat_diag(pcols)
  real(r8) colrh_diag(pcols)

!---------------------------------------------------------------------
!
!  3. 3D Cape
!
!---------------------------------------------------------------------
  real(r8) conv_energy_diag(pcols)
  real(r8) cf_tar_diag(pcols)
  real(r8) lwp_tar_diag(pcols)
  real(r8) cape_mean_diag(pcols)

  real(r8) LTS_diag(pcols)
  real(r8) EIS_diag(pcols)
  real(r8) mse_rat_diag(pcols)

  real(r8) cape_1_diag(pcols)
  real(r8) cape_2_diag(pcols)
  real(r8) cape_3_diag(pcols)
  real(r8) cape_4_diag(pcols)

  real(r8) cfraction_1_diag(pcols)
  real(r8) cfraction_2_diag(pcols)
  real(r8) cfraction_3_diag(pcols)
  real(r8) cfraction_4_diag(pcols)

  real(r8) cape_1_1_diag(pcols)
  real(r8) cape_2_1_diag(pcols)
  real(r8) cape_3_1_diag(pcols)

  real(r8) amp_diag(pcols)
  real(r8) amp_om5_diag(pcols)
  real(r8) amp_rain_diag(pcols)

  real(r8) dp_dn_diag(pcols)

  real(r8) f_1_diag(pcols)
  real(r8) f_2_diag(pcols)
  real(r8) f_3_diag(pcols)
  real(r8) f_4_diag(pcols)

  real(r8) hmuprain_start_diag(pcols)
  real(r8) hmuprain_surf_diag(pcols)
  real(r8) hmansnow_start_diag(pcols)
  real(r8) hmansnow_surf_diag(pcols)
  real(r8) hmanrain_surf_diag(pcols)

!-----------------------------------------------------------------------
!
!  5. 3D Start variables
!
!-----------------------------------------------------------------------
  real(r8) mass_start_1_diag(pcols)
  real(r8) mass_start_2_diag(pcols)

!-----------------------------------------------------------------------
!
! 7. 3D nonlinearities; org variables
!
!-----------------------------------------------------------------------
  real(r8) om4_diag(pcols)
  real(r8) om5_diag(pcols)
  real(r8) om6_diag(pcols)

  real(r8) km_width_diag(pcols)

! Rain variables
  real(r8) md_rat_diag(pcols) 
  real(r8) f_an_diag(pcols) 
  real(r8) f_ice_diag(pcols) 

  real(r8) av_z_ri_diag(pcols)
  real(r8) av_z_ri_rh_diag(pcols)
  real(r8) av_z_ri_dt_diag(pcols)
  real(r8) zpeak_diag(pcols)
  real(r8) precip_org_diag(pcols)

!-----------------------------------------------------------------------
!
!  8. 3D updraft eff
!
!-----------------------------------------------------------------------
  real(r8) updraft_eff_eff(pcols)
  real(r8) updraft_eff_det(pcols)
  real(r8) updraft_eff_rlp(pcols)

!-----------------------------------------------------------------------
!
!  9. 4DF tend
!
!-----------------------------------------------------------------------
  real(r8) tptend_real(pcols,pver)
  real(r8) rhtend_real(pcols,pver)
  real(r8) qtend_real(pcols,pver)
  real(r8) kmtend_save(pcols,pver)
  real(r8) rvtend_save(pcols,pver) 
  real(r8) rltend_save(pcols,pver) 
  real(r8) rlnew_save(pcols,pver) 
  real(r8) ritend_save(pcols,pver) 
  real(r8) utend_save(pcols,pver) 
  real(r8) vtend_save(pcols,pver) 

!-----------------------------------------------------------------------
!
!  10. 4DI
!
!-----------------------------------------------------------------------
  real(r8) fmass_diag(pcols,pver+1)
  real(r8) fmass_dn_diag(pcols,pver+1)
  real(r8) phalf_diag(pcols,pver+1)

!-----------------------------------------------------------------------
!
!  11. 4DF 
!
!-----------------------------------------------------------------------
  real(r8) ent_diag(pcols,pver)
  real(r8) det_diag(pcols,pver)

  real(r8) updet_1_diag(pcols,pver)
  real(r8) updet_2_diag(pcols,pver)

  real(r8) rh_dn_diag(pcols,pver)
  real(r8) drv_dndet_diag(pcols,pver)
  real(r8) dt_dndet_diag(pcols,pver)
  real(r8) dndet_diag(pcols,pver)
  real(r8) dnent_diag(pcols,pver)
  real(r8) mtdet_diag(pcols,pver)
  real(r8) mtent_diag(pcols,pver)
  real(r8) drv_up_diag(pcols,pver)
  real(r8) drv_an_diag(pcols,pver)

!-----------------------------------------------------------------------
!
!  12. 4DF 
!
!-----------------------------------------------------------------------
  real(r8) rhif_store(pcols,pver) 
  real(r8) tif_store(pcols,pver) 
  real(r8) hm_store(pcols,pver) 
  real(r8) prob_random_store(pcols,pver) 
  real(r8) ri_half_store(pcols,pver) 
  real(r8) rl_half_store(pcols,pver) 

!-----------------------------------------------------------------------
!
!  13. 4DF Variables scaling with updraft mass flux (i.e. updraft properties)
!
!-----------------------------------------------------------------------
  real(r8) mp_save(pcols,pver) 
  real(r8) bp_cape_save(pcols,pver) 
  real(r8) bp1_str_save(pcols,pver) 
  real(r8) bp1_prp_save(pcols,pver) 
  real(r8) bp1_itr_save(pcols,pver) 
  real(r8) bp2_str_save(pcols,pver) 
  real(r8) bp2_prp_save(pcols,pver) 
  real(r8) bp2_itr_save(pcols,pver) 
  real(r8) f_ent_1_save(pcols,pver) 
  real(r8) f_det_1_save(pcols,pver) 
  real(r8) f_ent_2_save(pcols,pver) 
  real(r8) f_det_2_save(pcols,pver) 
  real(r8) rlp_save(pcols,pver) 
  real(r8) rip_save(pcols,pver) 
  real(r8) hmp_mean_save(pcols,pver)
  real(r8) hm_diff_save(pcols,pver)

!-----------------------------------------------------------------------
!
!  14. 4DF Variables scaling with updraft detrainment mass
!
!-----------------------------------------------------------------------
  real(r8) tdiff_detrain_up_save(pcols,pver) 
  real(r8) cape_detrain_up_save(pcols,pver) 
  real(r8) mass_detrain_up_save(pcols,pver) 
  real(r8) hm_diff_det_save(pcols,pver) 

!-----------------------------------------------------------------------
!
!  15. 4DF Evaporation variables
!
!-----------------------------------------------------------------------
  real(r8) uprainevap_save(pcols,pver)
  real(r8) anrainevap_save(pcols,pver)
  real(r8) ansnowsubl_save(pcols,pver)
  real(r8) ansnowmelt_save(pcols,pver)
  real(r8) upraindown_save(pcols,pver)
  real(r8) anraindown_save(pcols,pver)

!-----------------------------------------------------------------------
!
!  16. 4DF Precipitation Production
!
!-----------------------------------------------------------------------
  real(r8) uprain_rlp_save(pcols,pver)
  real(r8) ansnow_rlp_save(pcols,pver)
  real(r8) uprain_rv_save(pcols,pver)
  real(r8) ansnow_rv_save(pcols,pver)
  real(r8) rlp_rvv_save(pcols,pver)

!-----------------------------------------------------------------------
!
!  17. 4DF Cloud Production/Loss
!
!-----------------------------------------------------------------------
  real(r8) uprain_rl_tend_save(pcols,pver)
  real(r8) rl_evap_tend_save(pcols,pver)
  real(r8) ri_evap_tend_save(pcols,pver)
  real(r8) ri_prod_tend_save(pcols,pver)
  real(r8) rl_det_tend_save(pcols,pver)  
  real(r8) ri_det_tend_save(pcols,pver)  
  real(r8) rl_vert_tend_save(pcols,pver)  
  real(r8) ri_vert_tend_save(pcols,pver)  

  real(r8) rl_dt_save(pcols,pver)  
  real(r8) rl_rh_save(pcols,pver)  
  real(r8) ri_dt_save(pcols,pver)  
  real(r8) ri_rh_save(pcols,pver)  
  real(r8) ri_rs_save(pcols,pver)  

  real(r8) col_ri_dt_save(pcols)  
  real(r8) col_ri_rh_save(pcols)  
  real(r8) col_rl_dt_save(pcols)  
  real(r8) col_rl_rh_save(pcols)  

!-----------------------------------------------------------------------
!
!  Not saved to history
!
!-----------------------------------------------------------------------
  real(r8) hm_uprain(pcols)
  real(r8) hm_ansnow(pcols)

  real(r8) rhnew(pcols,pver) 
!-----------------------------------------------------------------------
!
!  To get latitude
!
!-----------------------------------------------------------------------
  real(r8) :: rlat(pcols),this_lat

  real(r8) rvnew(pver)
  real(r8) rlnew(pver)
  real(r8) rinew(pver)
  real(r8) tnew(pver)
  real(r8) kmnew(pver)

!-----------------------------------------------------------------------
!
! general work fields (local variables):
!
!  changed p -> pp to avoid conflict with 1D p (maybe simpler to just avoid using)
!
!  Why isn't delt written as in IN variable? Does it matter?
!
!-----------------------------------------------------------------------
  real(r8) pp(pcols,pver)             ! w  grid slice of ambient mid-layer pressure in mbs.
  real(r8) const(pver)    
  real(r8) delt                       ! length of model time-step in seconds.
  real(r8) qvnew(pver),qv(pver)
  real(r8) qlnew(pver),ql(pver)
  real(r8) qinew(pver),qi(pver)
  real(r8) col_rl,col_rv,col_ri
  real(r8) col_rl_old,col_rv_old,col_ri_old
  real(r8) err_wat_corr
  real(r8) err_wat_corr_max
  real(r8) col_km_old,col_km
  real(r8) err_km_corr,err_km,dp_km
  real(r8) km_end,km_init
  real(r8) lapse_1_4 
!----------------------------------------------------------------------
!
!  Calculating omega
!
!----------------------------------------------------------------------
  real(r8) ppp,ooo,ppp_top,ppp_strat,dp_tot,dpp,ttt,prr,prob
!----------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------
  real(r8) rtt,km_small,dkm,ee,sec_day
  real(r8) rnnew_small,rlnew_small,rvnew_small,dmdry,dmdry_old,dmdry_new
  real(r8) dp_precip,totprecip,err_wat,dp_wat,small,dp_precip_dry
!----------------------------------------------------------------------
!
!   save variables for debugging
!
!----------------------------------------------------------------------
  real(r8) tdiff,rvdiff,rldiff,ridiff,kmdiff,hmdiff,udiff,vdiff,pdiff 
  real(r8) t_save(pver)
  real(r8) rv_save(pver)
  real(r8) rl_save(pver)
  real(r8) ri_save(pver)
  real(r8) km_save(pver)
  real(r8) hm_save(pver)
  real(r8) uwind_save(pver)
  real(r8) vwind_save(pver)
  real(r8) p_save(pver)

  real(r8) t_surf_cold,temp_hot,temp,sss,rss,rain_in_mmday
  real(r8) col_wat,col_wat_old
  real(r8) esat_cam4,rsat_cam4,rh_cam4
!----------------------------------------------------------------------
!
!  Integer internal variables
!  - shouldn't have to declare nstep here since in if_conv_params.f90
!
!----------------------------------------------------------------------
  integer i,j,m,ixcldliq,ixcldice,temp_ok
  integer k, kk
  integer nprint_wat,nprint_en,nprint_surf,nprint
!----------------------------------------------------------------------
!
!  FFF
!
!----------------------------------------------------------------------
   real(r8)  clat(pcols),area(pcols),pi,radians
   real(r8)  dx_angle,dy_angle,rearth,dx_rad,dy_rad,dx,dy,phi_rad
   real(r8)  tdpp, rho, dzz, AA, BB, CC, DD, rel_err, abs_err, BB_t, CC_t, abs_err_t
   real(r8)  DD_t, xxx, ri_half, zzz, tstep
   real(r8)  tcc, denom, esat_wat, rsat_wat
!----------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------
integer :: count_conv_start, count_conv_end
integer :: count_col_start, count_col_end

!  integer :: nstep

!----------------------------------------------------------------------
!
!  ********  COUNTER *****************
!
!----------------------------------------------------------------------
! call system_clock(count_conv_start, count_rate, count_max)
! print*,'--------- Called system_clock --------------'
! print*,'count_conv_start = ',count_conv_start

!----------------------------------------------------------------------
!
!  Variables needed by module if_conv_tend.f90
!  They are communicated to if_conv_tend.f90 through module if_conv_params.f90
!
!  nstep(SET ZERO)
!  tstep:
! 
!  ************ How to define tstep? *******************
! 
!  In normal cam4: 
!                prec         ,.5_r8*ztodt   ,ptend%lq,          ,pcnst    , &
!  delt  is passed in as 0.5*ztodt and should equal 30 min I think
! 
! 
!  tstep should probably be equal to ztodt. In zm_conv.F90 delt seems to be
!  always multiplied by 2.  
!  - keep in mind: as bits of mass/moisture/tracer entrained /detrained at model
!    levels become comparable with mass in layers, convective tendencies become
!    non-linear in tstep, so the value you assign to tstep does matter
!  - however likely better to keep tstep equal to the real model timestep not
!    ztodt. 
!
!   CHANGED MY MIND: easy to generate negative values if use smaller
!    tstep. Too risky.
!  
!   Sept 27/2010: back to using tstep = delt.
!   - using a longer timestep makes it much easier to generate negative values
!     inside the convection, and to prevent these, would be forced to put a 
!     limit on rainfall per timestep of about 70 mm/day, pretty low, especially
!     by tropical standards.    
!   - This will mean that negative tracer values will be generated at the n+1
!     timestep, but hopefully my checks will ensure that not generated when find
!     tracer values at the n=0 time step ...
!   - Feb 2013: It doesn't seem to make much difference to the simulations whether
!     you use tstep = delt or tstep = 2*delt
!   tstep is given a value here and then passed to if_conv_tend.f90 by being declared
!     as a variable in if_conv_params.f90 (Nov 13 2010: seems to be 900 s in if_conv_tend.f90
!     for some reason).
!  
!----------------------------------------------------------------------
temp_hot = 400.
t_surf_cold = 180.  ! gives a warning
tstep = 2._r8*delt
!tstep = delt   
!print*,'in if_conv.F90 tstep delt = ',tstep,delt
nstep = get_nstep()

!----------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------
   rgrav  = 1.0_r8/gravit
!----------------------------------------------------------------------
!
! Initialize necessary arrays.
! zero out variables not used in cam
! qtnd: q tendency
!
!----------------------------------------------------------------------
  qtnd(:,:) = 0._r8
  heat(:,:) = 0._r8
  prec_if(:ncol) = 0._r8
  rain_if(:ncol) = 0._r8
  snow_if(:ncol) = 0._r8

!----------------------------------------------------------------------
!
! Initialize tendencies
! - removed since not neccessary; since now loop over all columns.
!
!----------------------------------------------------------------------
  if (nstep == 1) then
!  do i = 1,100
!   temp = 200. + float(i)*1.
!   esat_cam4 = estblf(temp)
!   print*, 'esat = ',temp,esat_cam4,esat(temp),(esat_cam4-esat(temp))/esat(temp)
!  end do 
  endif

!----------------------------------------------------------------------
!
!  Check
!
!----------------------------------------------------------------------
! print*,'in if_conv.F90 nstep = ',nstep
!  if (pcols /= ncol) then
!    print*,'pcols = ',pcols
!    print*,'ncol = ',ncol
!    stop 'pcols not equal to ncol in if_conv.F90'
!  endif

  call get_rlat_all_p(lchnk, ncol, clat)
!----------------------------------------------------------------------
!
!  *********************************************************************
!  ***************    Start Loop over columns   *************************
!  *********************************************************************
!
!----------------------------------------------------------------------
  do i = 1,ncol

!----------------------------------------------------------------------
!
!  Get latitudes
!  - clat gives angles in radians
!  - then assign areas
!  - area is defined here as the fractional area of the earth's surface
!
!----------------------------------------------------------------------
  pi = 3.1415697
  radians = 360./(2.0*pi)
  dx_angle = 2.50
  dy_angle = 1.90
  rearth = 6378100.
  this_lat = clat(i)*57.296_r8  
  phi_rad = clat(i)
  dx_rad = (dx_angle*pi)/180.
  dy_rad = (dy_angle*pi)/180.
  dx = rearth*cos(phi_rad)*dx_rad
  dy = rearth*dy_rad
  area(i) = dx*dy/(4.*pi*rearth*rearth)
!   print*,'------------------'
!   print*,'i clat(i) = ',i,clat(i)
!   print*,'dx = ',dx
!   print*,'dy = ',dy
!   print*,'dx_rad = ',dx_rad
!   print*,'dy_rad = ',dy_rad
!   print*,'area(i) = ',area(i)

!----------------------------------------------------------------------
!
!  ********  COUNTER *****************
!
!----------------------------------------------------------------------
! call system_clock(count_col_start, count_rate, count_max)
! print*,'--------- Called system_clock --------------'
! print*,'count_col_start = ',count_col_start

!   print*,'in subroutine if_convr starting column i nstep = ',i,nstep
!----------------------------------------------------------------------
!
! On first timestep, initialize prev variables
!
!----------------------------------------------------------------------
  if (nstep == 0) then
    preciporg(i) = 0.
    massprev(i) = 0.
    cm5prev(i) = 0.
!----------------------------------------------------------------------
!
!  Define nonlinear prev variables
!  - These are the INPUTS into if_conv_tend.f90
!
!----------------------------------------------------------------------
  else
    precip_org = preciporg(i)
    cm5_prev = cm5prev(i)
    mass_prev = massprev(i)
  endif
!----------------------------------------------------------------------
!
!
!----------------------------------------------------------------------
  om4_diag(i) = 24.*3600.*0.01*oomega(i,pver-3)
  om5_diag(i) = 24.*3600.*0.01*oomega(i,pver-4)
  om6_diag(i) = 24.*3600.*0.01*oomega(i,pver-5)
!  om5_diag(i) = 0.5*24.*3600.*0.01*(oomega(i,pver-4)+oomega(i,pver-5))
!----------------------------------------------------------------------
!
!  Define surface height and pressure.
!
!  zs : local height of surface
!  paph : interface levels, pretty sure is in Pa
!  zf  = zi + zs is the interface height, but I don't use.
!
!  zh,phalf are "my" 1D arrays, so first index is surface.
!
!----------------------------------------------------------------------
  zs(i) = geos(i)*rgrav
  zh(1) = zs(1)
  phalf(1) = paph(i,pver+1)
!  print*,'surface height i zh(1) zs(i) = ',i,zh(1),zs(i)
!  print*,'phalf(1) = ',phalf(1)
!----------------------------------------------------------------------
!
!  Determine cloud liquid/ice indices
!
!----------------------------------------------------------------------
  call cnst_get_ind('CLDLIQ', ixcldliq)   ! Is this from me or default?
  call cnst_get_ind('CLDICE', ixcldice)   ! Is this from me or default?
!----------------------------------------------------------------------
!
!  Define basic arrays needed by conv_tend.f90
!
!  phalf (HALF pressure): get from paph
!  p (FULL pressure): define directly from pap(i,k)
!  z (FULL height): get from zm/zs   z(i,k) = zm(i,k) + zs(i)
!    zm : FULL height above the local surface
!    z = zm + zs would then be the absolute full height; consistent with
!    def in buoyan hmn = cp*t + grav*z + rl*q, so this is also the height I
!    should use to define my hm.
!  zh (HALF height): define from zf
!  rv: get from qh(i,k,1)
!  t: get from tt(i,k) is the input full level temperature (changed from t)
!  spec:
!
!----------------------------------------------------------------------
  do k = 1,pver
    phalf(pver+2-k) = paph(i,k)
    zh(pver+2-k) = zi(i,k) + zs(i)
    p(pver+1-k) = pap(i,k)
    om(pver+1-k) = oomega(i,k)
    z(pver+1-k) = zm(i,k) + zs(i)
    qv(pver+1-k) = qh(i,k,1)
    ql(pver+1-k) = qh(i,k,ixcldliq)
    qi(pver+1-k) = qh(i,k,ixcldice)
    t(pver+1-k) = tt(i,k)
!   print*,'k pver+2-k pver+1-k = ',k,pver+2-k,pver+1-k
!   print*,'k phalf(pver+2-k) = ',k,phalf(pver+2-k)
!   print*,'p(pver+1-k) = ',p(pver+1-k)
!   print*,'k pver+2-k  zh(pver+2-k) = ',k,pver+2-k,zh(pver+2-k)
!   print*,'rv(pver+1-k) = ',rv(pver+1-k)
!   print*,'rl(pver+1-k) = ',rl(pver+1-k)
   temp_ok = 0
   if ((t(pver+1-k) > 30.).and.(t(pver+1-k) < temp_hot)) temp_ok = 1
   if (temp_ok == 0) then
     print*,'------------ UNREALISTIC TEMPERATURE ON ENTRY TO if_conv.F90 ------'
     print*,'initialization nstep i = ',nstep,i
     print*,'initialization t(pver+1-k) = ',t(pver+1-k)
     print*,'initialization k z(pver+1-k) = ',k,z(pver+1-k)
     print*,'landfrac(i) = ',landfrac(i)
     print*,'i clat(i) = ',i,clat(i)
     print*,'latitude in degrees = ',clat(i)*57.296_r8
     stop 'temperature error on entry to if_conv.F90 '
   endif
  end do
! print*,'i zi(i,pver+1) = ',i,zi(i,pver+1)
!----------------------------------------------------------------------
!
!  Check for extremely cold surface temperatures
!
!----------------------------------------------------------------------
  if (t(1) < t_surf_cold) then
    print*,'----- WARNING: extremely cold surface temperatures  ---'
    print*,'---- This is on entry to if_conv.F90 ------------'
    print*,'--------- t(1) less than t_surf_cold ------'
    print*,'t(1) = ',t(1)
    print*,'t_surf_cold = ',t_surf_cold
    print*,'p(1) = ',p(1)
    print*,'hmp(1) = ',hmp(1)
    print*,'i clat(i) = ',i,clat(i)
    print*,'landfrac(i) = ',landfrac(i)
    print*,'precip_org = ',precip_org
    print*,'----------- leaving cold surface temperature ---'
  endif
!----------------------------------------------------------------------
!
!  Define rv/rl from qv/ql:
!
!  rv = mv/md = mv/(m - mv - ml) = mv/(m - qv*m - ql*m) = mv/[m(1-qv-ql)]
!  rv = qv/(1-qv-ql)
!
!----------------------------------------------------------------------
  do k = 1,pver
    rv(k) = qv(k)/(1. - qv(k) - ql(k) - qi(k))
    rl(k) = ql(k)/(1. - qv(k) - ql(k) - qi(k))
    ri(k) = qi(k)/(1. - qv(k) - ql(k) - qi(k))
  end do
!----------------------------------------------------------------------
!
!  Define zh:
!  July 2015: stick to their definition.
!
!----------------------------------------------------------------------
! zh(1) = zs(i)
! do k = 2,pver+1
!   rtt = rv(k-1) + rl(k-1) + ri(k-1)
!   tdpp = t(k-1)*(1. + epsi*rv(k-1))/(1. + rtt)
!   rho = p(k-1)/(rdd*tdpp)
!   dpp = phalf(k-1) - phalf(k)
!   dzz = dpp/(rho*g)
!   zh(k) = zh(k-1) + dzz
!   z(k-1) = 0.5*(zh(k-1) + zh(k))
!   print*,'----- k = ',k
!   print*,'rtt = ',rtt
!   print*,'rho = ',rho
!   print*,'p(k-1) = ',p(k-1)
!   print*,'dpp = ',dpp
!   print*,'t(k-1) = ',t(k-1)
!   print*,'tdpp = ',tdpp
!   print*,'dzz = ',dzz
!   print*,'zh(k) = ',zh(k)
!   print*,'z(k-1) = ',z(k-1)
! end do
! stop
!----------------------------------------------------------------------
!
!  Define uwind/vwind:
!
!  before "call momtran", winds is defined from the state as:
!   winds(:ncol,:pver,1) = state1%u(:ncol,:pver)
!   winds(:ncol,:pver,2) = state1%v(:ncol,:pver)
!
!  winds then used in call momtran:
!  call momtran      momtran
!     2               ncnst
!   winds             real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Wind array
!
!  Here:
!  call if_convr      if_convr
!     2               mcnst
!  winds              real(r8), intent(in) :: wind(pcols,pver,mcnst)  ! Wind array
!
!  Note ncnst -> mcnst
!
!  Can't associate winds with q, since too easilly confused with
!   constituents.
!
!  uwind/vwind: per unit dry mass, so a bit bigger than actual wind
!  - will have to modify tendencies accordingly
!
!----------------------------------------------------------------------
! print*,'Starting momentum'
  do m = 1, mcnst   
    if (domomtran(m)) then
      do k = 1,pver
        const(pver+1-k) = wind(i,k,m)
        if (m == 1) uwind(pver+1-k) = const(pver+1-k)*(1.+rv(pver+1-k)+rl(pver+1-k))
        if (m == 2) vwind(pver+1-k) = const(pver+1-k)*(1.+rv(pver+1-k)+rl(pver+1-k))
      end do
    endif
  end do
!----------------------------------------------------------------------
!
!  Define secondary quantities:
!
!  ri: set zero
!  rt: define
!  rs
!  rh:
!  rhod:
!  td:  Needs to be defined for cape calculation in if_conv_tend.f90
!  km:
!  hm:
!  bdflow:
!
!  How to define background rh?
!
!  In cam_diagnostics, defined as:
!      if (hist_fld_active('RELHUM')) then
!         call aqsat (state%t    ,state%pmid  ,tem2    ,ftem    ,pcols   , &
!            ncol ,pver  ,1       ,pver    )
!         ftem(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100._r8
!         call outfld ('RELHUM  ',ftem    ,pcols   ,lchnk     )
!
!
!----------------------------------------------------------------------
!print*,'Starting RH enthalpy calculation'
  do k = 1,pver
    rt(k) = rv(k) + rl(k) + ri(k)
!----------------------------------------------------------------------
!  Define rh using my definition
!  - is rh_cam4 needed at all? remove?
!----------------------------------------------------------------------
    esat_cam4 = estblf(t(k))
    if (p(k) > 1000.) then
      rsat_cam4 = eps*esat_cam4/(p(k) - esat_cam4)
    else
      rsat_cam4 = eps*esat_cam4/p(k)
    endif
    rh_cam4 = rv(k)/rsat_cam4
!   rs(k) = rsat_cam4
!   rh(k) = rh_cam4
    rs(k) = rsat(t(k),p(k))
    rh(k) = rv(k)/rs(k) 
    tcc = t(k)-tkelvin
    denom = 243.12 + tcc
    esat_wat = 100.*6.112*exp(17.62*tcc/denom)
    rs_wat(k) = eps*esat_wat/(p(k) - esat_wat)
    rh_wat(k) = rv(k)/rs_wat(k)
!   print*,'---------------- k = ',k
!   print*,'p(k) = ',p(k)
!   print*,'t(k) = ',t(k)
!   print*,'tcc = ',tcc
!   print*,'esat_wat = ',esat_wat
!   print*,'rv(k) = ',rv(k)
!   print*,'rs(k) = ',rs(k)
!   print*,'rs_wat(k) = ',rs_wat(k)
!   print*,'rh(k) = ',rh(k)
!   print*,'rh_wat(k) = ',rh_wat(k)
!   print*,'----------------'
!   if (k == 10) print*,'from if_conv.F90 k rh(k) = ',k,rh(k)
!   write(99,*) k,t(k),esat_cam4,esat(t(k)),rh_cam4,rh(k),rsat_cam4-rs(k)
!   print*,k,t(k),esat_cam4,esat(t(k)),rh_cam4,rh(k),rsat_cam4-rs(k)
!----------------------------------------------------------------------
!  Define density rho and density temperature
!----------------------------------------------------------------------
    ee = rv(k)*p(k)/(eps + rv(k))
    rhod(k) = (p(k) - ee)/(rdd*t(k))
    td(k) = t(k)*(1. + epsi*rv(k))/(1. + rv(k) + rl(k) + ri(k))
!----------------------------------------------------------------------
!  Define enthalpy km from t,rt,rv
!----------------------------------------------------------------------
    call enthalpy(km(k),t(k),rt(k),rv(k),ri(k))
    if ((km(k) > km_bad).or.(t(k) > 380.)) then
      print*,'======= km or t too big on entry if_conv.F90 ===='
      print*,'k = ',k
      print*,'km(k) = ',km(k)
      print*,'km_bad = ',km_bad
      print*,'t(k) = ',t(k)
      print*,'p(k) = ',p(k)
      print*,'z(k) = ',z(k)
      print*,'rt(k) = ',rt(k)
      print*,'rv(k) = ',rv(k)
      print*,'rl(k) = ',rl(k)
      print*,'ri(k) = ',ri(k)
      print*,'qv(k) = ',qv(k)
      print*,'ql(k) = ',ql(k)
      print*,'qi(k) = ',qi(k)
      print*,'rh(k) = ',rh(k)
      print*,'clat(i) = ',clat(i)
      stop 'km too big'
    endif
    hm(k) = hmoist(t(k),rt(k),rv(k),ri(k),z(k))
    bdflow(k) = 0._r8
!    print*,'k rh(k) km(k) hm(k) = ',k,rh(k),km(k),hm(k)
!   if (k == 1) then
!     print*,'----------- INITIAL DEF of km before call to get_tend ----------------'
!     print*,'km(1) = ',km(1)
!     print*,'t(1) = ',t(1)
!     print*,'rv(1) = ',rv(1)
!     print*,'rl(1) = ',rl(1)
!     print*,'ri(1) = ',ri(1)
!     print*,'-------------------------------------------------------------'
!   endif
  end do
! print*,'---------------------------------'
! stop 'got rh profile'
!----------------------------------------------------------------------
!
!  Define other constituents (spec)
!  - conv_tend.f90 OK with numsp = 0?
!  - could define spec array using slots from qh with m= 4,5,.... 
!
!  In convtran, there is an allowance for transporting either "dry" or
!    "moist" tracers via a setting of: cnst_get_type_byind(m).eq.'dry'
!  All my tracers are transported as "dry", but to do both would also
!    have to use this variable.
!
!  Could apply convective transport to cloud ice/water:
!  - How its done using convtran:
!
!  before call convtran need call to cnst_get_in to obtain ixcldliq/ixcldice
!     call cnst_get_ind('CLDLIQ', ixcldliq)
!     call cnst_get_ind('CLDICE', ixcldice)
!  Then set:
!    ptend_loc%lq(ixcldice) = .true.
!    ptend_loc%lq(ixcldliq) = .true.
!  and associate these with doconvtran as input argument to convtran to do
!   convective transports.
!
!  integer :: ixcldice, ixcldliq, ixcldpre    ! constituent indices for cloud liquid and ice water.
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!   Define save variables: basic state variables on input to convection
!   - Need to make sure that any starting variable used to define a tendency
!     is unchanged
!
!----------------------------------------------------------------------
  do k = 1,pver
    t_save(k) = t(k)
    rv_save(k) = rv(k)
    rl_save(k) = rl(k)
    ri_save(k) = ri(k)
    km_save(k) = km(k)
    hm_save(k) = hm(k)
    uwind_save(k) = uwind(k)
    vwind_save(k) = vwind(k)
    p_save(k) = p(k)
  end do

!----------------------------------------------------------------------
!
!  Input profile defined
!  - can call conv_unstable to determine convective switches
!
!  Reasons for non-zero tendencies:
!  (i) non-zero updrafts 
!  (ii) rv <-> rn exchange
!  (iii) precipitation of anvil rl 
!
!  It is important to always call get_tend. It calculates several quantities,
!  such as the nonlinear variables, that are needed every timestep.
!
!----------------------------------------------------------------------
!  print*,'calling get_tend i nstep = ',i,nstep
  call get_tend(did_convect,landfrac(i),om4_diag(i),om5_diag(i),om6_diag(i),this_lat,tstep)

!----------------------------------------------------------------------
!
! To save time: Only do this section if did convection
! - Aug 2013: using this switch does not seem to result in any increase
!   in speed; however is likely a good thing to do anyway.
!
!----------------------------------------------------------------------
if (did_convect == 1) then
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!  Will have on output:
!
! spectend,rvtend,utend,vtend,ndraft 
! dpsurf_up,dpsurf_dn (need?) _ But don't actually use these; compute
! tendencies from "new" values
!
!  All these tendencies are /day.
!
!  Debugging: check for changes in state variables (should not happen)
!
!----------------------------------------------------------------------
  small = 1.0E-10
  do k = 1,pver
    tdiff = abs(t_save(k)-t(k))
    rvdiff = abs(rv_save(k)-rv(k))
    rldiff = abs(rl_save(k)-rl(k))
    ridiff = abs(ri_save(k)-ri(k))
    kmdiff = abs((km_save(k)-km(k))/km(k))
    hmdiff = abs((hm_save(k)-hm(k))/hm(k))
    udiff = abs(uwind_save(k)-uwind(k))
    vdiff = abs(vwind_save(k)-vwind(k))
    pdiff = abs((p_save(k)-p(k))/p(k))
    if (tdiff > small) then
      print*,'k tdiff = ',k,tdiff
      print*,'t_save(k) = ',t_save(k)
      print*,'t(k) = ',t(k)
      stop 'tdiff too big: state variable temperature changed by convection'
    endif
    if (rvdiff > small) then
      stop 'rvdiff too big: state variable rv changed by convection'
    endif
    if (rldiff > small) then
      stop 'rldiff too big: state variable rl changed by convection'
    endif
    if (ridiff > small) then
      stop 'ridiff too big: state variable ri changed by convection'
    endif
    if (kmdiff > small) then
      stop 'kmdiff too big: state variable km changed by convection'
    endif
    if (hmdiff > small) then
      stop 'hmdiff too big: state variable hm changed by convection'
    endif
    if (udiff > small) then
      stop 'udiff too big: state variable uwind changed by convection'
    endif
    if (vdiff > small) then
      stop 'vdiff too big: state variable vwind changed by convection'
    endif
    if (pdiff > small) then
      stop 'pdiff too big: state variable p changed by convection'
    endif
  end do
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!  Define "heat" from tnew: 'PCONVT'
!
!  heat(pcols,pver)
!
!  call zm_convr      zm_convr
!  ptend_loc%s        heat(pcols,pver) ! heating rate (dry static energy tendency, W/kg)
!
!  Later, ptend_loc%s converted to K/s via:
!  ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
!  and output as 'ZMDT' in K/s.
!
!  This procedure is taken from if_driver.f90
!
!  - In this expression for dt, if dkm is J/kg dry air/day and drv is in
!    kg vapor/kg dry air/day, then dt will be K/day
!  - divide by 24*3600 to convert to K/s.
!  - multiply by cpair to get J/s (rate of change of dry enthalpy). 
!  - cpair should be a function of moisture (rv) but doesn't seem to be.
!  - Treating cp as temperature independent.
!
!----------------------------------------------------------------------
  sec_day = 24._r8*3600._r8
  rlnew_small = -1.0E-10
  rnnew_small = -1.0E-10
  rvnew_small = -1.0E-10
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!   Start loop over k
!   - rvtend, etc are arrays defined in if_conv_tend.f90 and refer to total
!     tendencies produced by up parcels.
!   - I should be updating surface pressure. In if_driver.f90 I do this by:
!     phalf(1) = phalf(1) + dpsurf_up + dpsurf_dn
!   - rlnew here is undefined: always zero?
!
!----------------------------------------------------------------------
  do k = 1,pver
    kmnew(k) = km(k) + kmtend(k)*(tstep/sec_day)
    rvnew(k) = rv(k) + rvtend(k)*(tstep/sec_day)
    rlnew(k) = rl(k) + rltend(k)*(tstep/sec_day)
    rinew(k) = ri(k) + ritend(k)*(tstep/sec_day)
    if (rlnew(k) < rlnew_small) then
!     print*,'WARNING rlnew negative in if_conv.F90 rlnew(=k) = ',rlnew(k)
!     rlnew(k) = 0.
    endif
    if (rvnew(k) < rvnew_small) then
      print*,'WARNING rvnew negative in if_conv.F90 rvnew(k) = ',rvnew(k)
!      stop
    endif
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!  - Note that cpair is actually the dry air specific heat (likely OK in this context
!    since DSE of model seems to use cpd rather than real specific heat).
!
!----------------------------------------------------------------------
    get_t_call = 3
    call t_from_km(kmnew(k),p(k),rvnew(k),rlnew(k),rinew(k),tnew(k),get_t_call)
    heat(i,pver+1-k) = ((tnew(k) - t(k))*cpair)/tstep
!   print*,'f_conv.F90 k tnew(k) - t(k) = ',tnew(k)-t(k)
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!   Check for unrealistic T increases
!
!----------------------------------------------------------------------
    if ( ((tnew(k) - t(k)) > 20.).or.(tnew(k) > temp_hot)) then
      print*,'-----------------------------------'
      print*,'WARNING: BIG T CHANGE in if_conv.F90'
      print*,'landfrac(i) = ',landfrac(i)
      print*,'colrh = ',colrh
      print*,'cape_1 = ',cape_1
      print*,'cape_2 = ',cape_2
      print*,'cape_3 = ',cape_3
      print*,'cape_4 = ',cape_4
      print*,'cfraction_1 = ',cfraction_1
      print*,'cfraction_2 = ',cfraction_2
      print*,'cfraction_3 = ',cfraction_3
      print*,'cfraction_4 = ',cfraction_4
      print*,'rh(1) = ',rh(1)
      print*,'rh(2) = ',rh(2)
      print*,'rh(3) = ',rh(3)
      print*,'rh(4) = ',rh(4)
      print*,'t(1) = ',t(1)
      print*,'t(2) = ',t(3)
      print*,'t(3) = ',t(3)
      print*,'t(4) = ',t(4)
      print*,'i clat(i) = ',i,clat(i)
      print*,'k = ',k
      print*,'tnew(k)-t(k) = ',tnew(k)-t(k)
      print*,'p(k) = ',p(k)
      print*,'old rl(k) = ',rl(k)
      print*,'rlnew(k) = ',rlnew(k)
      print*,'old rv(k) = ',rv(k)
      print*,'rvnew(k) = ',rvnew(k)
      print*,'old km(k) = ',km(k)
      print*,'old t(k) = ',t(k)
      print*,'tnew(k) = ',tnew(k)
      print*,'in mm/day uprain_start_rlp = ',uprain_start_rlp*3600.*24.
      print*,'in mm/day uprain_surf_rlp = ',uprain_surf_rlp*3600.*24.
      print*,'in mm/day uprain_surf_rv = ',uprain_surf_rv*3600.*24.
      print*,'in mm/day starting ansnow_conv = ',ansnow_conv*3600.*24.
      print*,'in mm/day starting ansnow_strat = ',ansnow_strat*3600.*24.
      print*,'in mm/day anrain_surf = ',anrain_surf*3600.*24.
      print*,'in mm/day ansnow_surf = ',ansnow_surf*3600.*24.
      do kk = 1,pver
        print*,'kk av_mp(k) = ',kk,z(kk)*0.001,av_mp(kk)
      end do
      print*,'------------ WONKY --------------------'
!     stop 'WONKY temperature issues'
    endif
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!   Check for very cold surface temperatures
!
!----------------------------------------------------------------------
    if ((k == 1).and.(tnew(k) < t_surf_cold)) then
      print*,'-----------------------------------'
      print*,'WARNING: VERY COLD SURFACE TEMP if_conv.F90'
      print*,'tnew(1) = ',tnew(1)
      print*,'p(1) = ',p(1)
      print*,'landfrac(i) = ',landfrac(i)
      print*,'colrh = ',colrh
      print*,'cape_1 = ',cape_1
      print*,'cape_2 = ',cape_2
      print*,'cape_3 = ',cape_3
      print*,'cape_4 = ',cape_4
      print*,'cfraction_1 = ',cfraction_1
      print*,'cfraction_2 = ',cfraction_2
      print*,'cfraction_3 = ',cfraction_3
      print*,'cfraction_4 = ',cfraction_4
      print*,'i clat(i) = ',i,clat(i)
      print*,'tnew(k)-t(k) = ',tnew(k)-t(k)
      print*,'p(k) = ',p(k)
      print*,'old rl(k) = ',rl(k)
      print*,'rlnew(k) = ',rlnew(k)
      print*,'old rv(k) = ',rv(k)
      print*,'rvnew(k) = ',rvnew(k)
      print*,'old km(k) = ',km(k)
      print*,'old t(k) = ',t(k)
      print*,'tnew(k) = ',tnew(k)
      print*,'in mm/day uprain_start_rlp = ',uprain_start_rlp*3600.*24.
      print*,'in mm/day uprain_surf_rlp = ',uprain_surf_rlp*3600.*24.
      print*,'in mm/day uprain_surf_rv = ',uprain_surf_rv*3600.*24.
      print*,'in mm/day starting ansnow_conv = ',ansnow_conv*3600.*24.
      print*,'in mm/day starting ansnow_strat = ',ansnow_strat*3600.*24.
      print*,'in mm/day anrain_surf = ',anrain_surf*3600.*24.
      print*,'in mm/day ansnow_surf = ',ansnow_surf*3600.*24.
      do kk = 1,pver
        print*,'kk av_mp(k) = ',kk,z(kk)*0.001,av_mp(kk)
      end do
      print*,'------------ WONKY --------------------'
!     stop 'WONKY temperature issues'
    endif
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!----------------------------------------------------------------------
    nprint_surf = 0
    if (nprint_surf == 1) then
    if (k == 1) then
      print*,'-------- SURFACE ANALYSIS in if_conv.F90  ---------------------'
      print*,'p(1) = ',p(1)
      print*,'km(1) = ',km(1)
      print*,'rv(1) = ',rv(1)
      print*,'rvnew(1) = ',rvnew(1)
      print*,'change in rv = ',rvnew(1) - rv(1)
      print*,'enthalpy change due to change in rv = ',(rvnew(1) - rv(1))*latent(t(1))
      print*,'enthalpy change due to change in T = ',(tnew(1) - t(1))*cpd
      print*,'rl(1) = ',rl(1)
      print*,'rlnew(1) = ',rlnew(1)
      print*,'change in rl = ',rlnew(1) - rl(1)
      print*,'t(1) = ',t(1)
      print*,'tnew(1) = ',tnew(1)
      print*,'change in surface temperature in if_conv.F90 = ',tnew(1)-t(1)
      print*,'surface  DSE tendency in if_conv.F90 = ',heat(i,pver)
      print*,'-------------------------------------------------------'
    endif
    endif
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!  Define qtnd from rvtend: 'PCONVB'
!
!  call zm_convr        zm_convr
!  ptend_loc%q(:,:,1)   qtnd(pcols,pver): specific humidity tendency (kg/kg/s)
!
!  rv = q/(1-q)
!  q = rv/(1+rv) 
!
!  How to I get qnew from rvnew?
!  - have to assume some constraint
!  - I think relevant constraint is that total mass in layer m = md+mv+ml is fixed
!
!  rvnew = (mv + dmv)/mdnew 
!  qvnew = (mv + dmv)/(mdnew + mvnew + mlnew) 
!  Define m = mdnew + mvnew + mlnew = md + mv + ml
!  qvnew = (mv + dmv)/m = rvnew*mdnew/m
!
!  m/mdnew = (mdnew + mvnew + mlnew)/mdnew = 1 + rvnew + rlnew + rnnew
!
!  So: qvnew = rvnew/(1 + rvnew + rlnew +rnnew)
!
!  since md/m = 1/(1+rv+rl+rn)
!
!  Technically the stratiform mass should be included in these conversions.
!  - probably true actually, but then would have to be fed into conv param,
!    so knew where mass of pressure intervals came from
!
!----------------------------------------------------------------------
    qvnew(k) = rvnew(k)/(1.+rvnew(k)+rlnew(k)+rinew(k))
    qtnd(i,pver+1-k) = (qvnew(k) - qv(k))/tstep
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!  Do all other constituent tendencies:
!
! call convtran <-> in convtran
!    ptend%lq       doconvtran(m)
!    state%q        real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array including moisture
!    pcnst          ncnst
!    ptend%q        real(r8), intent(out) :: dqdt(pcols,pver,ncnst)
!
!  Not sure why change from pcnst to ncnst
!  Try here to be as consistent as possible with convtran, except replace
!    q with qh.
!
!     Where is doconvtran set? Fed in as ptend_loc%lq as a subroutine argument.
!     Where is dqdt initialized to zero?
!
!   Feb 2013: ql(k) should be zero here. Any rlnew used to produce CLDLIQ.
!
!----------------------------------------------------------------------
    do m = 2, ncnst
      if (doconvtran(m)) then
        if (m == ixcldliq) then
          qlnew(k) = rlnew(k)/(1. + rvnew(k) + rlnew(k) + rinew(k))
          dqdt(i,pver+1-k,m) = (qlnew(k) - ql(k))/tstep
        elseif (m == ixcldice) then
          qinew(k) = rinew(k)/(1. + rvnew(k) + rlnew(k) + rinew(k))
          dqdt(i,pver+1-k,m) = (qinew(k) - qi(k))/tstep
        else
!         stop 'tracer not defined in module if_conv'
        endif
      endif
    end do
!----------------------------------------------------------------------
!
!  End loop over k
!
!----------------------------------------------------------------------
  end do
!----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!  Update momentum tendencies: 'ZMMTU' and 'ZMMTV'
!
! call momtran                          <->  in momtran
!     2                                      ncnst
! real(r8) :: winds(pcols, pver, 2)          q 
! real(r8) :: wind_tends(pcols, pver, 2)     real(r8), intent(out) :: dqdt(pcols,pver,ncnst)
!
!  Here:
!  call if_convr      if_convr
!  winds              real(r8), intent(in) :: wind(pcols,pver,mcnst)  ! Wind array
! real(r8) :: wind_tends(pcols, pver, 2)     real(r8), intent(out) :: dwdt(pcols,pver,mcnst)
!
!  ncnst -> mcnst
!
!  In addfield, 'ZMMTU' has units kg/m2/s
!  In addfield, 'ZMMTV' has units kg/m2/s
!
!  my uwind (per kg dry mass) will always be larger than actual u.
!  Therefore my tendencies would be too large, so must reduce.
!
!----------------------------------------------------------------------
  do m = 1, mcnst   
    if (domomtran(m)) then
    if (do_conv_mom_transport == 1) then  ! switch set in if_conv_tend.f90
      do k = 1,pver
        rtt = rv(pver+1-k) + rl(pver+1-k) 
        if (m == 1) then
          dwdt(i,k,m) = utend(pver+1-k)/(1. + rtt)
          dwdt(i,k,m) = dwdt(i,k,m)/sec_day
!          print*,'updating u tendency = ',k,dwdt(i,k,m)/sec_day
        elseif (m == 2) then
          dwdt(i,k,m) = vtend(pver+1-k)/(1. + rtt)
          dwdt(i,k,m) = dwdt(i,k,m)/sec_day
!          print*,'updating v tendency = ',k,dwdt(i,k,m)/sec_day
        else
          stop 'm out of range in momentum tendency dwdt definition'
        endif
      end do
    endif
    endif
  end do
! print*,'Got tendencies'
!-----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!   OK TO GET RID OF THIS?????
!
!  Define jctop/jcbot:
!
!  call zm_convr         zm_convr
!     jctop               jctop
!
!  - for physics_update, need top and bottom indices over which tendencies non-zero.
!  - loops go over all i=1,ncol (don't need gathered).
!  - presumably top/bot indices where rv/km tendencies exceed some threshold
!  - scattered quantities
!  - in zm jcbot is maxg level (max mse)
!  - in zm jctop is jt (basically lel, undilute detrainmnet level).
!
!  I think: convection occurs between surface level k = pver (usually) up
!     to level k = limcnv. Above, msg = limcnv - 1, so all levels from 
!     k = limcnv - 1 = msg to k = 1, should be free of convection.
!  This is a loop over k = msg + 1,pver where convective mass fluxes are
!     non-zero. msg + 1 + limcnv
!
!  How should I define jctop?
!  Starting from top: first level at which ttend exceeds some threshold
!  If there is no convection, should have jctop = jcbot = pver
!
!  in zm, initialization is 
!      jctop(i) = pver
!      jcbot(i) = 1
!  so have retained here 
!
!-----------------------------------------------------------------------
  km_small = 1.0E-02
  jctop(i) = pver
  do k = pver,1,-1
    if (kmtend(k) > km_small) then
      jctop(i) = pver+1-k
    endif
  end do
  jcbot(i) = 1
  do k = 1,pver
    if (kmtend(k) > km_small) then
      jcbot(i) = pver+1-k
    endif
  end do
!-----------------------------------------------------------------------
!
! ************* Inside did_convect = 1 ********************************
!
!  Define total surface 'IF' rain: 'PRECZ'  (m/s)
!
!  Should also define surface 'IF' snow
!
!  output     call zm_convr  zm_convr
!  'PRECZ'    prec           real(r8), intent(out) :: prec(pcols)
!
!  There are other addfld precipitation outputs defined, but the others seem
!    to be output only after zm_conv_evap)
!  - what to do about snow? Just give total
!
!  There is a variable snow_zmc that comes out of zm_conv_evap. I should
!    probably define this as an remaining unmelted snow ...
!
!  Not clear how energy conservation differentiates between surface rain/snow.
!
!  How is prec used to test water conservation?
!  In subroutine tphysbc, prec is called prec_zmc and is used to calculate the flux
!    condition for condensed water flz_cnd.
!
! flx_cnd(:ncol) = prec_zmc(:ncol) + rliq(:ncol)
! call check_energy_chng(state, tend, "convect_deep", nstep, ztodt, zero, flx_cnd, snow_zmc, zero)
!
!
! call diag_conv(state, ztodt,    &
!   prec_zmc, snow_zmc, prec_cmf, snow_cmf, prec_sed, snow_sed, prec_pcw, snow_pcw)
!
! Also convert prec from kg/m2/s to m/s. 1 kg/m2 = 1 mm, so just divide by 1000.
!
!-----------------------------------------------------------------------
  rain_if(i) = uprain_surf_rv + uprain_surf_rlp + anrain_surf  ! kg/m2/s
  snow_if(i) = ansnow_surf
  rain_if(i) = 0.001_r8*rain_if(i)
  snow_if(i) = 0.001_r8*snow_if(i)
  prec_if(i) = rain_if(i) + snow_if(i)

!-----------------------------------------------------------------------
!
!  Else No convection (did_convect = 0)
!
!-----------------------------------------------------------------------
  else
!-----------------------------------------------------------------------
!
!  Zero all tendencies
!  - this is probably done elsewhere so likely not neccessary
!
!-----------------------------------------------------------------------
  do k = 1,pver
    heat(i,k) = 0.
    qtnd(i,k) = 0.
    do m = 1, ncnst
      dqdt(i,k,m) = 0.
    end do
    do m = 1, mcnst   
      dwdt(i,k,m) = 0.
    end do
    qvnew(k) = qv(k)
    rvnew(k) = rv(k)
    rlnew(k) = rl(k)
    rinew(k) = ri(k)
    tnew(k) = t(k)
    if ((k == 1).and.(i == 1)) then
      print*,'SETTING rlnew equal rl = ',rl(k)
    endif
  end do
  jctop(i) = pver
  jcbot(i) = 1
!-----------------------------------------------------------------------
!
!  Zero all rain rates
!
!-----------------------------------------------------------------------
  rain_if(i) = 0.
  snow_if(i) = 0.
  prec_if(i) = 0.

!-----------------------------------------------------------------------
!
!  endif for did_convect = 1
!
!-----------------------------------------------------------------------
endif

!----------------------------------------------------------------------
!
!  Diagnostics: 
!  - why is rlnew different from rl_tar at the surface?
!
!----------------------------------------------------------------------
nprint = 0
if (nprint == 1) then
if (prec_if(i) < 0.00001) then
if (abs(clat(i))*57.296_r8 < 30.) then
  print*,'------------------'
  print*,'latitude = ',clat(i)*57.296_r8
  print*,'g/kg rl(1) = ',1000.*rl(1)
  print*,'g/kg rlnew(1) = ',1000.*rlnew(1)
  print*,'g/kg rlnew(1) - rl(1) = ',1000.*(rlnew(1)-rl(1))
  print*,'g/kg rl_tar(1) = ',1000.*rl_tar(1)
  print*,'fractional deviation rlnew from rl_tar = ',(rlnew(1)-rl_tar(1))/rl_tar
  print*,'g/kg/day one timestep rltend(1) = ',1000.*(tstep/sec_day)*rltend(1)
  print*,'g/kg/day one timestep rl_vert_tend(1) = ',1000.*(tstep/sec_day)*rl_vert_tend(1)
  print*,'ce_fact = ',ce_fact
  print*,'rh(1) = ',rh(1)
endif
endif
endif

!-----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  NOW ON: Should be mostly diagnostics stuff that needs to be defined
!   whether did convection or not.
!
!  If did not do convection, all flow and rain variables should be zeroed
!   in init_zero call anyway.
!
!
!  Define 4D "diag" variables for history file
!  - convert to m/s (strange units but just to be consistent with cam4)
!  - flip height index to be consistent with h0 history (pver + 1 = surface)
!  - retain mass flux units as kg/m2/s 
!
!-----------------------------------------------------------------------
  uprain_1_diag(i) = 0.001_r8*uprain_1  ! now in m/s
  uprain_2_diag(i) = 0.001_r8*uprain_2  ! now in m/s
  uprain_evap_diag(i) = 0.001_r8*uprain_evap  ! now in m/s
  anrain_down_diag(i) = 0.001_r8*anrain_down  ! now in m/s
  anrain_evap_diag(i) = 0.001_r8*anrain_evap  ! now in m/s
  ansnow_subl_diag(i) = 0.001_r8*ansnow_subl  ! now in m/s
  ansnow_melt_diag(i) = 0.001_r8*ansnow_melt  ! now in m/s
  uprain_down_diag(i) = 0.001_r8*uprain_down  ! now in m/s
  uprain_surf_rlp_diag(i) = 0.001_r8*uprain_surf_rlp  ! now in m/s
  uprain_surf_rv_diag(i) = 0.001_r8*uprain_surf_rv  ! now in m/s
  uprain_start_rlp_diag(i) = 0.001_r8*uprain_start_rlp  ! now in m/s
  uprain_rlp_rvv_diag(i) = 0.001_r8*uprain_rlp_rvv  ! now in m/s

  ansnow_conv_diag(i) = 0.001_r8*ansnow_conv
  ansnow_strat_diag(i) = 0.001_r8*ansnow_strat
  ansnow_strat_rv_diag(i) = 0.001_r8*ansnow_strat_rv
  ansnow_strat_ri_diag(i) = 0.001_r8*ansnow_strat_ri
  ansnow_surf_diag(i) = 0.001_r8*ansnow_surf
  anrain_surf_diag(i) = 0.001_r8*anrain_surf

! print*,'--------'
  do k = 1,pver+1
    fmass_diag(i,k) = fmass(pver+2-k)
    fmass_dn_diag(i,k) = fmass_dn(pver+2-k)
    phalf_diag(i,k) = phalf(pver+2-k)
!   print*,'k fmass_dn(pver+2-k) = ',k,fmass_dn(pver+2-k)
  end do

  do k = 1,pver
    ent_diag(i,k) = dmflow(1,pver+1-k)
    det_diag(i,k) = dmflow(2,pver+1-k)
    if ((ent_diag(i,k) < 0.).or.(det_diag(i,k) < 0.)) then
      print*,'ent_diag(i,k) = ',ent_diag(i,k)
      print*,'det_diag(i,k) = ',det_diag(i,k)
      print*,'i = ',i
      print*,'k = ',k
      print*,'dmflow(1,pver+1-k) = ',dmflow(1,pver+1-k)
      print*,'dmflow(2,pver+1-k) = ',dmflow(2,pver+1-k)
      stop 'DET or ENT negative problem'
    endif
    updet_1_diag(i,k) = updet_1(pver+1-k)
    updet_2_diag(i,k) = updet_2(pver+1-k)
    rh_dn_diag(i,k) = rh_dn(pver+1-k)
    drv_dndet_diag(i,k) = drv_dndet(pver+1-k)
    dt_dndet_diag(i,k) = dt_dndet(pver+1-k)
    dndet_diag(i,k) = dndet(pver+1-k)
    dnent_diag(i,k) = dnent(pver+1-k)
    mtdet_diag(i,k) = mtdet(pver+1-k)
    mtent_diag(i,k) = mtent(pver+1-k)
    drv_up_diag(i,k) = drv_up(pver+1-k)
    drv_an_diag(i,k) = drv_an(pver+1-k)
  end do

  cape_mean_diag(i) = cape_mean

  LTS_diag(i) = LTS
  EIS_diag(i) = EIS
  mse_rat_diag(i) = mse_rat

  cape_1_diag(i) = cape_1
  cape_2_diag(i) = cape_2
  cape_3_diag(i) = cape_3
  cape_4_diag(i) = cape_4

  amp_diag(i) = amp
  amp_rain_diag(i) = amp_rain
  amp_om5_diag(i) = amp_om5

  dp_dn_diag(i) = dp_dn

  f_1_diag(i) = f_1
  f_2_diag(i) = f_2
  f_3_diag(i) = f_3
  f_4_diag(i) = f_4

  colwat_diag(i) = colwat
  colrh_diag(i) = colrh
  cf_tar_diag(i) = cf_tar 
  lwp_tar_diag(i) = lwp_tar
  conv_energy_diag(i) = 0.001*conv_energy  ! convert to kJ
  conv_energy_all(i) = 0.001*conv_energy  ! convert to kJ
  i_cf_all(i) = i_cf
  cf_tar_all(i) = cf_tar

if ((cf_tar > 0.5).and.(lwp_tar < 10.)) then
  print*,'------  WIERD results in if_conv.F90 ---'
  print*,'latitude = ',clat(i)*57.296_r8
  print*,'conv_energy = ',conv_energy*0.001
  print*,'cf_tar = ',cf_tar
  print*,'lwp_tar = ',lwp_tar
  stop 'lwp_tar problem'
endif

!-----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  Start variables
!-----------------------------------------------------------------------
  mass_start_1_diag(i) = mass_start_1
  mass_start_2_diag(i) = mass_start_2

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  cfraction_1_diag(i) = cfraction_1
  cfraction_2_diag(i) = cfraction_2
  cfraction_3_diag(i) = cfraction_3
  cfraction_4_diag(i) = cfraction_4

  hmuprain_start_diag(i) = hmuprain_start
  hmuprain_surf_diag(i) = hmuprain_surf
  hmansnow_start_diag(i) = hmansnow_start
  hmansnow_surf_diag(i) = hmansnow_surf
  hmanrain_surf_diag(i) = hmanrain_surf

!-----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
! Diagnostic arrays for nonlinearities
! - these should be defined in the call to get_tend whether convection done
!   or not.
!-----------------------------------------------------------------------
  km_width_diag(i) = km_width

  col_ri_dt_save(i) = col_ri_dt
  col_ri_rh_save(i) = col_ri_rh
  col_rl_dt_save(i) = col_rl_dt
  col_rl_rh_save(i) = col_rl_rh

  md_rat_diag(i) = md_rat
  f_an_diag(i) = f_an
  f_ice_diag(i) = f_ice

!-----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  For diagnostics:
!  These quantities should be calculated in if_conv_tend.f90
!  precip_org should be in units of mm/day
!  Not that storing the value of precip_org of the previous timestep that
!   is used in the current calculation on nonlinear variables.
!-----------------------------------------------------------------------
  av_z_ri_diag(i) = av_z_ri
  av_z_ri_rh_diag(i) = av_z_ri_rh
  av_z_ri_dt_diag(i) = av_z_ri_dt
  zpeak_diag(i) = zpeak
  precip_org_diag(i) = precip_org
  updraft_eff_eff(i) = updraft_eff
  if ((updraft_eff_eff(i) < -490.).and.(updraft_eff_eff(i) > -500.)) then
    print*,'i updraft_eff_eff(i) in if_conv.F90 = ',i,updraft_eff_eff(i)
  endif
  updraft_eff_det(i) = updraft_det
  updraft_eff_rlp(i) = updraft_rlp
! print*,'------------------ i = ',i
! print*,'rcmax_lower = ',rcmax_lower

!-----------------------------------------------------------------------
!
!  Define rhnew (diagnostics)
!  XXX
!
!-----------------------------------------------------------------------
  do k = 1,pver
    rss = rsat(tnew(pver+1-k),p(pver+1-k))
    rhnew(i,k) = rvnew(pver+1-k)/rss
  end do

!-----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  CFF
!  Switch used in cldwat.F90
!  - turns off precip production from cldice/cldwat
!  - turns off cme - cond/evap of cldice/cldwat
!
!-----------------------------------------------------------------------
! print*,'lchnk i this_lat = ',lchnk,i,this_lat

  ifconvection_activated(i) = ifconvection
  cape_meann(i) = cape_mean
! print*,'in if_conv.F90 i ifconvection_activated(i) = ',i,ifconvection_activated(i)

!-----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  JUST FOR DIAGNOSTICS
!
!  These should be all /day from if_conv_tend.f90: don't change
!   - tendencies passed up via qtnd, just for diagnostics
!   - the up tendencies are just for diagnostic purposes in if_conv_tend.f90
!
!-----------------------------------------------------------------------
  do k = 1,pver
    rhtend_real(i,k) = (rhnew(i,k)-rhif_store(i,k))*((3600.*24.)/tstep)
    tptend_real(i,k) = (tnew(pver+1-k) - t(pver+1-k))*((3600.*24.)/tstep)
    qtend_real(i,k) = (qvnew(pver+1-k) - qv(pver+1-k))*((3600.*24.)/tstep)
!-----------------------------------------------------------------------
!  Look for large values
!-----------------------------------------------------------------------
    if (abs(qtend_real(i,k)) > 0.4) then
      print*,'========  LARGE VALUE OF QTEND_REAL in if_conv.F90 ====='
      print*,'qtend_real(i,k) = ',qtend_real(i,k)
      print*,'did_convect = ',did_convect
      print*,'i k = ',i,k
      print*,'pver+1-k = ',pver+1-k
      print*,'p(pver+1-k) = ',p(pver+1-k)
      print*,'km(pver+1-k) = ',km(pver+1-k)
      print*,'qvnew(pver+1-k) = ',qvnew(pver+1-k)
      print*,'qv(pver+1-k) = ',qv(pver+1-k)
      print*,'tnew(pver+1-k) = ',tnew(pver+1-k)
      print*,'t(pver+1-k) = ',t(pver+1-k)
      print*,'rvnew(pver+1-k) = ',rvnew(pver+1-k)
      print*,'rv(pver+1-k) = ',rv(pver+1-k)
      print*,'rlnew(pver+1-k) = ',rlnew(pver+1-k)
      print*,'rl(pver+1-k) = ',rl(pver+1-k)
      print*,'uprain_start_rlp in mm per day = ',3600.*24.*uprain_start_rlp
      print*,'ansnow_conv in mm per day = ',3600.*24.*ansnow_conv
      print*,'ansnow_strat in mm per day = ',3600.*24.*ansnow_strat
      print*,'uprain_rlp(pver+1-k) in mm per day = ',3600.*24.*uprain_rlp(pver+1-k)
      print*,'uprain_rv(pver+1-k) in mm per day = ',3600.*24.*uprain_rv(pver+1-k)
      print*,'ansnow_rlp(pver+1-k) in mm per day = ',3600.*24.*ansnow_rlp(pver+1-k)
      print*,'ansnow_rv(pver+1-k) in mm per day = ',3600.*24.*ansnow_rv(pver+1-k)
      print*,'========== GO ON  ==='
    endif
!-----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!-----------------------------------------------------------------------
    kmtend_save(i,k) = kmtend(pver+1-k)
    rvtend_save(i,k) = rvtend(pver+1-k)
    rltend_save(i,k) = rltend(pver+1-k)
    rlnew_save(i,k) = rlnew(pver+1-k)
    ritend_save(i,k) = ritend(pver+1-k)
    utend_save(i,k) = utend(pver+1-k)
    vtend_save(i,k) = vtend(pver+1-k)
!-----------------------------------------------------------------------
!  variables scaling with updraft mass flux
!-----------------------------------------------------------------------
    mp_save(i,k) = av_mp(pver+1-k)
    bp_cape_save(i,k) = bp_cape(pver+1-k)
    bp1_str_save(i,k) = bp1_str(pver+1-k)
    bp1_prp_save(i,k) = bp1_prp(pver+1-k)
    bp1_itr_save(i,k) = bp1_itr(pver+1-k)
    bp2_str_save(i,k) = bp2_str(pver+1-k)
    bp2_prp_save(i,k) = bp2_prp(pver+1-k)
    bp2_itr_save(i,k) = bp2_itr(pver+1-k)
    f_ent_1_save(i,k) = f_ent_1(pver+1-k)
    f_det_1_save(i,k) = f_det_1(pver+1-k)
    f_ent_2_save(i,k) = f_ent_2(pver+1-k)
    f_det_2_save(i,k) = f_det_2(pver+1-k)
    rlp_save(i,k) = av_rlp(pver+1-k)
    rip_save(i,k) = av_rip(pver+1-k)
    if (hmp_mean(pver+1-k) /= bad) then 
      hmp_mean_save(i,k) = 0.001*hmp_mean(pver+1-k)
    else
      hmp_mean_save(i,k) = hmp_mean(pver+1-k)
    endif
    if (hm_diff(pver+1-k) /= bad) then 
      hm_diff_save(i,k) = 0.001*hm_diff(pver+1-k)
    else
      hm_diff_save(i,k) = hm_diff(pver+1-k)
    endif
!-----------------------------------------------------------------------
!  variables scaling with updraft detrainment mass
!-----------------------------------------------------------------------
    tdiff_detrain_up_save(i,k) = tdiff_detrain_up(pver+1-k)
    cape_detrain_up_save(i,k) = cape_detrain_up(pver+1-k)
    mass_detrain_up_save(i,k) = mass_detrain_up(pver+1-k)
    if (hm_diff_det(pver+1-k) /= bad) then 
      hm_diff_det_save(i,k) = 0.001*hm_diff_det(pver+1-k)
    else
      hm_diff_det_save(i,k) = hm_diff_det(pver+1-k)
    endif
!-----------------------------------------------------------------------
!  Evaporation variables 
!-----------------------------------------------------------------------
    ansnowsubl_save(i,k) = ansnowsubl(pver+1-k)
    anrainevap_save(i,k) = anrainevap(pver+1-k)
    uprainevap_save(i,k) = uprainevap(pver+1-k)
    ansnowmelt_save(i,k) = ansnowmelt(pver+1-k)
!-----------------------------------------------------------------------
!  Down variables
!-----------------------------------------------------------------------
    upraindown_save(i,k) = upraindown(pver+1-k)
    anraindown_save(i,k) = anraindown(pver+1-k)
!-----------------------------------------------------------------------
!  Precip sources
!-----------------------------------------------------------------------
    uprain_rlp_save(i,k) = uprain_rlp(pver+1-k)
    ansnow_rlp_save(i,k) = ansnow_rlp(pver+1-k)
    uprain_rv_save(i,k) = uprain_rv(pver+1-k)
    ansnow_rv_save(i,k) = ansnow_rv(pver+1-k)
    rlp_rvv_save(i,k) = rlp_rvv(pver+1-k)
!-----------------------------------------------------------------------
!  Cloud sources
!  - units should be kg/kg/day
!-----------------------------------------------------------------------
    uprain_rl_tend_save(i,k) = uprain_rl_tend(pver+1-k) 
    rl_evap_tend_save(i,k) = rl_evap_tend(pver+1-k) 
    ri_evap_tend_save(i,k) = ri_evap_tend(pver+1-k)  
    ri_prod_tend_save(i,k) = ri_prod_tend(pver+1-k)  
    rl_det_tend_save(i,k) = rl_det_tend(pver+1-k)  
    ri_det_tend_save(i,k) = ri_det_tend(pver+1-k)
    rl_vert_tend_save(i,k) = rl_vert_tend(pver+1-k)  
    ri_vert_tend_save(i,k) = ri_vert_tend(pver+1-k)  
    rl_dt_save(i,k) = rl_dt(pver+1-k)  
    rl_rh_save(i,k) = rl_rh(pver+1-k)  
    ri_dt_save(i,k) = ri_dt(pver+1-k)  
    ri_rh_save(i,k) = ri_rh(pver+1-k)  
    ri_rs_save(i,k) = ri(pver+1-k)/rs(pver+1-k)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    rhif_store(i,k) = rh(pver+1-k)   ! input RH to convection
    tif_store(i,k) = t(pver+1-k)
    if (hm(pver+1-k) /= bad) hm_store(i,k) = 0.001*hm(pver+1-k)
!-----------------------------------------------------------------------
!  Calc prob_random
!-----------------------------------------------------------------------
ttt = tnew(pver+1-k)
prr = sigmoidal(ttt,t_pr_half,t_pr_scale,pr_min,pr_add)
if (use_add_lower == 1) then
if (ttt > t_add_lower) then
   prr = prr + (ttt - t_add_lower)*pr_add_inc
endif
endif
dzz = zh(pver+2-k) - zh(pver+1-k)
prob = prr*(dzz/1000.)
prob_random_store(i,k) = prob

ri_half_store(i,k) = sigmoidal(ttt,t_ri_half,t_ri_scale,ri_min,ri_add)
rl_half_store(i,k) = rl_scale
if (ri_half_store(i,k) < 0.) then
    stop 'ri_half is negative' 
endif
end do

!-----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  Define column new and old column quantities for testing conservation.
!  - in the case of surface precipitation the surface pressure phalf(1)
!    should be reduced. I have not done this here, since CAM4 does not
!    change the surface pressure. But I calculate the corrected column
!    water error below in err_wat_corr.
!
!-----------------------------------------------------------------------
  col_rv = 0._r8
  col_rl = 0._r8
  col_ri = 0._r8
  col_km = 0._r8
  col_rv_old = 0._r8
  col_rl_old = 0._r8
  col_ri_old = 0._r8
  col_km_old = 0._r8
  do k = 1,pver
    dmdry_old = (phalf(k) - phalf(k+1))/(g*(1. + rv(k) + rl(k) + ri(k)))
    dmdry_new = (phalf(k) - phalf(k+1))/(g*(1. + rvnew(k) + rlnew(k) + rinew(k)))
    col_rv = col_rv + dmdry_new*rvnew(k)
    col_rl = col_rl + dmdry_new*rlnew(k)
    col_ri = col_ri + dmdry_new*rinew(k)
    col_rv_old = col_rv_old + dmdry_old*rv(k)
    col_rl_old = col_rl_old + dmdry_old*rl(k)
    col_ri_old = col_ri_old + dmdry_old*ri(k)
    col_km_old = col_km_old + dmdry_old*km(k)
    col_km = col_km + dmdry_new*kmnew(k)
  end do
  col_wat = col_rv + col_rl + col_ri
  col_wat_old = col_rv_old + col_rl_old + col_ri_old
  hm_uprain(i) = hmuprain_surf
  hm_ansnow(i) = hmansnow_surf

!----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  Water error analysis
!  - units here look like mm water (or kg/m2) in a time step
!  - rain_if and snow_if should be the actual preecip reaching the surface.
!
!----------------------------------------------------------------------
  totprecip = rain_if(i)*1000._r8*tstep + snow_if(i)*1000._r8*tstep  ! have to convert from m/s
  dp_precip = totprecip*g
  dp_precip_dry = dp_precip/(1. + rvnew(1) + rlnew(1) + rinew(1))
  dp_wat = (dp_precip_dry/g)*(rvnew(1)+rlnew(1)+rinew(1))
  err_wat = totprecip + col_wat - col_wat_old  ! no correction for change in ps
  err_wat_corr = err_wat - dp_wat   ! corrected for change in ps
! print*,'err_wat = ',err_wat

!----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  Assign max allowable water error (kg/m2 of water in a timestep)
!  - probably multiply by 48 to get error in mm/day
!----------------------------------------------------------------------
! err_wat_corr_max = 0.000001
  err_wat_corr_max = 0.000005   ! inc after change BL scheme neg rlnew
!----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!   ************   ALL COLUMNS   ******************************
!
!  Store nonlinear variables
!
!  - these new values (current timestep) are defined in if_conv_tend.f90
!  - The old values are being over-written in this operation
!  - They are converted to pointers higher up and used in the next timestep
!  - precip_org should be in mm/day
!
!----------------------------------------------------------------------
  preciporg(i) = precip_org_new
  massprev(i) = mass_start_tot
  cm5prev(i) = cm5_new
  mass_inc_diag(i) = mass_start_tot - mass_prev
  if (mass_start_tot > 0.001) then
    mass_rat_diag(i) = (mass_start_tot - mass_prev)/mass_start_tot
  else
    mass_rat_diag(i) = 0.
  endif

!----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  Energy error analysis
!  - should give J/m2
!
!----------------------------------------------------------------------
  err_km = col_km + (hmuprain_surf + hmansnow_surf + hmanrain_surf)*tstep - col_km_old
  dp_km = (dp_precip_dry/g)*kmnew(1)
  err_km_corr = (err_km - dp_km)/tstep !  convert to W/m2
! print*,'i ps corrected error in W per m2 err_km_corr = ',i,err_km_corr

!----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!   ************   ALL COLUMNS   ******************************
!
!  Write out stuff if enthalpy error exceeds tolerance
!
!----------------------------------------------------------------------
  nprint_en = 1
  if (nprint_en == 1) then
  if (abs(err_km_corr) > 0.001) then ! should be W/m2
    print*,'==== enthalpy error exceeds tolerance: COLUMN ENERGY BUDGET ANALYSIS in if_conv.F90 ==='
    print*,'total precipitation reaching the surface in mm per timestep totprecip = ',totprecip
    print*,'dp_precip_dry = ',dp_precip_dry
    print*,'dp_km = ',dp_km
    print*,'kmnew(1) = ',kmnew(1)
    print*,'enthalpy error before ps correction err_km (J/m2) = ',err_km
    print*,'corrected (real) column enthalpy error (W/m2) err_km_corr = ',err_km_corr
    print*,'col_km = ',col_km
    print*,'col_km_old = ',col_km_old
    print*,'hmuprain_surf = ',hmuprain_surf
    print*,'hmansnow_surf = ',hmansnow_surf
    print*,'hmanrain_surf = ',hmanrain_surf
    print*,'tstep = ',tstep
    print*,'========================================================'
  endif
  endif

!----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!   ************   ALL COLUMNS   ******************************
!
!  EEE  check_energy type analysis
!
!  This is a separate reality check tham enthalpy conservation. The
!    model could be moving enthalpy around in an unrealistic way and
!    still conserve column enthalpy.
!
!  Dropping KE, cam4 calculates column energy change from vertical 
!  integral of:
!   
!  TE = dse + (lv0+lf0)*wv + lf0*wl
! 
!  must be consistent with the fluxes (dropping sensible and latent heat):
!
!  TF = - (flx_cnd(i) - flx_ice(i))*latice
!     = - total_rain*latice
!
!   dTE/dt = TF
!   VERT INT{ddse/dt} + (lv0+lf0)*d(wv_col)/dt + lf0*d(wl_col) = -total_rain*latice
!          AA         +       BB               +     CC        = DD
!
!   In cam4, they use constant lv0 (0 C).
!
!    Why am I having non-zero AA when all precip zero.
!
!----------------------------------------------------------------------
AA = 0.
BB_t = 0.
CC_t = 0.
do k = 1,pver
  dmdry_old = (phalf(k)-phalf(k+1))/(g*(1.+rv(k)+rl(k)+ri(k)))
  AA = AA + heat(i,pver+1-k)*dmdry_old
  BB_t = BB_t + (latent(t(k)) + fusion(t(k)))*dmdry_old*(rvnew(k)-rv(k))
  CC_t = CC_t + fusion(t(k))*dmdry_old*(rlnew(k)-rl(k))
end do
BB = (lv0 + lf0)*(col_rv - col_rv_old)/tstep
CC = lf0*(col_rl - col_rl_old)/tstep
DD = -lf0*(uprain_surf_rv + uprain_surf_rlp + anrain_surf)
DD_t = -fusion(t(1))*(uprain_surf_rv + uprain_surf_rlp + anrain_surf)
BB_t = BB_t/tstep
CC_t = CC_t/tstep

if (abs(AA) > 0.0001) then
  rel_err = ((AA + BB + CC) - DD)/AA
else
  rel_err = bad
endif

abs_err = (AA + BB + CC) - DD
abs_err_t = (AA + BB_t + CC_t) - DD

abs_err_diag(i) = abs_err
rel_err_diag(i) = rel_err

if (abs(abs_err) > 100.) then
print*,'---------  CAM4 ENERGY ANALYSIS  ---- i = ',i
print*,'rel_err = ',rel_err
print*,'abs_err = ',abs_err
print*,'AA + BB + CC = ',AA+BB+CC
print*,'AA + BB_t + CC_t = ',AA+BB_t+CC_t
print*,'DD = ',DD
print*,'DD_t = ',DD_t
print*,'dse heating (watts) AA = ',AA
print*,'col rv red (watts) BB = ',BB
print*,'BB_t = ',BB_t
print*,'col rl red (watts) CC = ',CC
print*,'CC_t = ',CC_t
print*,'in mm/day uprain_start_rlp = ',uprain_start_rlp*3600.*24.
print*,'in mm/day uprain_surf_rlp = ',uprain_surf_rlp*3600.*24.
print*,'in mm/day uprain_surf_rv = ',uprain_surf_rv*3600.*24.
print*,'in mm/day anrain_surf = ',anrain_surf*3600.*24.
print*,'in mm/day ansnow_surf = ',ansnow_surf*3600.*24.
print*,'in mm/day ansnow_conv = ',ansnow_conv*3600.*24.
print*,'in mm/day ansnow_strat = ',ansnow_strat*3600.*24.
print*,'column rv change (mm/day) = ',(col_rv-col_rv_old)*48.
print*,'column rl change (mm/day) = ',(col_rl-col_rl_old)*48.
print*,'column ri change (mm/day) = ',(col_ri-col_ri_old)*48.
if (totprecip < 0.001) then
do k = 1,pver
  print*,'rvnew(k) - rv(k) = ',rvnew(k)-rv(k)
end do
endif
print*,'----------------------------------'
endif

!----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!   ************   ALL COLUMNS   ******************************
!
!  Write out stuff if water error exceeds tolerance
!  Numerous differences in the way column water is calculated with cam4:
!  (1) convective scheme calculates dry mass of pressure levels as if
!      cloud liquid/ice from stratiform scheme were zero. This will lead
!      to some errors
!  (2) cam4 appears to calculate column water as if precip does not change
!      the surface pressure
!
!   Estimate of water error from ignoring pressure correction:
!   ---------------------------------------------------------
!   - assume rain = 100 mm/day for 24 hours
!   - change in water mass of column = 100 kg/m2 
!   - suppose rv = 0.02, so absolute error in water in 2 kg/m2
!   - So expect 2 % errors.
!
!   Estimate of enthalpy error from ignoring pressure correction:
!   ---------------------------------------------------------
!   - Make above assumptions and km near surface = 300 kJ/kg
!   - Then for a 100 kg/m2 mass error in a day, would get enthalpy
!     error of 300*100 kJ/m2
!   - from 100 mm/day for 24 hours, you have 100 kg/m2 of rain.
!   - multiplying by lv = 2.5E+06 , you have 2.5E+08 J/m2 heat released.
!   - the enthalpy error relative to the heat released in rain is:
!     300*100/2.5E+05 = 300*100/250000 = 3/25 = 1.5%
!  
!   check_energy.F90:
!   conservation tested on a total energy:
!   se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)  
!   where se(i) is the column dry static energy
!   se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
!   Here pdel/gravit is the mass thickness of a layer.
!   This includes the GPE, since s = cpd*T + g*z
!   In convection, some of the energy released goes into heat, and
!   another part to GPE. I think I should conserve this also ...
!  
!  
!----------------------------------------------------------------------
  nprint_wat = 1
  if (nprint_wat == 1) then
  if (abs(err_wat_corr) > err_wat_corr_max) then
  print*,'== water error exceeds tolerance: WATER BUDGET ANALYSIS in if_conv.F90 ====='
  print*,'column index i = ',i
  print*,'ps corrected (real) water error err_wat_corr = ',err_wat_corr
  print*,'Units likely kg per m2 in one timestep same as mm water'
  print*,'nstep = ',nstep
  print*,'using rvnew and old phalf: mm water col_rv = ',col_rv
  print*,'using rlnew and old phalf: mm water col_rl = ',col_rl
  print*,'using rinew and old phalf: mm water col_ri = ',col_ri
  print*,'using rv and old phalf: mm water col_rv_old = ',col_rv_old
  print*,'using rl and old phalf: mm water col_rl_old = ',col_rl_old
  print*,'using ri and old phalf: mm water col_ri_old = ',col_ri_old
  print*,'change in column rv = ',col_rv - col_rv_old
  print*,'change in column rl = ',col_rl - col_rl_old
  print*,'change in column ri = ',col_rl - col_rl_old
  print*,'loss of column water through surface rain in mm = ',rain_if(i)*1000.*tstep
  print*,'loss of column water through surface snow in mm = ',snow_if(i)*1000.*tstep
  print*,'total precip loss: totprecip (mm) = ',totprecip
  print*,'estimate of surface pressure decrease due to precip (Pa) dp_precip = ',dp_precip
  print*,'estimate of water change due to ps decrease dp_wat =',dp_wat
  print*,'total initial water (mm) col_wat_old = ',col_wat_old
  print*,'total final water (mm) col_wat = ',col_wat
  print*,'Change in column water (mm) expected from precip totprecip = ',totprecip
  print*,'change in column water (mm) col_wat - col_wat_old = ',col_wat - col_wat_old
  print*,'ps fixed column water error in mm err_wat = ',err_wat
  print*,'dp_surf corrected (should be zero) col water error (mm): err_wat - dp_wat = ',err_wat_corr
  print*,'in m/s rain_if(i) = ',rain_if(i)
  print*,'in m/s snow_if(i) = ',snow_if(i)
  print*,'in mm/day surface rain_if(i) = ',rain_if(i)*1000.*3600.*24.
  print*,'in mm/day surface snow_if(i) = ',snow_if(i)*1000.*3600.*24.
  print*,'in mm/day starting uprain_start_rlp = ',uprain_start_rlp*3600.*24.
  print*,'in mm/day uprain_surf_rlp = ',uprain_surf_rlp*3600.*24.
  print*,'in mm/day uprain_surf_rv = ',uprain_surf_rv*3600.*24.
  print*,'in mm/day starting ansnow_conv = ',ansnow_conv*3600.*24.
  print*,'in mm/day starting ansnow_strat = ',ansnow_strat*3600.*24.
  print*,'in mm/day anrain_surf = ',anrain_surf*3600.*24.
  print*,'in mm/day ansnow_surf = ',ansnow_surf*3600.*24.
  print*,'g = ',g
  print*,'========================================================'
  do k = 1,pver
!   print*,'k rvnew(k) rv(k) rvnew(k)-rv(k) = ',k,rvnew(k),rv(k),rvnew(k)-rv(k)
  end do
  if (abs(err_wat_corr) > 3.) stop 'column error too big'
! stop 'column water error too big in if_conv.F90'
  endif
  endif
!----------------------------------------------------------------------
!
!   ************  IN LOOP OVER COLUMNS   ******************************
!
!  Ouput stuff
!
!----------------------------------------------------------------------
! do k = 1,pver
  do k = 1,1
!   print*,'k qlnew(k) = ',k,qlnew(k)
  end do
! do k = 1,pver
  do k = 1,1
!   print*,'k (phalf(k) - phalf(k+1)) = ',k,(phalf(k) - phalf(k+1))
  end do
!  do k = 1,pver
  do k = 1,1
!   print*,'k qvnew(k) = ',k,qvnew(k)
  end do
!  print*,'========================================================'
!  print*,'qtnd(1,pver) = ',qtnd(1,pver)
!  print*,'qv(1) = ',qv(1)
!  print*,'qvnew(1) = ',qvnew(1)
!  print*,'tstep = ',tstep
!  print*,'========================================================'

   rrain_if(i) = rain_if(i)
!   print*,'in if_conv.F90 rain_if(i) = ',rain_if(i)

!----------------------------------------------------------------------
!
!  End loop over columns (i=1,ncol)
!
!---------------------------------------------------------------------
   end do

!---------------------------------------------------------------------
!
!   ************   AFTER LOOP OVER COLUMNS   **************************
!
! Desireable to have outfld here so don't have to pass arrays
! up to higher subroutines. 
!
! outfld is in: 
! /home/folkins/cesm1_0/models/atm/cam/src/control/cam_history.F90
!
!
! subroutine outfld (fname, field, idim, c)
!
! Arguments
!
! character(len=*), intent(in) :: fname ! Field name--should be 8 chars long
! integer, intent(in) :: idim           ! Longitude dimension of field array
! integer, intent(in) :: c              ! chunk (physics) or latitude (dynamics) index
! real(r8), intent(in) :: field(idim,*) ! Array containing field values
!
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!
!   ************   AFTER LOOP OVER COLUMNS   **************************
!
!  0. No corresponding diag array
!
!---------------------------------------------------------------------

  call outfld('PRECIF  ', prec_if   , pcols, lchnk)

!---------------------------------------------------------------------
!
!  1. 3D Rain diagnostics
!
!---------------------------------------------------------------------

  call outfld('UPRAIN_SURF_RLP', uprain_surf_rlp_diag   , pcols, lchnk)
  call outfld('UPRAIN_SURF_RV', uprain_surf_rv_diag , pcols, lchnk)
  call outfld('UPRAIN_START_RLP', uprain_start_rlp_diag , pcols, lchnk)
  call outfld('UPRAIN_RLP_RVV', uprain_rlp_rvv_diag   , pcols, lchnk)
  call outfld('UPRAIN_EVAP', uprain_evap_diag   , pcols, lchnk)
  call outfld('ANRAIN_DOWN', anrain_down_diag   , pcols, lchnk)
  call outfld('ANRAIN_EVAP', anrain_evap_diag   , pcols, lchnk)
  call outfld('ANSNOW_SUBL', ansnow_subl_diag   , pcols, lchnk)
  call outfld('ANSNOW_MELT', ansnow_melt_diag   , pcols, lchnk)
  call outfld('UPRAIN_DOWN', uprain_down_diag   , pcols, lchnk)
  call outfld('UPRAIN_1', uprain_1_diag   , pcols, lchnk)
  call outfld('UPRAIN_2', uprain_2_diag   , pcols, lchnk)

  call outfld('ANSNOW_CONV  ', ansnow_conv_diag   , pcols, lchnk)
  call outfld('ANSNOW_STRAT ', ansnow_strat_diag   , pcols, lchnk)
  call outfld('ANSNOW_STRAT_RV ', ansnow_strat_rv_diag   , pcols, lchnk)
  call outfld('ANSNOW_STRAT_RI ', ansnow_strat_ri_diag   , pcols, lchnk)
  call outfld('ANSNOW_SURF', ansnow_surf_diag   , pcols, lchnk)
  call outfld('ANRAIN_SURF', anrain_surf_diag   , pcols, lchnk)

!---------------------------------------------------------------------
!
!  2. 3D Column Water
!
!---------------------------------------------------------------------

  call outfld('COLWAT', colwat_diag  , pcols, lchnk)
  call outfld('COLRH', colrh_diag  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  3. 3D Cape
!
!---------------------------------------------------------------------

  call outfld('CONV_ENERGY', conv_energy_diag  , pcols, lchnk)
  call outfld('CF_TAR', cf_tar_diag  , pcols, lchnk)
  call outfld('LWP_TAR', lwp_tar_diag  , pcols, lchnk)
  call outfld('CAPE_MEAN', cape_mean_diag  , pcols, lchnk)

  call outfld('AMP', amp_diag  , pcols, lchnk)
  call outfld('AMP_OM5', amp_om5_diag  , pcols, lchnk)
  call outfld('AMP_RAIN', amp_rain_diag  , pcols, lchnk)

  call outfld('DP_DN', dp_dn_diag  , pcols, lchnk)

  call outfld('F_1', f_1_diag  , pcols, lchnk)
  call outfld('F_2', f_2_diag  , pcols, lchnk)
  call outfld('F_3', f_3_diag  , pcols, lchnk)
  call outfld('F_4', f_4_diag  , pcols, lchnk)

  call outfld('LTS', LTS_diag  , pcols, lchnk)
  call outfld('EIS', EIS_diag  , pcols, lchnk)
  call outfld('MSE_RAT', mse_rat_diag  , pcols, lchnk)

  call outfld('CAPE_1', cape_1_diag  , pcols, lchnk)
  call outfld('CAPE_2', cape_2_diag  , pcols, lchnk)
  call outfld('CAPE_3', cape_3_diag  , pcols, lchnk)
  call outfld('CAPE_4', cape_4_diag  , pcols, lchnk)

  call outfld('CFRACTION_1', cfraction_1_diag  , pcols, lchnk)
  call outfld('CFRACTION_2', cfraction_2_diag  , pcols, lchnk)
  call outfld('CFRACTION_3', cfraction_3_diag  , pcols, lchnk)
  call outfld('CFRACTION_4', cfraction_4_diag  , pcols, lchnk)

  call outfld('HMUPRAIN_START', hmuprain_start_diag  , pcols, lchnk)
  call outfld('HMUPRAIN_SURF', hmuprain_surf_diag  , pcols, lchnk)
  call outfld('HMANSNOW_START', hmansnow_start_diag  , pcols, lchnk)
  call outfld('HMANSNOW_SURF', hmansnow_surf_diag  , pcols, lchnk)
  call outfld('HMANRAIN_SURF', hmanrain_surf_diag  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  5. 3D start
!
!---------------------------------------------------------------------

  call outfld('MASS_START_1', mass_start_1_diag  , pcols, lchnk)
  call outfld('MASS_START_2', mass_start_2_diag  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  7. 3D Nonlinearities
!
!---------------------------------------------------------------------

  call outfld('OM5', om5_diag  , pcols, lchnk)

  call outfld('KM_WIDTH', km_width_diag  , pcols, lchnk)

  call outfld('MD_RAT', md_rat_diag , pcols, lchnk)
  call outfld('F_AN', f_an_diag , pcols, lchnk)
  call outfld('F_ICE', f_ice_diag , pcols, lchnk)

  call outfld('ZPEAK', zpeak_diag  , pcols, lchnk)
  call outfld('AV_Z_RI', av_z_ri_diag  , pcols, lchnk)
  call outfld('AV_Z_RI_RH', av_z_ri_rh_diag  , pcols, lchnk)
  call outfld('AV_Z_RI_DT', av_z_ri_dt_diag  , pcols, lchnk)
  call outfld('PRECIP_ORG', precip_org_diag  , pcols, lchnk)

  call outfld('MASS_INC', mass_inc_diag  , pcols, lchnk)
  call outfld('MASS_RAT', mass_rat_diag  , pcols, lchnk)

  call outfld('ABS_ERR', abs_err_diag  , pcols, lchnk)
  call outfld('REL_ERR', rel_err_diag  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  8. 3D Efficiencies
!
!---------------------------------------------------------------------

  call outfld('UP_EFF', updraft_eff_eff  , pcols, lchnk)
  call outfld('UP_RLP', updraft_eff_rlp  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  9. 4DF tend
!
!---------------------------------------------------------------------

  call outfld('TPTEND_REAL', tptend_real  , pcols, lchnk)
  call outfld('RHTEND_REAL', rhtend_real  , pcols, lchnk)
  call outfld('QTEND_REAL', qtend_real  , pcols, lchnk)
  call outfld('KMTEND  ', kmtend_save  , pcols, lchnk)
  call outfld('RVTEND  ', rvtend_save  , pcols, lchnk)
  call outfld('RLTEND  ', rltend_save  , pcols, lchnk)
  call outfld('RLNEW  ', rlnew_save  , pcols, lchnk)
  call outfld('RITEND  ', ritend_save  , pcols, lchnk)
  call outfld('UTENDIF', utend_save  , pcols, lchnk)
  call outfld('VTENDIF', vtend_save  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  10. 4DI
!
!---------------------------------------------------------------------
  call outfld('FMASS ', fmass_diag  , pcols, lchnk)
  call outfld('FMASS_DN', fmass_dn_diag  , pcols, lchnk)
  call outfld('PHALF', phalf_diag  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  11. 4DF 
!
!---------------------------------------------------------------------

  call outfld('ENT', ent_diag  , pcols, lchnk)
  call outfld('DET', det_diag  , pcols, lchnk)

  call outfld('UPDET_1', updet_1_diag  , pcols, lchnk)
  call outfld('UPDET_2', updet_2_diag  , pcols, lchnk)

  call outfld('RH_DN', rh_dn_diag  , pcols, lchnk)
  call outfld('DRV_DNDET', drv_dndet_diag  , pcols, lchnk)
  call outfld('DT_DNDET', dt_dndet_diag  , pcols, lchnk)
  call outfld('DNDET', dndet_diag  , pcols, lchnk)
  call outfld('DNENT', dnent_diag  , pcols, lchnk)
  call outfld('MTDET', mtdet_diag  , pcols, lchnk)
  call outfld('MTENT', mtent_diag  , pcols, lchnk)
  call outfld('DRV_UP', drv_up_diag  , pcols, lchnk)
  call outfld('DRV_AN', drv_an_diag  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  12. 4DF 
!
!---------------------------------------------------------------------

  call outfld('RHIF', rhif_store  , pcols, lchnk)
  call outfld('RHNEW', rhnew , pcols, lchnk)
  call outfld('TIF', tif_store  , pcols, lchnk)
  call outfld('HM', hm_store  , pcols, lchnk)
  call outfld('PROB_RANDOM', prob_random_store  , pcols, lchnk)
  call outfld('RI_HALF', ri_half_store  , pcols, lchnk)
  call outfld('RL_HALF', rl_half_store  , pcols, lchnk)

!---------------------------------------------------------------------
!
!  13. 4DF For updraft properties
!
!---------------------------------------------------------------------

  call outfld('MP', mp_save  , pcols, lchnk)
  call outfld('BP_CAPE', bp_cape_save  , pcols, lchnk)   
  call outfld('BP1_STR', bp1_str_save  , pcols, lchnk)   
  call outfld('BP1_PRP', bp1_prp_save  , pcols, lchnk) 
  call outfld('BP1_ITR', bp1_itr_save  , pcols, lchnk) 
  call outfld('BP2_STR', bp2_str_save  , pcols, lchnk)   
  call outfld('BP2_PRP', bp2_prp_save  , pcols, lchnk) 
  call outfld('BP2_ITR', bp2_itr_save  , pcols, lchnk) 
  call outfld('F_ENT_1', f_ent_1_save  , pcols, lchnk) 
  call outfld('F_DET_1', f_det_1_save  , pcols, lchnk) 
  call outfld('F_ENT_2', f_ent_2_save  , pcols, lchnk) 
  call outfld('F_DET_2', f_det_2_save  , pcols, lchnk) 
  call outfld('RLP', rlp_save  , pcols, lchnk)  
  call outfld('RIP', rip_save  , pcols, lchnk)  
  call outfld('HMP_MEAN', hmp_mean_save  , pcols, lchnk)  
  call outfld('HM_DIFF', hm_diff_save  , pcols, lchnk)  

!---------------------------------------------------------------------
!
!  14. 4DF For detrain properties
!
!---------------------------------------------------------------------

  call outfld('TDIFF_DETRAIN_UP', tdiff_detrain_up_save  , pcols, lchnk)   
  call outfld('CAPE_DETRAIN_UP', cape_detrain_up_save  , pcols, lchnk)   
  call outfld('MASS_DETRAIN_UP', mass_detrain_up_save  , pcols, lchnk)   
  call outfld('HM_DIFF_DET', hm_diff_det_save  , pcols, lchnk)   

!---------------------------------------------------------------------
!
!  15. 4DF Evaporation
!
!---------------------------------------------------------------------

  call outfld('UPRAINEVAP', uprainevap_save  , pcols, lchnk)   
  call outfld('ANRAINEVAP', anrainevap_save  , pcols, lchnk)   
  call outfld('ANSNOWSUBL', ansnowsubl_save  , pcols, lchnk)   
  call outfld('ANSNOWMELT', ansnowmelt_save  , pcols, lchnk)   

  call outfld('UPRAINDOWN', upraindown_save  , pcols, lchnk)   
  call outfld('ANRAINDOWN', anraindown_save  , pcols, lchnk)   

!---------------------------------------------------------------------
!
!  16. 4DF Precipitation Sources
!
!---------------------------------------------------------------------

  call outfld('UPRAIN_RLP', uprain_rlp_save  , pcols, lchnk)   
  call outfld('ANSNOW_RLP', ansnow_rlp_save  , pcols, lchnk)   
  call outfld('UPRAIN_RV', uprain_rv_save  , pcols, lchnk)   
  call outfld('ANSNOW_RV', ansnow_rv_save  , pcols, lchnk)   

  call outfld('RLP_RVV', rlp_rvv_save  , pcols, lchnk)   

!---------------------------------------------------------------------
!
!  17. 4DF Cloud production/Loss
!
!---------------------------------------------------------------------

  call outfld('UPRAIN_RL_TEND', uprain_rl_tend_save  , pcols, lchnk)   
  call outfld('RL_EVAP_TEND', rl_evap_tend_save  , pcols, lchnk)   
  call outfld('RI_EVAP_TEND', ri_evap_tend_save  , pcols, lchnk)   
  call outfld('RI_PROD_TEND', ri_prod_tend_save  , pcols, lchnk)   
  call outfld('RL_DET_TEND', rl_det_tend_save  , pcols, lchnk)   
  call outfld('RI_DET_TEND', ri_det_tend_save  , pcols, lchnk)   
  call outfld('RL_VERT_TEND', rl_vert_tend_save  , pcols, lchnk)   
  call outfld('RI_VERT_TEND', ri_vert_tend_save  , pcols, lchnk)   

  call outfld('RL_DT', rl_dt_save  , pcols, lchnk)   
  call outfld('RL_RH', rl_rh_save  , pcols, lchnk)   
  call outfld('RI_DT', ri_dt_save  , pcols, lchnk)   
  call outfld('RI_RH', ri_rh_save  , pcols, lchnk)   
  call outfld('RI_RS', ri_rs_save  , pcols, lchnk)   

  call outfld('COL_RI_DT', col_ri_dt_save  , pcols, lchnk)   
  call outfld('COL_RI_RH', col_ri_rh_save  , pcols, lchnk)   
  call outfld('COL_RL_DT', col_rl_dt_save  , pcols, lchnk)   
  call outfld('COL_RL_RH', col_rl_rh_save  , pcols, lchnk)   

!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
!  print*,'leaving if_convr'
   return
!---------------------------------------------------------------------
!
!  end subroutine if_convr
!
!---------------------------------------------------------------------
end subroutine if_convr
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
end module if_conv
