module if_conv_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the IF deep convection scheme
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!
!  Use statements
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair                              
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use if_conv,      only: if_convr
   use cam_history,  only: outfld, addfld, add_default, phys_decomp
   use perf_mod
   use cam_logfile,  only: iulog

   implicit none

   save
! IF added ... I removed this so physics_update.F90 could see physics_for_update
! May screw things up.
   private                         ! Make default type private to the module
!-------------------------------------------------------------------------------
!
! Public methods
!  - renamed these -> "if"
!
!-------------------------------------------------------------------------------
   public ::&
      if_conv_register,           &! seems to do very little
      if_conv_init,               &! basically now just addfld commands
      if_conv_tend                 ! return tendencies + output some fields
!-------------------------------------------------------------------------------
!
!  Private module data: allocatable commands
!
!  I removed all of the allocatable commands since they were all with gathered
!    quantities. retained jt below just as an example of how allocatable is used.
!  - mu,eu,du,md,ed,dp,dsubcld,maxg,ideep,lengath
!
!  What does the allocatable command do?
!
!  from : http://www.star.le.ac.uk/~cgp/f90course/f90.html#tth_sEc5.2
!
!  Allocatable arrays are more generally useful as their size may be set at any 
!  point. Only the rank has to be declared in advance, with a colon marking the 
!  each dimension:
!
!  REAL, ALLOCATABLE :: vector(:), matrix(:,:), three_d(:,:,:)
!
!  Later on, one sets the actual dimension bounds may then be set anywhere in the 
!  executable code (the lower bound is 1 by default).
!
!  ALLOCATE(vector(12345), matrix(0:511,0:255))
!
!  This means that most of the corresponding allocate commands in subroutine
!    subroutine if_conv_init(hypi) will be invalid, presumably, so should remove 
!    them also.
!
!  Interesting:
!
!  in zm_conv, jt is specified:    integer jt(pcols)
!  - i.e. of different rank
!
!  jt is passed to zm_conv as "jt(:,lchnk)", i.e. only one chunk at a time,
!   and it is this single chunk that is called jt at the lower level; presumably
!   completely different arrays, and no ambiguity created by giving them the
!   same name.
!
!  I guess the number of chunks is determined by the number of processors, and
!    the "begchunk:endchunk" range is only known lower down where the allocate
!    command is used ..... but actually know since "use ppgrid" above so why
!    bother with allocatable?
!
!   Had for jt:
!
!   integer, allocatable, dimension(:,:) :: jt   !(pcols,begchunk:endchunk)
!    wg top  level index of deep cumulus convection.
!
!-------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
! Purpose: 
!
!  IF (Feb 2013): The old purpose was to register ANVLIQ. However I have 
!     dropped this, so probably does nothing. I have retained in case 
!     I want to add a tracer in the future.
!
!  called by "subroutine convect_deep_register", which is called from
!    "subroutine initindx", in early stages of model.
!
!  - register fields with the physics buffer
!  - by these use commands it makes pbuf_times, pbuf_add available to
!    other subroutines in the module?
!  - what does idx do?
!
!  - The zm version of this subroutine was extremely simple.
!  - But seemed to me that this was the appropriate place to register
!    anvil liquid and anvil age.
!  - Using "subroutine stratiform_register" in "module stratiform" as my
!    template
!
!
!------------------------------------------------------------------------
subroutine if_conv_register

  use phys_buffer, only: pbuf_times, pbuf_add
!------------------------------------------------------------------------
!
!  Added this use statement on analogy with stratiform_register
!
!------------------------------------------------------------------------
  use constituents, only: cnst_add, pcnst
  use physconst,    only: mwdry, cpair, mwh2o, cpliq

  implicit none

  integer idx

!------------------------------------------------------------------------
!
! Register anvil liquid water and age, and determine index.
!
!  Two curious things:
!
!  (1) inconsistent number of arguments for "subroutine cnst_add" (SOLVED)
!  (2) mwdry and cpair used in "cnst_add" call for 'CLDLIQ', 'CLDICE' in
!      "subroutine stratiform_register"
!
!  ------------------------------------
!
! (1) In "module constituents", "cnst_add" has 10 arguments:
!
! subroutine cnst_add (name, mwc, cpc, qminc, &
!                      ind, longname, readiv, mixtype, cam_outfld, fixed_ubc)
!
! However, when this subroutine is called, usually fewer arguments are given.
! Some of the arguments are lableled optional in subroutine cnst_add
! for optional arguments, appear to need to explicitly specify variable name.
!
! Slots:
! (1) name (in)
! (2) constituent molecular wt (in)
! (3) constituent specific heat (in)
! (4) minimum value of mixing ratio (in) (Set both to zero)
! (5) global consituent index (out)
! (6) OPTIONAL: long name (in)
! (7) OPTIONAL: readiv (Put as false: do NOT read in from file)
!     true => read initial values from initial file (default: true)
! (8) OPTIONAL: mixtype='dry' (ignore, set myself)
! (9) OPTIONAL: cam_outfld:  (ignore)
! (10) OPTIONAL: fixed_ubc: upper boundary condition (ignore)
!
!  For anvil age: used a molecular weight and specific heat of zero.
!
!  ------------------------------------
!
!  (2) mwh20 and cpair
!
! Presumably the main reason to care about how these in my case is
!   for checking energy conservation. 
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
! Do I have to request buffer space for q? 
! No, I don't think so: this should be handled by a general request
!  for the tracer, or state type. Done in subroutine physics_type_alloc
!  of module physics_types.
! Where is the buffer space for q asked for?
!
! Request physics buffer space for fields that persist across timesteps.
! (removed some commands here)
!
!------------------------------------------------------------------------
end subroutine if_conv_register
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  Called from "subroutine convect_deep_init".
!    
!
!------------------------------------------------------------------------
subroutine if_conv_init(hypi)
!----------------------------------------
!
! Purpose:  declare output fields, initialize variables needed by convection
!
! removed since don't use analagous program:
!   use if_conv,        only: if_convi   
!
!----------------------------------------
  use cam_history,    only: outfld, addfld, add_default, phys_decomp
  use phys_buffer,    only: pbuf_times, pbuf_add
  use ppgrid,         only: pcols, pver
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use error_messages, only: alloc_err	
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is

  implicit none
!------------------------------------------------------------------------
!
! used in definition of limcnv but could probably be removed.
!
!------------------------------------------------------------------------
  real(r8),intent(in) :: hypi(plevp)        ! reference pressures at interfaces

  integer k, istat
!------------------------------------------------------------------------
!  added in cam5
!------------------------------------------------------------------------
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.

!print*,'starting subroutine if_conv_init'
!------------------------------------------------------------------------
!
! Allocate space for arrays private to this module
!
!  Removed all of these allocate statements; appeared to be all 
!    gathered quantities.
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!  addfld commands I have removed:
!
!------------------------------------------------------------------------
!
! Register fields with the output buffer
!
!  ZM addfld calls:
!
!   call addfld ('PRECZ   ','m/s     ',1,    'A','total precipitation from ZM convection',        phys_decomp)            
!   call addfld ('ZMDT    ','K/s     ',pver, 'A','T tendency - Zhang-McFarlane moist convection', phys_decomp)         
!   call addfld ('ZMDQ    ','kg/kg/s ',pver, 'A','Q tendency - Zhang-McFarlane moist convection', phys_decomp)         
!   call addfld ('ZMDICE ','kg/kg/s ',pver, 'A','Cloud ice tendency - Zhang-McFarlane convection',phys_decomp)            
!   call addfld ('ZMDLIQ ','kg/kg/s ',pver, 'A','Cloud liq tendency - Zhang-McFarlane convection',phys_decomp)                 
!   call addfld ('EVAPTZM ','K/s     ',pver, 'A','T tendency - Evaporation/snow prod from Zhang convection',phys_decomp)          
!   call addfld ('FZSNTZM ','K/s     ',pver, 'A','T tendency - Rain to snow conversion from Zhang convection',phys_decomp) 
!   call addfld ('EVSNTZM ','K/s     ',pver, 'A','T tendency - Snow to rain prod from Zhang convection',phys_decomp)          
!   call addfld ('EVAPQZM ','kg/kg/s ',pver, 'A','Q tendency - Evaporation from Zhang-McFarlane moist convection',phys_decomp) 
!   call addfld ('ZMFLXPRC','kg/m2/s ',pverp, 'A','Flux of precipitation from ZM convection'       ,phys_decomp)  
!   call addfld ('ZMFLXSNW','kg/m2/s ',pverp, 'A','Flux of snow from ZM convection'                ,phys_decomp)   
!   call addfld ('ZMNTPRPD','kg/kg/s ',pver , 'A','Net precipitation production from ZM convection',phys_decomp) 
!   call addfld ('ZMNTSNPD','kg/kg/s ',pver , 'A','Net snow production from ZM convection'         ,phys_decomp) 
!   call addfld ('ZMEIHEAT','W/kg'    ,pver , 'A','Heating by ice and evaporation in ZM convection',phys_decomp)
!   call addfld ('CMFMCDZM','kg/m2/s ',pverp,'A','Convection mass flux from ZM deep ',phys_decomp)                           
!   call addfld ('PRECCDZM','m/s     ',1,    'A','Convective precipitation rate from ZM deep',phys_decomp)                      
!   call add_default ('CMFMCDZM', 1, ' ')
!   call add_default ('PRECCDZM', 1, ' ')
!   call addfld ('PCONVB','Pa'    ,1 , 'A','convection base pressure',phys_decomp)                              
!   call addfld ('PCONVT','Pa'    ,1 , 'A','convection top pressure',phys_decomp)                        
!   call add_default ('PCONVB', 1, ' ')
!   call add_default ('PCONVT', 1, ' ')
!   call addfld ('CAPE',   'J/kg',       1, 'A', 'Convectively available potential energy', phys_decomp)  
!   call addfld ('FREQZM ','fraction  ',1  ,'A', 'Fractional occurance of ZM convection',phys_decomp)           
!   call add_default ('FREQZM', 1, ' ')
!   call addfld ('ZMMTT ', 'K/s',     pver, 'A', 'T tendency - ZM convective momentum transport',phys_decomp)            
!   call addfld ('ZMMTU',  'm/s2',    pver, 'A', 'U tendency - ZM convective momentum transport',  phys_decomp)            
!   call addfld ('ZMMTV',  'm/s2',    pver, 'A', 'V tendency - ZM convective momentum transport',  phys_decomp)  
!   call addfld ('ZMMU',   'kg/m2/s', pver, 'A', 'ZM convection updraft mass flux',   phys_decomp)                 
!   call addfld ('ZMMD',   'kg/m2/s', pver, 'A', 'ZM convection downdraft mass flux', phys_decomp)                 
!   call addfld ('ZMUPGU', 'm/s2',    pver, 'A', 'zonal force from ZM updraft pressure gradient term',       phys_decomp) 
!   call addfld ('ZMUPGD', 'm/s2',    pver, 'A', 'zonal force from ZM downdraft pressure gradient term',     phys_decomp)  
!   call addfld ('ZMVPGU', 'm/s2',    pver, 'A', 'meridional force from ZM updraft pressure gradient term',  phys_decomp)  
!   call addfld ('ZMVPGD', 'm/s2',    pver, 'A', 'merdional force from ZM downdraft pressure gradient term', phys_decomp)   
!   call addfld ('ZMICUU', 'm/s',     pver, 'A', 'ZM in-cloud U updrafts', phys_decomp)         
!   call addfld ('ZMICUD', 'm/s',     pver, 'A', 'ZM in-cloud U downdrafts', phys_decomp)         
!   call addfld ('ZMICVU', 'm/s',     pver, 'A', 'ZM in-cloud V updrafts', phys_decomp)       
!   call addfld ('ZMICVD', 'm/s',     pver, 'A', 'ZM in-cloud V downdrafts', phys_decomp)                      
!   call addfld ('ZMMTT ', 'K/s',     pver, 'A', 'T tendency - ZM convective momentum transport',phys_decomp)  
!   call addfld ('ZMMTU',  'm/s2',    pver, 'A', 'U tendency - ZM convective momentum transport',  phys_decomp)                  
!   call addfld ('ZMMTV',  'm/s2',    pver, 'A', 'V tendency - ZM convective momentum transport',  phys_decomp)                   
!   call addfld ('ZMMU',   'kg/m2/s', pver, 'A', 'ZM convection updraft mass flux',   phys_decomp)           
!   call addfld ('ZMMD',   'kg/m2/s', pver, 'A', 'ZM convection downdraft mass flux', phys_decomp)            
!   call addfld ('ZMUPGU', 'm/s2',    pver, 'A', 'zonal force from ZM updraft pressure gradient term',       phys_decomp)  
!   call addfld ('ZMUPGD', 'm/s2',    pver, 'A', 'zonal force from ZM downdraft pressure gradient term',     phys_decomp) 
!   call addfld ('ZMVPGU', 'm/s2',    pver, 'A', 'meridional force from ZM updraft pressure gradient term',  phys_decomp) 
!   call addfld ('ZMVPGD', 'm/s2',    pver, 'A', 'merdional force from ZM downdraft pressure gradient term', phys_decomp) 
!   call addfld ('ZMICUU', 'm/s',     pver, 'A', 'ZM in-cloud U updrafts', phys_decomp)           
!   call addfld ('ZMICUD', 'm/s',     pver, 'A', 'ZM in-cloud U downdrafts', phys_decomp)           
!   call addfld ('ZMICVU', 'm/s',     pver, 'A', 'ZM in-cloud V updrafts', phys_decomp)          
!   call addfld ('ZMICVD', 'm/s',     pver, 'A', 'ZM in-cloud V downdrafts', phys_decomp)             
!
!     
!------------------------------------------------------------------------
! 
!  Register fields with the output buffer
!
!   - In my case need addfld to add variables to the History Master List
!   - The zm code makes addfld calls to variables that are already in the
!     History Master list:
!   - e.g. in /home/folkins/cesm1_0/models/atm/cam/src/physics/cam/zm_conv_intr.F90
!     there is an addfld call for ZMMTT.
!   - But this is already in:
!  http://www.cesm.ucar.edu/models/cesm1.0/cam/docs/users_guide/hist_flds_fv_cam4.html
!   - don't know whether addfld does other things, or called just in case.
!
!  - I have removed many calls to addfld here - so will miss many  diagnostics(??)                
!  - What is the consequence of removing these calls to addfield?
!  - add IF snow?; change names
!  - subroutine "addfld" is in:
!    /home/folkins/cesm1_0/models/atm/cam/src/control/cam_history.F90
! 
!  There are 9 inout arguments, 3 optional:
!
!  subroutine addfld (fname, units, numlev, avgflag, long_name, &
!                     decomp_type, flag_xyfill, flag_isccplev, sampling_seq)
!
!  (1) character(len=*), intent(in) :: fname      ! field name--should be "max_fieldname_len" 
!                                                 ! characters long or less             
!      - Field name: 8-character field name, left-justified, alphanumeric or
!      spaces only
!  (2) character(len=*), intent(in) :: units      ! units of fname--should be 8 chars  
!      - Field units: 8-character units description. See Table 3.1
!  (3) character(len=1), intent(in) :: avgflag    ! averaging flag
!  (4) character(len=*), intent(in) :: long_name  ! long name of field
!  (5) integer, intent(in) :: numlev              ! number of vertical levels
!                                                 ! (dimension and loop)
!      - Number of vertical levels in the field.
!  (6) integer, intent(in) :: decomp_type         ! decomposition type
!      - Parallel decomposition type (i.e. is this a physics or dynamics
!      variable) 
!  (7) logical, intent(in), optional :: flag_xyfill ! non-applicable xy points
!                                                   ! flagged with fillvalue
!  (8) logical, intent(in), optional :: flag_isccplev ! levels are ISCCP levels
!                                                     !  not vertical
!  (9) character(len=*), intent(in), optional :: sampling_seq ! sampling sequence
!                                                ! if not every timestep,
!                                                ! how often field is sampled:
!                                                ! every other; only during LW/SW         
!                                                ! radiation calcs, etc
!
!  Strange: why is the order above different from the order given in the
!  calls to addfld below?
!
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  Are output locally in this program
!  - But these appear to be OBS: I do not output them to History files
!
!------------------------------------------------------------------------
call addfld ('IFDT    ','K/s     ',pver, 'A','T tendency - IF moist convection', phys_decomp)
call addfld ('IFDQ    ','kg/kg/s ',pver, 'A','Q tendency - IF moist convection', phys_decomp)
call addfld ('IFMTU',  'm/s2',    pver, 'A', 'U tendency - IF convective momentum transport',  phys_decomp)
call addfld ('IFMTV',  'm/s2',    pver, 'A', 'V tendency - IF convective momentum transport',  phys_decomp)
call add_default('IFDT     ', 1, ' ')
call add_default('IFDQ     ', 1, ' ')
call add_default('IFMTU    ', 1, ' ')
call add_default('IFMTV    ', 1, ' ')
!------------------------------------------------------------------------
!
!  Have copied this from ICLMRCU in:
!  http://www.cesm.ucar.edu/models/cesm1.0/cam/docs/users_guide/hist_flds_fv_cam4.html
!
!  Need to call add_default for ANVLIQ
!  - I guess you use this if you want to output it, but program shouldn't
!    crash if not there.
!
!     subroutine add_default (name, tindex, flag)
!  Purpose: Add a field to the default "on" list for a given history file
!
!  Also in:
!  - subroutine "addfld" is in:
!    /home/folkins/cesm1_0/models/atm/cam/src/control/cam_history.F90
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  0. 3D Rain variables
!
!------------------------------------------------------------------------
call addfld ('PRECIF  ','m/s     ',1,    'A','total precipitation from IF convection', phys_decomp)
call add_default('PRECIF   ', 1, ' ')

!------------------------------------------------------------------------
!
!  1. 3D Rain diagnostics 
!
!------------------------------------------------------------------------
call addfld ('UPRAIN_SURF_RLP','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('UPRAIN_SURF_RV','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('UPRAIN_START_RLP','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('UPRAIN_RLP_RVV','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('UPRAIN_EVAP','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('ANRAIN_DOWN','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('ANRAIN_EVAP','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('ANSNOW_SUBL','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('ANSNOW_MELT','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('UPRAIN_DOWN','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('UPRAIN_1','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call addfld ('UPRAIN_2','m/s     ',1, 'A','updraft rain from IF convection', phys_decomp)
call add_default('UPRAIN_SURF_RLP', 1, ' ')
call add_default('UPRAIN_SURF_RV', 1, ' ')
call add_default('UPRAIN_START_RLP', 1, ' ')
call add_default('UPRAIN_RLP_RVV', 1, ' ')
call add_default('UPRAIN_EVAP', 1, ' ')
call add_default('ANRAIN_DOWN', 1, ' ')
call add_default('ANRAIN_EVAP', 1, ' ')
call add_default('ANSNOW_SUBL', 1, ' ')
call add_default('ANSNOW_MELT', 1, ' ')
call add_default('UPRAIN_DOWN', 1, ' ')
call add_default('UPRAIN_1', 1, ' ')
call add_default('UPRAIN_2', 1, ' ')

call addfld ('ANSNOW_CONV','m/s     ',1, 'A','anvil snow from IF convection', phys_decomp)
call addfld ('ANSNOW_STRAT','m/s     ',1, 'A','anvil snow from IF convection', phys_decomp)
call addfld ('ANSNOW_STRAT_RV','m/s     ',1, 'A','anvil snow from IF convection', phys_decomp)
call addfld ('ANSNOW_STRAT_RI','m/s     ',1, 'A','anvil snow from IF convection', phys_decomp)
call addfld ('ANSNOW_SURF','m/s     ',1, 'A','anvil snow from IF convection', phys_decomp)
call addfld ('ANRAIN_SURF','m/s     ',1, 'A','anvil snow from IF convection', phys_decomp)
call add_default('ANSNOW_CONV', 1, ' ')
call add_default('ANSNOW_STRAT', 1, ' ')
call add_default('ANSNOW_STRAT_RV', 1, ' ')
call add_default('ANSNOW_STRAT_RI', 1, ' ')
call add_default('ANSNOW_SURF', 1, ' ')
call add_default('ANRAIN_SURF', 1, ' ')

call addfld ('ABS_ERR','m/s     ',1, 'A','    ', phys_decomp)
call addfld ('REL_ERR','m/s     ',1, 'A','    ', phys_decomp)
call add_default('ABS_ERR', 1, ' ')
call add_default('REL_ERR', 1, ' ')

!------------------------------------------------------------------------
!
!  2. 3D Column Water
!
!------------------------------------------------------------------------
call addfld ('COLWAT ','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('COLRH','kg water     ',1,    'A','column water',        phys_decomp)
call add_default('COLWAT', 1, ' ')
call add_default('COLRH', 1, ' ')

!------------------------------------------------------------------------
!
!  3. 3D Cape
!
!------------------------------------------------------------------------
call addfld ('CF_TAR','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('LWP_TAR','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('CONV_ENERGY','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('CAPE_MEAN','J/kg     ',1,    'A','IF density ',        phys_decomp)
call add_default('CF_TAR', 1, ' ')
call add_default('LWP_TAR', 1, ' ')
call add_default('CONV_ENERGY', 1, ' ')
call add_default('CAPE_MEAN', 1, ' ')

call addfld ('AMP','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('AMP_RAIN','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('AMP_OM5','kg water     ',1,    'A','column water',        phys_decomp)
call add_default('AMP', 1, ' ')
call add_default('AMP_RAIN', 1, ' ')
call add_default('AMP_OM5', 1, ' ')

call addfld ('DP_DN','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('F_1','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('F_2','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('F_3','kg water     ',1,    'A','column water',        phys_decomp)
call addfld ('F_4','kg water     ',1,    'A','column water',        phys_decomp)
call add_default('DP_DN', 1, ' ')
call add_default('F_1', 1, ' ')
call add_default('F_2', 1, ' ')
call add_default('F_3', 1, ' ')
call add_default('F_4', 1, ' ')

call addfld ('LTS','J/kg     ',1,    'A','IF density ',        phys_decomp)
call addfld ('EIS','J/kg     ',1,    'A','IF density ',        phys_decomp)
call addfld ('MSE_RAT','J/kg     ',1,    'A','IF density ',        phys_decomp)
call add_default('LTS', 1, ' ')
call add_default('EIS', 1, ' ')
call add_default('MSE_RAT', 1, ' ')

call addfld ('CAPE_1','J/kg     ',1,    'A','IF density ',        phys_decomp)
call addfld ('CAPE_2','J/kg     ',1,    'A','IF density ',        phys_decomp)
call addfld ('CAPE_3','J/kg     ',1,    'A','IF density ',        phys_decomp)
call addfld ('CAPE_4','J/kg     ',1,    'A','IF density ',        phys_decomp)
call add_default('CAPE_1', 1, ' ')
call add_default('CAPE_2', 1, ' ')
call add_default('CAPE_3', 1, ' ')
call add_default('CAPE_4', 1, ' ')

call addfld ('CFRACTION_1',  'no dimensions', 1 , 'A', 'cfraction_1',  phys_decomp)
call addfld ('CFRACTION_2',  'no dimensions', 1 , 'A', 'cfraction_2',  phys_decomp)
call addfld ('CFRACTION_3',  'no dimensions', 1 , 'A', 'cfraction_3',  phys_decomp)
call addfld ('CFRACTION_4',  'no dimensions', 1 , 'A', 'cfraction_4',  phys_decomp)
call add_default('CFRACTION_1', 1, ' ')
call add_default('CFRACTION_2', 1, ' ')
call add_default('CFRACTION_3', 1, ' ')
call add_default('CFRACTION_4', 1, ' ')

call addfld ('HMUPRAIN_START',  'no dimensions', 1 , 'A', 'hmuprain_start',  phys_decomp)
call addfld ('HMUPRAIN_SURF',  'no dimensions', 1 , 'A', 'hmuprain_surf',  phys_decomp)
call addfld ('HMANSNOW_START',  'no dimensions', 1 , 'A', 'hmansnow_start',  phys_decomp)
call addfld ('HMANSNOW_SURF',  'no dimensions', 1 , 'A', 'hmansnow_surf',  phys_decomp)
call addfld ('HMANRAIN_SURF',  'no dimensions', 1 , 'A', 'hmanrain_surf',  phys_decomp)
call add_default('HMUPRAIN_START', 1, ' ')
call add_default('HMUPRAIN_SURF', 1, ' ')
call add_default('HMANSNOW_START', 1, ' ')
call add_default('HMANSNOW_SURF', 1, ' ')
call add_default('HMANRAIN_SURF', 1, ' ')

!-----------------------------------------------------------------------
!
!  5. 3D Start variables
!
!-----------------------------------------------------------------------
call addfld ('MASS_START_1','J/kg     ',1,    'A',' mass_start ',        phys_decomp)
call add_default('MASS_START_1', 1, ' ')
call addfld ('MASS_START_2','J/kg     ',1,    'A',' mass_start ',        phys_decomp)
call add_default('MASS_START_2', 1, ' ')

call addfld ('MASS_INC','J/kg     ',1,    'A',' crat ',        phys_decomp)
call add_default('MASS_INC', 1, ' ')
call addfld ('MASS_RAT','J/kg     ',1,    'A',' crat ',        phys_decomp)
call add_default('MASS_RAT', 1, ' ')

!-----------------------------------------------------------------------
!
! 7. 3D nonlinearities; org variables
!
!-----------------------------------------------------------------------
call addfld ('OM5',  'no dimensions', 1 , 'A', 'amplify',  phys_decomp)
call add_default('OM5', 1, ' ')

call addfld ('KM_WIDTH',  'no dimensions', 1 , 'A', 'km_width',  phys_decomp)
call add_default('KM_WIDTH', 1, ' ')

call addfld ('MD_RAT',  'no dimensions', 1 , 'A', 'md_rat',  phys_decomp)
call addfld ('F_AN',  'no dimensions', 1 , 'A', 'f_an',  phys_decomp)
call addfld ('F_ICE',  'no dimensions', 1 , 'A', 'f_ice',  phys_decomp)
call add_default('MD_RAT', 1, ' ')
call add_default('F_AN', 1, ' ')
call add_default('F_ICE', 1, ' ')

call addfld ('AV_Z_RI',  'no dimensions', 1 , 'A', 'av_z_ri',  phys_decomp)
call addfld ('AV_Z_RI_RH',  'no dimensions', 1 , 'A', 'av_z_ri_rh',  phys_decomp)
call addfld ('AV_Z_RI_DT',  'no dimensions', 1 , 'A', 'av_z_ri_dt',  phys_decomp)
call addfld ('COL_RI_DT',  'no dimensions', 1 , 'A', 'col_ri_dt',  phys_decomp)
call addfld ('COL_RI_RH',  'no dimensions', 1 , 'A', 'col_ri_rh',  phys_decomp)
call addfld ('COL_RL_DT',  'no dimensions', 1 , 'A', 'col_rl_dt',  phys_decomp)
call addfld ('COL_RL_RH',  'no dimensions', 1 , 'A', 'col_rl_rh',  phys_decomp)
call addfld ('ZPEAK',  'no dimensions', 1 , 'A', 'zpeak',  phys_decomp)
call addfld ('PRECIP_ORG',  'no dimensions', 1 , 'A', 'precip_org',  phys_decomp)
call addfld ('PRECIP_REC',  'no dimensions', 1 , 'A', 'precip_rec',  phys_decomp)
call add_default('AV_Z_RI', 1, ' ')
call add_default('AV_Z_RI_RH', 1, ' ')
call add_default('AV_Z_RI_DT', 1, ' ')
call add_default('COL_RI_DT', 1, ' ')
call add_default('COL_RI_RH', 1, ' ')
call add_default('COL_RL_DT', 1, ' ')
call add_default('COL_RL_RH', 1, ' ')
call add_default('ZPEAK', 1, ' ')
call add_default('PRECIP_ORG', 1, ' ')
call add_default('PRECIP_REC', 1, ' ')

!-----------------------------------------------------------------------
!
!  8. 3D updraft eff
!
!-----------------------------------------------------------------------
call addfld ('UP_EFF',  'no dimensions', 1 , 'A', 'up eff',  phys_decomp)
call addfld ('UP_RLP',  'no dimensions', 1 , 'A', 'up rlp',  phys_decomp)
call add_default('UP_EFF', 1, ' ')
call add_default('UP_RLP', 1, ' ')

!-----------------------------------------------------------------------
!
!  9. 4DF tend
!
!-----------------------------------------------------------------------
call addfld ('TPTEND_REAL','K/day     ',pver, 'A','actual T tendency IF convection', phys_decomp)
call addfld ('RHTEND_REAL','K/day     ',pver, 'A','actual T tendency IF convection', phys_decomp)
call addfld ('QTEND_REAL','K/day     ',pver, 'A','actual Q tendency IF convection', phys_decomp)
call addfld ('KMTEND  ','J/day ',pver, 'A','up+down km tendency IF convection', phys_decomp)
call addfld ('RVTEND  ','kg/kg/day ',pver, 'A','up+down rv tendency IF convection', phys_decomp)
call addfld ('RLTEND  ','kg/kg/day ',pver, 'A','up+down rl tendency IF convection', phys_decomp)
call addfld ('RLNEW  ','kg/kg/day ',pver, 'A','up+down rl tendency IF convection', phys_decomp)
call addfld ('RITEND  ','kg/kg/day ',pver, 'A','up+down ri tendency IF convection', phys_decomp)
call addfld ('UTENDIF ','m/s/day ',pver, 'A','up+down u tendency IF convection', phys_decomp)
call addfld ('VTENDIF ','m/s/day ',pver, 'A','up+down v tendency IF convection', phys_decomp)
call add_default('RHTEND_REAL', 1, ' ')
call add_default('TPTEND_REAL', 1, ' ')
call add_default('QTEND_REAL', 1, ' ')
call add_default('KMTEND', 1, ' ')
call add_default('RLTEND', 1, ' ')
call add_default('RLNEW', 1, ' ')
call add_default('RITEND', 1, ' ')
call add_default('UTENDIF', 1, ' ')
call add_default('VTENDIF', 1, ' ')

!-----------------------------------------------------------------------
!
!  10. 4DI
!
!-----------------------------------------------------------------------
call addfld ('FMASS','kg/m2/s  ',pver+1, 'A','updraft mass flux IF convection', phys_decomp)
call addfld ('FMASS_DN','kg/m2/s  ',pver+1, 'A','updraft mass flux IF convection', phys_decomp)
call addfld ('PHALF','Pa       ',pver+1, 'A','interface pressure levels', phys_decomp)
call add_default('FMASS', 1, ' ')
call add_default('FMASS_DN', 1, ' ')
call add_default('PHALF', 1, ' ')

!-----------------------------------------------------------------------
!
!  11. 4DF 
!
!-----------------------------------------------------------------------
call addfld ('ENT   ','kg/m2/s  ',pver, 'A','updraft entrainment', phys_decomp)
call addfld ('DET   ','kg/m2/s  ',pver, 'A','updraft detrainment', phys_decomp)
call add_default('ENT', 1, ' ')
call add_default('DET', 1, ' ')

call addfld ('UPDET_1  ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call add_default('UPDET_1', 1, ' ')
call addfld ('UPDET_2  ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call add_default('UPDET_2', 1, ' ')

call addfld ('RH_DN','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call addfld ('DRV_DNDET ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call addfld ('DT_DNDET ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call addfld ('DNDET ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call addfld ('DNENT ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call addfld ('MTDET ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call addfld ('MTENT ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call addfld ('DRV_UP ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call addfld ('DRV_AN ','kg/m2/s  ',pver, 'A','nt', phys_decomp)
call add_default('RH_DN', 1, ' ')
call add_default('DRV_DNDET', 1, ' ')
call add_default('DT_DNDET', 1, ' ')
call add_default('DNDET', 1, ' ')
call add_default('DNENT', 1, ' ')
call add_default('MTDET', 1, ' ')
call add_default('MTENT', 1, ' ')
call add_default('DRV_UP', 1, ' ')
call add_default('DRV_AN', 1, ' ')

!-----------------------------------------------------------------------
!
!  12. 4DF 
!
!-----------------------------------------------------------------------
call addfld ('RHIF    ','no dim',pver, 'A','IF RH       ', phys_decomp)
call addfld ('RHNEW    ','no dim',pver, 'A','IF RH       ', phys_decomp)
call addfld ('TIF    ','no dim',pver, 'A','IF T       ', phys_decomp)
call addfld ('HM    ','no dim',pver, 'A','IF HM       ', phys_decomp)
call addfld ('PROB_RANDOM','no dim',pver, 'A','IF PROB_RANDOM       ', phys_decomp)
call addfld ('RI_HALF','no dim',pver, 'A','IF RI_HALF       ', phys_decomp)
call addfld ('RL_HALF','no dim',pver, 'A','IF RL_HALF       ', phys_decomp)
call add_default('RHIF    ', 1, ' ')
call add_default('RHNEW    ', 1, ' ')
call add_default('TIF    ', 1, ' ')
call add_default('HM    ', 1, ' ')
call add_default('PROB_RANDOM ', 1, ' ')
call add_default('RI_HALF', 1, ' ')
call add_default('RL_HALF', 1, ' ')

!-----------------------------------------------------------------------
!
!  13. 4DF Variables scaling with updraft mass flux (i.e. updraft properties)
!
!-----------------------------------------------------------------------
call addfld ('MP      ','kg/m2',pver, 'A','UPDRAFT MASS', phys_decomp)
call addfld ('BP_CAPE ','m/s2 ',pver, 'A','UPDRAFT BPR ', phys_decomp)
call addfld ('BP1_STR ','m/s2 ',pver, 'A','UPDRAFT BPR ', phys_decomp)
call addfld ('BP1_PRP ','m/s2 ',pver, 'A','UPDRAFT BM2 ', phys_decomp)
call addfld ('BP1_ITR ','m/s2 ',pver, 'A','UPDRAFT BM3 ', phys_decomp)
call addfld ('BP2_STR ','m/s2 ',pver, 'A','UPDRAFT BPR ', phys_decomp)
call addfld ('BP2_PRP ','m/s2 ',pver, 'A','UPDRAFT BM2 ', phys_decomp)
call addfld ('BP2_ITR ','m/s2 ',pver, 'A','UPDRAFT BM3 ', phys_decomp)
call addfld ('F_ENT_1 ','m/s2 ',pver, 'A','UPDRAFT BM3 ', phys_decomp)
call addfld ('F_DET_1 ','m/s2 ',pver, 'A','UPDRAFT BM3 ', phys_decomp)
call addfld ('F_ENT_2 ','m/s2 ',pver, 'A','UPDRAFT BM3 ', phys_decomp)
call addfld ('F_DET_2 ','m/s2 ',pver, 'A','UPDRAFT BM3 ', phys_decomp)
call addfld ('RLP     ','kg/kg',pver, 'A','UPDRAFT RLP ', phys_decomp)
call addfld ('RIP     ','kg/kg',pver, 'A','UPDRAFT RLP ', phys_decomp)
call addfld ('HMP_MEAN ','kg/kg',pver, 'A','UPDRAFT HMP ', phys_decomp)
call addfld ('HM_DIFF ','kg/kg',pver, 'A','UPDRAFT HM_DIFF ', phys_decomp)
call add_default('MP      ', 1, ' ')
call add_default('BP_CAPE ', 1, ' ')
call add_default('BP1_STR ', 1, ' ')
call add_default('BP1_PRP ', 1, ' ')
call add_default('BP1_ITR ', 1, ' ')
call add_default('BP2_STR ', 1, ' ')
call add_default('BP2_PRP ', 1, ' ')
call add_default('BP2_ITR ', 1, ' ')
call add_default('F_ENT_1 ', 1, ' ')
call add_default('F_DET_1 ', 1, ' ')
call add_default('F_ENT_2 ', 1, ' ')
call add_default('F_DET_2 ', 1, ' ')
call add_default('RLP     ', 1, ' ')
call add_default('RIP     ', 1, ' ')
call add_default('HMP_MEAN', 1, ' ')
call add_default('HM_DIFF ', 1, ' ')

!-----------------------------------------------------------------------
!
!  14. 4DF Variables scaling with updraft detrainment mass
!
!-----------------------------------------------------------------------
call addfld ('TDIFF_DETRAIN_UP','kg/kg',pver, 'A','UPDRAFT TDIFF_UP ', phys_decomp)
call addfld ('CAPE_DETRAIN_UP','kg/kg',pver, 'A','UPDRAFT CAPE ', phys_decomp)
call addfld ('MASS_DETRAIN_UP','kg/kg',pver, 'A','UPDRAFT MASS_UP ', phys_decomp)
call addfld ('HM_DIFF_DET','kg/kg',pver, 'A','HM_DIFF_DET ', phys_decomp)
call add_default('TDIFF_DETRAIN_UP', 1, ' ')
call add_default('CAPE_DETRAIN_UP', 1, ' ')
call add_default('MASS_DETRAIN_UP', 1, ' ')
call add_default('HM_DIFF_DET', 1, ' ')

!-----------------------------------------------------------------------
!
!  15. 4DF Evaporation variables
!
!-----------------------------------------------------------------------
call addfld ('UPRAINEVAP','kg/kg',pver, 'A','UPRAINEVAP ', phys_decomp)
call addfld ('ANRAINEVAP','kg/kg',pver, 'A','ANRAINEVAP ', phys_decomp)
call addfld ('ANSNOWSUBL','kg/kg',pver, 'A','ANSNOWSUBL ', phys_decomp)
call addfld ('ANSNOWMELT','kg/kg',pver, 'A','ANSNOWMELT ', phys_decomp)
call add_default('UPRAINEVAP', 1, ' ')
call add_default('ANRAINEVAP', 1, ' ')
call add_default('ANSNOWSUBL', 1, ' ')
call add_default('ANSNOWMELT', 1, ' ')

call addfld ('UPRAINDOWN','kg/kg',pver, 'A','UPRAINDOWN ', phys_decomp)
call addfld ('ANRAINDOWN','kg/kg',pver, 'A','ANRAINDOWN ', phys_decomp)
call add_default('UPRAINDOWN', 1, ' ')
call add_default('ANRAINDOWN', 1, ' ')

!-----------------------------------------------------------------------
!
!  16. 4DF Precipitation Production
!
!-----------------------------------------------------------------------
call addfld ('UPRAIN_RLP','kg/kg',pver, 'A','UPRAIN_RLP ', phys_decomp)
call addfld ('ANSNOW_RLP','kg/kg',pver, 'A','ANSNOW_RLP ', phys_decomp)
call addfld ('UPRAIN_RV','kg/kg',pver, 'A','UPRAIN_RV ', phys_decomp)
call addfld ('ANSNOW_RV','kg/kg',pver, 'A','ANSNOW_RV ', phys_decomp)
call addfld ('ANSNOW_RI','kg/kg',pver, 'A','ANSNOW_RI ', phys_decomp)
call addfld ('RLP_RVV','kg/kg',pver, 'A','RLP_RVV ', phys_decomp)
call add_default('UPRAIN_RLP', 1, ' ')
call add_default('ANSNOW_RLP', 1, ' ')
call add_default('UPRAIN_RV', 1, ' ')
call add_default('ANSNOW_RV', 1, ' ')
call add_default('ANSNOW_RI', 1, ' ')
call add_default('RLP_RVV', 1, ' ')

!-----------------------------------------------------------------------
!
!  17. 4DF Cloud Production/Loss
!
!-----------------------------------------------------------------------
call addfld ('UPRAIN_RL_TEND','kg/kg',pver, 'A','UPRAIN_RL_TEND', phys_decomp)
call addfld ('RL_EVAP_TEND','kg/kg',pver, 'A','RL_EVAP_TEND', phys_decomp)
call addfld ('RI_EVAP_TEND','kg/kg',pver, 'A','RI_EVAP_TEND', phys_decomp)
call addfld ('RI_PROD_TEND','kg/kg',pver, 'A','RI_PROD_TEND', phys_decomp)
call addfld ('RL_DET_TEND','kg/kg',pver, 'A','RL_DET_TEND', phys_decomp)
call addfld ('RI_DET_TEND','kg/kg',pver, 'A','RI_DET_TEND', phys_decomp)
call addfld ('RL_VERT_TEND','kg/kg',pver, 'A','RL_VERT_TEND', phys_decomp)
call addfld ('RI_VERT_TEND','kg/kg',pver, 'A','RI_VERT_TEND', phys_decomp)
call addfld ('RL_RH','kg/kg',pver, 'A','RL_RH', phys_decomp)
call addfld ('RL_DT','kg/kg',pver, 'A','RL_DT', phys_decomp)
call addfld ('RI_RH','kg/kg',pver, 'A','RI_RH', phys_decomp)
call addfld ('RI_RS','kg/kg',pver, 'A','RI_RS', phys_decomp)
call addfld ('RI_DT','kg/kg',pver, 'A','RI_DT', phys_decomp)
call add_default('UPRAIN_RL_TEND', 1, ' ')
call add_default('RL_EVAP_TEND', 1, ' ')
call add_default('RI_EVAP_TEND', 1, ' ')
call add_default('RI_PROD_TEND', 1, ' ')
call add_default('RL_DET_TEND', 1, ' ')
call add_default('RI_DET_TEND', 1, ' ')
call add_default('RL_VERT_TEND', 1, ' ')
call add_default('RI_VERT_TEND', 1, ' ')
call add_default('RL_RH', 1, ' ')
call add_default('RL_DT', 1, ' ')
call add_default('RI_RH', 1, ' ')
call add_default('RI_DT', 1, ' ')
call add_default('RI_RS', 1, ' ')

!------------------------------------------------------------------------
!
!  added in cam5:
!
!  IF: this is ZM stuff and could remove ...
!
!------------------------------------------------------------------------
    call phys_getopts(history_budget_out = history_budget)
    if ( history_budget ) then
       call add_default('EVAPTZM  ', 1, ' ')
       call add_default('EVAPQZM  ', 1, ' ')
       call add_default('ZMDT     ', 1, ' ')
       call add_default('ZMDQ     ', 1, ' ')
       call add_default('ZMDLIQ   ', 1, ' ')
       call add_default('ZMDICE   ', 1, ' ')

       if( cam_physpkg_is('cam5') ) then
          call add_default('ZMMTT    ', 1, ' ')
       end if

    end if
!------------------------------------------------------------------------
!
!  Removed: definition of limcnv
!
!  Limit deep convection to regions below 40 mb
!  Note this calculation is repeated in the shallow convection interface
!
!  hypi(1) would be the reference interface pressure at the top model level.
!  No reason why limcnv should not be retained I guess, though hard to see why
!   it should be needed.
!
!------------------------------------------------------------------------
end subroutine if_conv_init
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  Main Program
!  - called from subroutine convect_deep_tend
!
!  The variables pblht,tpert,qpert,fsns, fsnt,flns,flnt, go all the way to the top
!   of the module physpkg. How can I escape them????
!  - pblht and tpert are initialized in "subroutine phys_inidat" of module physpkg,
!    and are then inputs only downward to "subroutine tphysbc". The zm scheme
!    checks that the convection goes above the pblh, and stops if not, presumably
!    to prevent some coupling between pbl and zm schemes. I don't think I should
!    have to worry about.
!
!   - pblht and tpert are inputs from module physpkg to subroutine tphysbc
!
!
!
! ********************  Starting from the TOP ************************
!
!  subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf, cam_in, cam_out)
!
!    CONTAINS
!
! call tphysbc (ztodt, pblht(1,c), tpert(1,c), qpert(1,1,c),               &  (IN module physpkg)
!               fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), phys_state(c), &
!               phys_tend(c), pbuf,  fsds(1,c), landm(1,c),                &
!               cam_out(c), cam_in(c) )
!
!  BE CAREFUL WITH ANY CHANGES TO THESE ARGUMENTS:
!
!  subroutine tphysbc (ztodt,   pblht,   tpert,   qpert,   & (IN subroutine tphysbc)
!                   fsns,    fsnt,    flns,    flnt,    state,   &
!                   tend,    pbuf,    fsds,    landm,            &
!                   cam_out, cam_in )
!
!  FREELY GET RID OF ZM STUFF BELOW THIS POINT
!
!  Note changes here: dlf,pflx,rliq,snow_zmc not being passed up, but stop
!  in tphysbc: why are they needed in this module? I think mainly because
!  all the deep mass fluxes are used as inputs in "call convect_shallow_tend",
!  and rliq is needed for energy conservation, and variables like snow_zmc
!  are output in a "call diag_conv".
!  - To retain generality, I should keep passing them up ....
!
!  
!  call convect_deep_tend(  prec_zmc,   &   (IN subroutine tphybc)
!       pblht,    cmfmc,      cmfcme,             &
!       tpert,    dlf,        pflx,    zdu,       &
!       rliq,    &
!       ztodt,    snow_zmc,  &
!       state,   ptend, pbuf )
!
!  subroutine convect_deep_tend(prec    , &   (IN module convect_deep)
!    pblh    ,mcon    ,cme     ,          &
!    tpert   ,dlf     ,pflx    ,zdu      , &
!    rliq    , &
!    ztodt   ,snow    ,&
!    state   ,ptend   ,pbuf  )
!
! ** HERE IS WHERE I SHOULD START REMOVING ARGUMENTS SINCE I AM SPECIALIZING TO IF **
!
!    call zm_conv_tend(prec     , &    (IN module convect_deep)
!         pblh    ,mcon    ,cme     ,          &
!         tpert   ,dlf     ,pflx    ,zdu      , &
!         rliq    , &
!         ztodt   ,snow    ,&
!         jctop, jcbot , &
!         state   ,ptend   ,pbuf  )
!
! subroutine zm_conv_tend(prec    , &  (IN module zm_conv_intr)
!     pblh    ,mcon    ,cme     ,          &
!     tpert   ,dlf     ,pflx    ,zdu      , &
!     rliq    , &
!     ztodt   ,snow    ,&
!     jctop   ,jcbot , &
!     state   ,ptend_all   ,pbuf  )
!
!  The variables in the original "call zm_conv_tend" are identical to the
!    variables of the original "call convect_deep_tend" from tphysbc:
!
!
!  So, the variables are just being passed up the line from zm_convr:
!
! "subroutine tphysbc"  "module convect_deep"  "module zm_conv_intr"
!   prec_zmc                         
!
!   REMOVE arguments that have no relation to IF param (consider specialized)
!
!  I removed: pblh,mcon,cme,tpert,pflx,zdu,rliq,snow
!  in cam5 they added landfrac
!
!  Tendencies (current interpretation Feb 2013):
!  - ptend_loc is defined in the call to if_conv. The various "ptend_loc" are
!  then collected (maybe just one in my case) into "ptend_all", which is passed 
!  up as an argument here to
!  convect_deep_tend, which is then passed up as "ptend" to tphysbc.F90, 
!  where it is then applied to the "real", as opposed to local, physics 
!  state vector.
!
!------------------------------------------------------------------------
subroutine if_conv_tend(rain_if        ,snow_if                         ,&
   ztodt               ,jctop          ,jcbot                ,state     ,&
   ptend_all           ,landfrac       ,pbuf                            ,&
   preciporg, massprev, cm5prev)
!------------------------------------------------------------------------
!
!  Use statements
!
!------------------------------------------------------------------------
   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend, physics_tend
   use physics_types, only: physics_ptend_init,  physics_tend_init,physics_update
   use physics_types, only: physics_state_copy
   use physics_types, only: physics_ptend_sum

   use phys_grid,     only: get_lat_p, get_lon_p
   use time_manager,  only: get_nstep, is_first_step
   use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use constituents,  only: pcnst, cnst_get_ind
   use check_energy,  only: check_energy_chng
   use physconst,     only: gravit
!------------------------------------------------------------------------
!
!  IN/OUT Arguments
!   - got rid of zm stuff (mcon,dlf,pflx,cme,zdu)
!
!------------------------------------------------------------------------
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend_all      ! individual parameterization tendencies
   real(r8), intent(in) :: landfrac(pcols)            ! RBN Landfrac

   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(out) :: rain_if(pcols)   ! total precipitation from IF convection
   real(r8), intent(out) :: snow_if(pcols)   ! snow from IF convection
!------------------------------------------------------------------------
!
!  These are pointers higher up but arrays here
!
!------------------------------------------------------------------------
  real(r8), intent(inout) :: preciporg(pcols)  ! 
  real(r8), intent(inout) :: massprev(pcols)  ! 
  real(r8), intent(inout) :: cm5prev(pcols)  ! 
!------------------------------------------------------------------------
!
!  OBS
!
!------------------------------------------------------------------------
   real(r8) :: prec_if(pcols)
!------------------------------------------------------------------------
!
! Integer local variables
!
!------------------------------------------------------------------------
   integer :: i,k,m
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: nstep
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns
   integer :: itim, ifld              ! for physics buffer fields 
!------------------------------------------------------------------------
!
! Real local variables
!
!------------------------------------------------------------------------
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
!------------------------------------------------------------------------
!
! physics types
!
! Complex Issue:
! ==============
! What is state1? 
!   - In ZM, the q/T updraft tendencies, the evaporation tendencies,
!   and the tracer tendencies are all separately applied to the state vector. 
!   - A new state state1
!   is re-calculated after each of these calls, and this new state is the input
!   for the next tendency calculation. I am not sure why this is done.
!   - I don't think I need state1 since I do
!   everything together. The distinction between ptend_all and ptend_loc may no
!   longer be needed. The ptend_loc that comes back from "call if_convr" is
!   the same as ptend_all; no need to add the various ptend_loc from the different
!   calls. There is just 1 "call physics_ptend_sum(ptend_loc,ptend_all, state)"
!   here.
!
!------------------------------------------------------------------------
  type(physics_state) :: state1     ! locally modify for evaporation to use, not returned
  type(physics_tend ) :: tend       ! Physics tendencies (empty, needed for physics_update call)
  type(physics_ptend) :: ptend_loc  ! package tendencies
!------------------------------------------------------------------------
!
! physics buffer fields 
! - why are these needed? (see discussion below)
! - needed across timesteps; not sure why
!
!  cam5 adds evapcdp
!
!------------------------------------------------------------------------
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation  
!------------------------------------------------------------------------
!
!  Other arrays:
!  removed: pcont,pconb
!
!------------------------------------------------------------------------
   real(r8) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
!------------------------------------------------------------------------
!
! History output fields
! - removed cape
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
! Used in momentum transport calculation
!
!------------------------------------------------------------------------
   real(r8) :: winds(pcols, pver, 2)
   real(r8) :: wind_tends(pcols, pver, 2)
   logical  :: l_windt(2)   ! momentum switch for u and v
   integer  :: ii,ixcldliq,ixcldice

!------------------------------------------------------------------------
!
! Initialize 
!
!------------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem = 0._r8   
   wind_tends(:ncol,:pver,:) = 0.0_r8
!------------------------------------------------------------------------
!
!  What is state1?
!  type(physics_state) :: state1     ! locally modify for evaporation to use, not returned
!
!    subroutine physics_state_copy(state_in, state_out)
!    - this routine just transfers the contents of state_in to state_out
!    - in this case state to state1
!    - think you need this as the most recent "reference" for conservation tests.
!
!------------------------------------------------------------------------
!   print*,'calling physics_state_copy at start of if_conv_tend: state1 defined'
   call physics_state_copy(state,state1)   ! copy state to local state1.
!   print*,'calling physics_ptend_init for ptend_loc'
   call physics_ptend_init(ptend_loc)  ! initialize local ptend type
!   print*,'calling physics_ptend_init for ptend_all'
   call physics_ptend_init(ptend_all)  ! initialize output ptend type
!   print*,'calling physics_ptend_init for tend'
   call physics_tend_init(tend)        ! tend type here is a null place holder
!   print*,'finished calling physics_ptend_init for tend'
!------------------------------------------------------------------------
!
! Associate pointers with physics buffer fields
!
! Pointer cld:
! ------------
! real(r8), pointer, dimension(:,:) :: cld
! Pointer cld was later used as an argument in "call zm_conv_evap", where its 
!   called "cldfrc", an "in" argument
! real(r8), intent(in   ) :: cldfrc(pcols,pver) ! cloud fraction
! It seems like the cloud frcation helps determine how much of the falling
! precipitation falls through cloud and therefore does not evaporate.
! - I don't use this in my evaporation of convective precip (maybe later should)
! - probably safer to leave; no harm having the pointer associated with this
!   part of the buffer.
!
! Pointer ql:
! -----------
! Pointer ql was later used as an argument in "call zm_convr", also called ql
! !  wg * ql       grid slice of cloud liquid water.
! in zm_convr "ql" i used to determine the gathered version qlg, where its
!   used as an argument in "call cldprp"
! - I don't use but again probably not doing any harm
! - not sure if stratiform water or convective updraft water
!
! Pointer rprd:
! ------------
! Pointer rprd was later used as an output argument in "call zm_convr",
!  where also called "rprd", (rain production rate) and
!  an "in" argument "call zm_conv_evap" (where called "prdprec"), where
!  called precipitation production
! - not sure why this would have to be retained across timesteps.
!
! Pointer fracis:
! --------------
! - used as an argument in "call convtran"
! real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
! real(r8), intent(in) :: fracis(pcols,pver,ncnst) ! fraction of tracer that is insoluble
! - in convtran, fracis is an "in" argument, 
! - how is it determined? (presumably does not have to be 0 or 1, but somehow
!   parameterized to be intermediate values)
! - don't think I have to worry about this for now, only if using more chemistry.
!
! LEAVE ALL THESE FOR NOW; Do I need?
!
!------------------------------------------------------------------------
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('ICWMRDP')
   ql => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

   ifld = pbuf_get_fld_idx('RPRDDP')
   rprd => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

   ifld = pbuf_get_fld_idx('FRACIS')
   fracis  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1:pcnst)

!------------------------------------------------------------------------
!  cam5 adds:
!------------------------------------------------------------------------
   ifld = pbuf_get_fld_idx('NEVAPR_DPCU')
   evapcdp => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
!------------------------------------------------------------------------
!
! Begin with IAF convection parameterization
!
!------------------------------------------------------------------------
!  print*,'calling t_startf'
  call t_startf ('if_convr')
!  print*,'finished calling t_startf'
!------------------------------------------------------------------------
!
!  IAF: Added this so winds defined on entry into if_convr
!   - why using state1 here? (don't think it matters: copy routine above
!     should mean state and state1 are the same)
!
!------------------------------------------------------------------------
  winds(:ncol,:pver,1) = state1%u(:ncol,:pver)
  winds(:ncol,:pver,2) = state1%v(:ncol,:pver)
  l_windt(1) = .true.  ! turn on u mom transport
  l_windt(2) = .true.  ! turn on v mom transport
!------------------------------------------------------------------------
!
!  IAF: added 
!  - logical lq fed into if_conv as "doconvtran", and tells you to return
!    tendencies. To obey conservation
!    properties has to be a consistency between these logical flags and what
!    convection is actually transporting. 
!
!------------------------------------------------------------------------
  call cnst_get_ind('CLDLIQ', ixcldliq)
  ptend_loc%lq(ixcldliq) = .TRUE.
  call cnst_get_ind('CLDICE', ixcldice)
  ptend_loc%lq(ixcldice) = .TRUE.
!------------------------------------------------------------------------
!
!  New call to if_convr:
!  - added snow. This is computed by IAF param, but previously came out of
!    "call zm_conv_evap"
!
!  - why is there a 0.5 in fromt of ztodt????
!  - I re-define ptend_loc%lq after anyway, why call it here ...
!
!  Note that this call define ptend_loc, not ptend itself
!
!------------------------------------------------------------------------
!print*,'calling if_convr'

 call if_convr( lchnk        ,ncol          ,state%omega        ,              &
                state%t      ,state%q       ,l_windt            ,winds    ,    &
                2            ,state%pmid    ,state%pint         ,state%zm ,    &
                state%phis   ,state%zi      ,ptend_loc%q(:,:,1) ,ptend_loc%q,  &
                ptend_loc%s  ,wind_tends    ,jctop              ,jcbot    ,    &
                prec_if      ,.5_r8*ztodt   ,ptend_loc%lq       ,pcnst    ,    &
                snow_if      ,rain_if       ,preciporg          ,massprev, cm5prev, landfrac )

!print*,'finished if_convr'
!------------------------------------------------------------------------
!
!  removed a bunch of ZM outfld calls:
!
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  Removed:
!  call outfld('CAPE', cape, pcols, lchnk)        ! RBN - CAPE output
!
!  Fractional occurrence of ZM convection
!  Also don't bother outputting 'FREQZM' right now. 
!
!  freqzm(:) = 0._r8
!  do i = 1,lengath(lchnk)
!     freqzm(ideep(i,lchnk)) = 1.0
!  end do
!  call outfld('FREQZM  ',freqzm          ,pcols   ,lchnk   )
!
!  'ZMMU','ZMMD' output
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!  Changed tendencies so that do everything in one call.
!
!  Heating rate output:
!  - my calculation of heat should be the same as zm, so should not have 
!    to change much
!  - in ptend_loc%name  = 'if_convr' changed name to 'if_convr'
!
!  Momentum Transport
!
!  See section 3.3.3 to understand physics_ptend type.
!  Or look for type physics_ptend in module physics_types
!  'zm_convr' is the name of the parameterization associated with the tendency
!  lq(1) refers to specific humidity Q , l for logical
!  TRUE means, e.g. "true if dsdt is returned"; program tries
!    to update this quantity when next call physics_update
!
!------------------------------------------------------------------------
  ptend_loc%name  = 'if_convr'
  ptend_loc%ls = .TRUE.     ! adjust Temperature
  ptend_loc%lu = .TRUE.     ! adjust u
  ptend_loc%lv = .TRUE.     ! adjust v (non-zero?)
  ptend_loc%lq(1) = .TRUE.  ! adjust q - CLDLIQ set to .TRUE. above
                            ! CLDLIQ set to .TRUE. above but could be done here
                            ! CLDICE set to .TRUE. above but could be done here
!------------------------------------------------------------------------
!
!  Output T tendency 
!   - obtain from dse tendency by dividing by cpair (OK if z fixed).
!   - Dimensions of s are (J/kg/s)
!   - 'ZMDT' -> 'IFDT'
!
!------------------------------------------------------------------------
  ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
  call outfld('IFDT    ',ftem           ,pcols   ,lchnk   )
!------------------------------------------------------------------------
!
!  Output the Q tendency
!  - q array has dimension(pcols,pver,pcnst), so last 1 is for Q, but why other 2?
!  - Just outputting Q tendency at one spot?
!
!  Maybe add output of rl/ag tendencies later?
!
!------------------------------------------------------------------------
  call outfld('IFDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
  call t_stopf ('if_convr')
!------------------------------------------------------------------------
!
!  Define ptend for winds from wind_tends.
!
!  Output 'IF' Wind tendencies
!  How can I get away with "wind_tends(1,1,1)" here. Do I?
!  Even if I just pick first column, what about height dependence?
!  Maybe outfld is smart enough just to take the starting point, and loop
!   over pcols and lchnk ...
!
!------------------------------------------------------------------------
  ptend_loc%u(:ncol,:pver) = wind_tends(:ncol,:pver,1)
  ptend_loc%v(:ncol,:pver) = wind_tends(:ncol,:pver,2)
  call outfld('IFMTU', wind_tends(1,1,1), pcols, lchnk)
  call outfld('IFMTV', wind_tends(1,1,2), pcols, lchnk)
!------------------------------------------------------------------------
!
! Add tendency from 'if_convr' to tendencies from other processes
!  - again, this distinction between ptend_all and ptend_loc is likely redundant
!    here since I do all convection tendencies in one fell swoop
!  - ptend_all is what is passed up higher
!  - What about change in surface pressure? When is this adjusted in response
!    to change in water vapor?
!
!------------------------------------------------------------------------
!  print*,'calling physics_ptend_sum'
  call physics_ptend_sum(ptend_loc,ptend_all, state)
!  print*,'called physics_ptend_sum to define ptend_all'
!   print*,'ptend_all%q(1,16,4) = ',ptend_all%q(1,16,4)
!   print*,'ptend_all%q(1,16,5) = ',ptend_all%q(1,16,5)
!------------------------------------------------------------------------
!
!  Update physics state type state1 with ptend_loc 
!  - why state1?
!  - the real state itself is updated higher up in tphysbc 
!  - maybe don't need this?
!
!  Note that ptend_loc used here!!
!  ptend use din call from tphysbc.F90
!
!------------------------------------------------------------------------
! print*,'calling physics_update using ptend_loc for state1 from if_conv_intr.F90'
  call physics_update(state1, tend, ptend_loc, ztodt)
! print*,'finished calling physics_update using ptend_loc for state1 from if_conv_intr.F90'
!------------------------------------------------------------------------
!
!  Initialize ptend for next process
!  - zeroes all tendencies and set all logical variables to .false.
!
!------------------------------------------------------------------------
  call physics_ptend_init(ptend_loc)
!------------------------------------------------------------------------
!
!   Removed zm_conv_evap
!   Removed call momtran
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!  Removed: 
!  - Output apparent force from  pressure gradient
!  - Output in-cloud winds
!
!  Removed call convtran to transport gridscale cloud water and ice 
!  - Grid scale cloud water/ice are subject to convective transport in ZM
!  - now doing this in my scheme also.
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!  What does this do?
!
!  Note the comment below: changes to cloud water and ice will be applied
!    later, I believe in the call physics_update in subroutine tphysbc
!
! ptend_all will be applied to original state on return to tphysbc
! This name triggers a special case in physics_types.F90:physics_update()
!
!------------------------------------------------------------------------
  ptend_all%name = 'convect_deep'
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
end subroutine if_conv_tend
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  Removed: zm_conv_tend_2
!
!------------------------------------------------------------------------
end module if_conv_intr
