!------------------------------------------------------------------------
!
!   *****************    WHAT TO PUT HERE?  ************************
!
!   The only reason to put an array or array dimension here is if it is used 
!    in BOTH if_driver.f90 and if_conv_tend.f90. Otherwise put where used.
!   If not sure, safe thing is to put here.
!
! How to run (Sept 2012)
! - ssh -l folkins cluster.mathstat.dal.ca
! - f0lkins
! - qrsh
! - cd /users/clouds/folkins/BUBBLE
! - ifort -C shr_kind_mod.f90 if_conv_solvers.f90 if_conv_params.f90 if_conv_tend.f90 if_driver.f90 rad.conv.o 
! - may have to do: ifort -c rad.conv.f
! - a.out
!
!   ==============  OBSOLETE ================
!
!   ssh -p 44747 compute.mathstat.dal.ca  (can no longer log on - Sept 2012)
!   ?!.LlMm$1
!   cd BUBBLE
!
!   How to compile (April 2010)  - OBSOLETE
!   ---------------------------------------
! f77 -c rad.conv.f   (Have to recompile if change)
! f90 -C shr_kind_mod.f90 if_conv_solvers.f90 if_conv_params.f90 if_conv_tend.f90 if_driver.f90 rad.conv.o 
!
! -C   Check array references for out of range subscripts and
!          conformance.
! 
! -depend   Analyze loops for data dependence and restructuring.
!
! -fpover  Detect floating-point overflow in formatted input
!
!   Order is important: lowest modules used by other modules go first.
!   In cam4 compilation order is alphabetical, so changed naming to consistent with this.
!
!  Issues
!  -------
!
!  Break downdraft routine into subroutines
!   - current loop structure is confusing
!   - one possibility is to have a routine which takes a parcel at a particular
!     height, successively moves it down, and detrains. i.e. put almost verything
!     inside the innermost loop in a subroutine.
!  - start by having special subroutines for determining melting/evaporation
!
!  Find ways of peeling off diagnostics at end of if_conv_tend.f90 into subroutines
!
!  lp -> evap parameterization
!   - any errors
!   - sensitivity to changes? 
!   - what are the physical assumptions?
!   - justified to use the same param for rain + snow? Others out there?
!   - Think about imposing a variety of particle size spectra at the base of the
!  anvil, and a variety of lengths, and then have the particle size spectra
!  interact with the evaporation. Right now the Kessler param overides interactions
!  which might be interesting, and may make it harder for rain to reach the surface.
!  In reality, reasonably sized rain drops reach the surface quite easilly. One of
!  the reasons I don't have much stratiform rain reaching the surface is because
!  the mean state is dry relative to the enhanced convective state, so maybe not
!  worry about low stratiform precip efficiency.
!
!  Generalized vertical advection
!   - cam4 tracer transport uses net entrainment/detrainment mass fluxes at each level
!   - what is error introduced by not using mixing ent/det of individual parcels/plumes?,
!     (i.e. in assuming that the mixing ratio of every plume/parcel is the same).
!   - If adapt the current convtran, must be modified, since all downdrafts assumed to
!     detrain at the surface. This does seem the easiest approach. But use the in house
!     online rl/rv tendencies obtained from the cloud model.
!   - Could also write my own scheme in which use mass dmflow from individual parcels.
!     Would have to reduce number of updraft/downdraft parcels to avoid prohibitive memory.
!
!  z being defined correctly?  km <-> hm
!    In the model, the half pressure layers are fixed and fundamental 
!    (except the surface pressure which floats) since they define the mass of each layer. 
!    Variables of a layer should refer to the mass weighted average of the layer. 
!    The half level heights are fundamental, since the surface height is zero,
!    and all other heights can be defined in terms of the pressure levels,
!    the average density of the layer (from rv,rl,km), and the hydrostatic 
!    approximation. Are the full level heights z(i) neccessary? I think they
!    are really only needed to go between hm and km. I have used two ways of going
!    between hm and km. One way is to use mix, and km/rt, to solve for t(i),rv, and then
!    use these, plus z(i), in a call to hmoist to find hm. The other way is to just
!    set hm(i) = km(i) + (1+rt)*g*z. This latter way is preferable. In fact, hmoist
!    probably should not be neccessary. But how to define the full level heights?
!    Defining as z(i) = 0.5*(zhalf(i) + zhalf(i+1)) introduces error. This consistently 
!    overestimates the true mass weighted z since there should be more mass
!    in the layer below z(i)-zhalf(i) layer than the zhalf(i+1)-z(i) layer.
!    This overestimate may not be a problem provided the bias is consistent, since
!    the program only cares about differences in hm/km values.
!    April/2010: get non-zero values in kmtend due to numerical rounding. At last call
!      to hydro in previous timestep get new hm(i) from km(i),z(i),rtt. In next timestep
!      in conv_tendencies, find kmnew from hmnew,rtnew,z(i). Finding that for same values
!      of hmnew and hm(i),rtnew and rtt, and z(i), km and kmnew are not consistent, despite
!      using the same inter-relationship between hm/km/z/rtt. It's a small error, maybe
!      not worry about, but odd to have non-zero kmtend in stratosphere due to convection
!      (or just rounding error).
!
!   Onset of Penetrative Downdrafts: 
!     The best way to characterize the onset of penetrative downdrafts is to
!     say that they function as a source/sink of moist enthalpy to the background
!     atmosphere. Static downdrafts do this only a little, via removing moist enthalpy
!     from the background atmosphere, adding mass to layers, and inducing vertical mass
!     transports across pressure surfaces via hydrostatic adjustment.
!
!   Conservation of hm during vertical motion assumes all TKE dissipated
!   - remember to credit work during 70's (Richard Reed I think)
!
!   Melting at the Surface: 
!     - now get melting at the surface due to downdraft starting at
!     level ip just below melting level being transported all the way to the surface after
!     a few passes, then sitting there melting snow. Since has high ip, is one of first layers
!     new snow packets encounter. 
!     - also causes problems in generating negative remaining melted rain
!
!   Move to a more parcel based way of describing rain: e.g. with current formulation, can
!     have hmansnow_surf > 0 by all ansnow gone. Need a formulation in which this
!     could not happen and where km of rain more tightly coupled with rain T/mass
!
!  Vertical Velocities
!    - have vertical velocities in 3d gone down due to decrease in bmix?
!    - ultimately would like to have believable updraft velocities for microphysics
!      and to allow for spectrum of +-B as exists in real weak cumulus.
!    - role of dissipation?
!    - have added effect of mixing?
!
!  Precipitation Loading
!    - how much difference does this make?
! 
!  Role of km_width in rainfall variance
!   - likely that the homogeneity of km in the BL increases variance
!   - How strong is this effect
!   - since BL responds more uniformally to changes aloft? 
!
!   Debugging advice:
!   -----------------
!
! - It is probably best to start debugging by:
!   turn off conversion of updraft parcel condensate to rain (nprecip = 2)
!   turn off updraft mixing (mix_option = 0) 
!   turn off production of anvil condensate (anvil_option = 0)
!   turn off precip of anvil condensate (anvil_precip_option = 0)
!   turn off downdrafts (ndown = 0)
!   turn off penetrative downdrafts by increasing bdthresh
!   (In this case, updrafts remove moisture from precip by evap, so add mass to layers, 
!   so still trigger small vertical motions and non-zero kmtend).
!   turn off BD circulation (bdmassflux = 0)
!
!   This is the simplest possible version of the model.
!  
!   Check conservation properties. Is column enthalpy conserved?
!   If not:
!
!   (i) Check that hmflow + hmprecip = 0. Non-zero values here
!   indicate problems in the definitions of hmflow/hmansnow/hmuprain.
!   If hmflow + hmprecip = 0, but enthalpy is not conserved, the problem is
!   likely not in the definition of hmprecip or hmflow.
!
!   (ii) Check that sum_hmvert = 0 for updrafts and downdrafts. If non-zero, indicates
!     problems in the definition of hmvert.
!
!   (iii) Check that sum_dmvert = 0 for updrafts and downdrafts. If non-zero, indicates
!     problems in the definition of hmvert. (Similarly for sum_rlvert and sum_rvvert).
!
!   (iv) If hmflow + hmprecip = sum_hmvert = sum_dmvert = 0, check calculation of hmnew.
!
!   (v) Check that sum_kmdiff = 0 in if_conv_tend.f90 (should =0 for no precip)
!
!   Progressively add in processes and try and determine the source of the problem.
!
!   Column hm will not be conserved. Due to heat release
!   in conversion of rv to rl, the column will be heated and
!   the gravitational potential energy of the column increased also.
!   This should equal the change in column hm. 
!
!   Sources of column enthalpy errors:
!   ----------------------------------
!
! (1) Incorrect definition of hmflow (exchanges with background atmosphere)
! (2) Incorrect definitions of hmuprain, hmansnow
! (3) Incorrect calculation of hmvert
! (4) Incorrect calculation of hmnew from hmvert/hmflow
! (5) Incorrect calculation of kmtend from hmnew
! (6) Incorrect calculation of km(i) from kmtend 
! (7) Incorrect treatment of pressure levels, especially surface pressure
!
!   Sensitivities:
!   --------------
!
!   Sensitivity to changes in vertical resolution?
!    - prove numerical convergence
!
!   Sensitivity to changes in nlaunch and numdz
!   - important for numerical efficiency in cam4
!
!   Sensitivity of BL height and lapse rate to low-level Convective Downdrafts?
!   - do downdrafts from convective rain play a role in increasing BL height as 
!     approach ITCZ?
!   - often explained as due to weakening of subsidence but subsidence stronger 
!     near ITCZ.
!
!   hmflow -> kmflow
!    - This makes more sense: you are remiving and adding enthalpy at various 
!      levels, not hm; hm useful for vertical motions
!
!   Overall structure
!   =================
!
!   if_driver.f90 
!     - 1D driver
!     - defines state, loops over times, gets convective tendencies, calls radiation
!   module if_conv_solvers.f90 
!     - thermodynamic solvers that don't need arrays
!     - should be able to avoid chunks
!   module if_conv_tend.f90 
!     - analogous to zm_conv
!     - gets all convective tendencies from given state
!     - leaves state unchanged
!     - do this without looping over chunks?
!   compiled f77 radiation code
!
!------------------------------------------------------------------------
!
!  TO RUN IN CAM:
!   (i) comment out this module when running in CAM (include in 1D)
!   (ii) set usecam4 = 1
!   (iii) change nsteps  (TWO things to change)
!   (iv) set bdmassflux non-zero if needed (in if_driver.f90)
!
!------------------------------------------------------------------------
!module ppgrid
!integer, parameter :: pver = 37       !  number of full levels (must be consistent with nd)
!integer, parameter :: pcols = 1       !  dummy
!integer, parameter :: pverp = 1       !  dummy
!end module ppgrid
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
module if_conv_params
!------------------------------------------------------------------------
!
!  Use statements:
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
implicit none
!------------------------------------------------------------------------
!
!  Array Parameters
!  - have to put these array parameters here since needed to define av arrays
!  - numrv, numrl, etc refer to entrain/detrain and are all the same. Could
!    probably use same parameter
!
!------------------------------------------------------------------------
!integer, parameter :: nd = 37      !  number of full levels (redundant if using pver = nd) 
integer, parameter :: nd = pver     !  number of full levels 
integer, parameter :: numrv =  2    !  numrv: rvflow(numrv,nd)
integer, parameter :: numrl =  2    !  numrl: rlflow(numrl,nd)
integer, parameter :: numri =  2    !  numri: riflow(numri,nd)  May be zero
integer, parameter :: numdd = 2     !  numdd: dmflow(numdd,nd)
integer, parameter :: numhm = 2     !  hmflow(numhm,nd)
!------------------------------------------------------------------------
!
!   Convective Switches
!   nupdrafts = 1 implies some near surface layer is convective
!
!------------------------------------------------------------------------
integer :: nupdrafts            ! = 1 if non-zero updraft mass flux in the current timestep
!------------------------------------------------------------------------
!
!  maxlev needs to be here for nunstable
! 
!------------------------------------------------------------------------
integer, parameter :: maxlev = 4    !  Bottom number of layers checked for cape > 0
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
integer, parameter :: nlaunch = 2 !  
!------------------------------------------------------------------------
!
!  Aug 30, 2013: decreased to 6 from 8. Small increase in speed, perhaps 5%.
!  Sept, 2015: decreased to 4 from 5. Little change in model.
!  Feb, 2016: using 3 for BL,CG,DP
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!  Diagnostics:
!
!------------------------------------------------------------------------
integer, parameter :: num_pdf = 12
!------------------------------------------------------------------------
!
!  usecam4
!
!------------------------------------------------------------------------
!integer, parameter :: usecam4 = 0     !   1D
integer, parameter :: usecam4 = 1     !  shuts off diagnostics if in cam4
!------------------------------------------------------------------------
!
!  Array Parameters for tend arrays
!
!------------------------------------------------------------------------
integer, parameter :: numsptend = 3 !  1.updraft + ddraft advection 
                                    !  2.chemical production (ppbv/day) 
                                    !  3: chemical loss timescale (days)
!------------------------------------------------------------------------
!
!  arrays of background state 
!  - should probably be passed to if_conv_tend module through a use if_conv
!    statement
!  - basic state variables are p,km,rt (rv+rl+ri), plus anvil consensate
!  - km conservation is enforced in calls to hydro. This is because p levels are
!    considered fixed. New z and new hm are calculated in hydro. Therefore z and hm
!    are not basic state.
!  - However, hmp conservation is enforced during vertical motions (up/down parcels
!    and during induced subsidence where tendencies are calculated)
!
!------------------------------------------------------------------------
real(r8) :: p(nd)      ! BASIC
real(r8) :: km(nd)     ! BASIC
real(r8) :: rv(nd)     ! BASIC
real(r8) :: rl(nd)     ! BASIC   
real(r8) :: ri(nd)     ! BASIC 
!------------------------------------------------------------------------
!
!  Not Basic
!
!------------------------------------------------------------------------
real(r8) :: rt(nd)     ! not basic state
real(r8) :: t(nd)      ! not a basic state variable (since determined from p/km/rv/rl/ri/rp)
real(r8) :: om(nd)     ! not a basic state variable
real(r8) :: td(nd)     ! not a basic state variable
real(r8) :: z(nd)      ! not a basic state variable  (calculated using hydrostatic balance in hydro)
real(r8) :: hm(nd)     ! not a basic state variable  (from km + z)
real(r8) :: rs(nd)     ! not a basic state variable  (from T/p)
real(r8) :: rs_wat(nd) ! not a basic state variable  (from T/p)
real(r8) :: rh(nd)     ! not a basic state variable  (from rv/rs)
real(r8) :: rh_wat(nd) ! not a basic state variable  (from rv/rs)
real(r8) :: rhod(nd)   ! not a basic state variable  (from pd,T)
real(r8) :: uwind(nd)  
real(r8) :: vwind(nd)
real(r8) :: zh(nd+1)
real(r8) :: phalf(nd+1)
real(r8) :: dpsurf
real(r8) :: dpsurf_sf  ! change in surface pressure after applying surface vapor fluxes (positive)
!------------------------------------------------------------------------
!
!  Diagnostics: Averaging arrays for updrafts or downdrafts
!
!  - I think this is the logical place for arrays which are defined in
!    if_conv_tend.f90 (updraft or downdraft routine), but one wants to
!    access both from the climate model and if_driver.f90
!
!------------------------------------------------------------------------
real(r8) :: av_mp(nd)
real(r8) :: av_rip(nd),av_rlp(nd)
real(r8) :: av_updz,av_upmp
!------------------------------------------------------------------------
!
!   "tend" arrays
!
!   These arrays are calculated in subroutine conv_tendencies of module if_conv_tend
!   They are the main way the driver comminates with the convection. Better to define
!   here.
!
!  rvtend(nd)             Full
!      Water Vapor Mixing Ratio tendency due to advection + convection
!
!  rltend(nd)             Full
!      Background liquid condensate tendency due to advection, convection, precip, evap
!
!  ritend(nd)             Full
!      Background ice condensate tendency due to advection, convection, precip, evap
!
!  Need a flow array for every basic state variable (except height)
!
!------------------------------------------------------------------------
real(r8) :: utend(nd),vtend(nd)
real(r8) :: rvtend(nd)
real(r8) :: rltend(nd)  
real(r8) :: ritend(nd)  
real(r8) :: kmtend(nd)  
!------------------------------------------------------------------------
!
!  For determination of convective instability
!
!------------------------------------------------------------------------
integer :: nunstable(maxlev,nlaunch)   !  set equal to 1 if layer is unstable
!------------------------------------------------------------------------
!
!  For imposed Brewer-Dobson circulation
!  - should be set to zero unless doing longer run, 1D simulations
!
!------------------------------------------------------------------------
real(r8) :: bdflow(nd)
!------------------------------------------------------------------------
!
!  Miscellaneous (needed by convection)
!  days: used by if_conv_tend.f90 for diagnostics
!  nsteps - need this for cam4 use
!
!------------------------------------------------------------------------
integer :: nstep           ! here since needed in other modules - DON'T CHANGE
integer :: nsteps          ! here since needed in other modules - comment out if 1D
!------------------------------------------------------------------------
!
!  Have to be here since used for diagnostic purposes in if_conv_tend.f90
!  - for cam4 comment this out
!
!------------------------------------------------------------------------
!integer, parameter :: nsteps = 20  ! include if 1D; comment out if cam4
integer, parameter :: nav = 10
!integer, parameter :: nav = 4
!------------------------------------------------------------------------
!
!  npoption: choice for vertical resolution
!
!  this is only directly accessed by the driver, but is more convenient to keep
!    where nd defined (shouldn't do any harm to cam4)
!
!    to change vertical resolution must:
!    (i) change nd in conv_array_dim
!    (ii) change nv in mypara.file
!    (iii) change npoption in if_conv module
!    npoption = 1: nd = nv = 37 (default, low resolution)
!    npoption = 2: nd = nv = 72 (medium resolution)
!    npoption = 3: nd = nv = 100 (high resolution)
!
!------------------------------------------------------------------------
integer, parameter :: npoption = 1 
!  integer, parameter :: npoption = 2
!  integer, parameter :: npoption = 3 
!------------------------------------------------------------------------
!
!  end module
!
!------------------------------------------------------------------------
end module if_conv_params
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
