!--------------------------------------------------------------------------
!
!
!  Main Problem: Melting Cooling Depends on Background Atmosphere RH
!  -----------------------------------------------------------------
!  At SPARC sites colrh is consistently too high JJA: too much moisture
!    convergence in NH. Probable origin of monsoon problems
!  Try to get cooling under heating more. Is this a problem with rain
!    not being in phase with up motion in the lower troposphere?
!  l6:
!  - should look at obs event anoms in DJF/MAM, etc to see if differences
!  - h1.2.ps: T anomaly way too cold
!  - h1.5.ps: T anomaly way might be OK: Best? Rain profile really good
!  - h1.6.ps: T anomaly way too small (RH anomaly OK)
!  - size of colrh anomaly does not depend on RH
!  - Not sure of backround T problem or not.
!  - There is also clearly a lag in T anomaly. removing would help
!  - Is this related to monoon problem? Looks like stationary rain max; are rain
!    events not moving since strat cooling too weak?
!  - maybe I should decrease the role of evap in anvil mode but increase f_an so
!    get more melting cooling. This would guarantee stronger 600 hPa cooling ind
!    of colrh.
!  - Perhaps should move h1 domain to monsoon area.
!  - Related problem is that the buoyancies too high in h1.5.ps, when atm backround
!    has higher RH. Unlikely to be true.
!
!
!  Reduce cf_max_ri from 0.4: will this reduce cloud rad heating at high rain?
!
!  NEGATIVE RVNEW is likely most common error message
!
!  LWP Estimate: Fig 10 Woods paper
!  --------------------------------
!  avg of 0.5 g/m3 * 400 m = 200 g/m2 * 0.2 cloud cover = 40 g/m2
!
!  Strength of ML Cooling
!  ----------------------
!  - It is extremely strong (too cold) when background atmosphere is dry, 
!    and gets weak when atmosphere has high RH. 
!  - Probably need some mechanism to increase f_an at higher colrh.
!  - It is also possible that these anomaly patterns depend on the 
!    background RH.
!
!  Constant LWC (suggested by Rosenfield)
!  --------------------------------------
!  - could not get this to work
!  - However I clearly have too high LWC at low levels, and this
!    probably makes it too hard to produce warm rain.
!  - A constant value produced more disperse precipitation over land.
!
!  Problem with use cam_history in if_conv.F90 (Feb 2017)
!  ------------------------------------------------------
!  - I did a rm -r caml1 to totally remove an old run
!  - then when tried to build l1 again, got a number of problems
!  - May have been some compiler changes since last ran model.
!  - Shawn suggests:
!  - go to bld directory where routines are compiled.
!    e.g. /home/folkins/scratch/caml1/bld
!    Are cam_history.mod and cam_history.o there?
!    If not copy from somewhere else.
!    Go to another bld directory eg cams1/bld and copy over
!    cp cam_history.o ../../caml1/bld/cam_history.o
!    cp cam_history.mod ../../caml1/bld/cam_history.mod
!    cp histFileMod.o ../../caml1/bld/histFileMod.o (likely not needed)
!    Try again
!
!
!  Try a situation where IWC does not decrease as go toward the melting level
!
!  10741/10741: rv_mass_avail
!  8253/8239: rlp_rv_max = rs(ii)*(dmdry/mpp)*(rh_max_cond - rh(ii))
!  10806/10806: fraction_evap = mass_anrain_evap/rv_mass_avail
!  use_om5_a
!
!  Changing SST
!  ------------
!
!   Default sst file is:
!  /home/folkins/cesm1_0/inputdata/atm/cam/sst/sst_HadOIBl_bc_1.9x2.5_clim_c061031.nc
!  is 2.6 MB
!  I renamed this file default.nc and replaced with AMIP_4K.nc
!
!
!
!  To Do:
!  ------
!  - output KMTEND in h0.f90 and compare with diff downdrafts
!  - Lagrangian perspective on hm: if no shallow convection or
!    downdrafts hm would be a minimum above the BL, as rad cooling
!    decreased it.
!
!
!   Possible Origins of Monsoon/Hot Upper Trop Issues
!   -------------------------------------------------
!   - Upper Trop Too Hot JJA and too much rain in SE Asia
!   (1) Radiation: Cloud Radiative Feedbacks too strong
!   (2) Activation Dynamics: Lack of small spatial scale dynamical variability
!       - that might arise from small scale SST gradients as opposed to
!         highly smoothed
!   (3) Convection: Water Vapor Feedbacks too strong 
!     - colrh too high
!     - from h1.ps it looks like melting cold anomaly is too weak
!     - amp problem
!     - RH based detrainment problem
!     - mode 1/mode 2 problem
!
!  Setting zeros in hb_diff.F90: kvm and kvh to free values; cgh cgs = 0.
!  - but dtv_rain still not zero
!    cgh: counter-gradient term for heat [J/kg/m
!    cgs: counter-gradient star (cg/flux)
!    kvm: eddy diffusivity for momentum [m2/s]
!    kvh: eddy diffusivity for heat [m2/s]
!    kvq: eddy diffusivity for constituents [m2/s]
! 
!   Other outputs of hb_diff are: 
!    tpert 
!    qpert
!    ustar
!    ktopbl
!    ktopblmn
!    tke
! 
!  The routine where the diffusivities are used to calculate the
!     updated T/q is subroutine compute_vdiff in diffusion_solver.F90
! 
!  Removing HB BL in convective regions:
!  - likely the easiest place to do is physics_types.F90 
!  - Look for the diffusion call.
!  - have the latitude. Do not zero tendencies at the surface just at the levels.
! 
!  Why is RH at 16 km so high?  Solar heating at 16 km too high?
!  - have to limit detrainment
!
!  Looks like BL cooling caused by BL scheme - dtv_colrh: reduce?
!  - or do I need a BL mode?
!
!  - Getting incredibly strong buoyancies below 800 hPa: need more 
!    at low rain rates, and more mass flux. 
!  - remove layers with higher RH?
!
!  Increasing Variance/Sensitivity to RH
!  -------------------------------------
!  - decrease amp_low or amp_add
!  - increase f_det_amp, rh_f_det
!
!  Why two Modes?
!  --------------
!  - Keeping the lower troposphere warm enough seems to require 
!    disperse weak shallow convection (with detraining hm larger
!    than background).
!  - To maintain the positive buoyancy of the deep mode, and prevent 
!    excessive detrainment in the BL or mid-troposphere, this heating
!    should be spatially disjoint from deep convection, and the 
!    moistening of the atmosphere should mainly come from evaporation
!    rather than detrainment. For disperse shallow convection, the
!    main moisture source should likely be detrainment.
!  - To maintain the humidity of the mid-troposphere, it also seems 
!    neccessary to have some deep convection at low rain rates.
!  - Likely two two modes to achieve this: easier to define modes
!    with different properties, and modify relative amplitudes than
!    to have rain dependent properties.
!
!  Mass Flux closure - more natural approach?
!  -------------------------------------------
!  - by enforcing a constraint between rain and mass flux, effectively enforcing a 
!    "rainfall efficiency = net surface rain/starting mass flux". The efficiency
!    tends to decrease at high rain rates. This is likely because the mass flux
!    at high rain rates is being pushed to a "too large" value, i.e. The increase
!    in rain rate is being driven by the larger mass flux, rather than the rain
!    rate increasing because conditions (i.e. colrh) become more favorable.
! - is there a more naturalistic or less constrained way of doing this?
! - one method might be too impose an efficiency as a function of rain
!
!  2 km Stability Max not strong enough
!  ------------------------------------
!  - Dec 2016: got strong 2 km LR peak by increasing rate_down_uprain 
!    from 0.05 to 0.10.
!  - problem came from weighting rain to larger values
!  - Does not come from downdrafts; I think simply comes
!    from not having enough shallow heating. So if rain rates
!    higher shallow heating weakened, and 2 km max disappears.
!
!  Leaky Pipe Model
!  ----------------
!  - colrh sensitivity through mixing of the entire plume does not work.
!    To get the required colrh sensitivity, you have to entrain a lot; the
!    hm of the plume goes down too much, and the lower troposphere becomes
!    too cold.
!  - However, some plume mixing is required for updraft buoyancies to retain 
!    realistic values.
!  - I should actually calculate the additional heating that comes from having
!    the mean detraining hm at some positive bias with respect to the background.
!  - To get the required sensitivity of congestus to background RH, assume
!    that mixing on the edges of the cloud creates negative, or neutrally
!    buoyant mixtures, and that these then detrain.
! 
! 
!   Prescribed Mixing (from Buoyancy) vs Prescribed Bouyancy 
!   --------------------------------------------------------
!   - benefit of prescribed buoyancy: get better Deep Outflow Mode
!     This was a consistent problem with other kinds of mixing.
!     It must be closer to reality: plumes incrementally mix until the
!     buoyancy gets weaker, and mixing stops. Prescribing a particular
!     mixing value from a given starting buoyancy must be bad.
! 
!   Forced Detrainment for positive buoyancy?
!   ----------------------------------------
!   Possible that shallow detrainment detrains with a positive buoyancy, 
!     that it is not able to overcome dissipation, or drive a circulation.
! 
!  Average downward movement 
!  Int(p(det-ent))/Int(ent)
! 
!  Check radsw.F90, radiation.F90, pkg_cldoptics.F90, param_cldoptics.F90, 
!    cloud_fraction.F90 , cldwat.F90, radheat.F90
!
!   Lack of Dynamics (Aug 2016):
!   ---------------------------- 
!   - Currently model has poor dynamics: no MJO, Kelvin, etc
!   - This is hard to understand: my model is working very well in almost
!     every other respect.
!   - In Western Eq Pacific, my sea breeze circulations are too strong. Too
!     much rain over land. There is more rain over land, but this is too strong
!     in the model. It is possible that this is making the dynamics worse, but
!     not likely the only reason. It is probably a reflection of the rain to
!     not move enough.
!   - I had thought that the reason was because my cloud ice was too focused
!     over rainy regions, and replaced ri_rh with ri_om (coupled to omega).
!     But omega itself is quite localized with the rain so this does not help
!     spread out the cloud radiative heating. In any case the cloud heating 
!     at 227 hPa does not itself seem that localized, and in h2.ps my TOA LW
!     vs rain relationships are in reasonable agreement with CERES.
!   - Reducing a3_add from 6 to 3 increases rainfall variance. Concurrently,
!     appears also helpful to increase f_an_half, so that melting cooling does
!     not start reducing colrh too early. Also likely want lots of convective
!     downdrafts to drive down surface MSE, so convection forced to move 
!     somwhere else.
!   - It may help to have more radiative heating at the base of the ML, as a way
!     of triggering more stratiform waves. Stratiform waves may also be triggered
!     by penetrative downdrafts that extend the depth of the troposphere below
!     the ML.
! 
!   Maybe:
!   - Convective downdrafts reduce subsidence, so cool and moisten
!     the column, and increase colrh. Reduce downward advection of
!     low hm air.
!   - ML downdrafts generate the mid-level cooling max which contributes
!     to the mid-level GPHT anomaly which creates inward PGA and downward
!     motion which lowers colrh, and suppresses convective activity.
!
!  July 2016
!  ---------
!  - ECMWF: LTS correlates much better with LO CF than CE. Can't base
!    cloud scheme on CE. Should convert what I have to LTS or EIS.
!
!  Triangular Relationship
!  -----------------------
!  - should be able to model if this relationship is correctly
!    simulated.
!    CF - CE
!      \ /
!      LWP
!  - I am enforcing two of these connections: is third OK?
!
!
!  Strategy for semi-empirical model
!  ---------------------------------
! (1) Use OBS correlation with CE to get a "target" CF: cf_tar
! (2) Use CF:LWP relation to get target LWP: lwp_tar
! (3) Distribute lwp_tar in layer with highest RH
! (4) Assign cf_tar
! (5) Makes RH profile more realistic?
!  - what are model cm4:rh relations?
!
!  Advantages of Empirical Model/Publishable Arguments
!  ---------------------------------------------------
!  (i) uses two simple relations strongly constrained by observations
!  (ii) Does putting rl in one layer improve BL structure?
!   - Doing this is strongly supported by observations, since clouds
!     only about 200 - 300 m thick
!   - This might be what default CAM4 did already
!   - Likely helps generate correct relationship between albedo and LWP,
!     since if rl placed throughout BL, would likely generate low small
!     LO CF for given LWP and too small albedo
!  (iii) generates correct relationship with SST (fix low SST LWP problem?)
!      argue: SST dependence is a response to modulation of CAPE by SST
!      rather than an "instrinsic" effect of SST.
!  (iv) Look at Taylor diagram - correlation and deviation (individual months)
!  (v) Predicts high CE/low CE "split" in ALB vs LWP relationship
!      (Where does this come from? - regional variation?)
!      A bit hard to understand in the absence of a split in CF vs LWP
!  (vi) What is effect of cloud droplet diameter? Get better agreement
!     without any special dependence? Where does small degree of MODIS 
!     variability come from? Simplest assumption consistent with MODIS?
!     Does CAM4 use of a smaller value compensate for too smooth distribution
!     of rl in the vertical?
!
!
!  Important relationships:
!  ------------------------
!  22: CF vs CE
!  25: ALB scatterplot
!  27: CF scatterplot
!  29: LWP scatterplot
!  34: CF vs SST
!  36: LWP vs SST (now bad)
!  64: CF vs LWP
!  67: LWP vs CE (now bad)
!  73: ALB:ALB Correl
!  74: CF:CF Correl
!  75: LWP:LWP Correl
!  78: Mean ALB:ALB
!  79: Mean CF:CF
!  80: Mean LWP:LWP
!  95: ALB ERR
!  96: CF ERR
!  97: LWP ERR
!
!
!  BL Mixing Routines
!  diffusion_solver.F90: don't see any changes here no IF or IAF
!  hb_diff.F90
!
!  Reducing latent heat flux over the ocean:
!  /home/folkins/cesm1_0/models/csm_share/shr/shr_flux_mod.F90
!
!  excess LH flux problem: qneg4.F90 is where the error is written out
!  I have introduced changes to qdiff and to the BL mixing scheme,
!  The qdiff_factor was used in 
!    /home/folkins/cesm1_0/models/csm_share/shr/shr_flux_mod.F90
!  I have removed my changes (qdiff)
!  I think the water flux is called evap in shr_flux_mod.F90
!  In the "excess" calculation, uses qflx(i,1)
!  But for atm/ocn presumably the same as called "evap" in:
!    /home/folkins/cesm1_0/models/csm_share/shr/shr_flux_mod.F90
!
!  Ideas for seg fault
!  -------------------
!  - Check if March 7 had seg fault?
!  - Introduction of cloud continuum functions: used in March 7 archive
!  - Check for OK ts in:
!    vi /home/folkins/cesm1_0/models/csm_share/shr/shr_flux_mod.F90 DID
!  - Look for wierd surface temperatures/pressures
!    surface T called cam_in%ts
!    in stratiform_tend called ts(pcols) 
!  - Look for negative temperatures: doing this for air temperatures (DID)
!  - They suggest running with "state_debug_checks = .true." in the namelist (DID)
!  - Or run in DEBUG mode  (DID)
!  - restore old versions of files associated with BL mixing/radiation
!  - It seems to be associated with warnings about Max possible LH flux
!    exceeded.
!  - In radiation.F90 : added FUS and FDS by oct 2015; can't be the problem.
!    commented out all ozone stuff
!  - Identify better in the dynamics where occurring  HARD
!  - I have called a function a subroutine?
!  - drv_in problem?
!  - link with -g -traceback also
!  - adding -heap-arrays to ifort gets a run killed.
!
!  Good Discussion Seg Fault: 
!  https://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors
!  (3) user coding error
!  (4) exceed array bound
!  (5) calling function as subroutine
!
!  Error:
! (t_initf) Read in prof_inparm namelist from: drv_in
! utr_mpiwtime: not enabled
! GPTLinitialize: failure initializing MPI_Wtime: reverting underlying timer to gettimeofday
! forrtl: No such file or directory
! forrtl: severe (29): file not found, unit 99, file /home/folkins/scratch/cams2/bld/drv_in
!
!
!  Pointers:
!  ---------
!  preciporg 
! call pbuf_add('PRECIPORG' , 'global', 1,1,     1, idx)
!
!  if_conv_intr.F90
!  convect_deep.F90
!
!  removed precip_for_update, NDEEPORG, PRECIPPREC
!
!  Segmentation Fault
!  ------------------
!  Look at web site:
!  tp_core an issue with modifying radiation.F90 (to modify cloud ...
!  https://bb.cgd.ucar.edu/issue-modifying-radiationf90-modify-cloud-radi...
!  Jul 6, 2014 - The model configuration is CAM4 coupled with an aqua-planet 
!  slab ocean ..... The crash is in tp_core, part of the FV dycore, so 
!  your changes to ...
!
!  Big Clue:
!  IAF comment entering subroutine dyn_run in fv dyn_comp.F90
!  forrtl: severe (408): fort: (2): Subscript #2 of the array Q has value 101 
!  which is greater than the upper bound of 98
!  j = 81
!  in if
!  forrtl: severe (174): SIGSEGV, segmentation fault occurred
!  IAF comment entering subroutine xtpv in tp_core.F90
!  j =           81
!  in else
!  in abs
!  forrtl: severe (408): fort: (3): Subscript #1 of the array QTMPV has 
!    value -60 which is less than the lower bound of -48
!
!  Wed: 
!  recursive search : grep -Rin "call wshist" *
!  recursive search : grep -Rin "call cam_run4" *
!  recursive search : grep -Rin "subroutine cam_run1" *
!  recursive search : grep -Rin "subroutine stepon_run1" *
!  recursive search : grep -Rin "subroutine stepon_run1" *
!  recursive search : grep -Rin "subroutine dyn_run" *
!  IAF comment calling tp2c in sw_core.F90
!  vi /home/folkins/cesm1_0/models/atm/cam/src/dynamics/fv/tp_core.F90
!  vi /home/folkins/cesm1_0/models/atm/cam/src/dynamics/fv/sw_core.F90
!  vi /home/folkins/cesm1_0/models/atm/cam/src/dynamics/fv/cd_core.F90
!  vi /home/folkins/cesm1_0/models/atm/cam/src/dynamics/fv/dyn_comp.F90
!  vi /home/folkins/cesm1_0/models/atm/cam/src/dynamics/fv/stepon.F90
!  vi /home/folkins/cesm1_0/models/atm/cam/src/control/cam_comp.F90
!  2 places with "call cam_run4"
!  /home/folkins/cesm1_0/models/atm/cam/src/cpl_mct/atm_comp_mct.F90
!  /home/folkins/cesm1_0/models/atm/cam
!  4 places with "call wshist"
! (i) vi /home/folkins/cesm1_0/models/atm/cam/src/control/cam_history.F90
! (ii)  vi /home/folkins/cesm1_0/models/atm/cam/src/control/cam_comp.F90
! (iii) cd /home/folkins/cesm1_0/models/atm/cam/tools/rad_driver/.svn/text-base
!       vi ccsm_driver.F90.svn-base (But read only file)
! (iv) vi /home/folkins/cesm1_0/models/atm/cam/tools/rad_driver/ccsm_driver.F90
!  added print*,'exiting subroutine wshist in cam_history.F90' to cam_history.F90
!
!  Check for functions called as procedures:
!  add "-gen-interfaces -warn interfaces" to compile
!
!  Coding error:
!  Compile and link using the ifort driver and these options:
!  -g traceback
!
!  Using gdb
!  ---------
!  cd  /home/folkins/scratch/camxx/bld
!  gdb cam
!  run
!  bt
!  - problem is no stack
!  - likely go to run directory and use gdb.
!  - run directory = $wrkdir/$case
!    /home/$LOGNAME/scratch/camxx
!    /home/folkins/scratch/camxx
!  - added debug to xx.csh
!  Helpful websites:
!  http://unix.stackexchange.com/questions/132192/running-application-ends-with-segmentation-fault
!
!  BalanceCheck: soil balance error nstep =      1834 point =  7274 
!  imbalance =   -0.000005 W/m2
!  clm2: completed timestep         1834
!  QNEG4 WARNING from TPHYSAC  Max possible LH flx exceeded at    
!  1 points. , Worst excess =  -1.4067E-03, lchnk = 697, i =    8, 
!  same as indices lat =  74, lon =  97
!  QNEG4 WARNING from TPHYSAC  Max possible LH flx exceeded at    1 
!  points. , Worst excess =  -3.2901E-04, lchnk = 698, i =    8, 
!  same as indices lat =  74, lon =  98
!  WSHIST: writing time sample 834 to h-file 4 DATE=2000/02/08 NCSEC= 18000
!  forrtl: severe (174): SIGSEGV, segmentation fault occurred
!
!  - look in  radsw.F90 radheat.F90
!
!
!  TOA SW over high rain regions (April 2016)
!  -------------------------------------------
!  - How to generate curvature in TOA albedo versus TOA LW relationship? 
!    - crucial for generating better TOA SW vs TOA LW ratios, due to satauration
!      in TOA LW
!    - likely need very high IWP, and smaller ri_half, large cloud fractions
!      at high rain rates, and more frequent occurrence of high rain rates.
!  - plot time dependence of TOA SW and TOA LW in h1.gri.
!
!  MAIN PROBLEM (Aug 2016): 
!  - Need more moisture resistance: which mode is reducing moisture resistance?
!  - Add back in:
!  - It is possible that that having more penetrative downdrafts (relative to 
!    evaporation) increases the temperature of the lower troposphere.
!
!  2 km Stability Max From Penetrative Downdrafts
!  ----------------------------------------------
!  - Was hard to escape from idea that that the stability max was generated
!    by subsidence + BL convection. But this can't be true because the lapse
!    rate be low 2 km sharply increases, and isn't close to adiabatic at all,
!    as it is in the stratocumulous regime. Also hard to escape from the idea
!    that there wasn't some kind of continuity or common origin of the lapse
!    rate maxima in the weakly convective (1 km) and convective regimes (2 km),
!    that it is some kind of response to the increase in SST and more vigorous
!    BL convection. 
!
!  Cold Lower Troposphere Problem
!  -----------------------------
!  - lower troposphere appears to warm up as use more penetrative downdrafts
!    (and reduce or elimate pure evaporation). This is counterintuitive.
!    Possible: explanation warm by reducing subsidence (since convection scheme
!    is required to keep mass between pressure levels constant.)
!    This directly cools and moistens, but the net is increase in hm, so
!    when convection detrains, cooling from condensate evaporation is reduced.
!    Effectively, it allows one to tune the model to achieve a correct RH
!    profile by reducing evaporation.
!  - Make sure that parcels detrain with a larger positive bias in HM. I.e.
!    force to detrain with a positive buoyancy.
!
!  HOW TO WARM THE TROPOSPHERE
!  ----------------------------
!  - allow rain events to develop faster. This was a benefit of my
!    extrapolated mass flux closure, and I was able to generate a good
!    temperature profile. However, I had trouble maintaining
!    intermediate values of rainfall. And a cape dependence.
!
!   2 km LR max
!   ------------
!   - Very important for T of lower troposphere
!   - Associated with max in downdraft mass flux. If downdraft cooling
!     scales with mass flux, then cooling max would be location of max
!     stability.
!   - Consistent with minimum in TPTEND_REAL at 2 km (i.e. minimum 
!     warming)
!   - But what causes this max downdraft mass flux? I think mainly that
!     Melting downdrafts are only allowed to descend two levels.
!   - Can be lifted up by a decrease in t_melt.  
!   - Is the obs LR starting from the ML a saturated downward neutral 
!     buoyancy adiabat?
!
! Knobs for Moving Rain to higher colrh:
! --------------------------------------
! (1) km_width? NO
! (2) amplitudes: e.g. reduce a3_add, or increase a3_half
! (3) Lower b_tar
! (4) f_an: increase (seemed to have reached limit here)
! (5) increase rh_cond (but this may cool the lower trop)
! (6) Increase downdrafts 
!  rate_down
!  rh_down_max
!  b_precip_down - may be imp for pushing colrh to high values?
! (7) Increase melting
!   b_precip_melt
!
!  -------------
!
!  Why do I have no MJO? (Aug 2016)
!  --------------------------------
!  (1) Add momentum transport (not likely)
!  (2) Too strong rain over land of Maritime continent
!  (3) Upper trop LW heating not distributed evenly.
!      (make ice proportional to upward motion and RH?)
!  (4) If change lat_cloud: changes tropical dynamics?
!
!
!  Problems:
!  ---------
!  (1) a bit more colrh sensitivity
!  (2) TOA LW over oceans too low
!  (3) upper troposphere too warm
!  (4) 3-5 km too stable
!  (5) Too much rain over Borneo?
!  (6) 2 km LR max too high? 
!  (7) Upper trop RH too high?
!  (8) West Eq Pacific does not have low enough TOA LW
!
!  To Try:
!  -------
!  (1) km_add = 0.
!  (2) Constant reduction in rh_min_trop, etc
!  (3) ri_dt more sensitive to RH? (not just detrainment)
!
!  TOA LW too low (upper trop too warm):
!  ------------------------------------
!  - may have to decrease cloud fraction (increase ri_add
!    or lower t_ri_half = 205.0)
!    and increase ri_base_o = 9.0*0.001
!                 ri_base_l = 9.0*0.001
!
!  ------------------------------------------------
!
!  Current Problems: 
!  -----------------
!  - Only allow downdraft to go down 1 level at a time. Do a search on
!    i_lowest.
!  - increase amp_low = 0.1: better MJO
!  - Think about enhanced det at level 4 (seems to help)
!  - perhaps increase rh_ansnow/rh_uprain from 0.89 to 0.90
!  - can modify lat_rh_wat
!  - Why am I getting a peak at 55 mm/day? Maybe too easy to grow fast.
!    Need more RH resistance?  f_det_rh? reduce b_tar
!  - Melting level RH likely too big: reduce rh_melt_target?
!  - possibly increase md_rat_uprain to 800
!  - Is upper trop too hot?
!  - RH too high in JJA?
!  - or reduce rh_cond_min and rh_cond_max
!  - Set rate_evap_surf = 0.10 TRY LATER: DOES NOT WORK NOW
!    It is possible I need to increase surface RH at high rain to reduce
!    latent heat flux at high rain and give a negative feedback to reduce
!    intensity of monsoon. Cloud rad feedbacks another possibility.
!    (But watch doesn't lose colrh variance during events)
!  - Seem to be good things to do: try later
!    But don't solve max rain at 45 mm/day
!    a1_add = -0.95
!    mass_inc_max = 4.0 
!  - inc f_an_add = 1.2
!
!  ----------------
!
!
!  s4+l4:  
!   use_f_an_om = 1
!
!  s1+l1:  
!   f_evap_melt_max = 0.0 (Otherwise same as s4)
!   t_melt = 270. (from 269.)
!   b_precip_down = -0.02
!   rh_melt_target = 0.84 (DEF = 0.90)
!   f_an_half = 60 (DEF = 75)
!
!
!  l2+s2:  l6 + DOES NOT TURN OFF DOWNDRAFTS   SHAWN LIKES
!     b_precip_down = 0.20  KILLED at 28 and 19 MON
!    Doing this removes the 2 km stability max.
!    Tends to give create inc temp from 800 hPa to surface; too warm surface
!    Dynamics very different: in Fig3.Asym.IF_RAIN.pdf don't get the big MRG
!      peak at -3, Asym better, but have zero MJO.
!    Better correlations over land, better monsoon.
!
!  l6+s6:  b_tar = 0   SENT TO SHAWN MARCH 15  NEW DEF
!    b_tar = 0 (both from 0.008)
!    b_det(1) = 0.01 (DEF)
!    b_precip_down = -0.02 (DEF)
!    b_precip_melt = -0.02 (DEF)
!    - h1.1.ps better than h1.1.ps of l3
!    - h1.2.ps much better than h1.2.ps of l3
!
!  l3+s3: lwc_max = 3.0  DONE  DEFAULT: sent to Shawn
!    rh_min_ttl = 0.85  (from 0.90)
!    ri_base_o = 3.0*0.001
!    ri_base_l = 4.0*0.001
!
!  -------------------------
!
!  l4+s4: lwc_max = 1.0  LIKELY BEST Run to 34/50
!    a1_om5_scale = 30.0  
!    a1_om5_half = 100.0
!
!  l2+s2: l4+ reduce ice  NEW DEF (could be a bit cold)
!   ri_base_o = 3.0*0.001  (from l4 value 5)
!   ri_base_l = 4.0*0.001  (from l4 value 6)
!   amp_low = 0.0  (l4 value)
!   amp_add = 2.5  (l4 value)
!   f_an_add = 1.0 (l4 value)
!   rh_melt_target = 0.90 (l4 value)
!   rh_down_max = 0.95 (l4 value)
!   rate_down_anrain = 0.10 (l4 value)
!   a1_om5_half = 100 (l4 value)
!   a1_om5_scale = 30.0 (l4 value)
!
!  l4: KILLED at 35 months. NEW DEF. Lots more Kelvin waves
!   use_om5_a = 1
!
!  -------  no om5 --------------
!
!  l3+s3: NEW DEF?  KILLED at 14000 and NUM = 9
!    storing OM5
!    use_km_width_land = 1  (ONLY DIFF from l2)
!    rh_min_ttl = 0.85 (DEF)
!    rh_max_ttl = 1.00  (DEF)
!    outflow_power = 1.0
!    ri_base_o = 5.0*0.001 
!    ri_base_l = 6.0*0.001  
!    iwp_add = 0.0
!
!  s2(6000)+l2: NEW DEF 
!   outflow_power = 1.0
!   ri_base_o = 5.0*0.001  
!   ri_base_l = 6.0*0.001
!   ri_add = 3.5E-05
!   rl_rh_min = 0.85
!   rl_rh_base = 0.0005
!
!  ----------------
!
!  s1(6000)+l1: h1 looks good  NEW DEF
!   mass_inc_max = 5.0
!    outflow_power = 1.3
!    ri_base_o = 7.0*0.001
!    ri_base_l = 7.0*0.001
!
!  s2(6000)+l2: NEW DEF
!    rh_cond_min_ft = 0.40 (from 0.50)
!    rh_cond_max_ft = 0.65 (from 0.75)
!    rh_cond_min_bl = 0.55   ! i <= 5
!    rh_cond_max_bl = 0.80   ! i <= 5
!    rate_evap_surf = 0.0
!    A2
!
!   -------------
!
!  s3(6000)+l3: l1 + 
!   rate_evap_surf = 0.10
!
!  s1(6000)+l1: NEW DEF
!     MJO????
!     LOOKS LIKE RH IS TOO HIGH SO COULD GET WARMER ATM
!     BY REDUCING rh_cond_max; 
!   A2 from A1
!   mass_inc_max = 5.0
!   mass_rat_max = 0.5
!   md_rat_uprain = 1000.
!   rate_down_uprain = 0.05
!
!
!  =========================
!
!
!  l2: rl_tar = 0 for i > 5.
!      rl_o_base = 0.0
!      i_cf >= 3
! 
!  l1: lwp_amp = 30 from 40
!      cff_mult = 1.2
!
!  l1: new BL param (search qqq)
!  - too much LWP?
!  - CF too small
!
!  ---------
!
!  l1: 
!  - inc rl_rh_base (have to)
!  - ce_add = -0.80 (from -0.75)
!
!  --------------------
!
!    (i) Need to have a way of making SRF-500 LIQ DEC with rain rate
!    (ii) Need to have a way of increasing the reflectivity of this LIQ
!  Increasing cloud fractions would help a bit, but have to solve (i) first.
!  Increasing LIQ Water higher than CERES:
!    Blossey has 0.01 g/kg for 900-500 hPa: since this layer is 4000 kg/m2,
!    this corresponds to 40 g/m2 water. This is WAY more than CERES.
!  One reason for high water contents at low rain is that high liq water
!    likes more stable lapse rates.
!
!  TOA SW too small + too much rain over land
!  ------------------------------------------
!  - reduce ice crystal size by rei_reduce: not much effect
!  - increase ri_half and increase ice
!  - gci_IAF = 0.78: not much effect
!  - increase height of ice so get smaller particles: z_peak: not much effect
!
!  TOO LITTLE VARIANCE/Weak MJO
!  ----------------------------
!  - colrh distribution is OK
!  - But not enough rain at high rain rates, and rain should be shifted to higher colrh
!  - too weak stratiform heating. 
!  - May be related to too much LW heating at low rain rates.
!  - iwp may increase too slowly with rain.
!
!  RH OF LEVEL 4 TOO LOW
!  ----------------------
!
!  LESSON: MAIN CONVECTIVE OUTFLOW MUST BE > 12 km
!  -----------------------------------------------
!  - Otherwise, the transition to a balance between LS Ascent and Rad 
!    heating occurs too low and you will not get a MALR.
!
!  LESSON: SURFACE RH SHOULD BE CLOSE TO 0.80 EVERYWHERE IN THE TROPICS.
!  ---------------------------------------------------------------------
!  - otherwise LHFLX is too big, and cape at low rain rates is too small.
!  - But playing with the BL scheme is dangerous ....
!
!  Increasing Surface RH at low Rain Rate (< 1 mm/day):
!  ---------------------------------------------------
!  - this is a regime where convective drying has no effect on surface RH,
!    so nothing in the convective scheme itself can change.
!
!  LESSON: DO NOT PAY TOO MUCH ATTENTION TO CERES
!  ----------------------------------------------
!  - pay more attention to physical requirements in getting a good T profile.
!  - CERES IWP is way too high.
!
!  To Do:
!  - try and find out what determines rad_lwsw
!  - try again changes to asymmetry fraction?
!  - need to increase SW over land?
!  - need to increase TOA LW over ocean?
!  - plot RH of layers 1-3
!  - correlations with cloud fraction TOA LW/SW
!
!  How often is field archived?  
!  Where is rad_lwsw defined?
!  in   /home/folkins/cesm1_0/models/atm/cam/src/utils/cam_pio_utils.F90
!  type, public :: field_info
!    character(len=max_chars) :: sampling_seq     
! sampling sequence - if not every timestep, how often field is sampled
! (i.e., how often "outfld" is called):  every other; only during LW/SW
! radiation calcs; etc.
!
!  -----------
!
!  =================
!
!  ICE Paper:
!  ----------
!  Comparisons of global cloud ice from MLS, CloudSat, and correlative data sets
!  Figure 4: largest IWP around 30 g/m2.
!  "CloudSat mean IWC is generally 3–5 times greater than MLS ones, but both are 
!   greater than the ECMWF mean."
!
!  qdiff_factor = 0.85 
!  -------------------
!  - PRECT over OCEAN stable at bit less than 4
!  - CAPE_1 at WP sites collapsing to zero due to dec in BL RH and inc 
!    T of atm. Aloft becoming extremely dry. Rainfall collapses to equator.
!  - I have to look at trends in TOA LW and SW and compare with obs
!  - I am sustantially reducing the heating of the atmosphere by rain: if
!    my atm T does not go down as much as it should, I must be heating more
!    and more by radiation? Atmosphere losing ability to cool?
!
!
!  The high LHFLX problem is a  LOW rain rates/LOW COLRH problem. That is
!    when LHFLX is high (> 200), and where most of the area is.
!  Here the BL is too unstable and RH at surface too dry - moisture is
!    mixed off the surface too quickly.
!
!  Problems (July 9, 2015):
!  (i) TTL too cold
!  (ii) Too high BL in Caribean  SOLVED
!  (iii) Too high rainfall over land
!  (iv) ALB_MD is too small
!  (v) Too much cloud fraction at low rain 
!    - strange problem with ML cloud fraction at low cldice
!  (vi) Strange behaviour in 7sw_dn 
!    - reflects problem in rain diurnal cycle where rain 
!      occurs too close to solar noon.
!  (vii) Monthly LWP not peak at low value
!  (viii) Strongly negative LW cooling near 14 km: maybe
!     need more variation in the height of this peak, with
!     higher values at higher rain rates.
!
!
!  HI LW cooling too low 
!  ---------------------
!    How can I reduce my HI LW heating?
!    Mainly I think by reducing 200 hPa UP LW
!    500 hPa UP LW is good, so need more LW absorption in 500-200
!    Higher Cloud Fraction in this layer
!
!
!  ===============
!
!  EXTREME (extreme): rainfall becomes locked in place esp over land/mountains
!  ------------------------------------------------------------------------
!  - caused by too warm T, caused mainly by LHFLX
!    too high over the ocean. Need to increase surface rv.
!  - probably caused a bit also by albedo not high enough over land, esp LOW
!    albedo.
!  - As soon as rainfall over land has a positive bias, it produces a positive
!    bias in colrh, and any nonlinearity in mixing amplifies this bias.
!
!  Temperature Issues
!  ------------------
!  "Hot atmosphere problem" 
!   - occurs as increase IWP, and decrease atmospheric cooling.
!   The atmosphere heats up, and due to constant
!   SST's, convection occurs preferentially over the land. To combat, have to
!   find some way of decreasing the latent heat flux, which means increasing
!   the RH (inc rv) over the ocean surface layer.
!  Lower Atmospheric Temperatures (700 - 500 hPa)
!   - tend to be too warm by 0.5 - 1 K
!   - may be due to cold radiative heating: LW heating heating likely increased
!     due to overhead clouds.
!
!  rapid increase in layer 1 cape over land at sunrise
!  ---------------------------------------------
!  - originates from sudden increase in QSURF at 9:00 of 0.001 kg/kg
!  - at this time LHFLX = 200 W/m2 = 8E-05 kg/m2/s = 0.288 kg/m2/hour water vapor
!  - the bottom layer is 20 hPa or 200 kg/m2
!  - To get increase in this layer of 0.001 kg/kg, need to add .2 kg/m2 of 
!    water vapor
!  - these numbers are consistent: the modeled LHFLX could cause the spike in
!    QSURF. However, since the QSURF spike likely does not occur, water vapor must
!    get mixed up into the BL. Why is this not happening? BL mixing not vigorous
!    enough?
!
!
!  Surface Heat Fluxes (surface sensible heat fluxes)
!  --------------------------------------------------
!  - DON'T BOTHER: NOT THE PROBLEM
!  /home/folkins/cesm1_0/models/csm_share/shr/shr_flux_mod.F90
!  - This is where I use qdiff_factor
!  - However fluxes calculated here do not seem to be applied to state vector
!  - look for the call to shr_flux_atmOcn: this would be the subroutine where
!    the ocean/atm heat fluxes calculated
!  - recursive search to all lower directories: grep -Rin shr_flux_atmOcn *
!  - shr_flux_atmocn is called from seq_flux_mct.F90 (but can't figure out where
!    this file is).
!  - actually this just seems like atm/ocn and atm/ice
!  - The land fluxes are much more complicated: 
!    Somewhere in: /home/folkins/cesm1_0/models/lnd
!    Might be: array eflx_sh_tot(:) in 
!    /home/folkins/cesm1_0/models/lnd/clm/src/main/clm_atmlnd.F90
!    maybe "call CanopyFluxes" or "BareGroundFluxes"
!  I think a lot of the heat fluxes over land are calculated in:
!  /home/folkins/cesm1_0/models/lnd/clm/src/biogeophys/CanopyFluxesMod.F90
!  - From the description, SH and LH calculated from Monin-Obukhov similarity
!    theory (starting page 157 of the description).
!
!  How to Modify Sensible Heat Fluxes Over Land
!  --------------------------------------------
!  - DON'T BOTHER: NOT THE PROBLEM
!  - natural place to do may be in subroutine smooth in :
!  /home/folkins/cesm1_0/models/atm/cam/src/physics/cam-if/flux_avg.F90
!  would need to know when shflx is being called, know lat, and ocnfrac
!
!  Need more MH Cloud fraction at low rain: need more ice_mh at low rain
!  Need more heating ner 14-15 km to remove cold anomaly
!  ARE RH AND T PROBLEMS RELATED?
!  Why is cape low at lower colrh: just a sampling problem?
!  My swa hi is still increasing too early ...
!  Need episodic injection at mid-levels to get rare high cla and ice ...
!
!  Height dep rh_ri_min
!
!  Higher MH Cloud Fractions (esp over lower precip regions)
!  -------------------------
!  - I need much higher MH Cloud Fractions. One way is to increase ice_mh, or
!    make cloud fraction inc faster with ice amount.
!  - look for reduction in outgoing LW
!  - May be coupled: lwa_hi is too high at high rain rates.
!    So need more LW cooling in HI level: lwa_go go down
!    But need more LW heating in MH level: lwa_md go up
!  - But why are diurnal mean lw OK in OCEAN area? Maybe since mainly need to
!      get more cloud over lower precip regions
!
!  Cloud Strategy (May 2015)
!  -------------------------
!  - produce clouds that make the best T profile, or the best compromise between
!    the T profile and realistic ice/water clouds. Believe ceres TOA LW and SW
!    more than heating rates.
!  - painful lesson: model will only work if cloud rad heating produces a 
!    good T profile; so this is the only option. If by getting a good T profile,
!    you violate some aspects of ceres, this is OK. Also, need strong cloud
!    heating in the upper troposphere to get larger LS ascent in convective
!    regions, and reasonable RH.
!
!  Relative Humidity (May 2015)
!  ----------------------------
!  - Rainfall evap does not directly control RH in convective
!    regions. RH is mainly controlled by cloud radiative heating, which increases
!    LS ascent, and gives rise to LS export of mass. Rainfall evap cools, decreases
!    positive GPHT anomalies in the upper troposphere, reduces LS ascent, and may
!    not increase RH regionally, as effectively, as cloud radiative heating.
!  - Cloud heating inhibits convection by reducing cape. But by increasing ascent it
!    will increase rainfall, if convection is sufficiently sensitive to humidity.
!  - Perhaps an important issue is the extent to which rainfall is inhibited by
!    mixing, or by rainfall evaporation. Likely important to get this balance right.
!    Maybe only when rain filtered by mixing is rainfall associated with LS ascent.
!  - However may be true that melting cooling in the lower troposphere always
!    dries the mid-troposphere.
!
!  Cloud Radiative Heating - Ice Heterogeneity
!  ===========================================
!  - It is hard to get SW and LW hi heating correct at the same time.
!  - if I have reasonable hi ice and cloud fraction, my LW tends to be
!    correct but my SW absorption starts increasing too early.
!  - This may be because I am in the linear absorption regime: in this case, 
!    if the actual cloud ice in a grid cell is distributed heterogeneously,
!    so that absorption saturates where cloud ice is big, I will overestimate
!    SW absorption.
!  - Poor mans solution?: move my altitude of peak ice artificially low. In this
!    case, if the LW abs is mostly in the saturated regime, I will have not much
!    effect, but the SW absorption will be moved lower. But it may screw up LW
!    flux from below.
!
!
!  Adding Freezing Heating to Updrafts (April 2015)
!  ================================================
!
!  4 Calls to get_rv:
!  ------------------
!  Call to get_rv from cape: 
!    - No change is needed; rip = 0 always
!  Call to get_rv from up: DID
!    - up called in a virtual and real sense
!    - redefine the entry rtp and kmp to remove the ice part
!    - another call to get_t to get a final "homogeneous" T.
!  Call to get_rv from entrain: XXX
!    - same as from up
!    - modify kp and rtp
!    - second call to get_t
!  Call to get_rv from precipitate: XXX
!
!  Lack of ice in get_rv
!  ---------------------
!  - subroutine get_rv uses get_t,enthalpy,dksdt,rsat, etc.
!  - The definition subroutines enthalpy/hmoist do include ice.
!  - However, the functions calculating derivatives like dksdt/dhmsdt 
!    do not handle ice.
!  - when calling get_rv, have to do so for non-ice part of parcel only,
!    and keep ri fixed at pre-defined value.
!    Modify the input total condensate and enthalpy values (subtract the ice part)
!  - Then have to call another get_t to get final T with ice.
!    This will lead to some small inconsistency between sat rv and T. The rv
!    will be lower than it should be.
!  - get_rv called from subroutines cape, up, entrain, and precipitate
!    enthalpy is kept fixed at the "old" temperature. So to get a slightly modified
!    temperature where all condensate is at the same T, would need a call to get_t.
!  - I don't think get_t has to be changed: condensate amounts there are all 
!    assumed fixed.
!
!
!  - The program calculates new hmnew/rvnew/rlnew/rinew from the horizontal
!    flow variables hmflow/rvflow/rlflow/riflow and the vertical advection flows
!    hmvert/rvvert/rlvert/rvert.
!  - In CAM4, tnew is then calculated based on the assumption that that the water
!    partitioning is fixed (rv/rl/ri) and the known hm (actually km) using the
!    routine t_from_km.
!  - The calculation of the vertical flows hmflow/rvflow/rlflow/riflow is fixed
!  - However the calculation of the horizontal flow variables can be changed
!    depending on the scenario.
!
!  Tuning Issues
!  -------------
!  - so far I had been trying to get the iwp and cla vs rain relationships right.
!  - In doing, this I required much smaller cla for a given iwp than observations.
!  - If I use a correct cla vs iwp relationship, and keep my iwp vs rain fixed, 
!    my cla vs rain becomes way too high, and I get very bad TOA LW.
!  - This procedure also gives me a very bad IWP pdf with a spike at 20-30 g/m2
!    and very few low iwp values.
!  - I think t is impossible to tune iwp and cla vs rain and cla vs iwp all togethor.
!  - Part of the problem seems to come from the fact that at a given rain value, there
!    is a huge variance in iwp. e.g. at rain = 1 mm/day, observed iwp_hi can go from 
!    0 to 1000 g/m2, with a mean of 70 g/m2: very skewed with an extremely long tail. 
!    However the radiative effects scale as ln(iwp) (so maybe I should be tuning to
!    get a median ice value correct instead?). If I tune to get 
!    the mean iwp value right at a given rain rate, I am being strongly influenced
!    by a few very large iwp values, at a given rain rate, and then to get the correct
!    cla, I have to tune my cla versus iwp relationship to get a too small cla. Another
!    problem, is that I have too many weak rain rate values, and so tuning to get a good
!    iwp vs rain rate gives me a too large mean iwp value. This causes problems with
!    lots of subsidence of ice to the lower troposphere where it is a huge source of
!    cloud water, and also gives an extremely unrealistic ice_hi pdf.
!  - I think the cla vs ice relationship should be regarded as physically more fundamental.
!    It likely has less scatter than the iwp vs rain relationship (check?). So constrain
!    the cla vs iwp to be good, and reduce ri at low rain rates accordingly to get a good
!    cla vs rain relationship. This should improve my ice pdf, and overall mean ice
!    value.
!  - one problem with having a lower than observed iwp at low rain rates is that the solar
!    reflectivity might be too low.
!
!
!
!  Inconsistencies (I think related to very high variances in iwp vs rain)
!  -----------------------------------------------------------------------
!  The observed CF, iwp, and rain seem inconsistent 
!  At 1 mm/day: mean observed iwp = 70 g/m2
!  At this iwp, in the CF vs iwp relationship CF = 0.6, but observed CF = 0.3
!    at 1 mm/day.
!  Look at rain bin 10 (rain rate = 1.179 mm/day): 
!  In iwp_hi_RAIN.1x1, mean iwp = 74 g/m2
!  In cla_hi_RAIN.1x1, mean cla = 0.320 cloud fraction
!  cla_hi_ice_hi, iwp = 75 corresponds to line 126 of 200 and cla = 0.74
!
!  rain  cinf
!   1    0.693/4.615 = 0.15
!  10    0.519
!  50    0.85
!
!
!  How to Calculate Layer Fraction. (Maximal overlap and random overlap)
!  --------------------------------------------------------------------
!  - in (3) of Raisanen, MWR, 1998 Note, gives a formula for random
!    overlapped clouds: Ntot = 1 - PROD(1 - Ni)
!  - There is a problem with the "Contiguous" Cloud method for determining
!    when random overlap occurs. At increase vertical resolution, will introduce
!    more "breaks" in a cloud, get more random overlap, and increase total cloud
!    fractions. I think I solve this by imposing a random probability 
!    per pressure interval, so do not need Petri's solution. I also get smoother
!    cloud radiative heating profiles.
!  - TO see how to calculate overlap of several regions, see: ACP 12, 9097, 2012
!    Oreopolous, 3a and 3b.
!  - I think my overlap assumption is similar to "exponential random", maybe look
!    at: Hogan and Illingworth (QJ 2000), Mace and Benson-Troth (2002. They define
!    a "cloud overlap parameter". For vertically continuous clouds, becomes random
!    after 4 km.
!  - Also if you look in Naud: 
!    http://www.arm.gov/publications/proceedings/conf16/extended_abs/naud_c.pdf?id=59
!    C(j+k) = max(Cj,Ck) for maxmim overlap
!    C(j+k) = Cj + Ck - Cj*Ck for random overlap
!           = 1 - 1 + Cj + Ck - Cj*Ck
!           = 1 - (1 - Cj - Ck + Cj*Ck)
!           = 1 - (1 - Cj)*(1 - Ck)
!
!  Accessing new server:
!  ---------------------
!  - folkins-compute1 (from chase)
!  - sftp 192.168.0.238
!  - The optimal compiler flags would be -xHOST -O3
!
!  Converting to new server (Jan 20, 2015)
!  ---------------------------------------
!  - one issue is that I think I made some changes to the h1.f90 files
!    after Balagapol transferred. So I put h1.f90, h2.f90, h3.f90 onto
!    compute1. Not sure if any source code changes
!    were made ...
!
!  Changes to make compile on new computer (Balagapol helped):
!  -----------------------------------------------------------
!  - fix discussed at: https://bb.cgd.ucar.edu/node/1001972
!  - /home/folkins/cesm1_0/models/lnd/clm/src/main/ncdio.F90
!    In the file ncdio.F90 you can fix this problem by replacing the variable nan by the 
!    variable bigint. I did this.
!
!  Possible Papers
!  ===============
!  (1) Congestus formation over Hawaii
!    - What is dynamical source of cold T/colrh anomalies?
!    - Start by looking at radial response of T at HIL to rainfall.
!  (2) Diurnal paper
!    - use acars data to look at diurnal variation of cape/T/RH is SE US
!    - compare with ceres LW and SW diurnal: what is the role of clouds
!      in suppressing the 9 am increase in rainfall? Do to BL resistance,
!      and BL cloud formation, or overall reduced solar at surface, does
!      this also reduce 9 am increase in cape?
!    - look at role of LS vertical motion in modifying RH variation
!    - look at diurnal variation of TS - some NOAAA dataset?
!  (3) Cloud Rad Heating Dependence on Rain
!    - discuss random
!    - propose idea, compare LW and SW with CERES variation with colrh/rain and
!      geographic correlation
!
!  700 hPa (layer 6) Cold Anomaly
!  ------------------------------
!  - the T anomaly is clearly occurring at lower colrh, so does not appear to
!    be directly caused by the convection scheme
!  - way to solve by increasing cloud rad heating at some level? Doubt it
!  - LWP too high?
!  - Comparison with CERES suggests that my LW cooling below 500 hPa is to strong
!    and my SW heating too small. These discrepancies occur for rain < 20 mm/day.
!    Discrepancy could originate from an RH bias or clouds.
!
!  Penetrative Downdrafts and high RH in lower Troposphere
!  -------------------------------------------------------
!  - Model has two kinds of precipitation evaporation: purely local, and purely
!    non-local evaporation that drives air parcels all the way to the surface. Local
!    evaporation is self limiting in that will stop as local RH approaches rh_evap
!    for that rain type. Non-local evaporation can generate very high column RH, as
!    long as the RH of some level is less than the rh_down value. Therefore better to
!    think of non-local evaporation as being driven by the mean rh below the height
!    at which the parcel is entrained from the background atmosphere, rather than the
!    local RH. Likely physically realistic in that downdrafts that reach the surface
!    likely driven by continuous exposure to precipitation throughout downward 
!    trajectory.
!  - downdrafts increase the hm of the lower troposphere and decrease the hm of the BL.
!    (or decrease hm difference between them). I think the only mechanism to do this. 
!    So one way to tell how much you need this mechanism is to see how this difference
!    varies with rain rate, or column RH.
!
!  Strange Things about MJO ....
!  ------------------------------
!  - why do WK plots look the same for 5000 timesteps as 10000?
!    something wrong ??
!  - why do tiny changes in parameters have so much effect?
!  - why does it appear to degrade in NH summer? In NH summer, I mostly have 
!    westward moving waves, Kelvin and MJO both decrease
!
!  Balagapol's changes to script after changes to drives:
!  ------------------------------------------------------
!  - changed location of netcdf:
!    set netcdf = /software    
!  - added "-l netcdff" to a script file. netcdf now divided into c and fortran parts
!    He had some difficulty tracing back to the original scripts here, rather than
!    generated. It was probably: cesm1_0/models/atm/cam/bld/Makefile.in
!  - I had a "word too long error" that was caused by the csh shell shell which 
!    by default used the old bsd style csh. It had a word length limitation
!    of 1024. I have replaced that with tcsh, which doesn't have such limitations.
!
!  Using more than 1 core (oct 2014):
!  ---------------------------------
!  - (I tried this and it didn't work. Program failed conservation and convergence
!     tests. Maybe because my part of the code is not vectorized.)
!  - These are comments from Balagapol.
!  Please see the attached pdf file. That shows the final linking of the binary. 
!  It shows no reference to -openmp, which enables usage of multiple cores.
!  Please see this link - http://www.cesm.ucar.edu/models/cesm1.0/cam/docs/ug5_0/ug.html 
!  and the section -
!  "Configuring CAM for parallel execution" Apparently this model works with either 
!  openmp (for an smp/multi core computer) or using MPI 
!  (for smp/multi core computer or a cluster) 
!  "Note: The use of the -nthreads argument to configure implies building for SMP. 
!  This means that the OpenMP directives will be compiled. Hence, the specification 
!  -nthreads 1 is not the same as building for serial execution which is done via 
!  the -nosmp option and does not require a compiler that supports OpenMP."
!  The above is from the user guide. 
!  This is from /home/folkins/cesm1_0/models/atm/cam/bld/configure :
!  " turn off compilation of OMP directives.  For pure OMP set "-nthreads N -nospmd"" 
!  So on an 8 core computer, it seems like it should be
!  -nthreads 8 -nospmd in configure to enable parallel execution. Thanks. 
!  Another email:
!  Looks like this could be the line to enable openmp in exp.csh  -
!  $cfgdir/configure -test -nlev 26 -hgrid $hgrid -phys cam4 -ocn $ocn -dyn $dyn 
!  -chem none -verbose -nadv 3 -fc "ifort -O3 -inline all" -cc icc || echo 
!   "configure failed" && exit 1
!  In addition to -hgrid, -ocn etc, it could say -nospmd -nthreads 8 according to the 
!  user guide here -
!  http://www.cesm.ucar.edu/models/cesm1.0/cam/docs/ug5_0/ug.html
!
!  strawberry: has 8 cores and 2 threads/core.
!  Am I running on 1 core now?
!
!
!  Considerations When buying a new Computer:
!  ------------------------------------------
!  From Balagapol: 
!  I see that one instance is running (I think he just used top here)
!  16824 folkins   20   0 2211120 1.987g  10528 R 100.0 17.0  44:12.10 cam
!  This is using only one core and ~2 GB of ram. Does the model start using more cores 
!  after a while or is it serial
!  in nature? If it would start using more cores during the course of the run, 
!  then when looking at a new computer in future, you should 
!  consider more cores per cpu and a reasonable clock speed like 2.5 Ghz. If it is 
!  completely serial in nature, the cpu consideration
!  should be less cores per cpu, but a high clock speed  of ~ 3.4 Ghz or so. 
!  The current one runs at 2.27 Ghz discounting turbo frequencies. Thanks. 
!  My comments:
!  - strawberry has 8 cores and 2 threads per core. Multi-threading is I think
!    when you are using more than one core. So since I now run with one core,
!    one job is using 1/16 = 6.25% of the machine. To see this, use the top 
!    command and look at %Cpu(s) under "id" for idle. With cam4 running, I get 
!    93.7% idle, which means 6.3% active.
!  - The reason strawberry slows down with more than 4 jobs is that each cam4
!    simulation takes 2 GB RAM, and it has only 12 GB, so 5 jobs is starting 
!    to push the limits of the memory. Since strawberry has 8 cores, cpu's can
!    easilly handle more jobs.
!  - There is probably something in my changes which prevent cam4 from running
!    with more than 1 core, and causes it to crash when I try (segmentation
!    faults).
! 
! Speed Issues: cores/threads/cpu
! -------------------------------
! Got this comment from web: Looks like you should adjust number of threads in configure file.
!
! If you can run with threading on, then the thing to try is just increasing the
! number of threads until it fails. And check that the answers are identical
! independent of the thread count. That check is important evidence that threading
! is working correctly. Of coarse the other important check is that the wall time
! to complete the run is scaling appropriately for the increased compute resource.
! Your nodes sound plenty big enough to run 8 threads. Sometimes though a per
! thread stack size limit can be exceeded and this can typically be overcome by
! setting an environment variable. Check your compiler documentation for
! information on this.
!
! One other point is that we generally find the best performance at small core
! counts to come from either pure mpi configurations, or from hybrid configurations
! with a small number of threads per task. If your system is dual socket quad core
! chips, then that's 8 cores per node, and with hyperthreading you can assign up to
! 16 processes to the node. I would recommend looking at the performance of pure
! mpi with 16 tasks first, then moving to hybrid configurations with 8 tasks and 2
! threads per task, then try 4 tasks with 4 threads per task. I think it's unlikely
! that using 8 threads per task will perform well on this system.
! 
!
!  Strange 0 value problem in netcdf files (oct 2014): 
!  ---------------------------------------------------
!  - This occurred since was trying to run a multi-core job.
!  - It occurs in both h1 and h3 files (likely h2).
!  Can see using the command: ncdump -v date camrun.cam2.h3.2000-01-01-00000.nc > cc
!  the date variable is assigned zeroes for about 20 values early on.
!  - Occurred only after changes associated with new drives
!  - occurs in a variety of variables: PS_LON_150e_to_170e_LAT_10s_to_10n for
!    example if h1 history files.
!  - occurs for about half a day from 16-50 timesteps.
!
!  The CAPE-Rainfall Mystery
!  -------------------------
!  - cape has decreases at highest rain rates. This presents a problem for
!    the cloud base mass flux. Possible solutions for cape-based closures
!  (i) Can have an extremely small cape_scale, but extremely strong moisture
!      sensitivity, to suppres most of the mass flux except at very high 
!      column RH
!  (ii) Org parameters which directly increase the cloud base mass flux
!       (my amp factor), or indirectly by increasing the width of the km
!       spectrum in the BL.
!
!  module if_conv_tend
!  subroutine get_tend
!  subroutine init_zero
!  subroutine get_conv_fraction
!  subroutine cape
!  subroutine get_start_spectrum (start updraft t/rv spectrum from mass/km spectrum)
!  subroutine up
!  subroutine entrain
!  subroutine updrafts
!  subroutine up_init
!  subroutine get_var (for sigmoidal org variables)
!  subroutine precipitate
!  subroutine detrain (detrains updraft or downdraft parcel into atm)
!  function sigmoidal
!  subroutine get_stuff
!  subroutine bl_mix
!  subroutine get_km_spectrum
!  subroutine get_cape_spectrum
!  subroutine rl_to_rv
!  subroutine ri_remove
!  subroutine rv_to_ansnow_or_ri
!  subroutine rv_to_uprain
!  subroutine define_uprain_start
!  subroutine define_ansnow
!  subroutine evap
!  subroutine calc_vert
!  subroutine conv_tendencies (Calculate tendencies from vert and flow arrays)
!  subroutine calc_eff
!  subroutine calc_org_new (defines precip_org)
!
!  ansnow_conv
!  ansnow_strat
!  ansnow_surf
!  anrain_surf
!
!
!  Faster (faster) Speeding (speeding) up Program
!  -----------------------------------------------
!  - Oct 4, 2014: pretty sure that the change which sped up the program
!    so much was the removal of the -CB option.
!  - Looked at subroutine start in if_conv_solvers.f90. However, this is
!  unlikely to be a problem, since number of iterations never seems to
!  exceed 4.
!  - Look at whether need to call get_rv after precipitate.
!  - monitor the number of calls to get_rv from each type.
!  - maybe can do one level of vectorization ....??
!  - Glen suggested an option called "profile" which will rank lines with most
!    time spent, and an option "inline" for functions, which will put each
!    function in the code.
!  - Sept 8 2014: switched from "CB" to "O3" in exp.csh, and seemed to increase
!    speed by 10% (years/comp day = 0.157 with three other model runs at same
!    time).
!  - see: http://www.personal.psu.edu/jhm/f90/lectures/37.html
!  - ifort: https://computing.llnl.gov/tutorials/linux_clusters/man/ifort.txt
!  - try "-inline all"
!  - could try "-parallel"
!  - Oct 3 2014: I am not sure what caused my ~2.5 times increase in speed
!    about three weeks ago. I removed O3, and seems to have little effect. I am
!    using "-inline all", so that could be it. It could also be due to removal of
!    -CB option, of checking array bounds.
!
!   Full Pressure Levels 
!   --------------------
!   - pressures quite different depending on whether average includes mountains
!    Trop mean      H1 File
!   1. 979.6 hPa    997.43   about 20 hPa layer
!   2. 957.9 hPa    975.31   about 30 hPa layer
!   3. 917.6 hPa    934.19   about 51 hPa layer
!   4. 856.0 hPa    871.37   about 70 hPa layer
!   5. 777.6 hPa    791.48
!   6.              700    
! 
!  Cloud budget:
!  =============
!  Ice: CMEICE = 0
!       DISED = 0 (since I set the ice sedimentation velocity to zero)
!       REPARTICE /= 0
!  liq: CMELIQ = 0
!       REPARTLIQ /= 0 (non-zero between 9 and 10 km)
!       DLSED /= 0
!
!  The repartitioning term has not yet been affected: is calculated in stratiform.F90
!  I could lower the temps of the ice to water repartitioning.
! 
!
!  Factors that could be screwing up MJO
!  -------------------------------------
!  - Think the key is suppressing convection at high surface cape values.
!  - too large amplify could be making large rain events too easilly, so not getting
!    good suppressed phase.
!  - not enough low level cloud rad heating
!  - Randall explanation of MJO (Thayer-Calder and Randall):
!    We can see why the interactions of convection and water vapor are so critical 
!    to the MJO by looking at Q2. Figures 5 and 11 show that the discharge period is a 
!    time of intense precipitation in a nearly saturated column. When convection occurs
!    under these conditions, it is possible that evaporation of precipitation is 
!    inhibited, downdrafts are weaker, and convection is less effective in reducing 
!    instability. With weaker stabilization, con- vection intensifies. The end result 
!    is stronger conden- sation, heavier rain rates, and more intense upper-level latent 
!    heating.
!
!
!   ***************** IMPORTANT PARAMETERS  *******************
!
!      
!  You Need to Get Reasonable Ice Distribution:
!  --------------------------------------------
!   - To reduce solar at ground during day. Otherwise get too much rainfall over
!     land, and a worse diurnal cycle. May also get too persistent rain over Maritime
!     continent, and reduce ability of MJO to propagate.
!   - You could still use constrained radiative heating. However constrained rad
!     heating seems to contribute to the too high lower trop RH problem at the 
!     beginning of rain events (though this may be solved by turning off conv
!     transport); by focusing the heating/upward motion over higher colrh regions
!     may lower the RH of the subtropics too much, and probably undermines the
!     ML divergence. Certainly, the ML divergence is very sensitive to changes
!     in heating near the ML. During the dev of rain events as colrh goes up,
!     rad heating and upward motion at mid-levels will increase. This may undermine
!     the need for dynamically induced upward motion to increase lower trop RH.
!
!
!   Phase Relationship Between LS and Convective T Tendencies 
!   ---------------------------------------------------------
!   - If the LS flow produces the strat T anomaly, you would expect the strat
!     T response to occur exactly half-way between the max strat LS ascent
!     and descent. But the low level cooling slightly lags t = 0. This lag
!     is presumably generated by the convective stratform heating, which
!     would scale closer to the rain peak.
!   - The amplitude of the stratiform waves that produce the LS mid-level
!     divergence dipole should be related to the extent to which
!     convective heating amplifies or damps the T anomalies associated
!     with the waves.
!   - Right now, the Stratiform T response is essentially in phase
!     with the convective stratiform heating. I think I need to get the
!     stratiform conv heating to occur a bit earlier, and more nearly
!     in phase with the LS T tendency.
!   - I think: if the convective strat heating is 1/2 wavelength out of
!     phase with the LS strat heating, they would cancel, and there would
!     not be any LS strat waves. If the conv heating lags by 1/4 wavelength,
!     it would be half in phase, half out of phase, and no amplification 
!     would occur. To get amplification and coupling, would need to have
!     the conv heating displaced a bit earlier than peak rain, to be less
!     than 1/4 wavelngth out of phase. Convective stratiform heating should
!     in fact be shifted earlier because downdraft cooling would be larger
!     prior to t = 0, when the lower trop RH is lower. If congestus heating
!     is displaced after t = 0, due to the increased ML stability anomaly
!     at t = 0, that would be in phase with the LS flow.
!     
!   Recipe: How to Produce Correct Temp Anomaly Patterns
!   ----------------------------------------------------
!   - need a large increase in the fraction of precip as ansnow_conv at high
!     rain rates
!   - At high mid-level RH, evap cooling becomes less efficient. So much of the
!     mid-level cold anomaly probably comes from melting of large amounts of snow,
!     not evap. Melting cooling does not saturate out at high RH, so the relative
!     importance would increase with RH. It also is a pure decrease in the hm of the
!     lower troposphere, so helps avoid the column water catastrophe.
!   - Need some mechanism to weaken evap cooling at high RH to prevent excessive
!     column RH
!   - The mid-level cooling should contribute to a negative GPHT anomaly at 
!     mid-levels and drive preferential inflow near 500 hPa, and then this info
!     should go down (how to force it to go down? Probably need cooling to be 
!     below ML, so LS vertical motion wants to remove the T anomaly. A cold
!     anomaly above the ML will tend to push the converging air up.)
!   - After t = 0, you need low-level divergence and locally high GPHT. This probably
!     explains the increase in surface pressure after t = 0. However, this increase
!     can be too large if the BL cooling is too strong, since BL cooling will lower
!     GPHT in the lower troposphere, and make the formation of positive GPHT 
!     anomalies more difficult.
!
!   Why would anvil precip fraction increase at higher rain rates?
!   --------------------------------------------------------------
!   - may be related to higher buoyancies at higher RH, and larger transport of 
!     condensate to the upper troposphere, or anvil structure in which ice falls
!     outside the cloud, rather than inside.
!   - could be a positive feedback where higher rain rates are associated with more
!     organization, larger updrafts, more shear, and decreased mixing.
!     (is the this the Moncrieff idea?).
!   - maybe plot vertically integrated mass weighted buoyancy about high rain events
!   - plot observed variation in shear about high rain events.
!
!   Constrained Precip Evaporation
!   ------------------------------
!   - are the rh versus colrh profiles "universal"? Presumably would be to the
!     extent that the rain versus colrh profiles are realistic.
!   - constrains about a third of the problem. Still have to produce a correct
!     mass outflow (convective heating and drying) profile given a particular
!     RH profile. Requires that cloud bas mass flux closure be realistic.
!   - It is not clear to what extent the RH versus log precip relationship
!     is due to:
!     (i) rain evaporation 
!     (ii) increased precip efficiency/conv mass flux at higher RH
!     (iii) increased low-level conv/lower trop upward motion.
!     Therefore: "forcing" the relationship via constrained evap may not
!     be entirely justified. However, getting correct ML div anomaly prior
!     to peak rain would reassure that the causality is appropriate.
!
!   Useful Papers:
!   --------------
!   - Sahany, JAS, 2012: looks at entrainment assumptions consistent with column
!     water precip onset.
!
!
!   "Reverse Engineering"
!   ---------------------
!   David Straub: "One approach in modelling is to tune the relevant coefficients 
!     (within reasonable limits) in a way that gives the best possible fit to observed 
!     data. Ultimately, however, this is unsatisfactory since it makes the science 
!     non-deductive."
!   However:
!   (1) Want to create a model that "works" so can be used as a tool to address 
!       other problems, 
!   (2) Parameters derived from physical processes at small scales do not have the
!       same physical meaning when used at larger scales. More realistic to treat as
!       "empirical parameters", rather than parameters that can be attributed to 
!       specific physical processes. Nature is too nonlinear to attempt to attribute
!       processes occurring on the grid scale to processes occurring on the cloud 
!       scale.
!
!  Inward or Outward Spiral
!  ------------------------
!  - Including more detailed (smaller scale) cloud microphysical parameterizations
!    may give rise to more, or less, realistic cloud water/ice distributions, and
!    give rise to cloud radiative feedbacks which increase, or decrease, the
!    realism of cloud radiative heating profiles. May be a tension between
!    smaller scale microphysics (being able to in principle simulate cloud 
!    indirect effects), and simulation of the MJO.
!
!   Downdrafts Kill Convection in Two Ways:
!   ---------------------------------------
!   (1) Penetrative Downdrafts reduce BL moist static energy
!   (2) Cooling below the ML generates a negative GPHT anomaly which generates
!       convergent inflow and subsidence in the lower troposphere, choking off
!       the supply of column water. This mechanism is probably more efficient.
!
!   Three RH/Regimes
!   ----------------
!   - can demarcate three regimes based on primary fate of rlp: 
!      (i) evap regime at low colrh
!      (ii) uprain regime 
!      (iii) ansnow regime (high buoyancy regime, so more condensate goes above ML?)
!   - There is a critical lower trop RH at which updrafts entrain less and produce 
!     rainfall more efficiently.
!   - There is also a critical lower trop RH at which downdrafts shut off.
!
!  Is the CERES Cloud Radiative Heating Symmetric About t = 0?
!  -----------------------------------------------------------
!  - the addition of cloud radiative heating may undermine the tilt by supplying
!    upper tropospheric radiative heating prior to t = 0. I maybe should use the
!    higher resolution CERES dataset: isolate grid points at the SPARC stations
!    and get a long time series. The cloud radiative heating asymmetry is likely 
!    small; however the required heating asymmetry is small also.
!
!  Problem: Upper Trop RH too high at low col_rh/Too High Cloud Rad Heating
!  ------------------------------------------------------------------------
!  - likely the self-lofting problem: too much cloud rad heating and upward motion
!    into the TTL.
!  - Probably specifically caused by a SW cloud heating peak near 100
!    hPa (17 km), which by causing lifting and cooling makes it impossible for
!    ice to get removed from the TTL. Maybe add some precip timescale removal
!    to ice or increase fall speed? Want to get ice more localized around deep
!    convection.
!
!  MJO Test runs
!  -------------
!  - turn off mom transport
!
!  More strongly Tilted Convective Heating?
!  ----------------------------------------
!  - In JAS Holloway and Neelin, 2010, they seperate the column water increase into
!    a slow synoptic component and a sharp mesoscale increase near t = 0, perhaps
!    caused be mesoscale downdrafts. From their Fig 7, it does seem that I have no
!    business having downdrafts before the rain peak.
!
!  How are CAM4 WATER/ICE Handled? (June 2013)
!  -------------------------------------------
!  - Water and ice from CAM4 should be handled the same way as 
!    much as possible
!  - Within the convective scheme, I 
!    transport water/ice as a tracer in the background atmosphere.
!    Otherwise they will get lofted higher and higher by CAM4 background flow,
!    since not being forced to subside by the induced convective subsidence,
!    and I will get strange large values in the TTL.
!  - I do not entrain water/ice into updrafts/downdrafts on formation.
!  - When mixing updrafts/downdrafts with the background atmosphere, I also assume
!    the entrained parcel has no condensate
!  - The only way this can be a mistake is if condensate entrainment is 
!    a significant sink of cloud condensate, or modifies cloud updraft
!    properties much. I Doubt it.
!
!  CAM4 - Stratiform microphysics
!  ------------------------------
!  - Feb 2014: there seems to be a source of ice to water conversion in the cam4
!    microphysics around 6-7 km. There are an extremely large number of tendency terms
!    listed in stratiform.F90, so it would be difficult to determine which it is.
!  - The following tendencies are zero: MPDLIQ, MACPDLIQ, MPDICE, MACPDICE
!  - These four appear to be the total cldice and cldwat total tendencies from the
!    Morrison and revised Morrison microphysics. They only appear to be output for
!    cam5 so not a surprise they are zero.
!  - the cam4 microphysics tendencies appear to be: (based on physics package chosen)
!
!  ZMDLF  - Detrained liquid water from ZM convection  (LIKELY ZERO)
!  CME    - Rate of cond-evap within the cloud
!  DQSED  - Water vapor tendency from cloud sedimentation
!  DISED  - Cloud ice tendency from sedimentation
!  DLSED  - Cloud liquid tendency from sedimentation
!  HSED   - Heating from cloud sediment evaporation   (DON'T NEED)
!  CMEICE - Rate of cond-evap of ice within the cloud
!  CMELIQ - Rate of cond-evap of liq within the cloud
!  LIQ2PR - Rate of conversion of liq to precip       (KNOW IS VERY SMALL)
!  ICE2PR - Rate of conversion of ice to precip       (KNOW IS VERY SMALL)
!  HCME   - Heating from cond-evap within the cloud     (NOT NEEDED)
!  HEVAP  - Heating from evaporation of falling precip  (NOT NEEDED)
!  HFREEZ - Heating rate due to freezing of precip      (NOT NEEDED)
!  HMELT     - Heating from snow melt                   (NOT NEEDED)
!  HREPART   - Heating from cloud ice/liquid repartitioning   (NOT NEEDED)
!  HPROGCLD  - Heating from prognostic clouds                 (NOT NEEDED)
!  REPARTLIQ - Cloud liq tendency from cloud ice/liquid repartitioning  
!  REPARTICE - Cloud ice tendency from cloud ice/liquid repartitioning
!
!  Setting CME = 0 in convective regions
!  -------------------------------------
!  - Setting cme = 0 in cldwat.F90 seems to screw up the relationship between
!    TOA SW and LW with column cldwat and cldice. I am not sure why, but the
!    radiative routine likely needs some cloud properties it doesn't get if set
!    cme equal to zero.
!  - From stratiform.F90,:
!        cmeheat(i,k) = qme(i,k) * ( latvap + latice*fice(i,k) )
!        cmeice (i,k) = qme(i,k) *   fice(i,k)
!        cmeliq (i,k) = qme(i,k) * ( 1._r8 - fice(i,k) )
!  - The most consistent way is to set qme = 0
!    qme: Net stratiform condensation rate
!    qme(:ncol,:pver) = cmeliq(:ncol,:pver) + cmeiout(:ncol,:pver)
!  - This would handle cmeheat, cmeice, and cmeliq is one setting. Hopefully, cmeheat
!    is the real heating which determines the temperature tendency, and is not just 
!    a diagnostic variable ...
!  - For RK, qme is calculated in pcond, which is in cldwat.F90
!  - However, inside pcond, qme is called cme. Set equal to zero there.
!
!   
!    
!
!  ========================================
!
!
!   This module should be as self-contained as possible
!   - any array used here should be defined here
!   - any array dimension of variable used by both this module and the driver
!     should be defined in if_conv_params.f90
!   - should make coupling with cam4 easier
!
!------------------------------------------------------------------------
module if_conv_tend
!------------------------------------------------------------------------
!
!   INCLUDED
!
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!  use statements
!
!  I specifically list some variables here, other is use if_conv statements
!    in individual subroutines. There is some redundancy below, where things
!    listed here, so globally available, are agains listed in subroutines.
!
!  It's better to explicitly list what I am taking from if_conv.
!  Will make coupling easier.
!
!------------------------------------------------------------------------
use if_conv_params,   only: nupdrafts
use if_conv_params,   only: nd,nstep,nsteps,nav
use if_conv_params,   only: numrv,numrl,numri,numdd,numhm
use if_conv_params,   only: vwind,uwind
use if_conv_params,   only: z,p,t,td,hm,km,rv,rl,ri,rt,rs,om,rs_wat
use if_conv_params,   only: av_rlp,av_rip
use if_conv_params,   only: av_mp,av_updz,av_upmp
use if_conv_params,   only: rh,rhod
use if_conv_params,   only: phalf,zh
use if_conv_params,   only: dpsurf
use if_conv_params,   only: bdflow
use if_conv_params,   only: nunstable,maxlev,nlaunch
use if_conv_params,   only: usecam4 ! should have usecam4 = 1 if in cam4
use if_conv_params,   only: utend,vtend,rvtend,rltend,ritend,kmtend
use if_conv_solvers,  only: tkelvin,bad,tkelvin
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none

!------------------------------------------------------------------------
!
!  Cloud Fraction (cloud fraction) 
!
!  - The convective cloud fractions scale with the ZM convective mass flux.
!    So when I got rid of the ZM mass flux, I got rid of an additive source of
!    convective cloud fraction, so my total convective cloud fraction became
!    biased low. The best solution to increasing cloud fraction would be to 
!    introduce an explicit dependence of the convective cloud fractions on 
!    cldice. I can output SH_FRAC and DP_FRAC and confirm they are zero
!
!  - in cam4, convective cloud fraction is tied to convective mass flux,
!    while, "layered cloud" fraction is tied to RH. probably, the motivation
!    behind not directly relating cloud fraction to cldice or cldliq is that
!    the model becomes too unstable.
!
!   parameters used in cloud_fraction.F90
!
!   call cldfrc in stratiform.F90
!   subroutine cldfrc in cloud_fraction.F90
!
!  time mean cldice is supposed (Blossey) to be 0.00002
!
!  My method of getting cloud fractions:
!  - used if ifconvection_activated(i) == 1
!  - implemented in cloud_fraction.F90
!
! CFF: ri_cf is ice cloud fraction
!
!  This is designed to give good agreement with ceres cloud fraction
!     versus ice relationships.
!
!  tt = [log(ri) - log(ri_half)]/log(ri_scale)
!  ri_cf = 1/(1 + exp(-tt))
!
!  ri_half
!  - sigmoidal parameters given below.
!
!  For very small ri, want ri_cf -> 0.
!  So need large negative tt.
!  log(ri) will be large and negative, so log(ri) - log(ri_half) should
!    be large and negative.
!  So need log(ri_scale) small and positive
!
!  ri_scale
!  - To broaden ice distribution: increase ri_scale 
!    Jan 2017: inc to 3.5 and dec ice does not solve hot upper 
!    trop problem in JJA
!  - By spreading out cloud rad heating, larger ri_scale may favour
!    rainfall clustering on larger scales.
!  - Dec 2016: increasing from 3.0 to 3.5 big results in big reduction in
!    high rain rates Asian monsoon, and high rates generally.
!  - DO NOT MAKE ONE! (since can't take log)
!
!------------------------------------------------------------------------
real, parameter :: ri_scale = 3.  ! USED in cloud_fraction.F90
real, parameter :: cloud_fraction_max = 0.95  ! USED in cloud_fraction.F90
real, parameter :: cf_ri_min = 0.01  ! USED in cloud_fraction.F90

!------------------------------------------------------------------------
!
!  ri_half:
!
!  - CONTROLS ICE CLOUD FRACTION
!  - This is temperature dependent (early slide in h3.ps)
!  - increase ri_half: decreases cloud fraction
!  - by decreasing cloud fraction in the TTL can reduce LW heating.
!  - increasing ri_half in levels 11-13: this appears to increase the upward
!    LW flux and heat the TTL, but appears to decrease TOA SW
!
!  ri_min: value of ri_half at coldest T
!
!  In general: 
!  - increasing ri_half decreases CF and cools
!  - increasing ri_min -> cools
!
!  ri_add
!  - in h2.ps, the plot of CF versus layer ice increases quite early, which
!    suggests a smaller value of ri_add. But not sure if this is correct.
!
!  t_ri_half
!
!------------------------------------------------------------------------

real(r8), parameter :: ri_min = 5.0E-06
real(r8), parameter :: ri_add = 3.5E-05
real(r8), parameter :: t_ri_scale = 4.0
real(r8), parameter :: t_ri_half = 205.0

!------------------------------------------------------------------------
!
!  BL CF: Getting CF from CE
!  - this uses a logit function (inverse sigmoidal)
!  - This is tuned to generate the correct monthly mean relationship.
!    However, it appears that due to nonlinearity of this relationship,
!    the relationship that is enforced at 30 min timescales, has to be
!    larger than the monthly mean, to get agreement with the monthly mean.
!  qqq
!
!------------------------------------------------------------------------
real, parameter :: ce_half = 20.
real, parameter :: ce_scale = 5.
real, parameter :: cff_min = 0.65
!real, parameter :: cff_max = -0.54
real, parameter :: cff_max = -0.65

real, parameter :: land_fraction_max = 0.1

!------------------------------------------------------------------------
!
!  BL LWP: Getting LWP from CF
!
!------------------------------------------------------------------------
real, parameter :: cf1 = 0.79
real, parameter :: cf0 = -0.1
real, parameter :: lwp_amp = 40.

!------------------------------------------------------------------------
!
!  Liquid Cloud Fraction
!
!  cf_max_rl
!  - MAX amplitude of CFF (must be less than 1)
!
!  rl_half 
!  - the rl scale over which cf shows largest increase versus rl
!  - reduce to increase liquid CF 
!  - increase rl_half: decreases cloud fraction
!
!  Roughly:
!   CERES LWP   CERES CF  rl (100 hPa thick = 1000 kg/m2)
!    10 g/m2     0.05       0.01 g/kg
!    100 g/m2    0.80       0.10 g/kg
!
!   0.5*log(1 + rl/r0)
!   If r0 = 0.02, 0.01 -> log(1+0.05) = 0.05
!   If r0 = 0.02, 0.10 -> 0.5*log(1+5) = 0.9
!
!   cf_rl = rl_base*log(1.0 + (cldliq(i,k)/rl_scale))
!
!------------------------------------------------------------------------
real, parameter :: rl_base = 0.70
real, parameter :: rl_scale = 5.0E-05

!------------------------------------------------------------------------
!
!  conv_energy_danger  
!  - At very high cape, best thing to do is accelerated transition to
!    deep convection. Otherwise, BL scheme, and shallow convection
!    can just heat up column too much.
!  - cape of 2000 J/kg for 200 hPa thick = 2000*2000 = 4000000
!  - this big probably unrealistic, and adjust precip_org
!
!------------------------------------------------------------------------
!real, parameter :: conv_energy_danger = 4000.*1000.
real, parameter :: conv_energy_danger = 8000.*1000.
real, parameter :: t_danger = 330.

!------------------------------------------------------------------------
!
!  Radiative (radiative) effective liquid diameter.
!  - calculated in pkg_cldoptics.F90
!    reltab(ncol, t, landfrac, landm, icefrac, rel, snowh)
!  
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  Radiative (radiative) effective ice diameter.
!  
!  For discussion of where the variability in rel/rei comes from, see
!  "Cloud Optical Properties" in the manual, p.123.
!  "Over the ocean, the cloud drop effective radius for liquid water clouds, 
!   rel, is specified to be 14μm." Complex T dependent formula over land.
!   (Don't think this is correct)
!  "An ice particle effective radius, rei , is also diagnosed by CAM 4.0. 
!   Following Kristj ́ansson and Kristiansen [2000], the effective radius for 
!   ice clouds is now a function only of temperature, as shown in Figure 4.2."
!  From the top of page 125, it is pretty clear that rei is used to determine
!   the optical properties of ice. For water clouds, it seems that an effective
!   radius is defined that is used to find the 4 optical properties. For ice,
!   there is no mention of a distribution of radii.
!  The program that calculates the T dependent rei is: "subroutine reitab"
!    in pkg_cldoptics.F90
!  The model also seems to use an effective cloud top radius and number.
!   These are output as ACTREL and ACTREI, and in the model are referred to as
!      effc and effi.  (stratiform.F90)
!  Problem: AREI and AREL are output as the rei and rel values. However, this is
!    only for the MG scheme. rei is really what the radiation scheme uses. Where
!   is this value set? Maybe since rei just determine by formula, don't bother
!   outputting. rei is pretty clearly set in call param_cldoptics_calc
!  Maybe EFFICE and EFFLIQ are better to archive?
!  EFFICE   ', 'Micron  ', pver, 'A', 'Prognostic ice effective radius
!  EFFLIQ   ', 'Micron  ', pver, 'A', 'Prognostic droplet effective radius
!   Are these just microphysical? or actually what the radiation scheme uses.
!   If just microphysical, likely have no relevance to me.
!  radsw.F90: this is the actual location where the optical properties of
!    ice are calculated from rei and rel for the RK scheme (I am using).
!  Finally: archive REI and REL directly. These are the radiatively relevant
!    quantities.
!
!  June 2015:
!  - I fix rei in subroutine reitab in pkg_cldoptics.F90 (search IAFF)
!  - up till June 26 fixed at 28 (microns)
!  - rei_reduce = 0.6: does increase TOA SW albedo, but also increases
!    Upper tropospheric heating.
!
!  Sept 2015:
!  - rei_reduce = 0.80 doesn't seem to have much effect
!
!  Garrett Param:
!  - Strongly suppresses deep convective outflow relative to default. (TRUE)
!
!  April 2016:
!  - I should probably stick with Garrett. It doesn't solve the lack of SW
!    reflectivity but definitely improves
!  April 13, 2016: 
!  Garrett TOA Tropical OCEAN LW bias = -11 W/m2
!  Garrett TOA Tropical OCEAN SW bias = -2 W/m2
!  Sum = -13
!  CAM4 TOA Tropical OCEAN LW bias = -4 W/m2
!  CAM4 TOA Tropical OCEAN SW bias = > -15 W/m2 (not on plot)
!  Sum > -19
!  I think with this comparison, I have to keep to Garrett. The problem with
!    the default cam4 is that the SW bias is so huge.
!
!  REI
!------------------------------------------------------------------------
!integer, parameter :: use_IAF_rei = 0  ! Default
integer, parameter :: use_IAF_rei = 2  ! Garrett param

!------------------------------------------------------------------------
!
!  Parameters for Garrett rei scheme:
!
!  re(i,k) = 5.*exp((t(i,k)+75.-273.15)/39.)
!  re(i,k) = rei_min*exp((t(i,k)+t_add.-273.15)/t_bot)
!
!  Jan 2016: reducing rei_min from 5 to 4 has very little effect on TOA SW
!  t_add: changes the effective temperature. Reducing t_add should reduce
!         the effective temperature, and reduce rei for at T.
!         Jan 2016: decreasing t_add from 75 to 60 increases TOA SW by about
!         2 W/m2, but also decreases TOA LW.
!
!------------------------------------------------------------------------
real, parameter :: rei_min = 5.0  ! DEF
real, parameter :: t_add = 75.  ! DEF
real, parameter :: t_bot = 39.0  ! DEF

!------------------------------------------------------------------------
!
!  RRR
!  Change water radius microns
!
!  subroutine reltab in pkg_cldoptics.F90
!  Sep 2015: reduced to 6 from 8 and didn't see much effect
!  Aug 2016: turned off setting in pkg_cldoptics.F90 which put rel equal to
!    14 microns everywhere (wanted larger albedo over land to reduce rain).
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  ASYMMETRY PARAMETER
!
!  Ebert and Curry 1992 
!  asymmetry parameter: gi = ei + fi*rei
!  forward scattering parameter: fi = gi*gi
!  radsw.F90
!
!  The aerosol asymmetry parameter (g) is defined as the cosine-weighted
!  average of the phase function, where the phase function is the 
!  probability of radiation being scattered in a given direction. 
!  Values of can range from -1 for 180-degrees backwards scattering
!  to +1 for complete forward scattering.
!
!  Atmos. Chem. Phys., 13, 3185-3203, 2013
!  Median asymmetry parameters retrieved by the RSP range from 0.76 to 
!  0.78, and are generally smaller than those currently assumed in most 
!  climate models and satellite retrievals. 
!
! "Variation of ice crystal size, shape, and asymmetry parameter in tops 
!  of tropical deep convective clouds"
!  - Results for strongly convective periods show that, with increasing 
!    cloud top temperature, the distortion parameter generally decreases, 
!    while the asymmetry parameter and effective radius increase. 
!
! The only place gtot appears to be used is in the calculation of gtot
!   in radsw.F90.
!
!
!  June 2015: 
!  - setting the asymmetry parameter gci = 0.78 seems to have no 
!  effect on the albedo. (Not sure why ...)
!  - set in radsw.F90
!
!  GCI
!------------------------------------------------------------------------
integer, parameter :: use_IAF_gci = 0
real, parameter :: gci_IAF = 0.78

!------------------------------------------------------------------------
!
!   INCREASES ICE SINGLE SCATTER ALBEDO (Closer to 0.99999)  OBS?
!
!   - to reduce ice absorption of solar energy
!   - would increase cloud SW reflectivity a bit
!   - Jan 2016: with Garrett param, seems to be neccessary to reduce solar
!     absorption to prevent positive T anomalies in high rain regions (perhaps
!     wci_IAF = 0.5). Jan 2017: yes, esp in JJA SPARC sites.
!
!   set in radsw.F90
!   tmp2i = 1._r8 - cbarii - dbarii*rei(i,k)
!   wcl(i,k) = min(tmp2l,.999999_r8)
!   
!    tmp2i = wci_IAF*tmp2i + (1.-wci_IAF)*.999999_r8
!    wci(i,k) = min(tmp2i,.999999_r8)
!
!  - By decreasing solar absorption, would reduce upper trop temps
!  - But does not seem to decrease upper trop temps that much. 
!    Mainly decreases TTL temps, though this can likely be compensated
!    by increasing TTL ice.
!  - Increasing ri_half seems a better option to decreasing upper trop
!    temperatures.
!
!   - wci_IAF = 1: default
!   - wci_IAF = 0.5: means half way between default and 1
!   - wci_IAF = 0.1 increases SS albedo by 90%; very little SW absorption by
!     ice. 
!   - LATER ON: figured out that can get agreement with CERES HI SW by
!     moving det ice lower
!   - April 2016: using this with wci_IAF has a very modest effect on TOA SW
!     but likely an improvement, so keep.
!
!  WCI
!------------------------------------------------------------------------
integer, parameter :: use_IAF_wci = 0
real, parameter :: wci_IAF = 1.00  ! NO CHANGE DEF
!real, parameter :: wci_IAF = 0.2

!------------------------------------------------------------------------
!
!  Modify cloud water droplet single scatter albedo  OBS
!
!  - is a way to reduce cold T bias in the lower troposphere
!  - droplet single scatter albedo det by tmp2l in radsw.F90
!  - This parameter can play a significant role in reducing TOA SW.
!    I think it does this mainly by reducing the RH of the BL and
!    reducing BL RH clouds. However, it may also do this by direct
!    absorption of SW energy. Anyway, if you use even a value
!    such as tmp2l_IAF = 0.99, it is likely you will have to increase
!    ce_land to 900 or so to get large enough cloud fractions near the 
!    surface.
!  April 2016: it is likely better to not modify this, even though it 
!    does add the heat to where it is needed.
!  - don't know what should be realistic values for CD albedo
!  - reduces BL RH to below realistic values.
!  - with tmp2l_IAF = 0.99 I had to use rh_rl_land = 0.70 to get enough
!    cldliq over land to get a low enough albedo (and ce_land = 900).
!    This value of 0.99 for CD albedo is clearly too high. But maybe 
!    0.995 might be a good choice.
!  
!
!  tmp2l_IAF = 0.99: significant effect in reducing 700 hPa cold anomaly.
!         but probably too low for above reasons.
!
!------------------------------------------------------------------------
integer, parameter :: use_IAF_wcl = 0
real, parameter :: tmp2l_IAF = 0.99

!------------------------------------------------------------------------
!
!  ************  Weakening Boundary Layer Mixing  *****************
!
!  BLL
!
!  Sept 2016: reducing mixing coefficients again. 
!    The too strong upward mixing is almost
!    certainly the source of the excessive RH/cape of layers 2 amd 3.
!    It also tends to make a lapse rate that is too unstable, and a cold
!    anomaly near 800-900 hPa at low colrh. I tried turning it off entirely
!    in convective regions, but it becomes a problem what to replace it
!    with. The HB scheme should almost certainly be retained over land.
!    Over the ocean, diffusional mixing of adjacent layers based on lapse 
!    rate fails (warms/dries the surface too much), I think because mixing in 
!    the convective boundary layer is via plumes, not diffusion. One issue is
!    whether the difficulty in convection going above the BL is due to 
!    lack of instability, or RH resistance. 
!    - Another reason for modifying the
!    HB scheme is that, otherwise, LHFLX seems to big (since surface RH too low),
!    and therefore, rain is likely too high (higher than GPCP over the ocean)
!    the atmosphere gets too hot. The only other way to prevent surface RH
!    from going to low is to avoid convective removal unless RH goes above
!    0.80 or so, but this seems artificial. 
!    - The above comments were written before I knew the result of turning
!    off mixing but are mostly confirmed. Reducing mixing definitely reduces
!    rainfall over land (as does lowering amp_add). 
!    Now it is too low. However, this is unintentional.
!    It seems that my mixing reduction still occurs over land even though
!    trying to use landfrac to prevent.
!    - It would be desireable for the convective scheme to have its own BL,
!      scheme, since then the two parameterizations wouldn't be fighting each
!      other.
!
!   An issue is what to reduce. It seems kvm, kvh, kvq should be reduced. But
!     what about cgh and cgs? The countergradient terms are responsible for the
!     non-local transport (i.e. from plumes whose size is comparable with the
!     BL). From 1993 HB paper Eq. 3.8, might be better to leave cgs terms 
!     unchanged. A uniform reduction in transport is achieved by simply changing
!     the K diffusion constants. 
!
!  Aug 2015: decided not to do this. I am playing with things I don't
!     understand. It is true that the scheme mixes moisture off the
!     surface too fast and reduces surface RH too much. But just use 
!     qdiff_factor to compensate.
!     Too complex and too many parameters. I just end up with mush.
!
!  - three kv in hb_diff.F90:
!    kvm(pcols,pverp)         ! eddy diffusivity for momentum [m2/s]
!    kvh(pcols,pverp)         ! eddy diffusivity for heat [m2/s]
!    kvq(pcols,pverp)         ! eddy diffusivity for constituents [m2/s]
!  - The default BL mixing scheme is eddy_scheme= HB for Holtslag lag and 
!    Boville
!  - hb_diff.F90 doesn't calculate the effect on the tracers, just
!    the diffusivities.
!  - kvh seems to be calculated in call austausch_pbl, and then kvq is 
!    set equal to kvh.
!  - the effect of this BL mixing is tested in the 1993 paper, so the problem
!    is probably not the HB scheme itself, but that combined with my mixing, the
!    mixing is too much.
!
!  Dec 2014: 
!    tphysbc.F90: 
!      - calls convection, stratiform, and radiation
!      - all calls done for one chunk
!    tphysac.F90: 
!      - calls vertical_diffusion_tend
!    In physpkg.F90 all chunks done for tphysbc.F90, then all chunks done
!      for tphysbc.F90, so harder to exchange chunk arrays between if_conv.F90
!      and the diffusion scheme.
!
!  How is DTV calculated?
!  ----------------------
!  - identified with ptend&s in vertical_diffusion.F90
!  - It appears that dse is updated in a "call compute_vdiff" 
!    (in vertical_diffusion.F90), which goes to "subroutine compute_vdiff" 
!    (in diffusion_solver.F90). The dissipation (dtk) is added; the counter gradient
!    fluxes, and the explicit surface fluxes.
!  - Why does this DTV warming increase at high rain rates?
!    
!   tphysac.F90
!   - call vertical_diffusion_tend
!   vertical_diffusion.F90
!   - subroutine vertical_diffusion_tend
!     call compute_hb_diff
!     call compute_vdiff
!   diffusion_solver.F90
!   - subroutine compute_vdiff
!
!  May 2015: 
!  - without this, the model tends to produce excessively steep lapse rates
!    near the surface at low rain rates. 
!    Turning this on also seems to reduce the upward flux of moisture
!    and LHFLX (by keeping rh of surface air higher), and 
!    reduce RH/cape of levels 3 and 4. It also seems to lead
!    to situations where can have very large temperature increments. Therefore
!    I reduced these lapse rate values a bit to get a bit more mixing, even if
!    tends to make the BL a bit more unstable than it should be.
!  - in the past these were rain dependent.
!
!  July 2015:
!  - this has been turned off in land regions to avoid generating very high
!    surface temperatures (deserts?). Implemented in if_conv.F90
!  - tends to introduce noise.
!
!  August 2015:
!  - turned off since reduced BL mixing in NH winter oceans and lead to
!    v high surface RH and very unstable LR in bottom layer.
!  - However, is something good about this: seems to improve rainfall
!    correlations in the tropics.
!    
!  kv_hb_mix: lowers the diffusion constants. lowering will
!    increase surface RH.    
!    
!  cg_hb_mix: lowering to 0.5 seems to have no effect.    
!
!------------------------------------------------------------------------
integer, parameter :: reduce_hb_mix = 1
!real, parameter :: kv_hb_mix = 0.20
real, parameter :: kv_hb_mix = 0.30
real, parameter :: cg_hb_mix = 1.0

!------------------------------------------------------------------------
!
!  Reducing Latent Heat Flux over the Ocean:  (OBS?)  BLL
!  - qdiff factor
!
!  Look in: 
!  /home/folkins/cesm1_0/models/csm_share/shr/shr_flux_mod.F90
!  /home/folkins/cesm1_0/models/csm_share/shr
!  In SUBROUTINE shr_flux_atmOcn
!  Could modify delq:
!  delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
!  Or modify input qbot.
!  - shr_flux_atmocn is called from:
!  /home/folkins/cesm1_0/models/drv/driver/seq_flux_mct.F90
!   call shr_flux_atmocn (nloc , zbot , ubot, vbot, thbot, &
!                         shum , dens , tbot, uocn, vocn , &
!                         tocn , mask , sen , lat , lwup , &
!                         evap , taux , tauy, tref, qref , &
!                         duu10n,ustar, re  , ssq )
!   shum is the input to qbot
!   This routine doesn't provide any insight into how shum is calculated,
!     just extracted from rAttr:o
!   shum(n) = a2x%rAttr(index_a2x_Sa_shum,n)
!
!   May reduce trade wind circulation and bring in better agreement with
!     observations.
!
!   Jan 2016:
!   - used qdiff_factor = 1 and got better results.
!   - My atmosphere had been too cold. 
!   - reason for its introduction before may have been because ice was 
!     too high, and so atmosphere warmed up too much. Now using ice around
!     30 g/m2 in tropical mean.
!
!   Jan 2017: 
!   - disables qdiff_factor in shr_flux_mod.F90 since causing compilation
!     problems for zm runs
!
!------------------------------------------------------------------------
real, parameter :: qdiff_factor = 1.00
!real, parameter :: qdiff_factor = 0.90
!real, parameter :: qdiff_factor = 0.85

!------------------------------------------------------------------------
!
!  Parameters in if_conv_solvers.f90
!  - now set in if_conv_solvers.f90
!  - should look at sensitivity
!
!------------------------------------------------------------------------
!integer, parameter :: esat_ice = 1
!integer, parameter :: esat_ice = 0  !  use esat for water all T

!------------------------------------------------------------------------
!
!  ****  Set Default CAM4 cloud prod/loss = 0 in convective regions  ****
!
!  modify_cldat
!
!  cape_mean_activate
!  - Shut off precl production if cape_mean exceeds this value.
!  - The test occurs in if_conv.F90
!  - It is implemented in the routine cldwat.F90, via setting psaut = 0,
!    and in other places.
!  - Check profiles of ICE2PR and LIQ2PR in h1.ps to make sure precl is zero.
!
!  Nov 2014: 
!   - conv.F90, if modify_cldwat = 1 and
!     cape exceeds cape_mean_activate OR
!     latitude less than lat_min,
!   then switch ifconvection_activated(i) = 1.
!   In cldwat.F90, this sets a large number of cloud tendencies to zero.
!
!  Note: if surface T < t_min_cape, cape is not calculated, so cape_mean = 0.
!  So convection and IF clouds not used in this case also.
!
!------------------------------------------------------------------------
real, parameter :: cape_mean_activate = 100.
!real, parameter :: lat_cloud = 30.
real, parameter :: lat_cloud = 38.
integer, parameter :: modify_cldwat = 1
!integer, parameter :: modify_cldwat = 0

!------------------------------------------------------------------------
!
!  Set the ice sediment fall velocity to zero 
!
!  - in pkg_cld_sediment.F90
!  - the model has tendency to over-estimate lwp
!  - part of the problem is the ice fall velocity and subsequent re-partitioning
!    into liquid
!  - The ice fall velocity is set in pkg_cld_sediment.F90
!  - maybe don't do this. Maybe need ice to sink to moisten middle troposphere
!    and push cloud rad heating lower. Address overproduction of liquid
!    from falling ice by narrowing T window at which ice repartitions.
!  - Nov 2014: currently reduce_ice_fallspeed = 0.0. However, there is still
!    a strong downward motion in ice due to subsidence in the convective 
!    parameterization.
!  - Also used to set liquid fall speed equal to zero.
!
!------------------------------------------------------------------------
real, parameter :: reduce_ice_fallspeed = 0.0
!real, parameter :: reduce_ice_fallspeed = 1.  ! for no change

!------------------------------------------------------------------------
!
!  Modify re-partitioning of ice to water (cldwat.F90)
!
!  - IMPORTANT: I HAVE MODIFIED parameters that APPLY GLOBALLY!
!  - This is called the "repartitioning of stratiform condensate" and is given
!    by the output variables REPARTICE and REPARTLIQ.
!  - control parameter is tmin_fice in cldwat.F90
!  - It seems that ice to is too easily repartitioned to water, even
!    after I have set the sedimentation to zero.
!  - the fraction of ice seems to vary linearly from 0 for T > tmax_fice,
!    to 1 for T < tmin_fice.
!  In cldwat.F90
!  The default settings are:
!  tmax_fice = tmelt - 10._r8 (presumably tmelt = 273.15 K?)
!  tmin_fice = tmax_fice - 30._r8 
!  So, you only get pure ice for T < -40 C !!!
!  March 2014: I currently have modified to:
!   tmax_fice = tmelt - 5._r8 
!   tmin_fice = tmax_fice - 5._r8
!  This would give a 5 C transition range, I think of -5 C to -10 C for ice. 
!    This is better. Leads to a thicker
!    ice layer ( down to about 8 km as opposes to 9 km).
!  Mar 2014: changed to:
!   tmax_fice = tmelt + 5.
!   tmin_fice = tmax_fice + 5._r8
!  Want to disconnect ice as a source of ML cloud water - unrealistic. 
!    Better to evap, or otherwise convert to ansnow.
!  Mar 2014: changed to:
!   tmax_fice = tmelt 
!   tmin_fice = tmax_fice - 15._r8
!  Sep 2015: tried -5 to -20: more albedo? LWP to higher altitudes?
!    
!  What these temperature settings do is impose a "desired" condensate fraction 
!  that is ice : fice, and then the REPARTLIQ and REPARTICE tendencies are 
!  calculated to achieve these desired fractions. The best thing to do is
!  to shut these off entirely, at least in convective regions. However,
!  tried and failed: see comments in stratiform.F90. Try to get rid of 
!  ice before it sinks to ML.
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  Unphysical situations
!
!  Whether to set rv = 0 if get rvnew < 0.
!  - likely a good thing to do. However, adds water to a column and
!    violates water conservation.
!
!------------------------------------------------------------------------
integer, parameter :: zero_rv = 1
real, parameter :: km_bad = 450000.
real, parameter :: t_bad = 400.

!------------------------------------------------------------------------
!
!  Parameters for flow/vert arrays
!
!------------------------------------------------------------------------
integer, parameter :: nmom = 2      !  nmom: uflow(nmom,nd) and vflow(nmom,nd)

!------------------------------------------------------------------------
!
!  Parameters for rain arrays
!
!------------------------------------------------------------------------
integer, parameter :: nrain = 26     ! number of levels for various rain arrays

!------------------------------------------------------------------------
!
!  CAPE calculation
!  July 2014: It is probably better to go with reversible CAPE. There is
!    some kind of "bug" in the pseudo cape calculation, in which the values
!    come out to large. It may be related to the way that, at each level, 
!    the rl of the parcel is set to zero before the tdp calculation, which 
!    makes the parcel lighter than it should be. This defect would be related
!    to the discretization. There is no way to remove rlp between levels.
!
!------------------------------------------------------------------------
integer, parameter :: cape_reversible = 0   !  pseudo
!integer, parameter :: cape_reversible = 1   !  retain rlp

!------------------------------------------------------------------------
!
!  **********    Baseline Mass Flux Closure   ***************************
!
!  cape_scale: 
!  - decreasing will increase mass flux of the mode.
!    (Keep same for all modes and change amplitudes)
!  - decreasing tends to warm the atmosphere.
!
!------------------------------------------------------------------------
real(r8) :: cape_scale = 500.

!------------------------------------------------------------------------
!
!  **********    Baseline Mass Flux Closure   ***************************
!
!  tscale_cape: timescale for removal of positive CAPE air from BL in (s)
!
!  Used to get parcel mass:
!  mp(it) = conv_fraction(it,nl)*dmdry
!  conv_fraction(it,nl) = amp*(tstep/tscale_cape)*capp(it,nl)/cape_scale(nl)
!
!------------------------------------------------------------------------
real(r8) :: tscale_cape = 30.0*3600.

!------------------------------------------------------------------------
!
!  t_min_cape:
!  - Only calculate bulk cape and capp of layers if surface T larger than 
!    this value.
!  - Otherwise can get errors.
!
!------------------------------------------------------------------------
real(r8), parameter :: t_min_cape = 273.15

!------------------------------------------------------------------------
!
!  ********** BLL RH MIXING  (DO) **********************
!
!  ONLY DONE IF BL HAS CAPE
!  BL MIXING IS NOT DIFFUSIVE: must involve non-local (convective) hm
!  transport.
!
!  Aug 2016: actually not a good mechanism to decrease RH of BL.
!
!  - Boundary  boundary
!  - exchanges air between neighboring levels if rh(i) > rh_max_bl
!  - Grid cell RH > 0.90 is likely unrealistic. This "excess" water vapor
!    in this case is likely to be used to mix air upward until stopped
!    by an inversion.
!  - In principle, this upward transport could be accomplished by the
!    convective scheme, but at low colrh, it seems like convection has
!    a very hard time getting above the BL. This should help convection
!    get through the BL, by generating a more realistic RH profile in the 
!    BL.
!
!  July 2015:
! - significantly improves the BL lapse rate RH profiles in convective regions.
!   Without this extra mixing, temperatures near the surface are too cold.
!   HOWEVER, by mixing down drier air, it increases the LHFLX, increases the
!   rain rate, increases the temperature of the lowest layer of the model, makes
!   the SHFLX small or negative in high SST regions, and makes atmospheric
!   temperatures too warm overall. The problem here is that you have lots
!   of evaporation over the ocean, but the atmosphere is too warm to allow
!   much convection, so end up with a giant sea breeze circulation, with too
!   high rain rates over land. 
! - In trying to avoid the LHFLX catastrophe, two things that don't work are:
!   (i) increasing tscale_cape
!   (ii) restricting this mixing to levels 3 and 4.
! - However, since this is such a good thing to do, best to find other ways
!   to reduce LHFLX.
!
!  cape_rh_mix
!  - minimum cape to allow this mechanism to occur
!  - Surface RH too low if don't have this.
!
!  rh_max_bl
!  - it is probably unwise to set lower than 0.88. During high rain events
!    rh = 0.88 occurs in the BL.
!
!
!               do_rh_mix=1    do_rh_mix=0
! Jan Land rain     6.0            5.1
! Jan Ocean rain    4.5            4.0
!
!   dmd = amp_rh_factor*(rh(i) - rh_max_bl)*dmdry
!
!  NOW DOES ALL 4 LOWEST LEVELS
!
!  fraction_max:
!  - should probably be less than 0.1 to be realistic. Large values
!    lead to very strong warming and drying of the surface. LR near
!    surface likely unstable since after shflux and dynamics, and
!    before stabilizing from convection.
!
!------------------------------------------------------------------------
!integer, parameter :: do_rh_mix = 0
integer, parameter :: do_rh_mix = 1
real(r8), parameter :: rh_max_bl = 0.88
real(r8), parameter :: amp_rh_factor = 10.0
real(r8), parameter :: fraction_max = 0.10

real(r8), parameter :: cape_rh_scale = 2000.

!------------------------------------------------------------------------
!
!  *********************  Mass Flux Closure  ****************************
!
!  mass_start_max: 
!
!  - safety device to avoid very high rain rates
!  - is TOTAL initial convective mass removed in kg/m2 from a BL layer
!    in a timestep (30 minute)
!  - typically seem to get problems with rain rates > 200 mm/day
!  - In a 30 minute timestep, this is 4 kg/m2 of rain
!  - Suppose rv = 0.025, then this would require 
!    mass_start_max = 4/0.025 = 160 kg/m2
!  - Both of these constraints are currently implemented
!  - Feb 2016: extremely high rain rates are sometimes generated by
!    very high cape situations over land, where the mass flux is 
!    used to generate strongly entraining plumes. 
!
!  mass_inc_max:
!  - going to 4 results in no obvious improvement
!
!------------------------------------------------------------------------
real(r8), parameter :: conv_fraction_max = 0.6
real(r8), parameter :: mass_start_max = 200.
real(r8), parameter :: mass_start_thresh = 1.0
real(r8), parameter :: mass_inc_max = 5.0
real(r8), parameter :: mass_rat_max = 0.5

!------------------------------------------------------------------------
!
!  ncape: nstep at which write out a mixed parcel profile
!
!------------------------------------------------------------------------
integer, parameter :: ncape = 1             !   timestep at which write out cape

!------------------------------------------------------------------------
!
!  Parameters for the Mass Flux Closure
!
!  capp_thresh: 
!   - applied to INDIVIDUAL parcel in spectrum 
!
!------------------------------------------------------------------------
real(r8) :: capp_thresh = 50.

!------------------------------------------------------------------------
!
!  A threshold fraction of mass removed of a parcel in the mass spectrum.
!  If less then this, don't bother with convection
!  Intended to save computer time.
!
!------------------------------------------------------------------------
real(r8), parameter :: conv_fraction_min = 0.0001

!------------------------------------------------------------------------
!
!  General sigmoidal discussion:
!
!  Calling the sigmoidal function requires 4 parameters:
!  function sigmoidal(x_in,x_half,x_scale,y_min,y_max)
!  x_scale: normalize input variable to produce x value
!  x_half: value of x where function achieves half its value
!  y_min: value of y at x = 0.
!  y_max: value of y at large value of x
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!   Jan 2016: DO NOT make a function of colrh
!   
!  "Isolated cumulonimbus"
!  - should there be some deep convection at low rain rates/lower colrh?
!  - The difficulty with having none is that the upper and mid troposphere
!    seem to dry out. Transport does not seem to advect moist air fast enough
!    from convective regions. Also Atmospheric temps seem to get too cold, 
!    and cape too high. However, there should clearly be substantial resistance
!    to organized convection, which builds up from below, and which is triggered
!    by temperature fluctuations above the BL. This is also the way to get a
!    good div outflow before high rain events.
!
!  Sept 12, 2014: the use of amplify appears to undermine variability
!    on MJO spatial scales. On the other hand, cape does not increase
!    during high rain events.
!
!  Oct 2014: tried using amplify only over the land but still seemed to
!    undermine MJO.
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  Boundary Layer Inhomegeneity km_width
!
!
!  PPP
!------------------------------------------------------------------------
!integer, parameter :: use_km_width_land = 0 ! Do not use km_width land
integer, parameter :: use_km_width_land = 1 ! use km_width land

real(r8), parameter :: km_add = 2000.0
real(r8), parameter :: km_min = 0.
real(r8), parameter :: km_scale = 15.
real(r8), parameter :: km_half = 50.

!integer, parameter :: nl_depend = 0
integer, parameter :: nl_depend = 1

!------------------------------------------------------------------------
!
!   Organization Timescales 
!   - lag for defining precip_org
!   - has a strong effect on the dynamics. Should likely by made as small 
!     as possible (not too noisy) so cold anomaly under warm as close to
!     observed as possible.
!   - making tscale_org larger (2 hours) was a way of imposing an inertial
!     smoother on the rate of change of org variables (esp amp). But is not
!     good since leads to f_an lagging precip, and cold anomaly not under
!     warm. But reducing tscale_org leads to too fast growth of rain, and 
!     so rain in intermediate range is smaller. Better to explicitly build
!     in some damping in rate of change of mass flux.
!
!------------------------------------------------------------------------
real(r8), parameter :: tscale_org = 1.0
real(r8), parameter :: tstep_cam4 = 0.5  ! 

!------------------------------------------------------------------------
!
!  Switches
!
!  use_new_kmnew: results very similar with either option (july 29, 2015)
!
!------------------------------------------------------------------------
integer, parameter :: nprint_energy = 0
integer, parameter :: use_new_kmnew = 0

!------------------------------------------------------------------------
!
!  ******************   MIXING  ******************************
!
!  Buoyancy target entrainment
!
!
!------------------------------------------------------------------------

real(r8) :: b_tar(nlaunch) = &
!           (/ 0.008,  &  ! 1.
!              0.008/)    ! 2.
            (/ 0.000,  &  ! 1.
               0.000/)    ! 2.

real(r8) :: b_det(nlaunch) = &
           (/ 0.010,  &  ! 1.  DEF
!          (/ 0.000,  &  ! 1.
              0.000/)    ! 2.

real(r8) :: f_ent_max(nlaunch) = &
            (/  0.20,  &  ! 1.
                0.08/)    ! 2. To reduce ML outflow at high colrh

!real(r8) :: f_ent_max(nlaunch) = &
!            (/  0.24,  &  ! 1.
!                0.18/)    ! 2. To reduce ML outflow at high colrh

real(r8), parameter :: ent_inc_km = 0.01

integer, parameter :: iter_max = 40

!------------------------------------------------------------------------
!
!  ******************  UPDRAFT PARAMETERS  ******************************
!
!  RH    (1-RH)  old    new
!  0.90   0.1    0.00
!  0.85   0.15   0.00
!  0.82   0.18   0.00
!  0.80   0.20   0.04
!  0.70   0.30   0.24
!  0.60   0.40   0.44
!  0.50   0.50   0.64
!
!  Have a sigmoidal function of (1 - RH)
!  Have different amplitudes for different modes
!
!------------------------------------------------------------------------
integer, parameter :: use_drh_low = 0
real(r8), parameter :: rh_f_det = 0.82

! DEF
real(r8) :: f_det_rh(nlaunch) = &
            (/  2.0,  &  ! 1.
                0.00/)   ! 2. 

!real(r8) :: f_det_rh(nlaunch) = &
!            (/  1.0,  &  ! 1.
!                0.00/)   ! 2. 

real(r8), parameter :: drh_low = 0.0
real(r8), parameter :: drh_add = 0.5
real(r8), parameter :: drh_half = 0.40
real(r8), parameter :: drh_scale = 0.08

!------------------------------------------------------------------------
!
!  ******************  UPDRAFT PARAMETERS  ******************************
!
!  Updraft Parameters for breaking CIN
!
!  tke_start: Turbulent Kinetic Energy Scale in the BL.
!    - J/kg give to updraft parcel to break CIN
!    - should be essentially the same as setting a max CIN allowed for
!      convection to occur.
!    - if parcel does not achieve positive buoyancy before this used up, detrains
!    - every parcel is automatically lifted up at least one level. tke_start
!      begins to be used up after that (could change this; makes mass flux
!      insensitive to stability of the bottom level).
!    - if set tke_start = 0, may run into problems with convective shutdown
!    - Oct 2010: have impression that if make tke_start too big then get BL
!      detrainment peak too high, above 2 km, and not enough separation from
!      congestus mode detrainment. Lose trimodal detrainment character.
!      Also, lose clear decrease in BL RH and LR peak near 2km,
!      with too large tke_start (had 50 J/kg before).
!    - With rv conversion to uprain in the BL, it is possible that could build
!      up large heating at the top of the BL, and large convective inhibition,
!      if the heat is not mixed. So do not set too low.
!    - October 2012: Changing between tke_max = 0 (constant tke) and tke_max = 60
!      has a very minor effect on the simulations.
!    - I don't think you can set tke_start = 0. This is because the second
!      detrainment test will then almost always fail.
!
!  Feb 2016:
!  - these have a significant impact on the simulations: Mean BL RH, rainfall
!    variability, event divergence patterns, even mean T of the upper troposphere
!  - by lowering tke_start of a mode, you force the B profile of the mode to be
!    more positive, and lower its det profile in the BL. It might be important
!    that the tke be different between modes, to get more variability in BL
!    detrainment.
!
!  Sept 2016: Need to retain tke to prevent lots of BL detrainment at high colrh.
!
!------------------------------------------------------------------------
real(r8) :: tke_start(nlaunch) = &
          (/  5.0,  &  ! 1.
              10.0/)    ! 2.

real(r8) :: dp_tke_start = 1000.*100.

real(r8), parameter :: tke_low = 1.0
real(r8), parameter :: tke_add = 6.0
real(r8), parameter :: tke_half = 30.0
real(r8), parameter :: tke_scale = 15.0

!------------------------------------------------------------------------
!
!  Determine the relative strengths of the modes
!
!  conv_fraction(it,nl) = amp*(tstep/tscale_cape)*capp(it,nl)/cape_scale
!
!  amp_low:
!   - inc amp_low from 0.0 to 0.1 does not cool JJA upper trop
!
!  amp_add:
!   Dec 2016: dec amp_add from 2.5 to 2.2 no obvious improvement
!
!------------------------------------------------------------------------
integer, parameter :: use_om5 = 0   ! To

! This is what is done for use_om5 = 0   DEF  amp det from rain only

real(r8), parameter :: amp_low = 0.0
real(r8), parameter :: amp_add = 2.5
!real(r8), parameter :: amp_low = 0.4
!real(r8), parameter :: amp_add = 1.0
real(r8), parameter :: amp_half = 45.0
real(r8), parameter :: amp_scale = 25.0

!------------------------------------------------------------------------
!
!  OM5 initiation
!
!  Convective Initiation by Vertical Motion
!  ----------------------------------------
!  - have a switch for mass flux option
!  - check that both positive
!  - amp = amp_om5 + amp_rain
!
!  om5 likely pa/s
!------------------------------------------------------------------------

! This is what is done for use_om5 = 1    NOT USED

!
!  Here using both rain and om5 to determine convective mass flux
!  amp = amp_om5 + amp_rain
!

real(r8), parameter :: amp_om5_low = 0.0
real(r8), parameter :: amp_om5_add = 0.5
real(r8), parameter :: amp_om5_half = 100.
real(r8), parameter :: amp_om5_scale = 30.0

!  Used to restrict amp_om5 to low rain only
!
real(r8), parameter :: fff_low = 1.0
real(r8), parameter :: fff_add = -1.0
real(r8), parameter :: fff_half = 5.0
real(r8), parameter :: fff_scale = 2.0

! determines amp_rain
!
real(r8), parameter :: amp_rain_low = 0.0
real(r8), parameter :: amp_rain_add = 2.5
real(r8), parameter :: amp_rain_half = 45.0
real(r8), parameter :: amp_rain_scale = 25.0

!------------------------------------------------------------------------
!
!  Mixing/Deep mode partitioning
!
!  a1: strengthening may weaken 2 km stability max, and not warm BL.
!
!------------------------------------------------------------------------
integer, parameter :: use_om5_a = 1

!  use rain only; use_om5_a = 0
!
real(r8), parameter :: a1_low = 1.0
real(r8), parameter :: a1_add = -1.0
real(r8), parameter :: a1_half = 25.0
real(r8), parameter :: a1_scale = 10.0  ! 12 doesn't seem good

!  use_om5_a = 1
!
real(r8), parameter :: a1_om5_low = 1.0
real(r8), parameter :: a1_om5_add = -1.0
real(r8), parameter :: a1_om5_scale = 30.0  
real(r8), parameter :: a1_om5_half = 100.0
!real(r8), parameter :: a1_om5_scale = 40.0  
!real(r8), parameter :: a1_om5_half = 130.0

!------------------------------------------------------------------------
!
!  ******************  UPDRAFT PARAMETERS  ******************************
!
!  mix_option = 0: shut off mixing (e.g. to test conservation properties)
!  mix_option = 1: Forward + Backward B gradient mixing
!
!------------------------------------------------------------------------
!integer, parameter ::  mix_option = 0         !   TESTING: no mixing
integer, parameter :: mix_option = 1         !   default 

!------------------------------------------------------------------------
!
!  ******************  UPDRAFT PARAMETERS  ******************************
!
!  Updraft Detrainment
!
!  safety_entrain: do not entrain if entrainment so far exceeds this fraction of
!    layer mass. 
!
!  mass_ratio_max
!    - don't let updraft parcel get larger than this multiple of original mass
!    Bottom levels in cam4 extremely thin, so large increases in the mass
!    should not cause problems. Better to turn off entrainment when entrained
!    mass is larger than some fraction of the mass of the layer ... would have to 
!    be a function of nlaunch, say once (dm > dm_layer*2/nlaunch) to avoid
!    numerical problems
!    - June 2015: lowered this to 3, to try and reduce rain rates over tropical
!      mountain areas with high colrh. Plumes entrain all the way up and get bigger
!      and bigger, and produce too much rain.
!    - Aug 2016: increased from 3 to 10 and made now difference. Better to keep 
!      large. Have high rain rates over mountain areas, and a value of 3 does not
!      solve.
!
!------------------------------------------------------------------------
real(r8), parameter :: safety_entrain = 0.3    
real(r8), parameter :: mass_ratio_max = 10.0    

!------------------------------------------------------------------------
!
!  ******************  UPDRAFT PARAMETERS  ******************************
!
!  Updraft Precipitation Treatment
!
!  nprecip: determines threshold condensate loading for formation of precipitation
!          
!  nprecip = 1: call precip
!  nprecip = 2: do not call precip  (for determining water conservation errors)
!
!------------------------------------------------------------------------
integer, parameter :: nprecip = 1   !  Default 
!integer, parameter :: nprecip = 2  !  turn off all precipitation, for TESTING
                                    !  water conservation 

!------------------------------------------------------------------------
!
! Feb 2013: I tried this and it seems to make very little difference,
!    so probably not use. Tried to take into account surface pressure
!    reduction due to rain reaching surface, but very hard to implement
!    properly. 
!    60 mm/day = 60 kg/m2/day = 600 Pa/day = 6 hPa/day
!    In west pacific, obs ps changes by .4 hPa over 8 hours, so may be 
!    significant, though small compared with the vertical motion.
!
!------------------------------------------------------------------------
integer, parameter :: do_ps_reduction = 0
integer, parameter :: do_ps_reduction_comments = 0

!------------------------------------------------------------------------
!
!  Convective Momentum transport (momentum):
!
!  - there is also a switch l_windt in if_conv_intr.F90 that has a default
!    setting of true.
!  - This inserted here just for my own convenience.
!  - applied in if_conv.F90
!  - Sept 2014: mom transport puts all MJO variance into antisymmetric 
!    spectrum. Does seem to increase westward MRG waves. Doesn't seem like
!    a good idea. May be a problem with implementation.
!
!------------------------------------------------------------------------
integer, parameter :: do_conv_mom_transport = 0
!integer, parameter :: do_conv_mom_transport = 1

!------------------------------------------------------------------------
!
!   *******  Boundary (boundary) Layer BLL Drizzle drizzle **************
!
!  April 2016:
!  - important not to set rh_rv_to_uprain to low since high RH around 
!    cold SST regions is needed to generate BL clouds.
!  
!------------------------------------------------------------------------
!integer, parameter :: switch_rv_uprain = 0    ! Default: allow rv to uprain
integer, parameter :: switch_rv_uprain = 1    ! Default: allow rv to uprain
real(r8), parameter :: press_rv_to_uprain = 750.*100.
real(r8), parameter :: rh_rv_to_uprain = 0.92
real(r8), parameter :: tscale_rv_to_uprain = 2.*3600.

!------------------------------------------------------------------------
!
!  ***********  START CONDENSATE SECTION  *******************************
!  
!    Overview  
!  
!  - ri_dt refers to condensate produced by detrainment 
!  - ri_rh refers to condensate produced in situ
!  
!  Probably need ri_rh to push ri peak (and heating peak) above 12 km  
!  And to generate more widespread ri
!
!  ri profile may not equal ri_tar:
!  - modeled ri profile tends to be shifted below ri_tar, resulting in heating
!    profile that is too low.
!  (1) LS advection pushes up ri profile. ri will tend to exceed ri_tar
!    near the TTL but be below ri_tar below 250 hPa.
!  (2) When convection detrains rl, ri_tar is not yet known. Therefore detrain
!      ri to not allow ri + ri_added to exceed some multiple of initial ri. But
!      due to effect above, this will favour more ri detrainment below 250 hPa.
!      Combat this by forcing pressure dependence in ri_fraction_max.
!  (3) After convection is finished and ri_tar is calculated, the program 
!      affects ri as follows:
!     (i) call ri_remove: excess ri determined from ri(i) + ri_added(i) - ri_tar(i)
!         excess allocated to ansnow or evap depending on local rh
!
!  To increase Ice Height (more ice near 13 - 14 km):
!  (this is the height where small particles likely mean that ice is 
!   suspended longer)
!  p dependence to ri_fraction_max
!  rh_max_ttl : increasing may result in more ice in the TTL
!
!  However: likely dangerous to have ice profile too tightly concentrated
!   at one height: gives rise to constant uplift at that height, and
!   blocks detrainment.
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!  ICE RI_DT DETRAINMENT ICE
!  - parameters used to determine ri_dt
!
!  ri_base
!  - sets the magnitude of detrainment ice 
!  - Aug 2015: if goes higher than 0.05, ocean mean TOA LW
!    tends to be too low, unless prob_random < 0.2
!  - Sep 2015: 
!    To get more rainfall variance: helps to set ri_scale larger
!  - If TOA LW too low may be better to increase ri_half. This way
!    can maybe increase TOA LW but have less impact on TOA SW.
!
!  ri_dt = ri_base*tot_outflow
!  - I had before a log relationship. But pretty sure this would be wrong.
!    It would reduce cloud radiative feedbacks significantly, since more
!    rain would result in only modestly more cldice. And radiative effects
!    tend to already increasingly saturate with cldice anyway.
!
!  tot_outflow_min:
!  - set ri_det = 0 if less than this value. Intended to get
!  rid of thin ice layers that decreased TOA LW. This parameter
!  does increase TOA LW, but also decreases TOA SW (too much).
!
!  AUG 2016: Solution to low albedo problem over oceans is to make lots of
!    ice. Had previously used MLS as target, but this is likely too low. Need
!    thick ice clouds.
!
!  outflow_power:
! - the idea here was to make ice a nonlinear function of outflow, and
!   get curvature in TOA SW vs TOW LW relation.
! - But appears to contribute to a spike in rain rates near 45 mm/day
! - dec 2016: has a strong impact on the diurnal cycle over land: reducing
!   from 1.3 to 1.0 decreases rain during day (esp over Africa).
!
!   ri_dt(i) = gauss*ri_base*(tot_outflow**outflow_power)
!-------------------------------------------------------
real(r8), parameter :: ri_base_o = 5.0*0.001
real(r8), parameter :: ri_base_l = 6.0*0.001
!real(r8), parameter :: ri_base_o = 4.0*0.001
!real(r8), parameter :: ri_base_l = 5.0*0.001

real(r8), parameter :: outflow_power = 1.0
real(r8), parameter :: tot_outflow_min = 0.001

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!  ICE RI_DT DETRAINMENT ICE
!
!  Gaussian Peak  
!  - seems like no other way to keep ri_dt to a high enough altitude
!  - tried ri_dt proportional to outflow at that level but doesn't work.
!
!  Reasons for increasing z_peak with rain (OBS)
!  - suggested by ceres data
!  - with reduced SW heating from SS albedo reduction, need a mechanism
!    to smooth out v. strong LW cooling at cloud top.
!  - hopefully increase albedo at high rain with smaller ice
!  - maybe inc TT temps
!  - July 10 2015: no obvious advantage in doing this.
!
!  t_ice_det = 263.
!  - condensate detrains as liquid at T > t_ice_det
!  - this parameter shouldn't affect ice much since the gaussian profile
!    will mainly cut off ice detrainment anyway.
!
!  z_peak_peak_o
!  - Should likely be >= 12 km; otherwise hard to get good TOA LW values
!    is western tropical Pacific.
!  - I think the argument below was a problem when cf_max_ri was 0.9 and
!    the cooling rates at cloud tops could be huge.
!  - should likely be 12 or below. If put more detrained ice higher where
!    r_half is smaller, will get sharp peak in cloud fractions near 14-15
!    km and very strong peak in LW heating there. Also very high and 
!    presumably unrealistic ri/rs ratios. In the grid mean, likely this
!    ratio should never be above 0.2 or 0.1.
!
!  z_peak_width_above_av = 3.00*1000.
!  - BE CAREFUL: do not increase above 3 km with z_peak = 11.5 km;
!    likely get very low TOA LW and warm anomaliues at high rain;
!    try another mechanism to increase clouds and T in the TTL that
!    is not concentrated around high rain
!
!------------------------------------------------------------------------
real(r8), parameter :: t_ice_det = 265.

!real(r8), parameter :: z_peak_add = -1.50*1000.
!real(r8), parameter :: z_peak_add = -2.00*1000.
real(r8), parameter :: z_peak_add = 0.0*1000.
real(r8), parameter :: z_peak_half = 35.
real(r8), parameter :: z_peak_scale = 15.

real(r8), parameter :: z_peak_o = 12.0*1000.
real(r8), parameter :: z_peak_l = 12.0*1000.

real(r8), parameter :: z_peak_width_above_av = 2.25*1000.
real(r8), parameter :: z_peak_width_below = 4.0*1000.

!------------------------------------------------------------------------
!
!  Seasonal variation of z_peak
!
!------------------------------------------------------------------------
!real(r8), parameter :: z_peak_width_above_amp = 0.25*1000.
real(r8), parameter :: z_peak_width_above_amp = 0.00*1000.

!real(r8), parameter :: year_frac_add = -0.15
real(r8), parameter :: year_frac_add = -0.20

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!   WATER RL_RT  DETRAINMENT WATER
!
!  June 2016 rl philosophy
!  - After convection, rl_dt and rl_rh are used to define rl_tar
!  - The subroutines rl_remove and rv_to_rl are then used to adjust
!    the tendencies to try and match rl_tar. There is no actual
!    "detrainment" of rl in the model (not needed).
!  - a problem here in trying to get exact agreement between rlnew 
!    and rl_tar is that the subsidence tendencies are only know at the
!    end. But this should only be a problem at high rain rates.
!  - Should adopt same philosophy for ri.
!
!------------------------------------------------------------------------
real(r8), parameter :: rl_o_base = 0.00400*0.001
real(r8), parameter :: rl_l_base = 0.00400*0.001

!----------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION: IN SITU RH ICE  ***********************
!
!  - used to determine ri_rh
!  - ceres iwp and iwp patterns suggest not that strong a sensitivity
!    to RH.
!
!  July 2015:
!  - introduced different land/ocean ri_rh parameters
!  - My TOA LW over ocean was too low, but my TOA LW over land was OK.
!    Seem to need diff land/ocean parameters to get this right (if TOA
!    LW wrong screws up atm T).
!
!  ri_rh = (rh(i) - rh_ri_min)*ri_amp*log(1. + (rs(i)/rs_scale))
!
!  use_rh_max = 1
!  - in calculating in situ ice rh, ice production does not go higher if
!    the local rh exceeds rh_max. Intended to prevent excessive
!    RH/heating/upward motion feedbacks. rh_max should be the normal
!    upper limit, physically, for RH anyway.
!   
!-----------------------------------------------------------------------
real(r8), parameter :: f_ri_max = 0.5
real(r8), parameter :: ri_amp = 2.0E-05
real(r8), parameter :: rs_scale = 2.0E-05
real(r8), parameter :: rh_ri_min = 0.80

!-----------------------------------------------------------------------
!
!  WATER IN SITU RH CLOUD
!
!  July 2015:
!  - In CERES, the downward LW at the surface is larger over tropical
!    land than tropical ocean. This suggests that low-level cloud is
!    more common over tropical land than ocean. So need different 
!    parameterizations. Important to get this right, since if downward
!    loss at the surface of LW is wrong, mean T will be wrong also.
!
!  p_rl_rh_min 
!  - in situ water clouds restricted to below this pressure level.
!
!  Aug 2016: currently rl_rh_min and rl_rh_base used over land only.
!     ocean rl determined from a lwp from ce. I need a fix over land
!     to make these bigger.
!
!-----------------------------------------------------------------------
!real(r8), parameter :: rl_rh_min = 0.80
real(r8), parameter :: rl_rh_min = 0.85

!
! Increasing this parameter delays onset of increase in rain over land
!
!real(r8), parameter :: rl_rh_base = 0.001 
real(r8), parameter :: rl_rh_base = 0.0005

real(r8), parameter :: p_rl_rh_min = 500.*100.  

!------------------------------------------------------------------------
!
!  RL_RH rl_rh
!
!  ce_amp:  
!  - CF and LWP extremely sensitive to changes in ce_amp, since affects
!    both and variables are coupled.
!
!------------------------------------------------------------------------
integer, parameter :: use_ce_land = 0  ! Do not reduce rl over land from ce

real, parameter :: ce_amp = 0.7
real, parameter :: ce_cf = 200.
real, parameter :: ce_lwp = 40.

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!  LIMITS to ICE DETRAINMENT 
!   - problem in keeping ri close to a "target" value is that even a 
!   a tiny fraction of detrained condensate (relative to that going to
!   the part going to precipitation) can generate ri far in excess of the
!   target value. 
!
!  ri_factor
!   - do not allow old ri + ri_added to exceed this fraction of current ri
!
!  ri_fraction
!   - the allowed fraction of ice to detrain
!   - increases rapidly as go up
!
!------------------------------------------------------------------------
real(r8), parameter :: ri_fraction_max = 0.01

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!  CONDENSATE DETRAINMENT 
!
!  These parameters are mainly a way of increasing cloud at LOWER precip
!    and colrh.
!
!  From looking at the ri prod/loss terms versus colrh, the way to
!    increase cldice at lower colrh would be through detrainment, rather
!    than condensation.
!
!  Increasing these parameters does result in more net detrainment
!    of water to the background atmosphere, and therefore does
!    affect the background humidity.
!
!  full_detrain_only_ice 
!  - equal 1: only detrain condensate at full detrainment.
!  - Jan 2015: value of one appears to help prevent excessive LW heating
!    near 5/6 km and reduce warm T bias in the mid-troposphere. 
!    However, it will decrease IWP and increase SW at the suraface, and
!    TS over land regions, and possible produce too much cape and rainfall.
!  - March 2015: set equal to one to try and get more mid-level ice/liq at
!    rain rates; won't make high rain rate problem worse
!
!------------------------------------------------------------------------
integer, parameter :: full_detrain_only_ice = 0

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!  ICE REMOVAL 
!  ri -> rv 
!  ri -> ansnow_ri 
!
!  ri evap or conversion to ansnow_ri
!  - done if ri exceeds ri_tar
!  - ri -> rv at lower rh and ri -> ansnow_ri at higher rh
!  - essentially identical to rl -> rv or uprain_rl
!  - use rh_max here to define upper range of rh
!
!  f_ri_remove_max
!  - max fraction of ri allowed to remove.
!  - Nov 2015: increased to 0.9 to try and remove thin ice layers
!    that I thought were biasing TOA LW low. (Better way is to set
!    cf_min)
!
!  tscale_ri_rem_trop 
!  - This removes "excess" ri, i.e. in excess of ri_tar
!  - Aug 2016: increased to 4 hours in a way to smooth out heating from
!    clouds, and possible give a larger low TOA LW signature in western
!    Pacific. Reducing to 4 hours does not solve 15 km high RH problem.
!
!------------------------------------------------------------------------
!real(r8), parameter :: tscale_ri_rem_ttl = 1.0*3600.
real(r8), parameter :: tscale_ri_rem_ttl = 4.0*3600.
real(r8), parameter :: tscale_ri_rem_trop = 4.0*3600.
!real(r8), parameter :: tscale_ri_rem_trop = 1.0*3600.
real(r8), parameter :: f_ri_remove_max = 0.20

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!   WATER REMOVAL
!   rl -> rv 
!   rl -> uprain 
!
!  rlp_start_entrain:
!     whether to entrain rl at initial entrainment
!
!  f_rl_rem
!  - do not remove more than this fraction of initial value.
!  - To avoid negative rl (issue is vertical advection is done at the
!    end of convection, and can generate negative rl from small rl.)
!
!------------------------------------------------------------------------
real(r8), parameter :: rh_rl_evap_max = 0.95
real(r8), parameter :: rh_rl_evap_min = 0.85
!real(r8), parameter :: f_rl_rem = 0.90
real(r8), parameter :: f_rl_rem = 0.85
integer, parameter :: rlp_start_entrain = 1

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!  CONDENSATE ADVECTION  
!  
!  March 2015:
!  - appears to eliminate negative ri problems.
!  May 2015:
!  - pushing down the ice at high rain rates appears to be required to
!    combat LS ascent, and keep ice peak av_z_ri reasonably close to 10.0 km
!
! turn_off_cond_adv = 0: do not turn off condensate vertical advection
! turn_off_cond_adv = 1: turns off ri adv
! turn_off_cond_adv = 2: turns off both ri and rl vertical advection
!
!------------------------------------------------------------------------
integer, parameter :: turn_off_cond_adv = 0
!integer, parameter :: turn_off_cond_adv = 2

!------------------------------------------------------------------------
!
!  ***********  CONDENSATE SECTION  *********************************
!  
!     TTL Modifications of ri_rh
!  
!  rv_rem_max:
!  - never remove more than this fraction of existing rv
!
!------------------------------------------------------------------------
real(r8), parameter :: rv_rem_max = 0.50

!------------------------------------------------------------------------
!
!   **********  END CONDENSATE SECTION   ********************************
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!   ******  rv -> ansnow OR ri **************************************
!
!  parameters used in subroutine rv_to_ansnow_or_ri 
!
!  rh_max_trop
!  -  threshold RH for upper troposphere rv-to-ansnow conversion
!  -  Applied only at and above i_strat (t_strat)
!  - the amount of rv that goes to ri is determined by the maximum amount
!    that can be allowed up to ri_max
!  - April 2015: NOISE
!    In situations where the upper troposphere has high RH, ansnow_strat
!    becomes very noisy. Probably, an increase in ansnow_strat reduces RH
!    via direct removal of water and subsidence, so ansnow_strat on the 
!    next timestep goes down. But upward motion pushes RH back up for the
!    next timestep. This noise source appears to be strongest in the 300-500
!    hPa interval. One solution would be to make ansnow prod from rv
!    dependent on rv from last two timesteps
!  
!  p_ttl
!  - allow higher TTL RH of rh_max_ttl above this pressure level  
!  - should be interpolated
!  
!  rh_max_ttl 
!
!   These parameters are used to define an rh_max
!
!  REMOVAL of ri to RI or ANSNOW
!
!  lat_rh_wat: latitude at which start using RH with respect to water
!    Oct 16: value of 20 undermined warm anomaly in the upper troposphere
!      Try later. Dec 2016: tried again and was bad
!
!------------------------------------------------------------------------
real(r8) :: rh_min_trop = 0.75
real(r8) :: rh_max_trop = 0.80
!real(r8) :: rh_min_ttl = 0.90
real(r8) :: rh_min_ttl = 0.85
real(r8) :: rh_max_ttl = 1.00

real(r8) :: lat_rh_wat = 90.0
!real(r8) :: lat_rh_wat = 20.0

real(r8) :: p_ttl_top = 90.*100.
real(r8) :: p_ttl_bot = 150.*100.

!------------------------------------------------------------------------
!
!   ******  STRATIFORM PRECIP: rv -> ansnow PARAMETERS ********************
!
!   ansnow_rv
!   This switch also shuts off most of ri/rl production
!   This process is based on an rh_max profile, obtained from 
!     rh_max_trop, etc.   
!
!------------------------------------------------------------------------
!integer, parameter :: switch_rv_ansnow = 0   ! turn off
integer, parameter :: switch_rv_ansnow = 1    ! Default: allow rv to ansnow and ri
real(r8), parameter :: rv_fraction_remove = 0.60 ! safety to avoid neg rvnew 

!------------------------------------------------------------------------
!
!   ******  STRATIFORM PRECIP: rv -> ansnow PARAMETERS ********************
!
!   t_strat
!   - the temperature at which start converting rv -> ansnow
!   - making this a bit colder than 273 K appears to moisten the RH near
!     the ML, and reduce oscillations during high rain events.
!
!------------------------------------------------------------------------
!real(r8), parameter :: t_strat = 260. ! 
!real(r8), parameter :: t_strat = 270. ! 
real(r8), parameter :: t_strat = 267.

!------------------------------------------------------------------------
!
!   ***********  Cloud PARAMETERS ***********************************
!
!  - Reduce to get larger rlp in the upper troposphere and more moistening
!    from rlp -> rv
!  - Making smaller results in more melting cooling below the ML.
!    
!  i_rem_min
!  - do not remove any condensate if i less than or equal to this value    
!    
!  PPP
!------------------------------------------------------------------------
real(r8) :: f_rem(nlaunch) = &
             (/ 0.30,   &  ! 1.
                0.24/)     ! 2.

integer, parameter :: i_rem_min = 1

real(r8), parameter :: rlp_thresh = 0.000
real(r8), parameter :: rip_thresh = 0.000

!real(r8), parameter :: lwc_max = 2.0   !  in g/m3
real(r8), parameter :: lwc_max = 3.0   !  in g/m3
!real(r8), parameter :: t_lwc_min = 250.0   !  T lwc starts to decrease
!real(r8), parameter :: t_lwc_max = 185.0   !  T lwc goes to zero
!------------------------------------------------------------------------
!
!  **********   UPDRAFT CONDENSATE PARTITIONING  **************************
!
! f_ice
!
! Does this heat the upper troposphere? Probably a little.
! - But makes little overall difference. 
!   (At least with t_ice_high = 250; maybe more impact if lowered)
!
!------------------------------------------------------------------------
!real(r8), parameter :: t_ice_low = 200.00
!real(r8), parameter :: t_ice_high = 210.00
real(r8), parameter :: t_ice_low = 170.0
real(r8), parameter :: t_ice_high = 180.00

!------------------------------------------------------------------------
!
!  **********   UPDRAFT CONDENSATE PARTITIONING  **************************
!
!  f_an_add:
!   Note: when rain = f_an_half, f_an will equal f_an_min + 0.5*f_an_scale
!   
!  f_an_half
!  - Should probably be the maximum typical rain rate of rain events, i.e.
!    the rain rate at which want convergent mid-level inflow to start driving
!    down colrh. Appears that if lower a3_add to get more variance, should
!    increase this parameter.
!
!------------------------------------------------------------------------
integer, parameter :: use_f_an_om = 0

real(r8), parameter :: f_an_min = -0.1
real(r8), parameter :: f_an_add = 1.0
real(r8), parameter :: f_an_scale = 15.0
real(r8), parameter :: f_an_half = 60.

real(r8), parameter :: f_an_om_min = -0.1
real(r8), parameter :: f_an_om_add = 1.0
real(r8), parameter :: f_an_om_scale = 30.
real(r8), parameter :: f_an_om_half = 200.

!------------------------------------------------------------------------
!
!  **************   MELTING EVAP PARAMETERS  **************************
!
!  Melting downdrafts also experience evaporation. The amount of evaporation
!  is simply given by: 
!    drv_evap = rh_melt_target*rs_down - rv_down
! 
!  Aug 2016: tried initiating melting downdrafts at two different levels but
!    didn't work.
! 
!  mass_melt_ratio:
!  - used here: mass_down =  mass_melt_ratio*mass_ansnow
!  - smaller ratio should mean more likely to go down, but smaller mass flux
!  - So presumably, usually get more downward motion and greater cooling, for
!    larger mass_melt_ratio (provided air is dry enough that evaporation is
!    also strong enough to ensure negative buoyancy).
!  - This directly controls the lapse rate max at the melting level. It is
!    unlikely that mass_melt_ratio needs to be larger than 200. 300 is too big.
!    This bend may have a tendency to get bigger with time.
!
!  rh_melt_target
!  - rh_diff = rh_melt_target - rh_down
!  - can be tuned to give roughly correct 600 hPa RH max at t = +2 hours
!
!  b_precip_melt
!  - buoyancy added to downdraft from precip to help go down
!  - Maintains penetrative melting downdrafts at high RH. If value
!    is more negative can drive ML RH to higher value at high rain. A value
!    of -0.02 seems to be required. This probably means a local rain rate of
!    about 500 mm/day, using rho = 0.5 kg/m3.
!
!  iter_evap
!  - the number of iterations in the evap loop.
!  - equals max number of layers a melting downdraft parcel can move down.
!  - iter_evap = 3: not good (Aug 2016)
!
!  MMM
!
!------------------------------------------------------------------------
!integer, parameter :: switch_melt = 0      ! turn off melting
integer, parameter :: switch_melt = 1     ! default: melting done in downdrafts
real(r8), parameter :: mass_melt_ratio = 200.
real(r8), parameter :: b_precip_melt = -0.02  ! 
!real(r8), parameter :: b_precip_melt = 0.20

!real(r8), parameter :: rh_melt_target = 0.90
real(r8), parameter :: rh_melt_target = 0.84
!real(r8), parameter :: f_evap_melt_max = 0.10
!real(r8), parameter :: f_evap_melt_max = 0.05
real(r8), parameter :: f_evap_melt_max = 0.0
integer, parameter :: iter_evap = 2

!------------------------------------------------------------------------
!
!  **************   MELTING PARAMETER  **************************
!
!  Temperature at which start to melt
!  - justified in setting t_melt colder than melting temperature, since
!    each grid cell spans a range of T (+- 3).
!
!------------------------------------------------------------------------
real(r8), parameter :: t_melt = 270.
!real(r8), parameter :: t_melt = 271.

!------------------------------------------------------------------------
!
!  **************   EVAPORATION PARAMETERS  **************************
!
!------------------------------------------------------------------------
!integer, parameter :: nprint_evap = 1
integer, parameter :: nprint_evap = 0
real(r8), parameter :: z_evap_max = 16000.

!------------------------------------------------------------------------
!
!  ************  PENETRATIVE DOWNDRAFT PARAMETERS  ***************
!
!  HARD TO GENERATE PENTRATIVE DOWNDRAFTS
!  - without supersaturating downdraft parcel and adding negative buoyancy
!
!  Explanation:
!  The most important equation is likely:
!     f_evap = rs_down*drh_down*dmdry*rate_down
!  f_evap is the fraction of the rain to be evaporated.
!  It is assumed that this is proportional to the mass of the layer, the rs of
!  the downdraft parcel, and the amount by which the parcel rh is lower than
!  some maximum, and a proprtionality constant rate_down
!
!
!  - Detrains at lowest level at which b_av < 0,
!
!  rate_down_uprain
!  - tried to increase RH of upper Bl (levels 3 and 4) by increasing but didn't 
!    work.  (Jan 2015)
!
!  rate_down_anrain:
!  - controls the amount of anrain allocated to downdrafts. It is exactly
!    analagous to rate_evap_anrain
!  - the ratio (rate_down_anrain/rate_evap_anrain) should roughly reflect the
!    fraction of anrain evaporation used in downdrafts; however evaporation
!    is done before downdrafts, so acts on a larger anrain value
!
!  rh_down_max
!   - if this is too high BL will be too moist during high rain events, will
!     likely have particular difficulty simulating the decrease in BL RH during
!     high rain events.
!   - probably unphysical if exceeds 0.90
!
!  f_mass_max
!    - constraint on the mass_down/dmdry ratio
!    - allowed max size of ratio depends on mass_start_tot (large value 100)
!    - max allowed mass ratio is fff = f_mass_max*mass_start_tot
!    - So f_mass_max should likely be less than 0.01
!    - f_mass = mass_down/dmdry
!     - if (f_mass > fff) then
!    - mass_start_tot is kg/m2*s, but 
!    - Suppose f_mass_max = 0.01. Layers near the surface are 20 hPa; layers
!      near 5 km are 100 hPa. Suppose that a layer is 100 hpa. Then not allowed
!      to remove more than 10 hPa in 30 min. This corresponds to a mass flux of
!      24.*10 = 240 hPa/day. Or a max entrainment rate of 0.01*1000 kg/m2/30 min,
!      which is equal to 0.0055 kg/m2/s.
!
!   fff_min 
!    - min max allowed value of mass_down/dmdry
!    - needed because mass_start_tot = 0 for pure stratiform ansnow, and
!      mass_down was getting set to zero
!
!  f_evap_down_press
!  - puts upper limit on max amount of evaporation that can occur at level
!  - per 100 hPa
!  - potentially important parameter in controlling relative strengths of
!    melting and rain downdrafts.
!
!  md_rat
!  - used here:  mass_down = md_rat*mass_water_evap
!  - smaller md_rat should mean more evaporation per mass, more
!    negative buoyancy, and more "penetrative" downdrafts
!  - decreasing seems to increase T of lower troposphere, but also lower
!    colrh during events.
!  - Oct 2016: increased from 600 to 800 and atm temps seemed to inc.
!
!  b_precip_down
!  - physically this would reflect the buoyancy drag of precipitation
!  - It would enable downdrafts to go downward more easilly, and likely
!    detrain at a lower RH.
!  - Suppose local instaneous rain rate is 500 mm/day (likely reasonable
!    for convective precipitation). This is 500 kg/m2/day, or 0.0058 kg/m2/s.
!    Let R = precipitation density (kg water/m3). Assume rain speed is 5 m/s.
!    Then 0.0058 = R*5, and R = 0.001 kg/m3. Assume rho = 1 kg/m3. Then
!    B from rain = -0.01 m/s2
!  - March 2017: increasing to 0.2 improves the monsoon and results in better
!    rainfall correlations over land. And gets rid of wierd MRG max in asym
!    WK spectrum. But tends to cause warm bias at the surface.
!
!  f_drv_max
!  - maximum allowed fractional increase in rv of a downdraft due to evap
!  - required that be small as a way of preventing downdraft supersaturation
!  - hard to avoid superstauration since T cools as add water.  
!  - If make smaller, may have to increase iter_down
!  - In general, a smaller value should dry out downdrafts, since program
!    more accurately finds minimum amount of evaporation to achieve negative
!    buoyancy.
!
!  DDD
!------------------------------------------------------------------------
real(r8), parameter :: rate_down_uprain = 0.05
real(r8), parameter :: rate_down_anrain = 0.10

real(r8), parameter :: f_mass_max = 0.002
real(r8), parameter :: fff_min = 0.001
!real(r8), parameter :: f_evap_down_press = 0.10
real(r8), parameter :: f_evap_down_press = 0.15

!real(r8), parameter :: rh_down_allow = 0.90  ! max Background RH
real(r8), parameter :: rh_down_allow = 0.88  ! max Background RH
real(r8), parameter :: rh_down_max = 0.95    ! max parcel RH
!real(r8), parameter :: rh_down_max = 0.98    ! max parcel RH
real(r8), parameter :: f_drv_max = 0.05   ! 

real(r8), parameter :: t_down_min = 273.15

!------------------------------------------------------------------------
!
!  bp_down  (b_precip_down)
!  - Oct 2016: including a rain dependence gives a better MJO, but appears
!    to somewhat degrade T profile.
!  - physically this would reflect the buoyancy drag of precipitation
!  - It would enable downdrafts to go downward more easilly, and likely
!    detrain at a lower RH.
!  - Suppose local instaneous rain rate is 500 mm/day (likely reasonable
!    for convective precipitation). This is 500 kg/m2/day, or 0.0058 kg/m2/s.
!    Let R = precipitation density (kg water/m3). Assume rain speed is 5 m/s.
!    Then 0.0058 = R*5, and R = 0.001 kg/m3. Assume rho = 1 kg/m3. Then
!    B from rain = -0.01 m/s2
!
!  iter_down
!  - number of downdraft iterations
!
!  md_rat_uprain
!  - 1200 not imp from 1000
!
!  md_rat_anrain
!  - 800 is not an improvement from 600.
!
!------------------------------------------------------------------------
real(r8), parameter :: b_precip_down = 6.0
!real(r8), parameter :: b_precip_down = 0.20
!real(r8), parameter :: b_precip_down = 0.10
!real(r8), parameter :: b_precip_down = 6.0

!real(r8), parameter :: bp_down_min = 0.0
!real(r8), parameter :: bp_down_add = -0.02
!real(r8), parameter :: bp_down_half = 35.0
!real(r8), parameter :: bp_down_scale = 10.0


integer, parameter :: iter_down = 25

!integer, parameter :: use_bg_rh = 1   ! use BG RH:
integer, parameter :: use_bg_rh = 0   ! use parcel RH

real(r8), parameter :: md_rat_uprain = 1000.
real(r8), parameter :: md_rat_anrain = 600.

!------------------------------------------------------------------------
!
!  drh_down_scale
!  - determines factor that is used to take RH of background atmosphere
!    into account.
!  - increasing makes more linear and shifts downdrafts to drier environment
!    and by making penetrative downdrafts in env in which hm dec with height
!    faster, may warm the atm.
!
!  rh_factor = drh_down_scale*log(1. + (drh_down/drh_down_scale))
!  For drh_down_scale = 0.10:
!  drh_down     rh_factor
!    0.05        0.040
!    0.10        0.069
!    0.20        0.11
!    0.30        0.139
!
!------------------------------------------------------------------------
integer, parameter :: use_rh_scale = 1
real(r8), parameter :: drh_down_scale = 0.10

!------------------------------------------------------------------------
!
!  ************  PENETRATIVE DOWNDRAFT PARAMETERS  ***************
!
!  i_min_down: 
!    - lowest level from which penetrative downdrafts are initiated
!    - whether 3 or 4 doesn't seem to make much difference
!
!------------------------------------------------------------------------
integer, parameter :: i_min_down = 1

!------------------------------------------------------------------------
!
!  ************  PENETRATIVE DOWNDRAFT PARAMETERS  ***************
! 
!  Don't do downdrafts if anrain less than this value
!
!------------------------------------------------------------------------
real(r8), parameter :: anrain_thresh = 0.1

!------------------------------------------------------------------------
!
!  ******************  EVAPORATION PARAMETERS  ******************************
!
!  Evaporation Philosophy (May 2015)
!  ---------------------------------
!  - Likely that intermittent and inefficient deep convection is
!    important in regions of lower colrh to maintain RH
!  - The evaporation scheme must make the fractional evaporation increase
!    with lower RH.
!  - the current scheme does this much better than the old scheme. However, it
!    decreases the mid-level precip around 10 mm/day too much.
!
!  These processes are calculated in the same way but with different 
!    parameter choices
!
!  ANSNOW and UPRAIN evaporation is handled in a similar way, except that
!    rates may be different.
!
!  Issue: is high colrh during events created by the dynamics or by the
!    microphysics? If these parameters are too high (too much evaporation
!    during high rain events), then the high rain heating profile will be
!    too top heavy, and the lower troposphere too cold.
!
!  rate_evap_surf:
!  - surface evaporation to increase surface RH at high rain
!
!------------------------------------------------------------------------
real(r8), parameter :: rate_evap_uprain = 0.00
real(r8), parameter :: rate_evap_ansnow = 0.30  ! to moisten the upper trop
!real(r8), parameter :: rate_evap_ansnow = 0.00  ! to moisten the upper trop
real(r8), parameter :: rate_evap_anrain = 0.00

!real(r8), parameter :: rate_evap_surf = 0.10
real(r8), parameter :: rate_evap_surf = 0.0

!------------------------------------------------------------------------
!
!  ************  EVAPORATION PARAMETERS  ***************
!
!  As the rain rate goes up, a greater fraction of the grid cell will
!    be covered by precip, and evaporation can drive the mean RH to a 
!    higher value.
!
!  June 2015: keep constant unless proven need for variable
!  Uprain evaporation proportional to rh_diff = rh_uprain - rh(i)
!
!  rh_ansnow: used for both anrain and ansnow (evap and downdrafts)
!
!------------------------------------------------------------------------
!real(r8), parameter :: rh_uprain = 0.89  ! used for downdrafts also
!real(r8), parameter :: rh_ansnow = 0.89  ! 
real(r8), parameter :: rh_uprain = 0.90  ! used for downdrafts also
real(r8), parameter :: rh_ansnow = 0.90  ! 

!------------------------------------------------------------------------
!
!  ******************  EVAPORATION PARAMETERS  ******************************
!
!  Max fraction of uprain allowed to be evaporated at a level
!  - this is now an important parameter, especially at low rain rates, where
!    a significant fraction of the precipitation can be evaporated at a level.
!
!------------------------------------------------------------------------
real(r8), parameter :: max_evap = 0.40

!------------------------------------------------------------------------
!
!  Cloud Overlap 
!  -------------
!
!  - modifications made in pkg_cldoptics.F90
!  - Originally modified to make every cloud layer have random overlap.
!    This lowered my solar radiative heating, which helped warm the upper
!    troposphere, and resolved the cold bias. However, completely random
!    overlap means that there often no entirely cloud free areas in a 
!    grid cell even at relatively low cloud fraction. As a result, the
!    upward LW TOA is way too low (about 180 W/m2) even at relatively
!    low cldice. The LW heating of the column is therefore likely too
!    large, and it becomes hard to increase the cloud fraction to get
!    reasonable solar reflectivities.
!
!  - From Manual:
!   "The column is divided into sets of adjacent layers, called regions, 
!   in which the clouds are maximally overlapped.  The clouds are
!   randomly overlapped between different regions.  The number of
!   regions in a column is set by nmxrgn, and the range of pressures
!   included in each region is set by pmxrgn." 
!
!  Feb 2015: I thought that decreasing prob_random would decrease rad
!    heating in the upper troposphere, but doesn't seem to have much 
!    effect.
!
!  March 2015: changed to per km.
!  
!  IF CHANGE HERE MUST CHANGE IN h2.f90
!
!   ARRAY USED in cam4: use reverse order
!
!  March 2015: hoped to decrease HI LW abs by decreasing prob_random in
!    TTL and a bit below but didn't really work. This would have reduced
!    UP LW hitting clouds from below. But maybe need mid-level clouds to
!    block up LW and reduce high LW.
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!   Modifying ozone   OBS: commented out in radiation.F90; no effect
!  How to modify ozone?
!  --------------------
!  - appears to be initialized in 
!    /home/folkins/cesm1_0/models/atm/cam/src/chemistry/utils/prescribed_ozone.F90
!  - But hard to see where the ozone is stored or how to modify
!  - In radiation.F90, there is a call:
!    call rad_cnst_get_clim_gas('O3',  state, pbuf, o3)
!    which seems to define o3 for both shortwave and longwave 
!    computations (which is a pointer)
!  - These parameters below are all used in radiation.F90
!  - Jan 2016: use_shadoz = 1 doesn't seem to make much difference.
!
!------------------------------------------------------------------------
integer, parameter :: nmonths = 12
integer, parameter :: use_shadoz = 0
real(r8) :: oz_mod_7(nmonths) =   &
       (/1.50,  &  ! 1. JAN
         1.65,  &  ! 2. FEB
         1.60,  &  ! 3. MAR
         1.45,  &  ! 4. APR
         1.20,  &  ! 5. MAY
         1.22,  &  ! 6. JUN
         1.35,  &  ! 7. JUL
         1.58,  &  ! 8. AUG
         1.63,  &  ! 9. SEP
         1.74,  &  ! 10. OCT
         1.72,  &  ! 11. NOV
         1.66/)    ! 12. DEC

real(r8) :: oz_mod_8(nmonths) =   &
       (/1.25,  &  ! 1. JAN
         1.30,  &  ! 2. FEB
         1.38,  &  ! 3. MAR
         1.38,  &  ! 4. APR
         1.20,  &  ! 5. MAY
         1.34,  &  ! 6. JUN
         1.60,  &  ! 7. JUL
         1.80,  &  ! 8. AUG
         1.87,  &  ! 9. SEP
         1.88,  &  ! 10. OCT
         1.60,  &  ! 11. NOV
         1.45/)    ! 12. DEC

real(r8) :: oz_mod_9(nmonths) =   &
       (/1.15,  &  ! 1. JAN
         1.14,  &  ! 2. FEB
         1.22,  &  ! 3. MAR
         1.21,  &  ! 4. APR
         1.20,  &  ! 5. MAY
         1.30,  &  ! 6. JUN
         1.42,  &  ! 7. JUL
         1.35,  &  ! 8. AUG
         1.50,  &  ! 9. SEP
         1.42,  &  ! 10. OCT
         1.25,  &  ! 11. NOV
         1.34/)    ! 12. DEC

real(r8) :: oz_mod_10(nmonths) =   &
       (/1.18,  &  ! 1. JAN
         1.19,  &  ! 2. FEB
         1.30,  &  ! 3. MAR
         1.30,  &  ! 4. APR
         1.25,  &  ! 5. MAY
         1.35,  &  ! 6. JUN
         1.30,  &  ! 7. JUL
         1.20,  &  ! 8. AUG
         1.29,  &  ! 9. SEP
         1.30,  &  ! 10. OCT
         1.25,  &  ! 11. NOV
         1.34/)    ! 12. DEC

!------------------------------------------------------------------------
!
!  ***********  DETRAINING CONDENSATE EVAPORATION  *********************
!  
!  Used to define rlp_rv and rip_rv 
!  - fraction of detraining condensate to evap or precip
!  - third column is minimum value of f_precip for values of local
!    RH at and below first column.
!  - important to not keep decreasing f_precip at lower RH so that
!    precip efficiency does not keep going down at lower RH. Need
!    this to avoid low colrh 700-900 hPa cold anomaly.
!
!  How this works:
!  - if local RH > RH of second column all condensate converted to precip
!  - if local RH < RH of first column all condensate evaporated
!  - otherwise linearly interpolate
!  - However, third column may give minimum fraction of condensate allocated
!    to precip
!  - Thinking here is to force more evap in the upper troposphere. This should
!    cool the upper troposphere where I have a warm bias, and effectively heat
!    the lower troposphere by reducing evaporation.
!  
!  Aug 2016: Careful about increasing these as a way of moistening the lower  
!    troposphere, since could create lower trop temp cold bias. 
!  
!  Sept 2016: Seems to be needed only in the upper troposphere, to maintain
!     realistic RH (in lower troposphere, condensate plays a less important
!     role in water vapor budget). Appears also to warm the upper troposphere. 
!     Maybe better to avoid: try to stick to detrainment moistening at low
!     rain rates
!
!  drh_max: appears to help prevent very high RH at 16 km.
!
!  rh_cond_min_bl: 
!  - Non obvious imp from inc from 0.55 to 0.60
!
!------------------------------------------------------------------------

integer, parameter :: use_rh_cond = 1
!integer, parameter :: use_rh_cond = 0

real(r8), parameter :: rh_cond_min_ft = 0.40
real(r8), parameter :: rh_cond_max_ft = 0.65
real(r8), parameter :: rh_cond_min_bl = 0.55   ! i <= 5
real(r8), parameter :: rh_cond_max_bl = 0.80   ! i <= 5

real(r8), parameter :: rh_cond_rain = 0.1

real(r8), parameter :: rh_max_cond = 1.20

!------------------------------------------------------------------------
!
!  prob_random
!
! - April 2016: higher values in the TTL important for increasing LW heating
!   and TTL temps. Increasing in lower troposphere is a way of increasing
!   TOA SW without getting TOA LW too low.
! - Increasing prob_random increases LW heating in the upper troposphere,
!    increases upward motion, and RH. 
! - Low values likely required to prevent large heating rates at 
!   "realistic" ice values, and to prevent downward LW at surface
!   from being too large.
! - LW heating in the mid-troposphere likely sensitive to small
!   changes in prob_random, since upward LW flux is quite big.
! - Jan 16 : reducing pr_min, with no increase in ice, reduced TOA SW
!   significantly.
!
!  calculated in subroutine cldovrlap in pkg_cldoptics.F90
!  subroutine cldefr(l
!      real(r8), intent(in) :: t(pcols,pver)        ! Temperature
!
!------------------------------------------------------------------------

real(r8), parameter :: t_pr_half = 200.0
real(r8), parameter :: t_pr_scale = 3.0
real(r8), parameter :: pr_min = 1.0
real(r8), parameter :: pr_add = -0.90

!------------------------------------------------------------------------
!
!  To reduce random overlap at high column cldice.   
!  - maybe a good think to flatten out cloud heating curve
!  - implemented in pkg_cldoptics.F90
!  - uses tgicewp(i) for IWP
!  - should probably have iwp_add about -0.90. 
!  - Would still want to keep random high in TTL, but likely wouldn't
!    matter much since such high ice is rare.
!  - March 2017: ice reaches a high of 500 g/m2 at high rain.
!    Tried but no benefit. Would have to increase ice to compensate reduced
!    cloud rad heating.
!
!------------------------------------------------------------------------
real(r8), parameter :: iwp_half = 200.0
real(r8), parameter :: iwp_scale = 70.0
real(r8), parameter :: iwp_min = 1.0
!real(r8), parameter :: iwp_add = -0.90
real(r8), parameter :: iwp_add = 0.0

!
!  This would be to make clouds below the ML more
!  random with respect to clouds above
!
integer, parameter :: use_add_lower = 0
real(r8), parameter :: pr_add_inc = 0.030
real(r8), parameter :: t_add_lower = 273.15

!------------------------------------------------------------------------
!
!   cf_max_ri
!   ARRAY USED in cam4: use reverse order
!   - max cloud fraction
!   - should not be larger than cloud_fraction_max (since will be 
!     overwritten anyway).
!   - formula is     
!     log_ri = log(cldice(i,k))
!     tt = (log_ri - log(ri_half))/log(ri_scale)
!     cf_ri = 1. + exp(-tt)
!     cf_ri = cf_max_ri(k)/cf_ri
!
!------------------------------------------------------------------------
real(r8) :: cf_max_ri(nd) =   &
        (/0.40,  &  ! 1. 
          0.40,  &  ! 2. 
          0.40,  &  ! 3. 
          0.40,  &  ! 4. 
          0.40,  &  ! 5. 
          0.40,  &  ! DEF 6. 53 hPa
          0.40,  &  ! 7. 70 hPa
          0.40,  &  ! 8. 85 hPa
          0.40,  &  ! 9. 100 hPa
          0.40,  &  ! 10. 118 hPa
          0.40,  &  ! 11. 139 hPa
          0.40,  &  ! 12. 164 hpa
          0.40,  &  ! 13. 193 hPa
          0.40,  &  ! 14. 227 hPa
          0.40,  &  ! 15. 268 hPa
          0.40,  &  ! 16. 315 hPa
          0.40,  &  ! 17. 370 hPa
          0.40,  &  ! 18. 436 hPa
          0.40,  &  ! 19. 513 hPa
          0.40,  &  ! 20. 604 hPa
          0.40,  &  ! 21. 700 hPa
          0.40,  &  ! 22. 792 hPa
          0.40,  &  ! 23. 872 hPa
          0.40,  &  ! 24. 935 hPa
          0.40,  &  ! 25. 976 hPa
          0.40/)    ! 26. 998 hPa

!------------------------------------------------------------------------
!
!   cf_max_rl
!   ARRAY USED in cam4: use reverse order
!   - should be less than 1, probably decrease with height
!
!------------------------------------------------------------------------
real(r8) :: cf_max_rl(nd) =   &
        (/0.80,  &  ! 1. 
          0.80,  &  ! 2. 
          0.80,  &  ! 3. 
          0.80,  &  ! 4. 
          0.80,  &  ! 5. 
          0.80,  &  ! 6. 53 hPa
          0.80,  &  ! 7. 70 hPa
          0.80,  &  ! 8. 85 hPa
          0.80,  &  ! 9. 100 hPa
          0.80,  &  ! 10. 118 hPa
          0.80,  &  ! 11. 139 hPa
          0.80,  &  ! 12. 164 hpa
          0.80,  &  ! 13. 193 hPa
          0.80,  &  ! 14. 227 hPa
          0.80,  &  ! 15. 268 hPa
          0.80,  &  ! 16. 315 hPa
          0.80,  &  ! 17. 370 hPa
          0.80,  &  ! 18. 436 hPa
          0.80,  &  ! 19. 513 hPa
          0.80,  &  ! 20. 604 hPa
          0.80,  &  ! 21. 700 hPa
          0.80,  &  ! 22. 792 hPa
          0.80,  &  ! 23. 872 hPa
          0.90,  &  ! 24. 935 hPa
          0.90,  &  ! 25. 976 hPa
          0.90/)    ! 26. 998 hPa

!------------------------------------------------------------------------
!
!   END OF PARAMETERS
!
!
!   ***********   START ARRAYS   *****************************************
!
!   Start of flow/vert/tend arrays (entrainment/detrainment)
!   - after some thought, decided it was better to define these arrays here
!     and let other higher subroutines reference as neccessary
!
!  "flow" arrays
!
! dmflow(numdd,nd)               F     [kg/m2/s]
!    Dry Mass entrainment/detrainment
!    numdd = 1: dry mass entrainment
!    numdd = 2: dry mass detrainment
!
!     WATER "FLOWS"
!
!  rvflow(numrv,nd)   F    [kg water/m2/s]
!     1. entrainment  (rv -> rvp)
!     2. detrainment  (rvp -> rv)
!
!  rlflow(numrl,nd)   F    [kg water/m2/s]   
!     1. condensate entrainment  (rl -> rlp)  (liquid water)  
!     2. condensate detrainment  (rlp -> rv)  (liquid water)
!
!  riflow(numri,nd)   F    [kg water/m2/s]   
!     1. condensate entrainment  (rl -> rlp)  (liquid water)  
!     2. condensate detrainment  (rlp -> rv)  (liquid water)
!
!  hmflow(numhm,nd)               F   [J/m2/s] = [W/m2]
!    numhm = 1: hm entrainment : includes a variety of removal mechanisms of hm from a layer
!    numhm = 2: hm detrainment (includes source of hm from uprain evap)
!
!   uflow(nmom,nd)
!   nmom ordering
!     1. u entrainment
!     2. u detrainment
!
!   vflow(nmom,nd)
!   nmom ordering
!     1. v entrainment
!     2. v detrainment
!
!------------------------------------------------------------------------
real(r8) :: dmflow(numdd,nd)
real(r8) :: rvflow(numrv,nd),rlflow(numrl,nd),riflow(numri,nd)
real(r8) :: hmflow(numhm,nd)
real(r8) :: uflow(nmom,nd),vflow(nmom,nd)

!------------------------------------------------------------------------
!
!  "vert" arrays
!
!  dmvert(nd)                     F  [kg/m2/s]
!    change in mass/unit area/time due to vertical advection
!    Calculated from fmass
!
!  rvvert(nd)                     F  [kg vapor/m2/s]
!    change in water vapor mass/unit area/time due to vertical advection
!
!  rlvert(nd)                     F  [kg vapor/m2/s]
!    change in water vapor mass/unit area/time due to vertical advection
!
!  rivert(nd)                     F  [kg vapor/m2/s]
!    change in water vapor mass/unit area/time due to vertical advection
!
!  hmvert(nd)                     F  [J/m2/s]
!    change in mse/unit area/time due to vertical advection
!
!  uvert(nd)
!    tendency in u momentum due to vertical advection
!
!  vvert(nd)
!    tendency in v momentum due to vertical advection
!
!
!------------------------------------------------------------------------
real(r8) :: dmvert(nd)
real(r8) :: rvvert(nd),rlvert(nd),rivert(nd)
real(r8) :: hmvert(nd)
real(r8) :: uvert(nd),vvert(nd)

!------------------------------------------------------------------------
!
!  For Diagnostics Only
!
!------------------------------------------------------------------------
!real(r8) :: rlp_dett(nd)
real(r8) :: updet_1(nd)
real(r8) :: updet_2(nd)
real(r8) :: drv_dndet(nd)
real(r8) :: dt_dndet(nd)
real(r8) :: rh_dn(nd)
real(r8) :: mass_rh_dn(nd)
real(r8) :: dndet(nd)
real(r8) :: dnent(nd)
real(r8) :: mtdet(nd)
real(r8) :: mtent(nd)
real(r8) :: drv_an(nd)
real(r8) :: drv_up(nd)

real(r8) :: dp_dn

!------------------------------------------------------------------------
!
!  Diagnostics
!  - vertical motions in background atmosphere should not create hm/dry mass
!
!------------------------------------------------------------------------
real(r8) :: sum_hmvert,sum_dmvert
real(r8) :: sum_rlvert,sum_rvvert,sum_rivert

!------------------------------------------------------------------------
!
!   "mass" arrays
!
!   fmass(nd+1)         Half
!
!   Net Updraft/Downdraft Dry Mass Fluxes        [kg/m2*s]
!
!------------------------------------------------------------------------
real(r8) :: fmass(nd+1)
real(r8) :: fmass_dn(nd+1)   ! purely diagnostic
!------------------------------------------------------------------------
!
!  Cloud Prod/Loss Variables - Full Level
!
!  rl_evap(nd) - rv production from rl evap [kg water/m2*s]
!  ri_evap(nd) - rv production from ri evap [kg water/m2*s]
!  ri_prod(nd) - ri production from rv cond [kg water/m2*s]
!  rl_prod(nd) - rl production from rv cond [kg water/m2*s]
!
!------------------------------------------------------------------------
real(r8) :: rl_evap(nd)
real(r8) :: ri_evap(nd)
real(r8) :: ri_prod(nd)
real(r8) :: rl_prod(nd)
!------------------------------------------------------------------------
!
!  Diagnostics Cloud Prod/Loss Variables 
!  - Not needed for calculational purposes
!  - units are kg/kg/day rl/ri prod/loss
!  - These should sum to rltend/ritend
!
!------------------------------------------------------------------------
real(r8) :: uprain_rl_tend(nd)  ! evap
real(r8) :: rl_evap_tend(nd)  ! evap
real(r8) :: ri_evap_tend(nd)  ! evap
real(r8) :: ri_prod_tend(nd)  ! condensation
real(r8) :: rl_prod_tend(nd)  ! condensation
real(r8) :: rl_det_tend(nd)   ! conv detrainment
real(r8) :: ri_det_tend(nd)   ! conv detrainment 
real(r8) :: rl_vert_tend(nd)  ! vertical advection
real(r8) :: ri_vert_tend(nd)  ! vertical advection

!------------------------------------------------------------------------
!
!  Precip Prod/Loss Variables - Full Level
!  ---------------------------------------
!
!
!  Four Rain Sources:
!  ------------------
!  uprain_rlp(nd) - [kg water/m2*s]
!  uprain_rv(nd) - [kg water/m2*s]
!  uprain_rl(nd) - [kg water/m2*s]
!  ansnow_rlp(nd) - ansnow production from updrafts [kg water/m2*s]
!  ansnow_rv(nd) - ansnow production from rv cond [kg water/m2*s]
!  ansnow_ri(nd) - ansnow production from ri [kg water/m2*s]
!
!  2 ansnow sinks
!  --------------
!  anrainevap(nd) - evaporation of melted ansnow at level n [kg water/m2*s]
!
!  ansnowmelting
!  -------------
!  ansnowmelt(nd) - melting of ansnow at level n [kg water/m2*s]
!
!------------------------------------------------------------------------
real(r8) :: uprain_rlp(nd)
real(r8) :: uprain_rv(nd)
real(r8) :: uprain_rl(nd)
real(r8) :: ansnow_rlp(nd)
real(r8) :: ansnow_rv(nd)
real(r8) :: ansnow_ri(nd)
real(r8) :: anrainevap(nd)
real(r8) :: ansnowsubl(nd)
real(r8) :: anraindown(nd)
real(r8) :: upraindown(nd)
real(r8) :: uprainevap(nd)
real(r8) :: ansnowmelt(nd)

real(r8) :: rlp_rvv(nd)
!------------------------------------------------------------------------
!
!  uprain - updraft rain [kg water/m2*s]
!
!------------------------------------------------------------------------
real(r8) :: uprain_start_rlp  ! uprain from rlp -> uprain and rl -> uprain
real(r8) :: uprain_surf_rv
real(r8) :: uprain_surf_rlp
real(r8) :: uprain_rlp_rvv
!------------------------------------------------------------------------
!
!  ansnow variables
!
!  ansnow_conv - no evap/melt stratiform snow [kg water/m2*s]
!                 vertical sum of ansnow_rlp(nd)
!                 is all snow even if produced at
!                 T > 273 K
!
!  ansnow_strat - initial (before melt/subl) stratiform snow [kg water/m2*s]
!                 vertical sum of ansnow_rv + ansnow_ri; 
!
!  ansnow_strat_rv - ansnow production from "excess" rv
!
!  ansnow_strat_ri - ansnow production from "excess" ri
!
!------------------------------------------------------------------------
real(r8) :: ansnow_conv
real(r8) :: ansnow_strat
real(r8) :: ansnow_strat_rv
real(r8) :: ansnow_strat_ri
real(r8) :: ansnow_start   ! = ansnow_conv + ansnow_strat
real(r8) :: ansnow_surf
real(r8) :: anrain_surf
!------------------------------------------------------------------------
!
!  Diagnostic variable
!
!------------------------------------------------------------------------
real(r8) :: rlp_evap   ! evaporation of rlp due to detrainment 
!------------------------------------------------------------------------
!
!  MSE of precipitation
!  --------------------
!
!  - These are the actual hm of the precipitation, i.e. continuously updated
!    as uprain,ansnow evaporate. 
!  - Due to the definition of hm used, it is actually usually a negative number
!
!   hmansnow_surf  [hmansnow] = [J/m2*s] = W/m2
!
!  hmansnow:
!    - produced during precipitation
!      hmansnow_prod(ii) = hmansnow_prod(ii) + (dhm_an*mpp/tstep)
!                          NEED hmansnow_start?
!    - During ri to ansnow: 
!      subroutine ri_remove
!      hm_of_ice = cl*t(i) - fusion(t(i)) + g*z(i)
!      hm_ice = hm_of_ice*dri_ansnow*dmdry/tstep  ! removal of hm from a layer
!      hmansnow_start = hmansnow_start + hm_ice   WHY?
!      hmansnow_prod(i) = hmansnow_prod(i) + hm_ice
!    - subroutine rv_to_ansnow_or_ri
!      hm_of_ice = cl*t(i) - fusion(t(i)) + g*z(i)
!      hm_ice = hm_of_ice*drv_ansnow*dmdry/tstep ! removal of hm from a layer
!      hmansnow_start = hmansnow_start + hm_ice   WRONG?
!      hmansnow_prod(i) = hmansnow_prod(i) + hm_ice
!    - subroutine evap
!      hmansnow defined
!
!  CHECK TO SEE THAT 3 subroutines done then evap
!  - define hmansnow_start only after 3 subroutines.
!
!------------------------------------------------------------------------
real(r8) :: hmuprain_start
real(r8) :: hmuprain_surf
real(r8) :: hmansnow_start 
real(r8) :: hmansnow_surf  
real(r8) :: hmanrain_surf  
real(r8) :: hmansnow_prod(nd)
real(r8) :: hmuprain_prod(nd)
!------------------------------------------------------------------------
!
!   Miscellaneous (Diagnostics)
!
!------------------------------------------------------------------------
real(r8) :: diffcon(nd)       !  change in condensate due to precip formation
!------------------------------------------------------------------------
!
!   Starting updraft parcel arrays
!
!   - could probably find some way to do this better since only small
!     number of routines accessing
!
!------------------------------------------------------------------------
real(r8) :: tstart(maxlev,nlaunch)
real(r8) :: rvstart(maxlev,nlaunch),kmstart(maxlev,nlaunch)
real(r8) :: conv_fraction(maxlev,nlaunch)
real(r8) :: lnb(maxlev,nlaunch)
!------------------------------------------------------------------------
!
!   Parcel variables
!
!   - defining here as global to if_conv_tend module; probably not the best.
!
!------------------------------------------------------------------------
real(r8) :: tp(nd),tvp(nd),tdp(nd),rvp(nd),rlp(nd),uwindp(nd),vwindp(nd)
real(r8) :: rip(nd),rsp(nd),hmp(nd),mp(nd),rhp(nd)
!------------------------------------------------------------------------
!
!  Updraft buoyancy arrays
!
!  - these parcel variables are used mainly in if_conv_tend, so really belong there.
!  - However they are used for a diagnostics routine that is taken out of
!    module if_conv_tend and put in in convect.f90, since it is not needed in
!    cam4. Sometime causes compile problems, however.
!
!------------------------------------------------------------------------
real(r8) :: bp_start(nd)  !  updraft parcel buoyancy
real(r8) :: bp_precip(nd)   !  updraft parcel buoyancy
real(r8) :: bp_iter(nd)  !  updraft parcel buoyancy
real(r8) :: bp_freeze(nd)   !  updraft parcel buoyancy
real(r8) :: f_ent_1(nd)
real(r8) :: f_det_1(nd)
real(r8) :: f_ent_2(nd)
real(r8) :: f_det_2(nd)
real(r8) :: num_f_ent_1(nd)
real(r8) :: num_f_det_1(nd)
real(r8) :: num_f_ent_2(nd)
real(r8) :: num_f_det_2(nd)
real(r8) :: bp1_str(nd),bp1_prp(nd),bp1_itr(nd)
real(r8) :: num_bp1_str(nd),num_bp1_prp(nd),num_bp1_itr(nd)
real(r8) :: bp2_str(nd),bp2_prp(nd),bp2_itr(nd)
real(r8) :: num_bp2_str(nd),num_bp2_prp(nd),num_bp2_itr(nd)
real(r8) :: bp_cape(nd)

real(r8) :: hm_diff(nd)           !  change in parcel hmp from hm_start
real(r8) :: hm_diff_calc(nd)      !  change in parcel hmp from hm_start
real(r8) :: hmp_mean(nd)          !  
real(r8) :: hmp_mean_calc(nd)          !  
!------------------------------------------------------------------------
!
!  Variables needed by a several subroutines so better here.
!
!------------------------------------------------------------------------
real(r8) :: updraft_eff,updraft_det,updraft_rlp
real(r8) :: col_ri,col_rl,colrh,av_fr
real(r8) :: bwork_down1,col_rl_new,col_rl_prod,colwat,mf_inc
real(r8) :: av_z_ri,dmdry,summ,summ_dt,summ_rh
real(r8) :: av_z_ri_rh,av_z_ri_dt
real(r8) :: precip_decay_ratio
real(r8) :: surf_rain,cape_mean,cape_dp,conv_energy,rh_fact
real(r8) :: capp(maxlev,nlaunch)
real(r8) :: cape_1, cape_2, cape_3, cape_4
real(r8) :: cin_reg,cf_rl
integer :: i_ml,i_strat,ns_got,num_inc,ifconvection
!------------------------------------------------------------------------
!
!  get_rv_call: helps determine where call to get_rv is coming from.
!
!------------------------------------------------------------------------
integer :: get_rv_call
integer :: get_t_call
!------------------------------------------------------------------------
!
!  Updraft Detrain Diagnostics
!
!------------------------------------------------------------------------
real(r8) :: tdiff_detrain_up(nd),mass_detrain_up(nd)
real(r8) :: cape_detrain_up(nd)
real(r8) :: cape_detrain,hm_diff_detrain
real(r8) :: hm_diff_det(nd)

!------------------------------------------------------------------------
!
!  Amounts of precip evaporated or to downdrafts
!  [kg water/m2/s]
!
!------------------------------------------------------------------------
real(r8) :: anrain_down, uprain_down, anrain_evap, ansnow_subl, ansnow_melt
real(r8) :: uprain_evap
real(r8) :: uprain_1    ! starting uprain
real(r8) :: uprain_2    ! starting uprain
!------------------------------------------------------------------------
!
!  Precipitation Diagnostics 
!  - from calls to precipitation both from updrafts and detrain
!
!------------------------------------------------------------------------
real(r8) :: mass_rip_rv
real(r8) :: mass_rlp_rv
real(r8) :: mass_rip_an
real(r8) :: mass_rlp_up
real(r8) :: mass_rlp_an
!------------------------------------------------------------------------
!
!  Variables needed to monitor mean starting MSE of updrafts
!
!------------------------------------------------------------------------
real(r8) :: mass_start_1    ! has units kg/m2
real(r8) :: mass_start_2    ! has units kg/m2
real(r8) :: mass_start_tot    ! has units kg/m2
real(r8) :: hm_start    ! actual hm
real(r8) :: km_inc
real(r8) :: cfraction_1,cfraction_2,cfraction_3,cfraction_4

!------------------------------------------------------------------------
!
!  ??
!
!------------------------------------------------------------------------
real(r8) :: f_1,f_2,f_3,f_4

!------------------------------------------------------------------------
!
!  variables
!  
!------------------------------------------------------------------------
real(r8) :: z_peak,z_peak_rain,ri_base
real(r8) :: z_peak_width
real(r8) :: z_peak_width_above
real(r8) :: km_width
real(r8) :: a(nlaunch)
real(r8) :: tke_amp


!------------------------------------------------------------------------
!
!  Rain variables
!  RRR
!------------------------------------------------------------------------
real(r8) :: zpeak
real(r8) :: dz_ri_bot
real(r8) :: f_ice
real(r8) :: f_an
real(r8) :: amp,amp_om5,amp_rain
real(r8) :: md_rat

!------------------------------------------------------------------------
!
!  Cloud Ice variables
!
!------------------------------------------------------------------------
real(r8) :: ri_dt(nd)
real(r8) :: ri_rh(nd)
real(r8) :: ri_tar(nd)
real(r8) :: rh_max(nd)
real(r8) :: rh_min(nd)
real(r8) :: ri_added(nd)
real(r8) :: tscale_ri_rem(nd)

real(r8) :: col_ri_rh
real(r8) :: col_ri_dt
real(r8) :: col_rl_rh
real(r8) :: col_rl_dt
real(r8) :: tot_dmdry

!------------------------------------------------------------------------
!
!  Cloud water variables
!
!------------------------------------------------------------------------
real(r8) :: rl_dt_o
real(r8) :: rl_dt_l
real(r8) :: rl_rh_o
real(r8) :: rl_rh_l
real(r8) :: rl_dt(nd)
real(r8) :: rl_rh(nd)
real(r8) :: rl_tar(nd)
real(r8) :: rl_added(nd)

!------------------------------------------------------------------------
!
!  Various
!
!------------------------------------------------------------------------
real(r8) :: tt1, tt2, AA, BB, CC, xxx, yyy, zzz, tot, cee, om_rh_factor
real(r8) :: ttt, aaa, bbb, ppp, lwp_tar, cf_tar, rh_high,lwp_err
real(r8) :: cinf, ri_100, ri_amp_100, theta_7
real(r8) :: LTS, EIS, theta_700, z_700, z_lcl, gamma_850, p1, p2, z1, z2
real(r8) :: t1, t2, t_700, theta_0, mse_rat, ce_fact
integer :: nprint_ri,i_cf

!------------------------------------------------------------------------
!
!  Nonlinear variables: 
!
!------------------------------------------------------------------------
real(r8) :: precip_org_new  ! new value defined here for use in next timestep
real(r8) :: cm5_new  ! new value defined here for use in next timestep

!------------------------------------------------------------------------
!
!  prev variables: 
!
!   - These are defined in cam4 in if_conv.F90 and if_convr.F90 using pointers
!     to retain value from previous timestep
!   - should also be defined if running in 1D (but don't think they are)
!
!------------------------------------------------------------------------
real(r8) :: precip_org
real(r8) :: mass_prev
real(r8) :: cm5_prev

!------------------------------------------------------------------------
!
!  Various
!
!------------------------------------------------------------------------
integer :: have_precip
!------------------------------------------------------------------------
!
!  Counters
!  - get_rv is called an incredible number of times per timestep, and probably
!    where most of the time spent. Check to see where calls are comping from.
!
!------------------------------------------------------------------------
integer :: call_get_rv_cape 
integer :: call_get_rv_entrain
integer :: call_get_rv_precipitate
integer :: call_get_rv_up

integer :: count_tend_start
integer :: count_tend_end
integer :: count_rate, count_max
!------------------------------------------------------------------------
!
!  contains
!
!------------------------------------------------------------------------
contains
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!   Driver for getting convective tendencies
!   Given the column variables:
!    (1) init_zero (initializes various convective arrays to zero)
!    (3) Computes updraft tendencies
!    (4) Checks if updraft tendencies too large
!    (5) Computes downdraft tendenecies
!
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
!
!  subroutine get_tend
!
! - Feb 2014: at the current time, get_tend is called from cam4 at every
!   grid point, irrespective of whether convective instability is present
!   or not
!
!------------------------------------------------------------------------
subroutine get_tend(did_convect,land_fraction,om4,om5,om6,this_lat,tstep)
!------------------------------------------------------------------------
!
!  What to use
!   - may want to specify what from if_conv what to use
!     (actually not; better to specify what to use in the module, I think)
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
use if_conv_solvers, only :  latent,fusion,g,rsat
implicit none
!------------------------------------------------------------------------
!
!  In/Out Variables
!
!------------------------------------------------------------------------
integer, intent(out) :: did_convect
real(r8), intent(in) :: land_fraction
real(r8), intent(in) :: om4
real(r8), intent(in) :: om5
real(r8), intent(in) :: om6
real(r8), intent(in) :: this_lat
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local variables
!
!------------------------------------------------------------------------
integer :: nl,ni,i,i_surf,ilnb,n,itt,it
real(r8) :: cin,dmdry,cape_reg,dm_div,tot_outflow,gauss
real(r8) :: dpp,fff,om_av
real(r8) :: rvs_ml,rh_start,xx,rh_use
integer :: counter_print 
real(r8) :: ri_fract_mean,ri_fract_amp,year_fraction,ri_fract_ttl,ri_fract

!------------------------------------------------------------------------
!
!   *************   START   **********************************
!
!
!  call init_zero 
!
!    - initializes various convective arrays/tendencies to zero
!    - must be done for all flow/vert arrays, so if no updrafts/downdrafts
!      these arrays are still defined
!    - same applies to any diagnostic
!    - sets have_precip = 0
!
!------------------------------------------------------------------------
! print*,'starting get_tend nstep t(1) = ',nstep,t(1)
! print*,'calling init_zero'
  call init_zero

!------------------------------------------------------------------------
!
!  Determines 
!    - i_ml,colrh 
!
!------------------------------------------------------------------------
! print*,'calling get_stuff'
  call get_stuff

!------------------------------------------------------------------------
!
!  Calculate bulk cape values, cape_1, cape_2, cape_mean, etc
!  - These are just for diagnostics
!  - should put in subroutine
!  - assumes km(i) is known
!
!------------------------------------------------------------------------
  cape_dp = 0.
  cape_mean = 0.
  conv_energy = 0.
  do i_surf = 1,maxlev
    if (t(i_surf) > t_min_cape) then
!     print*,'CALLING CAPE from MAIN i_surf = ',i_surf
      call cape(t(i_surf),rv(i_surf),i_surf,ilnb,cape_reg)
      if (i_surf == 1) cape_1 = cape_reg
      if (i_surf == 2) cape_2 = cape_reg
      if (i_surf == 3) cape_3 = cape_reg
      if (i_surf == 4) cape_4 = cape_reg
      if (cape_reg > 0.) then
        dpp = 0.01*(phalf(i_surf) - phalf(i_surf+1))
        cape_dp = cape_dp + dpp
        cape_mean = cape_mean + cape_reg*dpp
        conv_energy = conv_energy + (cape_reg*100.*dpp/g)
      endif
    endif
  end do

!rint*,'---------- GGG -----'
!print*,'cape_1 = ',cape_1
!print*,'cape_2 = ',cape_2
!print*,'cape_2 = ',cape_3
!print*,'cape_4 = ',cape_4

  if (cape_dp > 1.) then
    cape_mean = cape_mean/cape_dp
  endif

  did_convect = 1

!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
ifconvection = 0
if (modify_cldwat == 1) then
  if (cape_mean > cape_mean_activate) ifconvection = 1
  if (abs(this_lat) < lat_cloud) ifconvection = 1
! print*,'------ in if_conv_tend.f90  this_lat = ',this_lat
! print*,'cape_mean = ',cape_mean
! print*,'ifconvection = ',ifconvection
endif

if (ifconvection == 0) then
! print*,'Returning from if_conv.f90: should be no changes in clouds'
! print*,'------------------------------------------'
  return
else
! print*,'Continuing with convection in if_conv_tend.f90'
! print*,'------------------------------------------'
endif

!------------------------------------------------------------------------
!
!  What to do when conv_energy is very high?
!  - May not cause a problem.
!  - could convert all plumes to deep.
!  - shut off convection (not realistic).
!  - This problem is likely now avoided by having plumes in which entrainment
!    is limited. The usual problem is that plumes would be strongly entraining
!    and give rise to extremely high surface temperatures.
!
!------------------------------------------------------------------------
if (conv_energy > conv_energy_danger) then
  print*,'----- DANGEROUS conv_energy or dangerous near surface T ----------'
  print*,'nstep = ',nstep
  print*,'conv_energy = ',conv_energy
  print*,'conv_energy_danger = ',conv_energy_danger
  print*,'land_fraction = ',land_fraction
  print*,'cape_1 = ',cape_1
  print*,'cape_2 = ',cape_2
  print*,'cape_3 = ',cape_3
  print*,'cape_4 = ',cape_4
  print*,'t(1) = ',t(1)
  print*,'t(2) = ',t(2)
  print*,'t(3) = ',t(3)
  print*,'t(4) = ',t(4)
  print*,'rh(1) = ',rh(1)
  print*,'rh(2) = ',rh(2)
  print*,'rh(3) = ',rh(3)
  print*,'rh(4) = ',rh(4)
  print*,'colrh = ',colrh
! do nl = 1,nlaunch
!   a(nl) = 0.
! end do
  print*,'----------------------------------------'
endif

!------------------------------------------------------------------------
!
!  Define ORG variables
!  PPP
!
!------------------------------------------------------------------------
om_av = -0.33333*(om4+om5+om6)

if (use_f_an_om == 0) then
  f_an = sigmoidal(precip_org,f_an_half,f_an_scale,f_an_min,f_an_add)
else
  if (om_av < 0.) then
    f_an = 0.
  else
    f_an = sigmoidal(om_av,f_an_om_half,f_an_om_scale,f_an_om_min,f_an_om_add)
  endif
endif 
if (f_an < 0.) f_an = 0.
if (f_an > 1.) f_an = 1.0
!print*,'om_av = ',om_av
!print*,'f_an = ',f_an

km_width = sigmoidal(precip_org,km_half,km_scale,km_min,km_add)
if (use_km_width_land == 0) km_width = km_width*(1.-land_fraction)

tke_amp = sigmoidal(precip_org,tke_half,tke_scale,tke_low,tke_add)

if (use_om5 == 0) then
  amp = sigmoidal(precip_org,amp_half,amp_scale,amp_low,amp_add)
  if (amp < 0.0) amp = 0.
!-------------------------------------------------------
!  For xxx want net up.
!  - om5 is positive down, so use -1
!  - cm5_prev is negative down so keep as is.
!-------------------------------------------------------
else
  fff = sigmoidal(precip_org,fff_half,fff_scale,fff_low,fff_add)
  amp_om5 = fff*sigmoidal(-om5,amp_om5_half,amp_om5_scale,amp_om5_low,amp_om5_add)
  amp_rain = sigmoidal(precip_org,amp_rain_half,amp_rain_scale,amp_rain_low,amp_rain_add)
  if (amp_om5 < 0.0) amp_om5 = 0.
  if (amp_rain < 0.0) amp_rain = 0.
  amp = amp_om5 + amp_rain
  if (precip_org > 1.) then
!   print*,'-----------'
!   print*,'precip_org = ',precip_org
!   print*,'om5 = ',om5
!   print*,'amp_om5 = ',amp_om5
!   print*,'amp_rain = ',amp_rain
!   print*,'fff = ',fff
!   print*,'amp = ',amp
  endif
endif


if (use_om5_a == 0) then
  a(1) = sigmoidal(precip_org,a1_half,a1_scale,a1_low,a1_add)
else
  a(1) = sigmoidal(-om5,a1_om5_half,a1_om5_scale,a1_om5_low,a1_om5_add)
! print*,'---------------'
! print*,'-om5 = ',-om5
! print*,'a(1) = ',a(1)
endif

if (a(1) > 1.0) a(1) = 1.0
if (a(1) < 0.0) a(1) = 0.0

a(2) = 1.0 - a(1)

!------------------------------------------------------------------------
!
!  Impose Limits
!
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!
!   Diagnostics: LTS and EIS
!   LTS = theta_700 - theta_0    
!   EIS = LTS - (z_700 - z_LCL) Gamma_850
!   LTS, EIS, theta_700, z_700, z_LCL, gamma_850
!   gamma__850 is the lapse rate of a moist adiabat!!
!
!------------------------------------------------------------------------
do i = 1,nd-1
  p1 = p(i)
  p2 = p(i+1)
  if ((p1 > 70000.).and.(p2 < 70000.)) then
    cinf = (p1 - 70000.)/(p1 - p2)
    t1 = t(i)
    t2 = t(i+1)
    z1 = z(i)
    z2 = z(i+1)
    t_700 = (1.-cinf)*t1 + cinf*t2
    z_700 = (1.-cinf)*z1 + cinf*z2
    theta_700 = t_700*((10./7.)**0.286)
    theta_0 = t(1)*((100000./p(1))**0.286)
    LTS = theta_700 - theta_0
    EIS = LTS - (z_700 - z_lcl)*gamma_850
    mse_rat = hm(6)/hm(1)  ! LEV6 very close to 700 hPa
!   print*,'------------ EIS calculation  -----'
!   print*,'cape_1 = ',cape_1
!   print*,'z_lcl = ',z_lcl
!   print*,'z_700 = ',z_700
!   print*,'t_700 = ',t_700
!   print*,'theta_700 = ',theta_700
!   print*,'theta_0 = ',theta_0
!   print*,'gamma_850 = ',gamma_850
!   print*,'t(5) - t(4) = ',t(5) - t(4)
!   print*,'z(5) - z(4) = ',z(5) - z(4)
!   print*,'LTS = ',LTS
!   print*,'EIS = ',EIS
!   print*,'mse_rat = ',mse_rat
  endif
end do

!------------------------------------------------------------------------
!
!  Determine 
!  - km spectrum of lower layers
!  - assumes km(i) is known
!  - should be extremely fast: simply assigns km values.
!  - here is where km_width should be added
!
!------------------------------------------------------------------------
! print*,'calling get_km_spectrum'
  call get_km_spectrum

!------------------------------------------------------------------------
!
!  Determine 
!  - starting properties of air parcel not trivial due to km_width
!  - rvstart, tstart spectrum
!  - these are needed to find the cape spectrum
!  - this might be time consuming: most of the work here is done in a
!    call to "subroutine start" located in if_conv_solvers.f90. This
!    is a crude brute force program finding t/rv from given fixed RH/p/km.
!
!------------------------------------------------------------------------
! print*,'calling get_start_spectrum'
  call get_start_spectrum

!------------------------------------------------------------------------
!
!  Determine capp(it,nl)
!  - cape spectrum
!  - starting t = tstart(it,nl) and rv = rvstart(it,nl)
!
!------------------------------------------------------------------------
! print*,'calling get_cape_spectrum'
  call get_cape_spectrum

!------------------------------------------------------------------------
!
!  Determines conv_fraction
!  - the fraction of each of the layer for each convective mode that is 
!    removed via convection
!  - can be calculated as soon as capp is known
!  - determines nunstable and nupdrafts
!
!------------------------------------------------------------------------
  call get_conv_fraction(land_fraction,tstep)

!------------------------------------------------------------------------
!
!  Produce uprain_start from rv
!  - mainly to heat and dry the BL via drizzle formation
!
!------------------------------------------------------------------------
  if (switch_rv_uprain == 1) call rv_to_uprain(tstep)

!------------------------------------------------------------------------
!
!   ****************   UPDRAFT IF   *************************************
!
!  nupdrafts = 1 implies there is at least one convective parcel 
!  Call updrafts:
!     - calculates "flow" variables for updrafts
!     - Pretty sure no reason to call updrafts if nupdrafts = 1. No 
!       initialization needed here.
!
!------------------------------------------------------------------------
! print*,'nupdrafts = ',nupdrafts
  if (nupdrafts == 1) then
!   print*,'calling updrafts t(1) = ',t(1)
    call updrafts(land_fraction,tstep)
  endif

!------------------------------------------------------------------------
!
!   *********  START CONDENSATE SECTION  *********************
!
!  - Must be after detrainment profile defined
!  - affects ansnow so must be before evap
!
!  Define ri_dt(i) and rl_dt(i) from detrainment profile
!  - Note that rl_dt is no longer generated from direct detrainment but
!    but from rv -> conversion in rv_to_rl, but make sure rl_dt known first 
!
!------------------------------------------------------------------------
!print*,'Start condensate section'
tot_outflow = 0.
do i = 1,nd
  if (t(i) < t_ice_det) then 
    tot_outflow = tot_outflow + dmflow(2,i)
  endif
end do

if (tot_outflow < tot_outflow_min) tot_outflow = 0.

year_fraction = float(nstep)/(365.*48.)
z_peak_width_above = z_peak_width_above_av +   &
          z_peak_width_above_amp*cos(2.*3.141569*(year_fraction+year_frac_add))

z_peak = land_fraction*z_peak_l + (1.-land_fraction)*z_peak_o
ri_base = land_fraction*ri_base_l + (1.-land_fraction)*ri_base_o

z_peak_rain = sigmoidal(precip_org,z_peak_half,z_peak_scale,z_peak,z_peak_add)

!print*,'z_peak = ',z_peak
!print*,'z_peak_rain = ',z_peak_rain
!print*,'precip_org = ',precip_org

do i = 1,nd
  ri_dt(i) = 0.
  rl_dt_o = 0.
  rl_dt_l = 0.
  if (t(i) < t_ice_det) then 
    if (z(i) > z_peak_rain) then 
      z_peak_width = z_peak_width_above
!   Turn off seasonal variation for levels 7(20) and 8(19)
!   Likely not much effect at these two levels
!     if ((i == 20).or.(i == 19)) then 
!       z_peak_width = z_peak_width_above_av
!     endif
    else
      z_peak_width = z_peak_width_below
    endif
    if (tot_outflow > 1.0E-06) then
      gauss = (z(i) - z_peak_rain)/z_peak_width
      gauss = exp(-gauss*gauss)
      ri_dt(i) = gauss*ri_base*(tot_outflow**outflow_power)
!     print*,'--------'
!     print*,'i z(i) = ',i,z(i)
!     print*,'gauss = ',gauss
!     print*,'i = ',i
!     print*,'z(i) = ',z(i)
!     print*,'----'
    else
      ri_dt(i) = 0.
    endif
  else
    ri_dt(i) = 0.
    dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
    dm_div = (dmflow(2,i)/dmdry)*3600.*24.
    rl_dt_o = rl_o_base*dm_div
    rl_dt_l = rl_l_base*dm_div
  endif
  rl_dt(i) = land_fraction*rl_dt_l + (1.-land_fraction)*rl_dt_o
end do

!------------------------------------------------------------------------
!
!  *************  Condensate Section  *****************
!
!   Define ri_rh and ri_tar
!   - note here that cinf can be > 1
!
!------------------------------------------------------------------------
do i = 1,nd
  if (t(i) > t_ice_det) then
    ri_tar(i) = 0.
    ri_rh(i) = 0.
  else
!------------------------------------------------------------------------
!  t< t_ice_det
!------------------------------------------------------------------------
    if (rh(i) > rh_ri_min) then
      rh_use = rh(i)
      ri_rh(i) = (rh_use - rh_ri_min)*ri_amp*log(1. + (rs(i)/rs_scale))
      if (ri_rh(i) > f_ri_max*rs(i)) ri_rh(i) = f_ri_max*rs(i)
    else
      ri_rh(i) = 0.
    endif
  endif
!------------------------------------------------------------------------
!   Define ri_tar   
!------------------------------------------------------------------------
  ri_tar(i) = ri_rh(i) + ri_dt(i)
end do

!------------------------------------------------------------------------
!
!  BL Cloud
!  - first find cf_tar: needed in cloud_fraction.F90
!  - then find lwp_tar
!  qqq
!
!------------------------------------------------------------------------
cee = conv_energy*0.001
ttt = (log(1.+cee) - log(1.+ce_half))/log(1.+ce_scale)
cf_tar = cff_min + (cff_max/(1. + exp(-ttt)))
if (cf_tar < 0.) cf_tar = 0.
if (cf_tar > 1.) cf_tar = 1.

!print*,cee,cf_tar

aaa = 1./(cf1 - cf0)
bbb = -aaa*cf0
ppp = aaa*cf_tar + bbb
lwp_tar = lwp_amp*log(1. + (ppp/(1.-ppp)))
if (lwp_tar < 0.) lwp_tar = 0.

if ((cf_tar > 0.5).and.(lwp_tar < 10.)) then
  print*,'------  WIERD results  ---'
  print*,'cee = ',cee
  print*,'cf_tar = ',cf_tar
  print*,'aaa = ',aaa
  print*,'bbb = ',bbb
  print*,'ppp = ',ppp
  print*,'ppp/(1.-ppp) = ',ppp/(1.-ppp)
  print*,'log(1. + (ppp/(1.-ppp))) = ',log(1. + (ppp/(1.-ppp)))
  print*,'lwp_amp = ',lwp_amp
  print*,'lwp_tar = ',lwp_tar
  stop 'lwp_tar problem'
endif

!-------------------------------
!  put i_cf one higher that highest RH
!  (makes physical sense)
!-------------------------------
i_cf = 0.
rh_high = 0.0
do i = 1,5
  if (rh(i) > rh_high) then
    rh_high = rh(i)
    i_cf = i 
  endif
end do

if (t(i_cf+1) > t_ice_det) i_cf = i_cf + 1

!print*,'------------'
!print*,'cee = ',cee
!print*,'cf_tar = ',cf_tar
!print*,'lwp_tar = ',lwp_tar
!print*,'rh_high = ',rh_high
!print*,'i_cf = ',i_cf

if (i_cf == 0) stop 'i_cf is zero'

!------------------------------------------------------------------------
!
!  *************  Condensate Section  *****************
!
!   Use rl_max and rh to define rl_tar
!
!------------------------------------------------------------------------

do i = 1,nd
  if (t(i) > t_ice_det) then  
!------------------------------------------------------------------------
!  OCEAN
!------------------------------------------------------------------------
    if (land_fraction < land_fraction_max) then
      if (i == i_cf) then
        dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
        rl_rh(i) = lwp_tar*0.001/dmdry
      else
        rl_rh(i) = 0.
      endif
    else
      rh_fact = rl_rh_base*(rh(i) - rl_rh_min)
      if (rh_fact < 0.) rh_fact = 0.
      ce_fact = ce_amp/(1.+log(1.+(cee/ce_lwp)))
      if (use_ce_land == 0) ce_fact = 1.0
      rl_rh(i) = ce_fact*rh_fact
    endif
!------------------------------------------------------------------------
!  Define 
!------------------------------------------------------------------------
    rl_tar(i) = rl_rh(i) + rl_dt(i)
!------------------------------------------------------------------------
!  t < t_ice_det
!------------------------------------------------------------------------
  else
    rl_rh(i) = 0.
    rl_dt(i) = 0.
    rl_tar(i) = 0.
  endif
end do

!-------------------------------------------------------------------------
!
!  *************  Condensate Section  *****************
!
!   Check if rl_rh and lwp_tar are consistent
!
!-------------------------------------------------------------------------
xxx = 0.
do i = 1,nd
  dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
  xxx = xxx + dmdry*rl_rh(i)*1000.
end do
lwp_err = xxx - lwp_tar

if (t(1) > 285.) then
if (abs(lwp_err) > 1.) then
if (land_fraction < land_fraction_max) then
  print*,'lwp_err = ',lwp_err
  print*,'lwp_tar = ',lwp_tar
  print*,'i_cf = ',i_cf
  print*,'xxx = ',xxx
  print*,'this_lat = ',this_lat
  print*,'land_fraction = ',land_fraction
  do i = 1,nd
    print*,'i rl_rh t = ',i,rl_rh(i),t(i)
  end do
  print*,'----------  moving on --------'
! stop 'lwp_err too big in if-conv_tend.f90'
endif
endif
endif


!-------------------------------------------------------------------------
!
!  *************  Condensate Section  *****************
!
!  Define rh_max,rh_min
!
!-------------------------------------------------------------------------
  do i = 1,nd
    if (p(i) > p_ttl_bot) then
      rh_max(i) = rh_max_trop
      rh_min(i) = rh_min_trop
      tscale_ri_rem(i) = tscale_ri_rem_trop
    elseif (p(i) < p_ttl_top) then
      rh_max(i) = rh_max_ttl
      rh_min(i) = rh_min_ttl
      tscale_ri_rem(i) = tscale_ri_rem_ttl
    else
      cinf = (p(i) - p_ttl_top)/(p_ttl_bot - p_ttl_top)
      rh_max(i) = cinf*rh_max_trop + (1.-cinf)*rh_max_ttl
      rh_min(i) = cinf*rh_min_trop + (1.-cinf)*rh_min_ttl
      tscale_ri_rem(i) = cinf*tscale_ri_rem_trop + (1.-cinf)*tscale_ri_rem_ttl
    endif
  end do

!-------------------------------------------------------------------------
!
!  *************  Condensate Section  *****************
!
!  Define average height of ri entry (before tendencies applied)
!
!-------------------------------------------------------------------------
av_z_ri = 0.
av_z_ri_rh = 0.
av_z_ri_dt = 0.

col_ri_rh = 0.
col_ri_dt = 0.
col_rl_rh = 0.
col_rl_dt = 0.

summ = 0.
summ_dt = 0.
summ_rh = 0.

do i = 1,nd
  dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
!
  av_z_ri = av_z_ri + z(i)*dmdry*ri(i)
  summ = summ + dmdry*ri(i)
!
  av_z_ri_dt = av_z_ri_dt + z(i)*dmdry*ri_dt(i)
  summ_dt = summ_dt + dmdry*ri_dt(i)
!
  av_z_ri_rh = av_z_ri_rh + z(i)*dmdry*ri_rh(i)
  summ_rh = summ_rh + dmdry*ri_rh(i)
!
  col_ri_dt = col_ri_dt + 1000.*dmdry*ri_dt(i)  ! convert to g/m2
  col_ri_rh = col_ri_rh + 1000.*dmdry*ri_rh(i)  ! convert to g/m2
  col_rl_dt = col_rl_dt + 1000.*dmdry*rl_dt(i)  ! convert to g/m2
  col_rl_rh = col_rl_rh + 1000.*dmdry*rl_rh(i)  ! convert to g/m2
end do

!print*,'summ = ',summ
!print*,'summ_dt = ',summ_dt
!print*,'summ_rh = ',summ_rh

if (summ > 0.001) then
  av_z_ri = av_z_ri/summ
else
  av_z_ri = bad
endif

if (summ_dt > 0.001) then
  av_z_ri_dt = av_z_ri_dt/summ_dt
else
  av_z_ri_dt = bad
endif

if (summ_rh > 0.001) then
  av_z_ri_rh = av_z_ri_rh/summ_rh
else
  av_z_ri_rh = bad
endif

!print*,'av_z_ri = ',av_z_ri
!print*,'av_z_ri_dt = ',av_z_ri_dt
!print*,'av_z_ri_rh = ',av_z_ri_rh

!------------------------------------------------------------------------
!
!  Convert rv to ansnow or rv
!
!  - other sources of ansnow calculated later
!  - if ri+ri_added < ri_tar then rv remove is decreased by an appropriate
!    amount, so the TTL does not become starved of water to produce condensate.
!  - need to know ri_tar at this point.
!  - should ALWAYS call this, so can produce ri in non-convective regions
!
!------------------------------------------------------------------------
  if (switch_rv_ansnow == 1) call rv_to_ansnow_or_ri(tstep,this_lat)

!------------------------------------------------------------------------
!
!  *************  Condensate Section  *****************
!
!  Remove excess condensate.
!  - using "target ri/rl profiles
!  - if this process is too efficient could dry out the TTL.
!
!------------------------------------------------------------------------
  call ri_remove(tstep)
  call rl_remove(tstep)
  call rv_to_rl(tstep)

!------------------------------------------------------------------------
!
!  Define uprain_start
!  - sets have_precip = 1 if nonzero uprain
!  - must be called after condensate section since uprain can be produced
!    from precip of excess rl
!
!------------------------------------------------------------------------
  call define_uprain_start

!------------------------------------------------------------------------
!
!  Define ansnow
!  - computes have_precip
!  - sets have_precip = 1 if nonzero ansnow
!
!------------------------------------------------------------------------
  call define_ansnow
! print*,'have_precip now defined = ',have_precip

!------------------------------------------------------------------------
!
!  Precip Evaporation and downdrafts
!
!------------------------------------------------------------------------
  call evap(tstep)

!------------------------------------------------------------------------
!
!  Enhanced BL Mixing (BLL)
!  - has to be after cape calculation and after rain calculation
!
!------------------------------------------------------------------------
  call bl_mix(land_fraction,tstep)

!------------------------------------------------------------------------
!
!  call calc_vert
!
!  - calculates the vert arrays defining vertical advection for various 
!     quantities.
!  - vert arrays needed to find tendencies.
!  - any process that adds or removes mass from a layer should generate
!    vertical motions. Currently, only if have_precip = 1 or nupdrafts = 1.
!  - all vert arrays should be initialized to zero previously
!  - this program does not have to be called for non-zero rl/ri evap only.
!
!------------------------------------------------------------------------
  call calc_vert(tstep)

!------------------------------------------------------------------------
!
!  call conv_tendencies
!
!   - This should be called every time a column passes the minimum threshold
!     for convective instability.
!   - Calculate updraft/downdraft/cloud tendencies
!   - these should be all initialized as zero, so don't have to call this
!     routine if no convection.
!     
!------------------------------------------------------------------------
! print*,'calling conv_tendencies have_precip = ',have_precip
  call conv_tendencies(land_fraction,tstep)

!------------------------------------------------------------------------
!
!  Calculate updraft_eff
!
!------------------------------------------------------------------------
  call calc_eff
!------------------------------------------------------------------------
!
!  Calculate preciporg
!  - always do this whether have convection or not.
!
!------------------------------------------------------------------------
  call calc_org_new(om5)
!------------------------------------------------------------------------
!
!  Print counters
!
!------------------------------------------------------------------------

counter_print = 0
if (counter_print == 1) then
 print*,'--------------------  COUNTERS ---------------------'
 print*,'call_get_rv_cape = ',call_get_rv_cape
 print*,'call_get_rv_entrain = ',call_get_rv_entrain
 print*,'call_get_rv_precipitate = ',call_get_rv_precipitate
 print*,'call_get_rv_up = ',call_get_rv_up
 print*,'cape_1 = ',cape_1
 print*,'cape_2 = ',cape_2
 print*,'cape_3 = ',cape_3
 print*,'cape_4 = ',cape_4
 print*,'conv_energy = ',conv_energy
 print*,'------------------------------------------------------'
endif

!------------------------------------------------------------------------
!
!   end
!
!------------------------------------------------------------------------
! print*,'exiting get_tend t(1) = ',t(1)
end subroutine get_tend
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!                subroutine init_zero
!
!   Sets various arrays to zero that are re-defined every timestep
!   Should be called at the start of every convective timestep
!   Does updrafts and downdrafts
!
!------------------------------------------------------------------------
subroutine init_zero
!------------------------------------------------------------------------
!
!  What to use
!  - should change this : do I really need it?
!
!------------------------------------------------------------------------
use if_conv_solvers
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!  In/Out Variables
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!   Local variables
!
!------------------------------------------------------------------------
integer :: i,it,nl,n,ns,nm,nums,nf,iz,ndz,ipass
real(r8) :: dz_int,dz_min,dz_max

!------------------------------------------------------------------------
!
! calc_vert not always called
!
!------------------------------------------------------------------------
  dpsurf = 0.
  have_precip = 0
  bwork_down1 = 0.
!------------------------------------------------------------------------
!
!  set counters to zero
!
!------------------------------------------------------------------------
  call_get_rv_cape = 0
  call_get_rv_entrain = 0
  call_get_rv_precipitate = 0
  call_get_rv_up = 0
!------------------------------------------------------------------------
!
!  These variables are used to calculate mass weighted initial MSE of
!  updraft parcels.
!
!------------------------------------------------------------------------
  mass_start_1 = 0.
  mass_start_2 = 0.
  mass_start_tot = 0.
  hm_start = 0.
  cin_reg = 0.
  cape_dp = 0.
  updraft_eff = 0.
  updraft_rlp = 0.
  cfraction_1 = 0.
  cfraction_2 = 0.
  cfraction_3 = 0.
  cfraction_4 = 0.
  cape_1 = 0.
  cape_2 = 0.
  cape_3 = 0.
  cape_4 = 0.
  cape_mean = 0.
  conv_energy = 0.

!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
  mass_rip_rv = 0.
  mass_rlp_rv = 0.
  mass_rip_an = 0.
  mass_rlp_up = 0.
  mass_rlp_an = 0.
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
  anrain_down = 0.
  uprain_down = 0.
  uprain_1 = 0.
  uprain_2 = 0.
  anrain_evap = 0.
  uprain_evap = 0.
  ansnow_subl = 0.
  ansnow_melt = 0.

!------------------------------------------------------------------------
!
!  ********** FLOW VARIABLES TO ZERO ****************
!
!   Set all flow variables to zero
!
!------------------------------------------------------------------------
  do i = 1,nd
!------------------------------------------------------------------------
!
!  flow arrays with a "ndraft" and nd arguments
!
!------------------------------------------------------------------------
      do n = 1,numrv
        rvflow(n,i) = 0.
      end do
      do n = 1,numrl
        rlflow(n,i) = 0.
      end do
      do n = 1,numri
        riflow(n,i) = 0.
      end do
      do n = 1,numdd
        dmflow(n,i) = 0.
      end do
      do n = 1,numhm
        hmflow(n,i) = 0.
      end do
      do nm = 1,nmom
        uflow(nm,i) = 0.
        vflow(nm,i) = 0.
      end do
!------------------------------------------------------------------------
!
!  end loop over nd
!
!  ********** END FLOW VARIABLES TO ZERO ****************
!
!------------------------------------------------------------------------
  end do
!------------------------------------------------------------------------
!
!   Set all vert arrays to zero
!
!------------------------------------------------------------------------
  do i = 1,nd
    dmvert(i) = 0.
    rvvert(i) = 0.
    rlvert(i) = 0.
    rivert(i) = 0.
    hmvert(i) = 0.
    uvert(i) = 0.
    vvert(i) = 0.
!   rlp_dett(i) = 0.
    updet_1(i) = 0.
    updet_2(i) = 0.
    dndet(i) = 0.
    dnent(i) = 0.
    mtdet(i) = 0.
    mtent(i) = 0.
    drv_up(i) = 0.
    drv_an(i) = 0.
    drv_dndet(i) = 0.
    dt_dndet(i) = 0.
    rh_dn(i) = 0.
    mass_rh_dn(i) = 0.
  end do
  dp_dn = 0.
!------------------------------------------------------------------------
!
!  Set fmass to zero
!
!------------------------------------------------------------------------
  do i = 1,nd+1
    fmass(i) = 0.
    fmass_dn(i) = 0.
  end do
!------------------------------------------------------------------------
!
!  Set hm of precip equal to zero.
!
!------------------------------------------------------------------------
  hmuprain_start = 0.
  hmuprain_surf = 0.
  hmansnow_start = 0.
  hmanrain_surf = 0.
  hmansnow_surf = 0.
!------------------------------------------------------------------------
!
!  Set uprain to zero.
!
!------------------------------------------------------------------------
  uprain_surf_rlp = 0.
  uprain_surf_rv = 0.
  uprain_start_rlp = 0.
  uprain_rlp_rvv = 0.
!------------------------------------------------------------------------
!
!  Set ansnow to zero.
!
!------------------------------------------------------------------------
  ansnow_conv = 0.
  ansnow_strat = 0.
  ansnow_strat_rv = 0.
  ansnow_strat_ri = 0.
  ansnow_start = 0.
  anrain_surf = 0.
  ansnow_surf = 0.
!------------------------------------------------------------------------
!
!  Vertically integrated rlp evap
!
!------------------------------------------------------------------------
  rlp_evap = 0.
!------------------------------------------------------------------------
!
!  Set precipitation sources/sinks equal to zero (full levels)
!
!------------------------------------------------------------------------
  do i = 1,nd
    uprain_rlp(i) = 0.
    uprain_rv(i) = 0.
    uprain_rl(i) = 0.
    ansnow_rlp(i) = 0.
    ansnow_rv(i) = 0.
    ansnow_ri(i) = 0.
    anraindown(i) = 0.
    upraindown(i) = 0.
    uprainevap(i) = 0.
    anrainevap(i) = 0.
    ansnowsubl(i) = 0.
    ansnowmelt(i) = 0.
    hmuprain_prod(i) = 0.
    hmansnow_prod(i) = 0.
    rlp_rvv(i) = 0.
  end do

!------------------------------------------------------------------------
!
!  Set Cloud sources/sinks equal to zero (full levels)
!
!------------------------------------------------------------------------
  do i = 1,nd
    rl_evap(i) = 0.
    ri_evap(i) = 0.
    ri_prod(i) = 0.
    rl_prod(i) = 0.
    uprain_rl_tend(i) = 0.
    rl_evap_tend(i) = 0.
    ri_evap_tend(i) = 0.
    ri_prod_tend(i) = 0.
    rl_prod_tend(i) = 0.
    rl_det_tend(i) = 0.
    ri_det_tend(i) = 0.
    rl_vert_tend(i) = 0.
    ri_vert_tend(i) = 0.
    ri_added(i) = 0.
    rl_added(i) = 0.
    ri_rh(i) = 0.
    ri_dt(i) = 0.
    ri_tar(i) = 0.
    rl_rh(i) = 0.
    rl_dt(i) = 0.
    rl_tar(i) = 0.
  end do

  lwp_tar = 0.
  cf_tar = 0.

!------------------------------------------------------------------------
!
!  "tend" arrays
!
!------------------------------------------------------------------------
  do i = 1,nd
    kmtend(i) = 0.
    utend(i) = 0.
    vtend(i) = 0.
    rvtend(i) = 0.
    rltend(i) = 0.
    ritend(i) = 0.
  end do

!------------------------------------------------------------------------
!
!  Diagnostics
!
!------------------------------------------------------------------------
  do i = 1,nd

!------------------------------------------------------------------------
! variables which scale with detraining updraft mass flux
!------------------------------------------------------------------------
    tdiff_detrain_up(i) = 0.
    cape_detrain_up(i) = 0.
    mass_detrain_up(i) = 0.
!
    hm_diff_det(i) = 0.

!------------------------------------------------------------------------
! variables which scale with updraft mass flux
!------------------------------------------------------------------------
    av_rlp(i) = 0.
    av_rip(i) = 0.
    av_mp(i) = 0.
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
    hm_diff(i) = bad
    hm_diff_calc(i) = 0.
    hmp_mean(i) = bad
    hmp_mean_calc(i) = 0.

    f_ent_1(i) = 0.
    num_f_ent_1(i) = 0.
    f_det_1(i) = 0.
    num_f_det_1(i) = 0.

    f_ent_2(i) = 0.
    num_f_ent_2(i) = 0.
    f_det_2(i) = 0.
    num_f_det_2(i) = 0.

    bp1_str(i) = 0.
    bp1_prp(i) = 0.
    bp1_itr(i) = 0.
    bp2_str(i) = 0.
    bp2_prp(i) = 0.
    bp2_itr(i) = 0.

    num_bp1_str(i) = 0.
    num_bp1_prp(i) = 0.
    num_bp1_itr(i) = 0.
    num_bp2_str(i) = 0.
    num_bp2_prp(i) = 0.
    num_bp2_itr(i) = 0.

    bp_cape(i) = 0.

  end do

!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
  do it = 1,maxlev
  do nl = 1,nlaunch
    nunstable(it,nl) = 1
    capp(it,nl) = 0.
    lnb(it,nl) = 0.
  end do
  end do

!------------------------------------------------------------------------
!
!   end
!
!------------------------------------------------------------------------
end subroutine init_zero
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine get_conv_fraction
!
!  Defines 
!  - updraft mass spectrum
!  - nunstable(it,nl) should be defined on input
!
!------------------------------------------------------------------------
subroutine get_conv_fraction(land_fraction,tstep)
!------------------------------------------------------------------------
!
!   Use
!
!------------------------------------------------------------------------
use if_conv_solvers, only : rsat,g,hmerr_max,enthalpy
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: land_fraction
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: xxx,tmass,xx,dmdry,error,cin,slope,cinf
real(r8) :: mass_fraction,totalmass,rtt
real(r8) :: small,ratio, rh_diff, mass_new, mass_min, mass_max
integer :: it,i,nl,ilnb,ns

nupdrafts = 0

!------------------------------------------------------------------------
!
!  Loop over it 
!
!------------------------------------------------------------------------
  do it = 1,maxlev
  do nl = 1,nlaunch
!------------------------------------------------------------------------
!
!   Determine conv_fraction 
!   - fractional convective removal of mass of a PARCEL of a layer.
!   - Modify conv_fraction of bottom two levels if high RH.
!
!------------------------------------------------------------------------
if (capp(it,nl) > 0.) then

conv_fraction(it,nl) = amp*a(nl)*(tstep/tscale_cape)*capp(it,nl)/cape_scale

else

conv_fraction(it,nl) = 0.

!------------------------------------------------------------------------
!
! endif for capp > 0.
!
!------------------------------------------------------------------------
endif

!------------------------------------------------------------------------
!
!  end loop over nl
!  end loop over it (from 1 to maxlev)
!
!------------------------------------------------------------------------
  end do
  end do

!------------------------------------------------------------------------
!
!  Check for negative conv_fraction
!
!------------------------------------------------------------------------
do it = 1,maxlev
do nl = 1,nlaunch
  if (conv_fraction(it,nl) < 0.) then
    print*,'conv_fraction(it,nl) = ',conv_fraction(it,nl)
    print*,'capp(it,nl) = ',capp(it,nl)
    stop 'conv_fraction negative'
  endif
end do
end do

!------------------------------------------------------------------------
!
!  To save computation
!
!------------------------------------------------------------------------
do it = 1,maxlev
do nl = 1,nlaunch
  if (conv_fraction(it,nl) < conv_fraction_min) then
    nunstable(it,nl) = 0
    conv_fraction(it,nl) = 0.
  else
    nunstable(it,nl) = 1
    nupdrafts = 1
  endif
end do
end do

!------------------------------------------------------------------------
!
!  Do not let conv_fraction(it,nl) exceed conv_fraction_max
!
!------------------------------------------------------------------------
  do it = 1,maxlev
  do nl = 1,nlaunch
    if (conv_fraction(it,nl) > conv_fraction_max) then
      conv_fraction(it,nl) = conv_fraction_max
    endif
  end do
  end do

!------------------------------------------------------------------------
!
!  Calculate mass_start_tot
!
!------------------------------------------------------------------------
  mass_start_tot = 0.
  do it = 1,maxlev
  do nl = 1,nlaunch
    dmdry = (phalf(it)-phalf(it+1))/(g*(1.+rv(it)+rl(it)+ri(it)))
    mass_start_tot = mass_start_tot + conv_fraction(it,nl)*dmdry
  end do
  end do

!print*,'----------'
!print*,'mass_start_tot = ',mass_start_tot
!print*,'mass_prev = ',mass_prev

!------------------------------------------------------------------------
!
!   Start of three constraints
!   - order of implementation of constraints matters
!
!   1. Do not let mass_start_tot > mass_start_max
!
!------------------------------------------------------------------------
  mass_new = mass_start_tot
  if (mass_start_tot > mass_start_max) then
    mass_new = mass_start_max
    print*,'------- mass_start_tot exceeds limit -------------------------'
    print*,'nstep = ',nstep
    print*,'mass_start_tot = ',mass_start_tot
    print*,'mass_start_max = ',mass_start_max
    print*,'land_fraction = ',land_fraction
    print*,'cape_1 = ',cape_1
    print*,'cape_2 = ',cape_2
    print*,'cape_3 = ',cape_3
    print*,'cape_4 = ',cape_4
    print*,'conv_energy = ',conv_energy
    print*,'t(1) = ',t(1)
    print*,'t(2) = ',t(2)
    print*,'t(3) = ',t(3)
    print*,'t(4) = ',t(4)
    print*,'rh(1) = ',rh(1)
    print*,'rh(2) = ',rh(2)
    print*,'rh(3) = ',rh(3)
    print*,'rh(4) = ',rh(4)
    do it = 1,maxlev
    do nl = 1,nlaunch
     print*,'it conv_fraction(it,nl) = ',it,conv_fraction(it,nl)
    end do
    end do
    print*,'reducing all conv_fraction by ratio = ',ratio
    print*,'--------------------------------'
  endif

!------------------------------------------------------------------------
!
!  2. Do not let allow excessive absolute change in starting mass
!
!------------------------------------------------------------------------
  mass_min = mass_prev - mass_inc_max
  mass_max = mass_prev + mass_inc_max
  if (mass_start_tot < mass_min ) mass_new = mass_min
  if (mass_start_tot > mass_max ) mass_new = mass_max

!  print*,'After ABSOLUTE'
!  print*,'mass_min = ',mass_min
!  print*,'mass_max = ',mass_max
!  print*,'mass_new = ',mass_new

!------------------------------------------------------------------------
!
!  3. Do not let allow excessive relative change in starting mass
!
!------------------------------------------------------------------------
if (mass_prev > mass_start_thresh) then
  mass_min = (1. - mass_rat_max)*mass_prev
  mass_max = (1. + mass_rat_max)*mass_prev
  if (mass_start_tot < mass_min ) mass_new = mass_min
  if (mass_start_tot > mass_max ) mass_new = mass_max
! print*,'After RELATIVE'
! print*,'mass_min = ',mass_min
! print*,'mass_max = ',mass_max
! print*,'mass_new = ',mass_new
endif

!------------------------------------------------------------------------
!
!  Define new conv_fraction and mass_start_tot
!
!------------------------------------------------------------------------
if (mass_start_tot > mass_start_thresh) then
  ratio = mass_new/mass_start_tot
! print*,'mass_new = ',mass_new
! print*,'mass_start_tot = ',mass_start_tot
! print*,'ratio = ',ratio
  if (ratio < 0.) stop 'ratio negative'
  mass_start_tot = mass_new
  do it = 1,maxlev
  do nl = 1,nlaunch
    conv_fraction(it,nl) = ratio*conv_fraction(it,nl)
  end do
  end do
endif

!------------------------------------------------------------------------
!
!   END   (subroutine get_conv_fraction)
!
!------------------------------------------------------------------------
end subroutine get_conv_fraction
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine cape
!
!  Inputs
!  ------
!  - tempstart: starting temperature
!  - rvvstart: starting rv
!  - istart: starting height index in pre-defined vertical grid
!
!  Outputs
!  -------
!  - cape_out: 
!  - ilnb
!
!------------------------------------------------------------------------
subroutine cape(tempstart,rvvstart,istart,ilnb,cape_out)
!------------------------------------------------------------------------
!
!   Use
!
!------------------------------------------------------------------------
use if_conv_solvers, only : rdd,rsat,hmoist,epsi,g,get_rv
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: tempstart
real(r8), intent(in) :: rvvstart
integer, intent(in) :: istart
integer, intent(out) :: ilnb
real(r8), intent(out) :: cape_out
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: xxx,ff,dt,dlnp,small,diff
real(r8) :: dtd,rtp,kmp,t_in,bpp
real(r8) :: rsatstart,riistart,rtstart,z_top,rh_start
integer :: iwant,ngot,i,nprint,allow_sat,i_top,ii,ngot_lcl
!------------------------------------------------------------------------
!
!  Initialize
!
!------------------------------------------------------------------------
cape_out = 0.
ilnb = 0
!------------------------------------------------------------------------
!
!  Determine i_top
!  - This is required, because will generate extremely cold temperatures
!    in the stratosphere when try to calculate cape in Arctic regions, 
!    which will violate imposed "reasonable" temperature range in 
!    subroutine get_rv.
!
!------------------------------------------------------------------------
  z_top = 16000.
  i_top = 0
  ngot = 0
  do i = 1,nd
    if ((ngot == 0).and.(z(i) > z_top)) then
      i_top = i
      ngot = 1
    endif
  end do
  if (ngot == 0) stop 'ngot equals zero in cape'

!------------------------------------------------------------------------
!
!   Make sure rvvstart > 0
!
!------------------------------------------------------------------------
  if (rvvstart < 0.) then
    print*,'negative rvvstart in subroutine cape'
    stop 'negative rvvstart in subroutine cape'
  endif
  if (istart <= 0) stop 'istart incorrect in cape'
!------------------------------------------------------------------------
!
!   Calculate rsatstart
!   Check if rvvstart < rsatstart
!    - probably doesn't matter
!
!------------------------------------------------------------------------
  rsatstart = rsat(tempstart,p(istart))
  small = 1.0E-08
  if ((rvvstart-rsatstart) > small) then
!   print*,'rvvstart = ',rvvstart
!   print*,'rsatstart = ',rsatstart
!   print*,'diff = ',rsatstart-rvvstart
!   print*,'WARNING: supersaturated rvvstart in subroutine cape'
  endif
!------------------------------------------------------------------------
!
!  Initialize parcel properties at level of origin (i = istart)
!  You have to give an arbitrary starting mass for the air parcel.
!
!------------------------------------------------------------------------
  tp(istart) = tempstart
  tdp(istart) = tempstart*(1. + epsi*rvvstart)/(1. + rvvstart)
  tvp(istart) = tempstart*(1. + epsi*rvvstart)/(1. + rvvstart)
  rvp(istart) = rvvstart
  rsp(istart) = rsatstart
  rlp(istart) = 0.        ! Assume no initial condensate
  riistart = 0.
  rtstart = rvvstart
  rip(istart) = 0.
  hmp(istart) = hmoist(tempstart,rtstart,rvvstart,riistart,z(istart))
!------------------------------------------------------------------------
!
!   Diagnostics
!
!------------------------------------------------------------------------
  nprint = 0
  rh_start = rvvstart/rsatstart
  if (nprint == 1) then
  if (tempstart > 300.) then
  if (rh_start > 0.8) then
    print*,'-----  cape starting info ---------------'
    print*,'istart = ',istart
    print*,'tempstart = ',tempstart
    print*,'rh_start = ',rh_start
    print*,'rtstart = ',rtstart
    print*,'rvvstart = ',rvvstart
    print*,'riistart = ',riistart
    print*,'z(istart) = ',z(istart)
    print*,'-----------------------------------------'
  endif
  endif
  endif
!------------------------------------------------------------------------
!
!  Only enter loop if i_top is above start level
!
!------------------------------------------------------------------------
if (i_top >= (istart+1)) then
!------------------------------------------------------------------------
!
!  Determine tdp: parcel density temperature profile
!
!------------------------------------------------------------------------
  do i = istart+1,i_top
!------------------------------------------------------------------------------
! Define quantities for call to get_rv.
!------------------------------------------------------------------------------
    rvp(i) = rvp(i-1)
    rlp(i) = rlp(i-1)
    rip(i) = rip(i-1)
    hmp(i) = hmp(i-1)
    rtp = rvp(i) + rlp(i) + rip(i)
    kmp = hmp(i) - (1. + rtp)*g*z(i)
!------------------------------------------------------------------------------
!  Check for negative kmp
!------------------------------------------------------------------------------
    if ((kmp < 0.).or.(kmp > km_bad)) then
      print*,'======== kmp NEGATIVE or too BIG PROBLEM in cape'
      print*,'kmp = ',kmp
      print*,'i = ',i
      print*,'z(i) = ',z(i)
      print*,'rvp(i) = ',rvp(i)
      print*,'hmp(i) = ',hmp(i)
      print*,'i_top = ',i_top
      print*,'z(i_top) = ',z(i_top)
      print*,'----- Start quantities:-----' 
      print*,'istart = ',istart
      print*,'p(istart) = ',p(istart)
      print*,'z(istart) = ',z(istart)
      print*,'t(istart) = ',t(istart)
      print*,'tp(istart) = ',tp(istart)
      print*,'tdp(istart) = ',tdp(istart)
      print*,'tvp(istart) = ',tvp(istart)
      print*,'rvp(istart) = ',rvp(istart)
      print*,'rsp(istart) = ',rsp(istart)
      print*,'rlp(istart) = ',rlp(istart)
      print*,'hmp(istart) = ',hmp(istart)
      print*,'--- Background Temperature profile: ----'
      do ii = 1,i
        print*,'ii = ',ii
        print*,'z(ii) = ',z(ii)
        print*,'t(ii) = ',t(ii)
        print*,'tp(ii) - t(ii) = ',tp(ii)-t(ii)
        print*,'hmp(ii) - hm(ii) = ',hmp(ii) - hm(ii)
        print*,'rvp(ii) = ',rvp(ii)
        print*,'hm(ii) = ',hm(ii)
      end do 
      stop 'kmp problem'
    endif
!------------------------------------------------------------------------------
!
!  Defining t_in:
!  - get_rv needs an initial guess for a best T
!  - this is used to determine of the unsaturated solution works.
!  - This guess assumes the temperature at the lower level, plus the same lapse
!    rate as the background atm.
!
!------------------------------------------------------------------------------
   t_in = tp(i-1) + t(i) - t(i-1) 
!------------------------------------------------------------------------------
!  Diagnostics
!------------------------------------------------------------------------------
    if (nprint == 1) then
    if (tempstart > 300.) then
    if (rh_start > 0.8) then
      print*,'-----IN cape BEFORE call to get_rv INPUT VALUES'
      print*,'i = ',i
      print*,'z(i) = ',z(i)
      print*,'p(i) = ',p(i)
      if (i > 1) print*,'tp(i-1) = ',tp(i-1)
      print*,'t(i) = ',t(i)
      if (i > 1) print*,'t(i-1) = ',t(i-1)
      print*,'t_in = ',t_in
      print*,'kmp = ',kmp
      print*,'hmp(i) = ',hmp(i)
      print*,'hm(i) = ',hm(i)
      print*,'hmp(i) - hm(i) = ',hmp(i)-hm(i)
      print*,'rvp(i) = ',rvp(i)
      print*,'rlp(i) = ',rlp(i)
      print*,'rip(i) = ',rip(i)
      print*,'----------------------'
    endif
    endif
    endif
!------------------------------------------------------------------------------
!  Call get_rv from cape to get tp(i)
!  get_rv keeps ri fixed at the input value
!  allow_sat = 0: since do not want parcel rh > 1
!------------------------------------------------------------------------------
    get_rv_call = 1
    allow_sat = 0
    call_get_rv_cape = call_get_rv_cape + 1
    call get_rv(kmp,  &   ! pure in
                rtp,  &   ! pure in
               p(i),  &   ! pure in  
             rvp(i),  &   ! pure out
             rlp(i),  &   ! inout
             rip(i),  &   ! in only (fixed)
               t_in,  &   ! inout
        get_rv_call,  &   ! pure in
          allow_sat)     ! pure in
!------------------------------------------------------------------------
!  Make changes to rlp depending on reversible or pseudo.
!------------------------------------------------------------------------
    tp(i) = t_in
    rip(i) = 0.
    if (cape_reversible == 0) then
      rlp(i) = 0.
      rtp = rvp(i) + rlp(i) + rip(i)
      hmp(i) = hmoist(tp(i),rtp,rvp(i),rip(i),z(i)) 
    endif
!------------------------------------------------------------------------
!  Define tdp(i) From tp(i)
!------------------------------------------------------------------------
    tdp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rvp(i) + rlp(i) + rip(i))
    tvp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rvp(i))
    rhp(i) = rvp(i)/rsat(tp(i),p(i))
!------------------------------------------------------------------------------
!  Diagnostics
!------------------------------------------------------------------------------
    if (nprint == 1) then
    if (rh_start > 0.8) then
    if (tempstart > 300.) then
      print*,'----- IN cape AFTER call to get_rv'
      print*,'i = ',i
      print*,'z(i) = ',z(i)
      print*,'p(i) = ',p(i)
      print*,'kmp = ',kmp
      print*,'hmp(i) = ',hmp(i)
      print*,'output tp(i) = ',tp(i)
      print*,'rvp(i) = ',rvp(i)
      print*,'rlp(i) = ',rlp(i)
      print*,'rip(i) = ',rip(i)
      print*,'Saturated rs = ',rsat(tp(i),p(i))
      print*,'RH on output from get_rv = ',rvp(i)/rsat(tp(i),p(i))
      print*,'----------------------'
    endif
    endif
    endif
!------------------------------------------------------------------------
!  End loop over heights.
!------------------------------------------------------------------------
  end do
!------------------------------------------------------------------------
!
!   Determine LNB
!    ilnb: start at top and count down to first level with dtd > 0
!
!------------------------------------------------------------------------
  ilnb = 0
  ngot = 0
  do i = i_top,istart,-1
    dtd = tdp(i) - td(i)
    if ((dtd > 0.0001).and.(ngot == 0)) then
      ilnb = i
      ngot = 1
    endif
  end do
!------------------------------------------------------------------------
!
!   Determine LCL and gamma_850
!
!------------------------------------------------------------------------
if (istart == 1) then
  z_lcl = bad
  ngot_lcl = 0
  do i = istart,nd
    if ((rhp(i) > 0.98).and.(ngot_lcl == 0)) then
      ngot_lcl = 1
      z_lcl = z(i)
    endif
  end do
  gamma_850 = -(tp(5)-tp(4))/(z(5)-z(4))
endif
!------------------------------------------------------------------------
!    
!   Calculate CAPE
!   - Use tdp and ilnb to calculate cape
!
!   CAPE = integral Rd(tdp - td)dlnp
!         dlnp = log(100.*p(i+1))-log(100.*p(i))
!
!   LFC look for first height at which bpp > 0.
!
!------------------------------------------------------------------------
! print*,'Starting cape calc  -----'
  if (ilnb /= 0) then
    do i = istart,ilnb
      bpp = g*(tdp(i) - td(i))/td(i)
      cape_out = cape_out + bpp*(z(i+1) - z(i))
      if (i == nd) stop 'NOT ALLOWED in if_conv_tend.f90'
! CCC
      if (istart == 1) then
        bp_cape(i) = bpp
      endif
!     print*,'z(i) = ',z(i)*0.001
!     print*,'bpp = ',bpp
!     print*,'cape_out = ',cape_out
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
    end do
  endif
! print*,'FINAL CAPE_OUT = ',cape_out
  if (cape_out < 0.) cape_out = 0.
!------------------------------------------------------------------------
!
!   Diagnostics
!
!------------------------------------------------------------------------
  if (nprint == 1) then
  if (tempstart > 300.) then
  if (rh_start > 0.8) then
    print*,'cape_out = ',cape_out
    print*,'-------  leaving cape -----'
  endif
  endif
  endif
!------------------------------------------------------------------------
! endif for i_top > istart+1
!------------------------------------------------------------------------
endif
!------------------------------------------------------------------------
!
!  end
!
!------------------------------------------------------------------------
end subroutine cape
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine get_start_spectrum
!
!  Determines a starting km/rv/t spectrum
!  - rh of all parcels is exactly equal to background rh at that level
!    (implemented as a constraint in subroutine start)
!  - parcels do not have any condensate (even if background rl non zero) 
!  - starting parcel km spectrum should be very close but not exactly equal to
!    initial spectrum kmstart: kmstart is modified here
!  - Desireable to assign air parcels starting at the same
!    level the same RH, since this is probably most realistic. Avoids
!    strange starting buoyancies, excessive entrainment, and strange high
!    LCL's for low RH parcels
!  - in previous versions was trying to get average km/rv of spectrum exactly
!    the same as background. But this is not needed for conservation, since
!    kmflow,rvflow are based on actual final solutions, as determined here.
!
!  Inputs:
!
!  istart: starting level of updraft parcel spectrum
!
!------------------------------------------------------------------------
subroutine get_start_spectrum
!------------------------------------------------------------------------
!
!   Use
!
!------------------------------------------------------------------------
use if_conv_solvers, only : hmoist,start,enthalpy,rsat
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: km_got,rh_got,t_got,rv_got
real(r8) :: tb,pb,rvb,rlb,rsb,rhb,km_diff,rh_diff,t_diff
integer :: nl,it
!------------------------------------------------------------------------
!
!  Loop over it
!
!------------------------------------------------------------------------
do it = 1,maxlev
!------------------------------------------------------------------------
!
!   Check rv > 0 (not really needed)
!
!------------------------------------------------------------------------
  if (rv(it) < 0.) then
    print*,'it z(it) = ',it,z(it)
    stop 'rv(it) negative at start of get_start_spectrum'
  endif
!------------------------------------------------------------------------
!
!  Define background quantities
!
!------------------------------------------------------------------------
  tb = t(it)
  pb = p(it)
  rvb = rv(it)
  rlb = rl(it)
  rsb = rsat(tb,pb)
  rhb = rvb/rsb
!------------------------------------------------------------------------
!
!  Loop over number of updraft parcels at a level
!  - determines km_got,t_got,rv_got
!  - km_got should be the same as kmstart(it,nl) to numerical accuracy
!  - rhb and rh_got should be identical
!  - "subroutine start" is in if_conv_solvers.f90
!
!------------------------------------------------------------------------
  do nl = 1,nlaunch
!   print*,'kmstart(it,nl) = ',kmstart(it,nl)
    call start(kmstart(it,nl),km_got,rhb,rlb,rh_got,tb,t_got,rv_got,pb)
!   print*,'back from start'
    km_diff = abs(kmstart(it,nl) - km_got)/km_got
    rh_diff = abs(rhb - rh_got)/rh_got
    t_diff = abs(tb - t_got)/t_got
!------------------------------------------------------------------------
!
!   km discrepancy
!   - a km_diff isn't a conservation problem so not much to worry about
!
!------------------------------------------------------------------------
    if (km_diff > 1.0E-05) then
     print*,'********  KM DIFFERENCE in subroutine start **********'
     print*,'km_diff = ',km_diff
     print*,'rh_diff = ',rh_diff
     print*,'t_diff = ',t_diff
     print*,'km_got = ',km_got
     print*,'rhb = ',rhb
     print*,'rh_got = ',rh_got
     print*,'tb = ',tb
     print*,'t_got = ',t_got
     print*,'rv_got = ',rv_got
     print*,'pb = ',pb
     stop 'km_diff problem'
    endif
!------------------------------------------------------------------------
!
!   rh discrepancy
!   - again, having an rh discrepancy is not a big problem. Doesn't
!     cause conservation problems.
!
!------------------------------------------------------------------------
    if (rh_diff > 1.0E-04) then
     print*,'********  RH DIFFERENCE in subroutine start **********'
     print*,'km_diff = ',km_diff
     print*,'rh_diff = ',rh_diff
     print*,'t_diff = ',t_diff
     print*,'km_got = ',km_got
     print*,'it kmstart(it,nl) = ',it,nl,kmstart(it,nl)
     print*,'rhb = ',rhb
     print*,'rh_got = ',rh_got
     print*,'tb = ',tb
     print*,'t_got = ',t_got
     print*,'rv_got = ',rv_got
     print*,'pb = ',pb
!    stop
    endif
!------------------------------------------------------------------------
!
!   Set "start" variables equal to "got" variables
!
!------------------------------------------------------------------------
    kmstart(it,nl) = km_got
    tstart(it,nl) = t_got
    rvstart(it,nl) = rv_got
!------------------------------------------------------------------------
!
!  end loop over nlaunch
!
!------------------------------------------------------------------------
  end do
!------------------------------------------------------------------------
!
!  End loop over it
!
!------------------------------------------------------------------------
end do
!------------------------------------------------------------------------
!
!  end
!
!------------------------------------------------------------------------
end subroutine get_start_spectrum
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Subroutine up
!
!   Inputs:
!   i: current model level; moving parcel to level i+1
!
!
!   Comments:
!   - the current value of hmp(i) must be defined
!   - the sum rvp(i)+rlp(i)+rip(i) must be up to date.
!   - hmp(i) and total water are assumed to be conserved;
!   - z(i) and z(i+1) must be up to date
!   - mp(i+1),tp(i+1),rvp(i+1),ri(i+1),rlp(i+1),rsp(i+1) are all updated
!
!    From input hmp(i) and z(i+1), define kmp at level i+1
!    use get_rv to solve for t/rv/rl/ri at level i+1
!
!------------------------------------------------------------------------
 subroutine up(i)
!------------------------------------------------------------------------
!
!   Use
!
!------------------------------------------------------------------------
use if_conv_solvers, only : get_rv,epsi,g,rsat,hmoist
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
integer, intent(in) :: i
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: rtp,kmp
real(r8) :: hm_check, hm_error
integer :: nums, allow_sat, nprint_upp
!------------------------------------------------------------------------
!
!  Define total parcel water and parcel enthalpy at model level i+1
!  hm = km + (1. + rt)*g*z
!  - looks a bit strange but basically assuming that hmp(i+1) = hmp(i),
!    and rtp conserved also.
!
!------------------------------------------------------------------------
  nprint_upp = 0
  if (i == nd) stop 'i too high in subroutine up of if_conv_tend.f90'
  rtp = rvp(i) + rlp(i) + rip(i)
  kmp = hmp(i) - (1. + rtp)*g*z(i+1)
  rip(i+1) = rip(i)
  tp(i+1) = tp(i)   ! Need initial estimate of t for input to get_rv
!------------------------------------------------------------------------
!
!  Check for negative kmp
!
!------------------------------------------------------------------------
  if (kmp < 0.) then
    print*,'=========== NEGATIVE kmp or ======================'
    print*,'i kmp = ',i,kmp
    print*,'hmp(i) = ',hmp(i)
    print*,'rtp = ',rtp
    print*,'rvp(i) = ',rvp(i)
    print*,'rlp(i) = ',rlp(i)
    print*,'rip(i) = ',rip(i)
    print*,'z(i+1) = ',z(i+1)
    stop 'stopping for negative kmp on entry to subroutine up'
  endif
!------------------------------------------------------------------------
!
!  Call get_rv to determine parcel properties at level i+1
!  rvp(i+1),rlp(i+1),rip(i+1) undefined on entry to get_rv
!  get_rv should exactly conserve rtp and kmp
!
!  Calling get_rv from subroutine up
!
!------------------------------------------------------------------------
! print*,'calling get_rv from subroutine up i kmp = ',i,kmp
  allow_sat = 0
  get_rv_call = 2
  call_get_rv_up = call_get_rv_up + 1
  if (nprint_upp == 1) print*,'Calling get_rv from subroutine up'
  call get_rv(kmp,  &  ! pure in
              rtp,  &  ! pure in
           p(i+1),  &  ! pure in
         rvp(i+1),  &  ! pure out
         rlp(i+1),  &  ! pure out
         rip(i+1),  &  ! pure in
          tp(i+1),  &  ! inout (need to estimate ice enthalpy)
      get_rv_call,  &  ! pure in
        allow_sat)
!------------------------------------------------------------------------
!
!  Check if tp in range
!
!------------------------------------------------------------------------
  if ((tp(i+1).ge.180.).or.(tp(i+1).le.320.)) then
  else
    print*,'tp(i+1) = ',tp(i+1)
    print*,'kmp = ',kmp
    print*,'rtp = ',kmp
    print*,'i p(i+1) = ',i,p(i+1)
    print*,'rlp(i+1) = ',i,rlp(i+1)
    print*,'rvp(i+1) = ',i,rvp(i+1)
    stop 'tp out of range after call to get_rv in subroutine up'
  endif
!------------------------------------------------------------------------
!
!   Update parcel mp, rsp, hmp
!   - this update only neccessary for actual move up
!
!------------------------------------------------------------------------
  if (i == nd) stop 'i will be out of range in subroutine up'
  mp(i+1) = mp(i)
  rsp(i+1) = rsat(tp(i+1),p(i+1))
  hmp(i+1) = hmp(i)
  tdp(i+1) = tp(i+1)*(1. + epsi*rvp(i+1))/(1. + rtp)
  tvp(i+1) = tp(i+1)*(1. + epsi*rvp(i+1))/(1. + rvp(i+1))
  uwindp(i+1) = uwindp(i)
  vwindp(i+1) = vwindp(i)
  if (mp(i) < 0.) stop 'mp(i) negative in subroutine up'

if (nprint_upp == 1) then
  print*,'At end of subroutine up current value hmp(i+1) = ',hmp(i+1)
  print*,'At end of subroutine up current value rip(i+1) = ',rip(i+1)
  print*,'At end of subroutine up current value rlp(i+1) = ',rlp(i+1)
  print*,'At end of subroutine up current value rvp(i+1) = ',rvp(i+1)
  print*,'At end of subroutine up current value tp(i+1) = ',tp(i+1)
endif

!------------------------------------------------------------------------
!
!  Check for hmp consistency
!
!------------------------------------------------------------------------
  rtp = rvp(i+1) + rlp(i+1) + rip(i+1)
  hm_check = hmoist(tp(i+1),rtp,rvp(i+1),rip(i+1),z(i+1)) 
  hm_error = abs(hm_check - hmp(i+1))/hmp(i+1)
  if (hm_error > 1.0E-06) then
    print*,'-------------------------------'
    print*,'hm inconsistency in subroutine up'
    print*,'hm_error = ',hm_error
    print*,'hm_check = ',hm_check
    print*,'hmp(i+1) = ',hmp(i+1)
    print*,'-------------------------------'
    stop 'hm_error problem'
  endif
!------------------------------------------------------------------------
!
!   END
!
!------------------------------------------------------------------------
end subroutine up
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  Subroutine entrain
!
!  August 2014: currently called only from updrafts
!
!  i  : level at which background air is being mixed with the parcel.
!  sigma : fractional entrainment of ambient air (sigma = dMd/dMd+Mdp,
!          where dMd is entrained dry mass and Mdp is the initial dry mass
!          of the parcel)
!
!------------------------------------------------------------------------
subroutine entrain(i,sigma,tstep)
!------------------------------------------------------------------------
!
!   Use
!
!------------------------------------------------------------------------
use if_conv_solvers, only : get_rv,hmoist,rsat,enthalpy,epsi,g
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
integer, intent(in) :: i
real(r8), intent(in) :: sigma
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: nums,allow_sat
real(r8) :: dmd,rtm,kmm,kp,kmb_eff,hmb_eff,rim
real(r8) :: rtb,rtp,rlb,rib,rvb,rhb,rsb,rspp
real(r8) :: hm_check, hm_error

!------------------------------------------------------------------------
!
!   Define background properties
!    - Note that do not entrain ri
!    - Feb 2015: entrain rl as a way of reducing rl at high rain rates, 
!      especially at the surface. Very likely a good thing.
!    - Note that must define an hmb_eff since background hm(i) includes rl,ri,
!      and may be entraining rv not equal to background.
!
!------------------------------------------------------------------------
  rsb = rsat(t(i),p(i))
  rspp = rsat(tp(i),p(i))
  rvb = rv(i) 
  if (rvb > rsb) rvb = rsb
  if (rvb > rspp) rvb = rspp
  rlb = rl(i)
  rtb = rvb + rlb
  rib = 0.
  hmb_eff = hmoist(t(i),rtb,rvb,rib,z(i))  ! careful: no identical arguments
  kmb_eff = hmb_eff - (1.0 + rtb)*g*z(i)

!------------------------------------------------------------------------
!
!  Determine parcel properties
!  Could have a check here that parcel and/or background variables
!    consistent with km def using function enthalpy
!  - rvp and rlp global variables
!  - ice not now included as parcel variable
!
!------------------------------------------------------------------------
  rtp = rvp(i) + rlp(i) + rip(i)
  kp = hmp(i) - (1.0000000+rtp)*g*z(i)

!------------------------------------------------------------------------
!
!   Determine total total enthalpy and total water of mixture (required
!      as inputs to get_rv)
!
!   kmb_eff,kp : ambient/parcel enthalpies
!   sigma is the ambient dry mixing fraction.
!   The expressions for km are based on per unit dry mass, so
!    the mixing fractions sigma must be based on dry mass ....
!   ri is purely diluted during mixing (note rib = 0.)  
!  
!------------------------------------------------------------------------
  kmm = sigma*kmb_eff + (1-sigma)*kp
  rtm = sigma*rtb + (1-sigma)*rtp
  rim = sigma*rib + (1-sigma)*rip(i)

!------------------------------------------------------------------------
!
!  Call get_rv to determine properties of mixture.
!
!  Inputs:
!  (i) kmm: total enthalpy of a mixture
!  (ii) rtm: total total water of a mixture
!  (iii) pb: background pressure
!  (iv) tg: guess temperature
!  (v) rim
!
!  Outputs:
!  (i) rvp(i)
!  (ii) rlp(i)
!  (iii) rip(i)
!  (iv) tp(i)
!
!  Calling get_rv from subroutine entrain
!
!------------------------------------------------------------------------
! print*,'kmm sigma kmb_eff kp = ',kmm,sigma,kmb_eff,kp
! print*,'i rtm rtb rtp = ',i,rtm,rtb,rtp
! print*,'p(i) rvp(i) rlp(i) rip(i) tp(i) = ',p(i),rvp(i),rlp(i),tp(i)
! print*,'calling get_rv from subroutine entrain kmm = ',kmm
  allow_sat = 0
  get_rv_call = 3
  call_get_rv_entrain = call_get_rv_entrain + 1
! print*,'Calling get_rv from entrain'
  call get_rv(kmm,  &
              rtm,  &
             p(i),  &
           rvp(i),  &
           rlp(i),  &
              rim,  &  ! pure in (fixed)
            tp(i),  &  ! inout: need estimate on input for ice enthalpy
      get_rv_call,  &  ! pure in
      allow_sat)
!  print*,'finished get_rv'

!------------------------------------------------------------------------
!
!  Set new parcel variables to mixed variables.
!  No need to define rtpi on exit.
!  April 21/2010: used kmm to re-define hmp, instead of hmoist. Should be
!    the same but avoid using moist unless absolutely neccessary.
!
!------------------------------------------------------------------------
  tdp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rtm)
  tvp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rvp(i))
  rsp(i) = rsat(tp(i),p(i))
  hmp(i) = kmm + (1.0+rtm)*g*z(i)
  rip(i) = rim

!------------------------------------------------------------------------
!
!  Check for hm consistency
!  - April 2015: A bit troubled I had to set hm_error to a relatively
!    high value. Not sure why is neccessary ... get_rv should exactly
!    conserved km, and just adding and subtracting gravitational part.
!
!------------------------------------------------------------------------
  hm_check = hmoist(tp(i),rtm,rvp(i),rip(i),z(i))
  hm_error = abs(hm_check - hmp(i))/hmp(i)
  if (hm_error > 1.0E-05) then
     print*,'--- hm inconsistency after call to get_rv from entrain ------'
     print*,'hm_error = ',hm_error
     print*,'hm_check = ',hm_check
     print*,'hmp(i) = ',hmp(i)
     print*,'z(i) = ',z(i)
     print*,'tp(i) = ',tp(i)
     print*,'rvp(i) = ',rvp(i)
     print*,'rlp(i) = ',rlp(i)
     print*,'rip(i) = ',rip(i)
     print*,'rtm = ',rtm
     print*,'rvp(i) + rlp(i) + rip(i) = ',rvp(i)+rlp(i)+rip(i)
     print*,'-----------------------------'
     stop 'hm_error after get_rv problem'
  endif

!------------------------------------------------------------------------
!
!   Update parcel mass:
!   dmd = sigma*mass/(1 - sigma)
!   new parcel mass: mass = mass + dmd
!   entrained vapor mass = rv(i)*dmd
!
!------------------------------------------------------------------------
  dmd = mp(i)*(sigma/(1.-sigma))
  mp(i) = mp(i) + dmd

  if (dmd < 0.) then
    print*,'----- dmd negative problem in subroutine entrain -----'
    print*,'dmd = ',dmd
    print*,'sigma = ',sigma
    stop 'dmd negative'
  endif
!------------------------------------------------------------------------
!
!   Modify flow variables (should be for updrafts only):
!
!   Note that this would have to modified for non-zero updraft ice
!
!------------------------------------------------------------------------
    dmflow(1,i) = dmflow(1,i) + (dmd/tstep)
    rvflow(1,i) = rvflow(1,i) + (dmd*rvb/tstep)
    rlflow(1,i) = rlflow(1,i) + (dmd*rlb/tstep) ! 
    riflow(1,i) = riflow(1,i) + (dmd*rib/tstep) ! But no entrain since rib = 0

!   if (i.eq.3) then
!   print*,'---------------------------------------------------'
!   print*,'in subroutine entrain for i = ',i
!   print*,'kg vapor removed from layer i by entrainment = ',dmd*rv(i)
!   print*,'(dmd*rv(i)/tstep) = ',(dmd*rv(i)/tstep)
!   print*,'dmd = ',dmd
!   print*,'sigma = ',sigma
!   print*,'mp(i) = ',mp(i)
!   print*,'---------------------------------------------------'
!   endif
!   endif
    hmflow(1,i) = hmflow(1,i) + (dmd*hmb_eff/tstep)
    uflow(1,i) = uflow(1,i) + (dmd*uwind(i)/tstep)
    vflow(1,i) = vflow(1,i) + (dmd*vwind(i)/tstep)

!------------------------------------------------------------------------
!
!   Update parcel momentum
!
!------------------------------------------------------------------------
    uwindp(i) = sigma*uwind(i) + (1. - sigma)*uwindp(i)
    vwindp(i) = sigma*vwind(i) + (1. - sigma)*vwindp(i)

!------------------------------------------------------------------------
!
!   End
!
!------------------------------------------------------------------------
end subroutine entrain
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  Subroutine updrafts
!
!------------------------------------------------------------------------
subroutine updrafts(land_fraction,tstep)
!------------------------------------------------------------------------
!
!   Use
!
!------------------------------------------------------------------------
use if_conv_solvers, only : epsi,g,hmoist,rsat,hmerr_max,enthalpy,latent,cl,get_t
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: land_fraction
real(r8), intent(in) :: tstep

!------------------------------------------------------------------------
!
!  Local Variables
!
!------------------------------------------------------------------------
integer :: i,full_detrain,nl,it,iter
integer :: jj,ipp,nm,nprint_up,from_detrain
real(r8) :: xxx,rtp,rtpi,rtpf,dpp,xx,pp,dzz,tot_f_ent
real(r8) :: hmdiag(nd,nlaunch),massdiag(nd,nlaunch)
real(r8) :: xmass,kmerror,kmwant,kkk
real(r8) :: dmdry,rtt,ripp,mass_ratio
real(r8) :: totmass, sigma, bp, drv_det
real(r8) :: tke_remain,b_av,ff_det,fff
real(r8) :: tot_cond, rlp_old, rip_old, tp_old, kmp, mp_old, rvp_old
real(r8) :: t_high, t_low, hm_error, hmpp, rtp_old, km_test, km_error
real(r8) :: rv_old, rl_old, ri_old, dmdry_it

!------------------------------------------------------------------------
!  zzz
!
!   **********  BEFORE ALL LOOPS   *******************
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!   Loop over launch levels
!
!------------------------------------------------------------------------
  do it = 1,maxlev
! print*,'it = ',it
!------------------------------------------------------------------------
!
!   Loop over numbers of parcels launched
!
!------------------------------------------------------------------------
  do nl = 1,nlaunch
! print*,'nl = ',nl


!------------------------------------------------------------------------
!
!  This is a switch for printing out updraft quantities for diagnostics.
!  It generates a huge amount of output, so only turn on if running the
!  model for ONE timestep only!
!
!------------------------------------------------------------------------
nprint_up = 0

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!------------------------------------------------------------------------
  if (it == 1) cfraction_1 = cfraction_1 + conv_fraction(it,nl)
  if (it == 2) cfraction_2 = cfraction_2 + conv_fraction(it,nl)
  if (it == 3) cfraction_3 = cfraction_3 + conv_fraction(it,nl)
  if (it == 4) cfraction_4 = cfraction_4 + conv_fraction(it,nl)

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!  The criterion for instability is defined in a previous subroutine.
!  - parcel in updraft spectrum exceeds capp_thresh, etc
!
!------------------------------------------------------------------------
  if (nunstable(it,nl) == 1) then

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!   call up_init.
!   - initialization to zero of various quantities
!
!------------------------------------------------------------------------
! print*,'called up_init'
  call up_init

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!  Define mp(it) 
!  - mp should have units kg/m2
!
!------------------------------------------------------------------------
  dmdry_it = (phalf(it)-phalf(it+1))/(g*(1.+rv(it)+rl(it)+ri(it)))
  mp(it) = conv_fraction(it,nl)*dmdry_it

  if (mp(it) < 0.) then
     print*,'conv_fraction(it,nl) = ',conv_fraction(it,nl)
     print*,'dmdry_it = ',dmdry_it
     print*,'phalf(it) = ',phalf(it)
     print*,'phalf(it+1) = ',phalf(it+1)
     print*,'it = ',it
     stop 'mp(it) initialized negative'
  endif

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!  Initialize updraft parcel properties at starting level it
!    - put this into subroutine?
!    - initialization identical for every (it)
!    - re-initialize starting parcel for every pass
!    - momentum initialization: assume scales as fraction of layer mass
!
!  Feb 2015: important to include rl entrainment as a sink of rl near
!     the surface to avoid "fog" at high rain rates. Subsidence pushes
!     rl to the surface so need a way to "fight back", and evaporation
!     is not very efficient at high RH near the surface.
!  March 2015: it is desireable to entrain rl, but seems to contribute to
!     negative values at high rain rates.
!
!------------------------------------------------------------------------
  rvp(it) = rvstart(it,nl)
  tp(it) = tstart(it,nl)
  if (rlp_start_entrain == 1) then
    rlp(it) = rl(it)
  else
    rlp(it) = 0.
  endif
  tdp(it) = tp(it)*(1.+epsi*rvp(it))/(1.+rvp(it)+rlp(it)+rip(it))
  tvp(it) = tp(it)*(1.+epsi*rvp(it))/(1.+rvp(it))
  rsp(it) = rsat(tp(it),p(it))
  rip(it) = 0.
  hmp(it) = kmstart(it,nl) + (1. + rvp(it) + rlp(it) + rip(it))*g*z(it)
  uwindp(it) = uwind(it)
  vwindp(it) = vwind(it)

! print*,'it kmstart(it,nl) = ',it,kmstart(it,nl)
!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!   Diagnostics
!
!------------------------------------------------------------------------

if (nprint_up == 1) then
  print*,'------------ INITIAL CONDITIONS  --------'
  print*,'Starting it = ',it
  print*,'Initial rvp(it) = ',rvp(it)
  print*,'Initial tp(it) = ',tp(it)
  print*,'Initial tdp(it) = ',tdp(it)
  print*,'Initial hmp(it) = ',hmp(it)
  print*,'Background hm(it) = ',hm(it)
  print*,'Initial capp(it,nl) = ',capp(it,nl)
  print*,'Background t(it) = ',t(it)
  print*,'Background rl(it) = ',rl(it)
  print*,'Background td(it) = ',td(it)
  print*,'Parcel temperature offset = ',tp(it) - t(it)
  print*,'Parcel rv offset = ',rvp(it) - rv(it)
endif

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!  Set "flow" arrays
!  - in subroutine "updrafts"
!
!------------------------------------------------------------------------
  dmflow(1,it) = dmflow(1,it) + (mp(it)/tstep)
  rvflow(1,it) = rvflow(1,it) + (mp(it)*rvp(it)/tstep)
  rlflow(1,it) = rlflow(1,it) + (mp(it)*rlp(it)/tstep)
  riflow(1,it) = riflow(1,it) + (mp(it)*rip(it)/tstep)

! print*,'it starting parcel water = ',(rvp(it)+rlp(it))*mp(it)
! if (it.eq.3) then
!   print*,'------------------------------------------------'
!   print*,'mp(it) = ',mp(it)
!   print*,'rv(it) rvp(it) = ',rv(it),rvp(it)
!   print*,'kg of vapor removed for this mp(it)*rvp(it) = ',mp(it)*rvp(it)
!   print*,'------------------------------------------------'
! endif

  hmflow(1,it) = hmflow(1,it) + (mp(it)*hmp(it)/tstep)
  uflow(1,it) = uflow(1,it) + (mp(it)*uwindp(it)/tstep)
  vflow(1,it) = vflow(1,it) + (mp(it)*vwindp(it)/tstep)

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!  Diagnostics for calculating mean updraft starting MSE
!  - these variables should be mass weighted of all updraft parcels
!    from all initial levels.
!  - Make sure initialized at the before this routine is entered, or
!    here before all loops.
!  - units of mass_start and mp(i): kg/m2
!
!------------------------------------------------------------------------
  if (nl == 1) mass_start_1 = mass_start_1 + mp(it)
  if (nl == 2) mass_start_2 = mass_start_2 + mp(it)

  hm_start = hmp(it)

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!  Check for km inconsistency in starting parcel spectrum
!  - likely redundant but keep for now
!
!------------------------------------------------------------------------
  call enthalpy(kmwant,tp(it),rvp(it)+rlp(it)+rip(it),rvp(it),rip(it))
  kmerror = abs(kmwant-kmstart(it,nl))/kmwant
  if (kmerror > hmerr_max) then
    print*,'========== Problem in START PARCEL Spectrum ====='
    print*,'kmwant = ',kmwant
    print*,'kmstart(it,nl) = ',kmstart(it,nl)
    print*,'kmerror = ',kmerror
    print*,'tp(it) = ',tp(it)
    print*,'rvp(it) = ',rvp(it)
    print*,'rlp(it) = ',rlp(it)
    print*,'rip(it) = ',rip(it)
    print*,'rsp(it) = ',rsp(it)
    print*,'z(it) = ',z(it)
    print*,'it = ',it
    stop 'km inconsistency at start of updrafts'
 endif

!------------------------------------------------------------------------
!
!   **********  0. INITIALIZATION   *******************
!
!  Initialize tke_remain with appropriate value for cloud mode.
!
!------------------------------------------------------------------------
  tke_remain = tke_amp*tke_start(nl)

!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
  full_detrain = 0
  from_detrain = 0

!------------------------------------------------------------------------
!
!  **********************************************************************
!  ****************   FINISHED INITIALIZATION  **************************
!  **********************************************************************
!
!  Start loop over parcel heights i
!
!------------------------------------------------------------------------

  do i = it,nd
! print*,'i = ',i

!------------------------------------------------------------------------
!
!  **********  1. INITIAL CHECKS/DIAGNOSTICS  ***************************
!
!  Initialize some stuff to zero (so don't get run time complaints)
!  - all these different buoyancies for diagnostics
!
!------------------------------------------------------------------------
  mass_ratio = 0.

  bp_start(i) = 0.
  bp_precip(i) = 0.
  bp_iter(i) = 0.
  bp_freeze(i) = 0.

!------------------------------------------------------------------------
!
!   If for updraft parcel not yet detrained
!
!------------------------------------------------------------------------
  if (full_detrain == 0) then

!------------------------------------------------------------------------
!
!  **********  INITIAL CHECKS/DIAGNOSTICS  ***************************
!
!  Check to see if updraft parcel has reached top level
!
!------------------------------------------------------------------------
  if (i == nd) then 
    print*,'------------------------------------------------'
    print*,'i it z(i) = ',i,it,z(i)
    print*,'tp(i-1) tp(i) = ',tp(i-1),tp(i)
    stop 'strange situation: updraft parcel has reached top model level'
  endif

!------------------------------------------------------------------------
!
!  **********  INITIAL CHECKS/DIAGNOSTICS  ***************************
!
!   Check for negative parcel mass
!
!------------------------------------------------------------------------
  if (mp(i) < 0.) then
    print*,'i mp(i) = ',i,mp(i)
    stop 'mp(i) negative at start of updraft loop'
  endif

!------------------------------------------------------------------------
!
!  **********  INITIAL CHECKS/DIAGNOSTICS  ***************************
!
!   Define buoyancy at start of level
!   Make sure td(i) and tdp(i) are defined.
!   - tdp(i) should be defined from subroutine up, or from inititialization
!     at start level
!
!------------------------------------------------------------------------
  rtp = rvp(i) + rlp(i) + rip(i)
  tdp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rtp)
  tvp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rvp(i))
  bp_start(i) = g*(tdp(i) - td(i))/td(i)

! print*,'z(i) bp_start(i) = ',0.001*z(i),bp_start(i)

  if (nl == 1) then
    bp1_str(i) = bp1_str(i) + bp_start(i)
    num_bp1_str(i) = num_bp1_str(i) + 1.
  elseif (nl == 2) then
    bp2_str(i) = bp2_str(i) + bp_start(i)
    num_bp2_str(i) = num_bp2_str(i) + 1.
  endif

!------------------------------------------------------------------------
!
!  **********  INITIAL CHECKS/DIAGNOSTICS  ***************************
!
!------------------------------------------------------------------------
if (nprint_up == 1) then
  print*,'---------  AT start of model level ----------'
  print*,'i it = ',i,it
  print*,'starting buoyancy before precip = ',bp_start(i)
  print*,'tke_remain = ',tke_remain
  print*,'Initial parcel mass = ',mp(i)
  print*,'Parcel hmp(i) = ',hmp(i)
  print*,'Background hm(i) = ',hm(i)
  print*,'hmdiff = ',hmp(i)-hm(i)
  print*,'Background Temperature t(i) = ',t(i)
  print*,'Parcel Temperature tp(i) = ',tp(i)
  print*,'Background Density Temperature td(i) = ',td(i)
  print*,'Parcel Density Temperature tdp(i) = ',tdp(i)
  print*,'Parcel rvp(i) = ',rvp(i)
  print*,'Parcel rlp(i) = ',rlp(i)
  print*,'Parcel rip(i) = ',rip(i)
  print*,'Parcel condensate before precipitate rlp(i) = ',rlp(i) 
  print*,'-------'
endif

!---------------------------------------------------------------------
!
!  ********************  Update tke_remain   *********************
!
!---------------------------------------------------------------------

if (i > it) then
  b_av = 0.5*(bp_start(i) + bp_start(i-1))
  if (b_av < 0.) tke_remain = tke_remain + b_av*(zh(i)-zh(i-1))
endif
if ( (p(it) - p(i)) > dp_tke_start ) tke_remain = 0.

if (nprint_up == 1) then
print*,'-----  TKE -----------'
print*,'bp_start(i) = ',bp_start(i)
print*,'bp_start(i-1) = ',bp_start(i-1)
print*,'b_av = ',b_av
print*,'tke_remain = ',tke_remain
print*,'Starting tke = ',tke_amp*tke_start
print*,'-------------------------'
endif

!------------------------------------------------------------------------
!
!  ***********  2. Check for FULL DETRAIN   ****************************
!
!  Detrainment check
!  - Initialize full_detrain = 0. Then look for reasons to detrain.
!
!------------------------------------------------------------------------
full_detrain = 0

if (i == nd) full_detrain = 1


if ((i > it).and.(bp_start(i) < b_det(nl)).and.(tke_remain < 0.0001)) full_detrain = 1

if (nprint_up == 1) then
  print*,'-----------'
  print*,'Outcome of full_detrain calculation full_detrain = ',full_detrain
  print*,'i it = ',i,it
  print*,'bp_start(i) = ',bp_start(i)
  print*,'tke_remain = ',tke_remain
  print*,'-----------'
endif

!------------------------------------------------------------------------
!
!   Diagnostics
!
!   Detrainment test for starting level
!   - no way to test how much tke used up in going to higher level.
!   - currently always go up one level
!
!------------------------------------------------------------------------

if (nprint_up == 1) then
if ((full_detrain == 1).and.(z(i) < 3000.)) then
  print*,'DEEP detrain at low altitude'
  print*,'i z(i) = ',i,z(i)
  print*,'it = ',it
  print*,'cape_1 = ',cape_1
  print*,'cape_2 = ',cape_2
  print*,'cape_3 = ',cape_3
  print*,'cape_4 = ',cape_4
  print*,'cfraction_1 = ',cfraction_1
  print*,'cfraction_2 = ',cfraction_2
  print*,'cfraction_3 = ',cfraction_3
  print*,'cfraction_4 = ',cfraction_4
  do jj = 1,i
    print*,'jj z(jj) = ',jj,z(jj)
    print*,'rh(jj) = ',rh(jj)
    print*,'z(jj) = ',0.001*z(jj)
    print*,'rlp(jj) = ',rlp(jj)
  end do
!  stop 'DEEP detrain too low'
  print*,'--------  KEEP GOING ---------'
endif
endif

!------------------------------------------------------------------------
!
!  ***********  2. FULL DETRAIN   ****************************
!
!  - detrain returns with mp(i) = 0 if full_detrain = 1
!
!------------------------------------------------------------------------
if (full_detrain == 1) then

  ipp = i   !  dummy variable for parcel index
  cape_detrain = capp(it,nl)
  if (nprint_up == 1) print*,'Calling detrain for full detrain'
  call detrain(i,ipp,it,full_detrain,bad,nl,tstep)   ! updraft final detrain

    if (i == it) then
!     print*,'Full detrain '
!     print*,'Detraining at i it = ',i,it
!     stop 'should not detrain at starting level'
    endif

!------------------------------------------------------------------------
!  Diagnostics
!------------------------------------------------------------------------
if (nprint_up == 1) then
  print*,'---------- DID FULL DETRAIN ------------'
  print*,'hmp(i) - hm(i) = ',hmp(i) - hm(i)
endif

!------------------------------------------------------------------------
!
!  else for full_detrain = 0
!
!  entrain
!  - sigma should be positive
!
!------------------------------------------------------------------------
else

if (nprint_up == 1) print*,'Starting entrainment iteration'
bp = bp_start(i)
tot_f_ent = 0.

if (i > 1) then
  dzz = 0.5*(z(i+1)-z(i-1))*0.001
  if (i == nd) stop 'i too high in if_conv_tend.f90'
else
  dzz = (z(i+1)-z(i))*0.001
endif

!------------------------------------------------------------------------
!
!  ***************  LOOP  ***********************************
!
!  Start entrainment loop
!  tot_f_ent : total entrainment so far per km.
!
!------------------------------------------------------------------------
do iter = 1,iter_max

  if (nprint_up == 1) then 
    print*,'testing for entrainment iteration = ',iter
    print*,'bp = ',bp
    print*,'b_tar(nl) = ',b_tar(nl)
    print*,'tot_f_ent = ',tot_f_ent
    print*,'f_ent_max(nl) = ',f_ent_max(nl)
  endif

if (bp > b_tar(nl)) then
if (tot_f_ent < f_ent_max(nl)) then

!------------------------------------------------------------------------
!
!  Assign entrainment
!
!------------------------------------------------------------------------

  sigma = ent_inc_km*dzz

  if (nprint_up == 1) print*,'sigma = ',sigma

!------------------------------------------------------------------------
!
!  Parcel mass must not exceed some multiple of starting mass
!
!------------------------------------------------------------------------
   mass_ratio = mp(i)/(dmdry_it*conv_fraction(it,nl))

   if (mass_ratio > mass_ratio_max) then
!    print*,'*******  mass_ratio larger than mass_ratio_max = ',mass_ratio_max
!    print*,'DANGER: updraft parcel has grown much larger than initial mass'
!    print*,'AND WORSE: trying to get even bigger: prevent this nonsense'
     sigma = 0.
   endif

   tot_f_ent = tot_f_ent + (sigma/dzz)   ! Entrainment per km

   if (nprint_up == 1) print*,'updated tot_f_ent = ',tot_f_ent

!------------------------------------------------------------------------
!
!  ***************  ENTRAINMENT  *****************************
!
!  entrain
!
!------------------------------------------------------------------------

   if (nprint_up == 1) print*,'Calling entrain for sigma = ',sigma

   call entrain(i,sigma,tstep)

!------------------------------------------------------------------------
!
!  ***************  ENTRAINMENT  *****************************
!
!  recalculate buoyancy: bp
!
!------------------------------------------------------------------------

rtp = rvp(i) + rlp(i) + rip(i)
tdp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rtp)
bp = g*(tdp(i) - td(i))/td(i)

!print*,'iter bp b_tar(nl) = ',iter,bp,b_tar(nl)

if (nprint_up == 1) print*,'Buoyancy after entrain = ',bp

!------------------------------------------------------------------------
!
!  endif for tot_f_ent < f_ent_max
!  endif for bp > b_tar
!  enddo for iter = 1, iter_max
!
!------------------------------------------------------------------------
endif
endif
end do

!------------------------------------------------------------------------
!
!  ******   Detrainment  **********
!
!------------------------------------------------------------------------

fff = 0.

drv_det = rh_f_det - rh(i)
fff = f_det_rh(nl)*drv_det

if (use_drh_low == 1) then
  xxx = 1.-rh(i)
  if (xxx < 0.) xxx = 0.
  fff = sigmoidal(xxx,drh_half,drh_scale,drh_low,drh_add)
  if (fff < 0.) fff = 0.
  fff = fff*f_det_rh(nl)
endif

!print*,'xxx = ',xxx
!print*,'fff = ',fff

if ((drv_det > 0.).and.(i > it)) then
  ff_det = dzz*fff
  if (ff_det > 0.95) ff_det = 0.95
  full_detrain = 0.
  ipp = i   !  dummy variable for parcel index
  call detrain(i,ipp,it,full_detrain,-ff_det,nl,tstep)
  if (nprint_up == 1) then
    print*,'--------- in drv_det > 0 -------------'
    print*,'drv_det = ',drv_det
    print*,'rh(i) = ',rh(i)
    print*,'rv(i) = ',rv(i)
    print*,'fff = ',fff
    print*,'drv_det = ',drv_det
    print*,'dzz = ',dzz
    print*,'fff = ',fff
    print*,'ff_det = ',ff_det
    print*,'------------'
  endif
endif

if (nl == 1) then
  f_det_1(i) = f_det_1(i) + (fff/dzz)
  num_f_det_1(i) = num_f_det_1(i) + 1.
else
  f_det_2(i) = f_det_2(i) + (fff/dzz)
  num_f_det_2(i) = num_f_det_2(i) + 1.
endif

!print*,'i  it fff = ',i,it,fff
!print*,'(fff/dzz) = ',(fff/dzz)
!print*,'dzz = ',dzz

!------------------------------------------------------------------------
!
!  ******   Entrainment  **********
!
!  Buoyancy after final entrainment
!
!------------------------------------------------------------------------
rtp = rvp(i) + rlp(i) + rip(i)
tdp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rtp)
tvp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rvp(i))
bp_iter(i) = g*(tdp(i) - td(i))/td(i)

if (nprint_up == 1) then 
  print*,'Finished entrainment loop'
  print*,'z(i) bp_iter(i) = ',0.001*z(i),bp_iter(i)
  print*,'Total fractional entrainment tot_f_ent = ',tot_f_ent
endif

if (nprint_up == 1) then
if (bp_iter(i) < 0.) then
  print*,'------- bp_iter negative  ----'
  print*,'Buoyancy at end of iterations = ',bp_iter(i)
  print*,'same as bp_start(i) ? = ',bp_start(i)
  print*,'i it = ',i,it
  print*,'tke_remain = ',tke_remain
  print*,'-----------------'
endif
endif

!------------------------------------------------------------------------
!
!  ******   Entrainment Diagnostics  **********
!
!------------------------------------------------------------------------

   if (nl == 1) then
     f_ent_1(i) = f_ent_1(i) + tot_f_ent
     num_f_ent_1(i) = num_f_ent_1(i) + 1.
     bp1_itr(i) = bp1_itr(i) + bp_iter(i)
     num_bp1_itr(i) = num_bp1_itr(i) + 1.
   else
     f_ent_2(i) = f_ent_2(i) + tot_f_ent
     num_f_ent_2(i) = num_f_ent_2(i) + 1.
     bp2_itr(i) = bp2_itr(i) + bp_iter(i)
     num_bp2_itr(i) = num_bp2_itr(i) + 1.
   endif

!------------------------------------------------------------------------
!
!  ***************  FREEZING   *********************************
!
!  Where to put freezing?
!  - probably not that important
!  - If occurred at the start of a level, would increase local buoyancy
!    and increase entrainment during backward mixing.
!  - If occurred after backward mixing, would increase local buoyancy
!    and decrease entrainment during forward mixing.
!  - If occurred after forward mixing, would increase buoyancy and prevent 
!    detrainment
!  - To increase outflow height, likely best to put just before final 
!    detrainment
!
!------------------------------------------------------------------------
tot_cond = rlp(i) + rip(i)
rtp_old = rlp(i) + rvp(i) + rip(i)
rlp_old = rlp(i)
rip_old = rip(i)
rvp_old = rvp(i)
tp_old = tp(i)
mp_old = mp(i)

f_ice = 0.

if (tp(i) > t_ice_high) then
  rlp(i) = tot_cond
  rip(i) = 0.
  f_ice = 0.
elseif (tp(i) < t_ice_low) then
  rlp(i) = 0.
  rip(i) = tot_cond
  f_ice = 1.
else
  f_ice = (t_ice_high - tp(i))/(t_ice_high-t_ice_low)
  rlp(i) = (1.-f_ice)*tot_cond
  rip(i) = f_ice*tot_cond
  if (f_ice < 0.) stop 'f_ice negative'
  if (f_ice > 1.) stop 'f_ice larger than 1'
endif

!------------------------------------------------------------------------
!
!  ***************  FREEZING   *********************************
!
!  Update tp after (possibly) modified condensate partitioning
!
!------------------------------------------------------------------------
! print*,'AFTER FREEZING Updating parcel properties: calling get_t to get tp(i)'
! print*,'hmp(i) = ',hmp(i)

  rv_old = rvp(i)
  rl_old = rlp(i)
  ri_old = rip(i)
  

  rtp = rvp(i) + rlp(i) + rip(i)
  kmp = hmp(i) - (1. + rtp)*g*z(i)
  tp(i) = get_t(kmp,rtp,rvp(i),rip(i))
  tdp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rtp)
  tvp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rvp(i))

  call enthalpy(km_test,tp(i),rvp(i)+rlp(i)+rip(i),rvp(i),rip(i))
  km_error = abs(km_test - kmp)/kmp

  hmpp = hmoist(tp(i),rtp,rvp(i),rip(i),z(i))
  hm_error = abs(hmpp - hmp(i))/hmp(i)

  if ((km_error > 1.0E-06).or.(hm_error > 1.0E-06)) then
    print*,'==========  KM or HM INCONSISTENCY AFTER FREEZING  ========='
    print*,'from parcel properties hmpp  = ',hmpp
    print*,'actual hmp(i) = ',hmp(i)
    print*,'hm_error = ',hm_error
    print*,'kmp = ',kmp
    print*,'km_test = ',km_test
    print*,'km_error = ',km_error
    print*,'kmp - km_test = ',kmp - km_test
    print*,'z(i) = ',z(i)
    print*,'g = ',g
    print*,'tp(i) = ',tp(i)
    print*,'rtp = ',rtp
    print*,'rtp_old = ',rtp_old
    print*,'rtp_old - rtp = ',rtp_old-rtp
    print*,'another rtp diff = ',rtp-rvp(i)-rlp(i)-rip(i)
    print*,'rvp(i) = ',rvp(i)
    print*,'rlp(i) = ',rlp(i)
    print*,'rip(i) = ',rip(i)
    print*,'rvp(i) - rvp_old = ',rvp(i)-rv_old
    print*,'rvp(i) - rlp_old = ',rlp(i)-rl_old
    print*,'rip(i) - rip_old = ',rip(i)-ri_old
    print*,'tp(i) - tp_old = ',tp(i)-tp_old
    print*,'f_ice = ',f_ice
    print*,'======================'
    stop 'AFTER FREEZING '
  endif

if (nprint_up == 1) then
  print*,'------------------ i = ',i
  print*,'Temperature change after freezing = ',tp(i)-tp_old
  print*,'Initial temperature = ',tp_old
  print*,'Final temperature = ',tp(i)
  print*,'Initial rtp rtp_old = ',rtp_old
  print*,'current rtp = ',rtp
  print*,'Initial rvp rvp_old = ',rvp_old
  print*,'current rvp(i) = ',rvp(i)
  print*,'Initial rip rip_old = ',rip_old
  print*,'New rip rip(i) = ',rip(i)
  print*,'Initial rlp = ',rlp_old
  print*,'New rlp rlp(i) = ',rlp(i)
  print*,'kmp = ',kmp
  print*,'-----------'
endif

!------------------------------------------------------------------------
!
!  Buoyancy after feezing
!
!------------------------------------------------------------------------
  bp_freeze(i) = g*(tdp(i) - td(i))/td(i)


!------------------------------------------------------------------------
!
!  ***************  1. PRECIPITATION  ***********************************
!
!   - Condensate removal/precip formation
!   - updates parcel density temperature
!   - does not update rsp
!   - all parcel quantities at parcel index i must be defined
!   - not neccessary to have parcel variables as arguments since 
!     precipitate in module, but is convenient
!   - diffcon is just for diagnostics
!   - only do if not starting level
!
!   Calling precipitate from updrafts
!
!------------------------------------------------------------------------
  rsp(i) = rsat(tp(i),p(i))
  rtpi = rvp(i) + rlp(i) + rip(i)
! print*,'parcel water mass before precipitate = ',(rvp(i)+rlp(i))*mp(i)
  full_detrain = 0
  from_detrain = 0

! print*,'Calling precipitate from updrafts'
  if (mp(i) > 1.0E-08) then
    call precipitate(tp(i),      &
                     z(i),       &
                        i,       &
                   rvp(i),       &
                   rlp(i),       &
                   rip(i),       &
                   hmp(i),       &
                   mp(i),        &
                   tdp(i),       &
                   from_detrain, &
                   full_detrain, & 
                   nl,tstep)
   endif

!------------------------------------------------------------------------
!
!  ***************  2. PRECIPITATE  ***********************************
!
!  Diagnostics
!  - Check for unphysical increase in parcel water
!
!------------------------------------------------------------------------

! print*,'parcel water mass after precipitate = ',(rvp(i)+rlp(i))*mp(i)
  rtpf = rvp(i) + rlp(i) + rip(i)
  diffcon(i) = rtpi - rtpf
  if (diffcon(i) < -1.0E-10) then 
      print*,'------------ PROBLEM AFTER PRECIPITATE -------'
      print*,'i diffcon(i) = ',i,diffcon(i)
      print*,'rtpi  = ',rtpi
      print*,'rtpf  = ',rtpf
      print*,'rvp(i) = ',rvp(i)
      print*,'rlp(i) = ',rlp(i)
      print*,'rip(i) = ',rip(i)
      stop 'diffcon(i) negative'
  endif

!------------------------------------------------------------------------
!
!  ***************  1. PRECIPITATION  *****************************
!
!  recalculate buoyancy: bp_precip
!
!------------------------------------------------------------------------

rtp = rvp(i) + rlp(i) + rip(i)
tdp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rtp)
tvp(i) = tp(i)*(1. + epsi*rvp(i))/(1. + rvp(i))
bp_precip(i) = g*(tdp(i) - td(i))/td(i)

!print*,'z(i) bp_precip(i) = ',0.001*z(i),bp_precip(i)

 if (nl == 1) then
   bp1_prp(i) = bp1_prp(i) + bp_precip(i)
   num_bp1_prp(i) = num_bp1_prp(i) + 1.
 else
   bp2_prp(i) = bp2_prp(i) + bp_precip(i)
   num_bp2_prp(i) = num_bp2_prp(i) + 1.
 endif

!------------------------------------------------------------------------
!
!  ***********  6. MOVE UP   ****************************
!
!  Move up
!
!------------------------------------------------------------------------
call up(i)

if (nprint_up == 1) print*,'MOVED UP'

!------------------------------------------------------------------------
!  Diagnostics
!------------------------------------------------------------------------
if (nprint_up == 1) then
  print*,'----  FINISHED UPDRAFT LEVEL  --'
endif

!------------------------------------------------------------------------
!
!  **********   END OF CALCULATIONS FOR A LEVEL  ************************
!
!  For diagnostics:
!  - by using mp_old here should still pick up the values AT FULL DETRAINMENT
!    mp(i) comes back zero after full detrainment
!  - variables that scale with the updraft parcel mass mp(i)
!  - average arrays should be initialized to zero in up_init
!
!------------------------------------------------------------------------
  av_rlp(i) = av_rlp(i) + mp_old*rlp(i)
  av_rip(i) = av_rip(i) + mp_old*rip(i)

  hm_diff_calc(i) = hm_diff_calc(i) + mp_old*(hmp(i) - hm_start)
  hmp_mean_calc(i) = hmp_mean_calc(i) + mp_old*hmp(i)
  av_mp(i) = av_mp(i) + mp_old

!------------------------------------------------------------------------
!
!   inner endif for full_detrain = 0  (i.e. parcel not yet detrained)
!
!------------------------------------------------------------------------
endif

!------------------------------------------------------------------------
!
!   outer endif for full_detrain = 0  (i.e. parcel not yet detrained)
!
!------------------------------------------------------------------------
endif

!------------------------------------------------------------------------
!
!   End loop over updraft levels for a particular (it) updraft parcel.
!   (i.e. starting level and launch index)
!   i = 1,nd
!
!------------------------------------------------------------------------
end do

!------------------------------------------------------------------------
!
!   endif for launch of (it,nl) parcel with nunstable(it,nl) = 1
!
!------------------------------------------------------------------------
  endif

!------------------------------------------------------------------------
!
!  End loop over index nl of nl = 1,nlaunch (parcels).
!
!------------------------------------------------------------------------
  end do
  close(25)

!------------------------------------------------------------------------
!
!  End loop over starting launch levels it = 1,maxlev
!
!------------------------------------------------------------------------
  end do

!------------------------------------------------------------------------
!
!   ****************** END ALL UPDRAFT LOOPS *********************
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  Normalize mass weighted updraft quantities
!  - do not set av_mp to a bad number: keep value of small
!  - also do not set updarft variables to bad. They are weighted by
!    av_mp in h1.f90.
!
!------------------------------------------------------------------------
  do i = 1,nd
    if (av_mp(i) > 1.0E-10) then
      hmp_mean(i) = hmp_mean_calc(i)/av_mp(i) 
      hm_diff(i) = hm_diff_calc(i)/av_mp(i) 
 
      av_rlp(i) = av_rlp(i)/av_mp(i) 
      av_rip(i) = av_rip(i)/av_mp(i) 
    else
      hmp_mean(i) = bad
      hm_diff(i) = bad
!
      av_rlp(i) = 0.
      av_rip(i) = 0.
    endif

   xxx = num_f_ent_1(i)
   if (xxx > 0.5) then 
     f_ent_1(i) = f_ent_1(i)/xxx
   else
     f_ent_1(i) = bad
   endif

   xxx = num_f_ent_2(i)
   if (xxx > 0.5) then 
     f_ent_2(i) = f_ent_2(i)/xxx
   else
     f_ent_2(i) = bad
   endif

   xxx = num_f_det_1(i)
   if (xxx > 0.5) then 
     f_det_1(i) = f_det_1(i)/xxx
   else
     f_det_1(i) = bad
   endif

   xxx = num_f_det_2(i)
   if (xxx > 0.5) then 
     f_det_2(i) = f_det_2(i)/xxx
   else
     f_det_2(i) = bad
   endif

   xxx = num_bp1_str(i)
   if (xxx > 0.5) then 
     bp1_str(i) = bp1_str(i)/xxx
   else
     bp1_str(i) = bad
   endif

   xxx = num_bp1_itr(i)
   if (xxx > 0.5) then 
     bp1_itr(i) = bp1_itr(i)/xxx
   else
     bp1_itr(i) = bad
   endif

   xxx = num_bp1_prp(i)
   if (xxx > 0.5) then 
     bp1_prp(i) = bp1_prp(i)/xxx
   else
     bp1_prp(i) = bad
   endif

   xxx = num_bp2_str(i)
   if (xxx > 0.5) then 
     bp2_str(i) = bp2_str(i)/xxx
   else
     bp2_str(i) = bad
   endif

   xxx = num_bp2_itr(i)
   if (xxx > 0.5) then 
     bp2_itr(i) = bp2_itr(i)/xxx
   else
     bp2_itr(i) = bad
   endif

   xxx = num_bp2_prp(i)
   if (xxx > 0.5) then 
     bp2_prp(i) = bp2_prp(i)/xxx
   else
     bp2_prp(i) = bad
   endif

!------------------------------------------------------------------------
!  end loop over i
!------------------------------------------------------------------------
  end do

!------------------------------------------------------------------------
!  OBS?
!------------------------------------------------------------------------
  if (av_upmp > 1.0E-08) then
    av_updz = av_updz/av_upmp
  else
    av_updz = bad
  endif

!------------------------------------------------------------------------
!
!  Diagnostics 
!
!  - normalize quantities which scale with mass_detrain
!  - These variables are for a individual call to updrafts only. When averaging
!    over columns, need to carry both tdiff_detrain and mass_detrain to take
!    overall mass weighted average.
!
!------------------------------------------------------------------------
do i = 1,nd
  if (mass_detrain_up(i) > 1.0E-10) then
    tdiff_detrain_up(i) = tdiff_detrain_up(i)/mass_detrain_up(i)
    cape_detrain_up(i) = cape_detrain_up(i)/mass_detrain_up(i)
    hm_diff_det(i) = hm_diff_det(i)/mass_detrain_up(i)
  else
    tdiff_detrain_up(i) = 0.0
    cape_detrain_up(i) = 0.0
    hm_diff_det(i) = bad
  endif
end do

!------------------------------------------------------------------------
!
!  END   (subroutine updrafts)
!
!------------------------------------------------------------------------
end subroutine updrafts
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Sets various arrays to zero that are re-defined every up parcel
!     timestep.
!
!------------------------------------------------------------------------
subroutine up_init
!------------------------------------------------------------------------
!
!   Use
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: i
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
  do i = 1,nd
    diffcon(i) = bad
    rvp(i) = 0.
    rlp(i) = 0.
    rip(i) = 0.
    rsp(i) = 0.
    mp(i) = 0.
    tp(i) = 0.
    tdp(i) = 0.
    tvp(i) = 0.
    hmp(i) = 0.
  end do
!------------------------------------------------------------------------
!
!   end
!
!------------------------------------------------------------------------
end subroutine up_init
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!    subroutine precipitate
!
!   - Jan 2014: this program is now called from both updrafts and detrain,
!     so that what was done to the "excess condensate" in both cases would
!     be consistent.
!   - constructed to avoid reference to parcel arrays
!   - The excess condensate goes to uprain, ansnow_conv, or detraining rv.
!     In the case of detraining rv, it is assumed that the liquid condensate
!     evaporates immediately in the background atmosphere.
!
!   - Jan 2010: what to do about momentum of precipitate? Decided to add to
!     momentum of updraft parcel whch produced it. Precipitate assumed to
!     have zero momentum, and instead leaves momentum behind.
!
!------------------------------------------------------------------------
subroutine precipitate(tt,zz,ii,rvpp,rlp_in,rip_in,hmpp,mpp,tdpp,  &
                       from_detrain,full_detrain,nl,tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only :  rdd,latent,hmoist,epsi,cl,ci,fusion,g,hmerr_max,get_rv,tkelvin
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!  rvpp: 
!   - The incoming value here is actually not used. The outgoing value is
!     purely determined from the value coming out of the call to get_rv used
!     to determine the final parcel properties. However, it seems unlikely
!     that the incoming and outgoing values should be different, since this
!     program simply removes condensate from the parcel, and should not
!     increase supersaturation to decrease rvpp.
!
!------------------------------------------------------------------------
real(r8), intent(inout) :: tt     ! can increase to freezing warming of ansnow
real(r8), intent(in) :: zz        !
integer, intent(in) :: ii         ! height index
real(r8), intent(inout) :: rvpp   ! Should not change
real(r8), intent(inout) :: rlp_in ! liquid condensate of parcel. Is adjusted.
real(r8), intent(inout) :: rip_in ! ice condensate of parcel 
real(r8), intent(inout) :: hmpp   ! updated parcel hm
real(r8), intent(in) :: mpp       ! dry mass
real(r8), intent(inout) :: tdpp   ! updated density temperature
integer, intent(in) :: from_detrain ! 
integer, intent(in) :: full_detrain ! 
integer, intent(in) :: nl           !  launch index
real(r8), intent(in) :: tstep

!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: hmppp,hmerr,hminit,dhm_rv,dhm_up,dhm_an,hm_of_wat,hm_of_ice
real(r8) :: rtpp,dmdry,kmpp,diff,tiny,summ
integer :: i, allow_sat
real(r8) :: rip_an,rlp_up,rlp_an,rip_ri,rlp_rv,rip_rv
real(r8) :: rlp_precip,dp,error,rlp_add,rlp_excess,rip_excess
real(r8) :: rvf,ri_fraction, rho
real(r8) :: rh_max,rh_min,pp,drh
real(r8) :: rlp_rem,rlp_max, tt_in, rvpp_in
real(r8) :: rip_rem,rip_max,f_precip, rlp_rv_max
real(r8) :: rh_cond_min, rh_cond_max, lwc_use

!print*,'entering precipitate ii = ',ii
!------------------------------------------------------------------------
!
!  Check if tt in range
!
!------------------------------------------------------------------------
  if ((tt.ge.180.).or.(tt.le.320.)) then
  else
    print*,'tt = ',tt
    stop 'tt out of range on entry to precipitate'
  endif

!------------------------------------------------------------------------
!
!  Check for negative rlp_in on entry
!
!------------------------------------------------------------------------
  if ((rlp_in < 0.).or.(rip_in < 0.)) then
    print*,'----- rlp_in or rip_in NEGATIVE on entry to subroutine precipitate --------'
    print*,'rvpp = ',rvpp
    print*,'rlp_in = ',rlp_in
    print*,'rip_in = ',rip_in
    print*,'ii = ',ii
    print*,'zz = ',zz
    print*,'tt = ',tt
    print*,'hmpp = ',hmpp
    print*,'mpp = ',mpp
    print*,'--------  move on -----------'
!   stop 'rlp_in negative on entry to subroutine precipitate'
  endif

!------------------------------------------------------------------------
!
!   Define rtpp
!
!------------------------------------------------------------------------
  rtpp = rvpp + rlp_in + rip_in
  if (rtpp < 0.) then
    print*,'----- rtpp NEGATIVE on entry to subroutine precipitate --------'
    print*,'rvpp = ',rvpp
    print*,'rlp_in = ',rlp_in
    print*,'rip_in = ',rip_in
    print*,'ii = ',ii
    print*,'zz = ',zz
    print*,'tt = ',tt
    print*,'hmpp = ',hmpp
    print*,'mpp = ',mpp
    stop 'rtpp negative on entry to subroutine precipitate'
  endif

!------------------------------------------------------------------------
!
!  Check for hm consistency
!  in subroutine precipitate
!  For a function, all arguments are "in" variables so none
!    will be adjusted.
!
!------------------------------------------------------------------------
! print*,'calling hmoist zz = ',zz
  hmppp = hmoist(tt,rtpp,rvpp,rip_in,zz)
! print*,'Got hmpp = ',hmpp
  hmerr = abs(hmpp - hmppp)/hmpp
  if (hmerr > 1.0E-06) then
    print*,'============== Start of precipitate ================'
    print*,'Inconsistency between actual + input hm'
    print*,'relative error hmerr = ',hmerr
    print*,'hm calculated from input temp rt rvp = ',hmppp
    print*,'input hm = ',hmpp
    print*,'zz = ',zz
    print*,'tt = ',tt
    print*,'rtpp = ',rtpp
    print*,'rvpp = ',rvpp
    print*,'rlp_in = ',rlp_in
    print*,'rip_in = ',rip_in
    print*,'from_detrain = ',from_detrain
    print*,'full_detrain = ',full_detrain
    stop 'hmerr problem'
  endif

!------------------------------------------------------------------------
!
!  Define lwc_use
!
!------------------------------------------------------------------------
!  if (tt > t_lwc_min) then
!    lwc_use = lwc_max
!  elseif (tt > t_lwc_max) then
!    xxx = (t_lwc_min - tt)/(t_lwc_min - t_lwc_max)
!    lwc_use = (1.-xxx)*lwc_max 
!  else
!    lwc_use = 0.
!  endif

!------------------------------------------------------------------------
!
!   Determine rlp_max
!
!   - max allowed rlp of parcel
!   - set equal to zero for detraining air parcel
!
!------------------------------------------------------------------------
  if (from_detrain == 0) then  ! calling precipitation from updraft
!   print*,'finding rlp_max'
    if ((ii <= i_rem_min).or.(rlp_in < rlp_thresh)) then
      rlp_rem = 0.
    else
      rho = 0.5*(phalf(ii) + phalf(ii+1))/(rdd*tt)
      rlp_max = (lwc_max*0.001)/rho
      if (ii > 1) then
        rlp_rem = rlp_in*f_rem(nl)*(z(ii) - z(ii-1))*0.001
        if ((rlp_in-rlp_rem) > rlp_max) then 
          rlp_rem = rlp_in - rlp_max
        endif
      else
        rlp_rem = rlp_in*f_rem(nl)*(z(ii+1) - z(ii))*0.001
        if (ii == nd) stop 'ii too high in if_conv_tend.f90'
      endif
    endif
    rlp_max = rlp_in - rlp_rem
    if (rlp_max < 0.) rlp_max = 0.
  else
    rlp_max = 0.                ! for detraining air parcel
  endif

!  if (from_detrain == 0) then  ! calling precipitation from updraft
!    rho = 0.5*(phalf(ii) + phalf(ii+1))/(rdd*tt)
!    rlp_max = (lwc_use*0.001)/rho
!    if (rlp_max < 0.) rlp_max = 0.
!    if (rlp_max > rlp_in) rlp_max = rlp_in
!  else
!    rlp_max = 0.                ! for detraining air parcel
!  endif

!rip_max = rlp_max
!if (rip_max > rip_in) rip_max = rip_in

!------------------------------------------------------------------------
!
!   Determine rip_max
!
!   - max allowed rip of parcel
!   - set equal to zero for detraining air parcel
!   - Remove a fixed amount per km.
!
!------------------------------------------------------------------------
  if (from_detrain == 0) then  ! calling precipitation from updraft
    if ((ii <= i_rem_min).or.(rip_in < rip_thresh)) then
      rip_rem = 0.
    else
      if (ii > 1) then
        rip_rem = rip_in*f_rem(nl)*(z(ii) - z(ii-1))*0.001
      else
        rip_rem = rip_in*f_rem(nl)*(z(ii+1) - z(ii))*0.001
        if (ii == nd) stop 'ii too high in if_conv_tend.f90'
      endif
    endif
    rip_max = rip_in - rip_rem
    if (rip_max < 0.) rip_max = 0.
  else
    rip_max = 0.                ! for detraining air parcel
  endif

!------------------------------------------------------------------------
!
!  Define rlp_excess (condensate to be removed from parcel).
!
!------------------------------------------------------------------------
  rlp_excess = rlp_in - rlp_max
  rip_excess = rip_in - rip_max

!  if (rlp_excess < 0.) rlp_excess = 0.
!  if (rip_excess < 0.) rip_excess = 0.

!------------------------------------------------------------------------
!
!   **************   Part I: Calculate rlp_up, rip_an   ***********
!
!------------------------------------------------------------------------
  if ((rip_excess > 0.).or.(rlp_excess > 0.)) then

!------------------------------------------------------------------------
!
!   *********  Inside rlp_excess > 0  ***********************************
!
!  Update outgoing parcel condensate
!  - When precipitate is being called for detrain rtpp and rvpp should now
!    be equal
!
!------------------------------------------------------------------------
     rlp_in = rlp_max   ! Define new value for "out"
     rip_in = rip_max   ! Define new value for "out"

     rtpp = rvpp + rlp_in + rip_in  ! adjusted value of total parcel rt

!------------------------------------------------------------------------
!
!  Determine rl_fraction: fraction of rlp_excess to grid scale condensate
!    (June 2016: now set to zero)
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  Determine ri_fraction: fraction of rip_excess to condensate
!
!------------------------------------------------------------------------
     summ = rl(ii) + ri(ii) + ri_added(ii)
     if (summ > ri_tar(ii)) then 
       ri_fraction = 0.
     else
       ri_fraction = ri_fraction_max
     endif

!print*,'ri_fraction = ',ri_fraction
!------------------------------------------------------------------------
!
!  Sometimes set detrained condensate to zero if not full detrain
!
!------------------------------------------------------------------------
  if ((full_detrain_only_ice == 1).and.(full_detrain == 0)) then
    ri_fraction = 0.
  endif

  if (ri_fraction < 0.) ri_fraction = 0.

!------------------------------------------------------------------------
!
!   *********  Inside rlp_excess > 0  ***********************************
!
!  Define values allocated to condensate: rip_ri
!
!------------------------------------------------------------------------
  rip_ri = ri_fraction*rip_excess

!------------------------------------------------------------------------
!
!   *********  Inside rlp_excess > 0  ***********************************
!
!------------------------------------------------------------------------
if (ii <= 5) then
  rh_cond_min = rh_cond_min_bl
  rh_cond_max = rh_cond_max_bl
else
  rh_cond_min = rh_cond_min_ft
  rh_cond_max = rh_cond_max_ft
endif

!------------------------------------------------------------------------
!
!   *********  Inside rlp_excess > 0  ***********************************
!
!   Calculate f_precip
!   - fraction of excess rlp to precip
!
!------------------------------------------------------------------------
if (use_rh_cond == 1) then
  if (rh(ii) > rh_cond_max) then
    f_precip = 1.   
  elseif (rh(ii) > rh_cond_min) then
    xxx = (rh(ii) - rh_cond_min)/(rh_cond_max - rh_cond_min)
    f_precip = (1.-xxx)*rh_cond_rain + xxx
  else
    f_precip = rh_cond_rain
  endif
! if (z(ii) > 13000.) then
! print*,'--------------------'
! print*,'z(ii) = ',0.001*z(ii)
! print*,'rh(ii) = ',rh(ii)
! print*,'f_precip = ',f_precip
! endif
else
  f_precip = 1.0
endif
  
!------------------------------------------------------------------------
!
!   *********  Inside rlp_excess > 0  ***********************************
!
!  Balance of rlp_excess/rip_excess goes into uprain and ansnow
!  Check if excessive rlp detrainment (doesn't seem to help)
!  - mass of rlp detrained = rlp_rv*mpp
!  - enhancement in local rv = rlp_rv*mpp/dmdry
!  - enhancement in local RH: drh = (rlp_rv/rs(ii))*(mpp/dmdry)
!  - at max rlp_rv rh(ii)+drh = rh_max_cond
!  - rh(ii) + (rlp_rv_max/rs(ii))*(mpp/dmdry) = rh_max_cond
!  - (rlp_rv_max/rs(ii))*(mpp/dmdry) = rh_max_cond - rh(ii)
!  - rlp_rv_max/rs(ii) = (dmdry/mpp)*(rh_max_cond - rh(ii))
!  - rlp_rv_max = rs(ii)*(dmdry/mpp)*(rh_max_cond - rh(ii))
!
!  More complex exp involving rlp_dett (not used)
!  - accumulated mass of rlp detrained = rlp_dett + rlp_rv*mpp
!  - enhancement in local rv = (rlp_dett+rlp_rv*mpp)/dmdry
!  - enhancement in local RH: drh = (rlp_dett+rlp_rv*mpp)/(rs(ii)*dmdry)
!  - at max rlp_rv
!  - rh(ii) + drh = rh_max_cond
!  - rh(ii) + (rlp_dett+rlp_rv*mpp)/(rs(ii)*dmdry) = rh_max_cond
!  - rh(ii) + (rlp_dett/rs(ii)*dmdry) + (rlp_rv*mpp)/(rs(ii)*dmdry) = rh_max_cond
!  - (rlp_rv_max/rs(ii))*(mpp/dmdry) = rh_max_cond - rh(ii) - (rlp_dett/rs(ii)*dmdry)
!  - rlp_rv_max/rs(ii) = (dmdry/mpp)*[rh_max_cond - rh(ii) - (rlp_dett/rs(ii)*dmdry)]
!  - rlp_rv_max = rs(ii)*(dmdry/mpp)*[rh_max_cond - rh(ii) - (rlp_dett/rs(ii)*dmdry)]
!
!------------------------------------------------------------------------
  rlp_rv = (1.-f_precip)*rlp_excess
  dmdry = (phalf(ii) - phalf(ii+1))/(g*(1.+rv(ii)+rl(ii)+ri(ii)))
  rlp_rv_max = rs(ii)*(dmdry/mpp)*(rh_max_cond - rh(ii))

! More complex expression for rlp_rv_max involving rlp_dett
!  rlp_rv_max = rh_max_cond - rh(ii) - (rlp_dett(ii)/(rs(ii)*dmdry))
!  rlp_rv_max = rs(ii)*(dmdry/mpp)*rlp_rv_max

  if (rlp_rv_max < 0.) then
    rlp_rv = 0.
  elseif (rlp_rv > rlp_rv_max) then
    rlp_rv = rlp_rv_max
  endif

!  Update rlp_dett
!  rlp_dett(ii) = rlp_dett(ii) + rlp_rv*mpp

!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
  rlp_precip = rlp_excess - rlp_rv
  if (rlp_precip < 0.) rlp_precip = 0.

  if (ii >= i_ml) then
    rlp_up = (1.-f_an)*rlp_precip
    rlp_an = f_an*rlp_precip
  else
    rlp_up = rlp_precip
    rlp_an = 0.
  endif

  rip_an = f_precip*(rip_excess - rip_ri)
  rip_rv = (1.-f_precip)*(rip_excess - rip_ri)

!------------------------------------------------------------------------
!
!  Check for negative
!
!------------------------------------------------------------------------
  if ((rlp_rv < 0.).or.  &
      (rlp_up < 0.).or.  &
      (rlp_an < 0.).or.  &
      (rip_rv < 0.).or.  &
      (rip_an < 0.)) then
    print*,'NEGATIVE rlp_up or rip_an in precipitate'
    print*,'rlp_max = ',rlp_max
    print*,'rip_max = ',rip_max
    print*,'rlp_excess = ',rlp_excess
    print*,'rip_excess = ',rip_excess
    print*,'rlp_in = ',rlp_in
    print*,'rip_in = ',rip_in
    print*,'rlp_rv = ',rlp_rv
    print*,'rlp_up = ',rlp_up
    print*,'rlp_an = ',rlp_an
    print*,'rip_rv = ',rip_rv
    print*,'rip_an = ',rip_an
  endif

!------------------------------------------------------------------------
!
!  *******   Both rlp_excess and rip_excess < 0  *******************
!
!  - set the 6 possibilities to zero.
!  - leave rlp_in unchanged.
!  - could likely have a return at this point since doing nothing.
!
!------------------------------------------------------------------------
  else
    rlp_up = 0.      ! rlp to updraft rain
    rlp_an = 0.      ! rlp to updraft rain
    rip_an = 0.      ! rlp to ansnow
    rip_ri = 0.
    rlp_rv = 0.
    rip_rv = 0.
  endif

!------------------------------------------------------------------------
!
!  Define precipitation quantities
!
!    - just for diagnostics.
!    - These are quantities which scale with the mass of precipitating parcels.
!      This includes updraft parcels that are still going up, and parcels that
!      are detraining.
!
!------------------------------------------------------------------------
  if ((rip_an + rlp_up + rlp_an) > 0.00000001) then
    rlp_add = rlp_up + rlp_an + rip_an + rip_ri + rlp_rv + rip_rv
    mass_rlp_rv = mass_rlp_rv + rlp_rv*mpp
    mass_rip_rv = mass_rip_rv + rip_rv*mpp
    mass_rip_an = mass_rip_an + rip_an*mpp
    mass_rlp_up = mass_rlp_up + rlp_up*mpp
    mass_rlp_an = mass_rlp_an + rlp_an*mpp
  endif

!------------------------------------------------------------------------
!
!  For testing: 
!  If nprecip = 2 : turn off ansnow and uprain production.
!    Add condensate to background rv source.
!    Leave other possibilities unchanged.
!    PROBLEM HERE: rlp not allocated (was to rlp_rv) So the nprecip = 2
!      option would likely not conserve water
!
!------------------------------------------------------------------------
   if (nprecip == 2) then
     rlp_rv = 0.
     rip_rv = 0.
     rip_an = 0.
     rlp_up = 0.
     rlp_an = 0.
     rip_ri = 0.
   endif

!------------------------------------------------------------------------
!
!  Set rlp's to zero if very small.
!
!------------------------------------------------------------------------
 if (abs(rlp_rv) < 1.0E-14) rlp_rv = 0.
 if (abs(rip_rv) < 1.0E-14) rip_rv = 0.
 if (abs(rlp_up) < 1.0E-14) rlp_up = 0.
 if (abs(rlp_an) < 1.0E-14) rlp_an = 0.
 if (abs(rip_an) < 1.0E-14) rip_an = 0.
 if (abs(rip_ri) < 1.0E-14) rip_ri = 0.

!------------------------------------------------------------------------
!
!  Check for negative values. 
!
!------------------------------------------------------------------------
 if ((rlp_up < 0.).or.  &
     (rlp_an < 0.).or.  &
     (rip_an < 0.).or.  &
     (rip_ri < 0.)) then
   print*,'===========  rlp problem  =============='
    print*,'rlp_max = ',rlp_max
    print*,'rip_max = ',rip_max
    print*,'rlp_excess = ',rlp_excess
    print*,'rip_excess = ',rip_excess
    print*,'rlp_in = ',rlp_in
    print*,'rip_in = ',rip_in
   print*,'rlp_up = ',rlp_up
   print*,'rlp_an = ',rlp_an
   print*,'rip_an = ',rip_an
   print*,'rip_ri = ',rip_ri
   print*,'f_precip = ',f_precip
   print*,'rip_excess = ',rip_excess
   print*,'rip_ri = ',rip_ri
   print*,'rip_an = f_precip*(rip_excess - rip_ri)'
   stop 'one of rlp negative in precipitate'
 endif

!------------------------------------------------------------------------
!
!   **************   PRECIPITATE PART II - hm CHANGES ***********************
!
!  The 7 rlp's are now determined (may all be zero)
!
!  Determine MSE changes
!
!   rl has units kg water/kg dry mass
!   precip added = (change in rl)*(dry parcel mass)
!                = (kg water/kg dry air)*(kg dry air/m2)
!                = kg water/m2
!
!   dhmuprain: reduction in parcel hm due to loss of condensate (per kg dry air)
!              to updraft rain
!   [dhmuprain] = [kg water/kg dry air]*[J/kg water*K]*[K]
!             = [/kg dry air]*[J]
!             = [J/kg dry air]
!
!   units of hmflow are: J/m2/s
!
!   hmflow(numhm,nd)
!   numhm = 2: detrainmnet
!
!   Adjust hmuprain:
!   ----------------
!   hmuprain is absolute MSE in J/m2*s, so need to multiply by mass
!     of updraft parcel. Also add at all levels below (evaporation in
!     downdrafts can reduce).
!
!   From above:
!   [dkmp_rain] = [J/kg dry air]
!   Units of hmuprain = [dkmp_rain]*[mpp]*[1/tstep]
!                   = [J/kg dry air]*[kg dry air/m2]*[1/s]
!                   = [J/m2*s]
!
!  dhm_rv 
!  - the hm that goes into the background atmosphere from
!    the total detraining condensate. It doesn't matter which component of the 
!    background atmosphere this hm is allocated to, since T is later
!    determined from fixed rl/ri (independently determined from rlflow/riflow
!    and rlvert/rivert and the previous values).
!  - What is happening here is that dhm_rv, dhm_up, and dhm_an are simply
!    DEFINED here using the calculated values of the 4 rip_ri, rlp_up,
!    rlp_an, and rip_an.
!  - These 5 terms are then SUBTRACTED from the PARCEL hm when hmpp is
!    defined below as hmpp = hmpp - dhm_up - dhm_rv - dhm_an. This is
!    how conservation of mse is implemented. For example, rlp_an would
!    represent a loss of mse from the parcel to ansnow. However, since hm 
!    of ice is less than hm of water, the removal of this ice hm would 
!    effectively warm the remaining air parcel, by decreasing hmpp less than
!    it would otherwise.
!
!  hm_of_ice:
!  Suppose T = 220 K, z = 1200 m, cl = 4188, lf0 = 3.337E+05
!  hm_of_wat = 921360 + 12000 = 933360
!  hm_of_ice = 599660
!  
!
!
! PPP
!------------------------------------------------------------------------
hm_of_wat = cl*tt + g*zz
hm_of_ice = hm_of_wat - fusion(tt)

! hm to background atm
dhm_rv = rlp_rv*hm_of_wat + (rip_ri+rip_rv)*hm_of_ice  

dhm_up = rlp_up*hm_of_wat              ! hm going to uprain
dhm_an = (rlp_an+rip_an)*hm_of_ice     ! hm going to ansnow

!------------------------------------------------------------------------
!
!   **************   PRECIPITATE PART II - hm CHANGES ***********************
!
!  Adjust outflows to background atmosphere
!
!  - this is to take into account the portion of the condensate that is
!    assumed to evaporate immediately upon detrainment. This should result
!    in a cooling of the background atmosphere, provided the condensate that
!    has been assigned to rlp_rv is retained as rv in the background
!    atmosphere.
!  - dhm_rv is calculated above as if the condensate is liquid water (was rlp)
!  - The middle 1 below is for updrafts, but is actually arbitray whether 
!    allocated to updrafts or downdrafts.
!
!------------------------------------------------------------------------
  hmflow(2,ii) = hmflow(2,ii) + (dhm_rv*mpp/tstep)  ! To background atm hm

  riflow(2,ii) = riflow(2,ii) + (rip_ri*mpp/tstep) 
  rvflow(2,ii) = rvflow(2,ii) + ((rlp_rv+rip_rv)*mpp/tstep)

  dmdry = (phalf(ii) - phalf(ii+1))/(g*(1.+rv(ii)+rl(ii)+ri(ii)))

!------------------------------------------------------------------------
!
!  Update the cumulative amount of condensate produced
!
!------------------------------------------------------------------------
  ri_added(ii) = ri_added(ii) + (rip_ri*mpp/dmdry)   ! precip+detrain

! print*,'dmdry = ',dmdry

!------------------------------------------------------------------------
!
!   **************   PRECIPITATE PART II - hm CHANGES ***********************
!
!  Update  uprain prod 
!  Do not adjust hmflow, since hm is going to hmuprain not the background 
!    atmosphere. The loss of hm to uprain/ansnow is taken into account through
!    a reduction in the parcel hm.
!
!------------------------------------------------------------------------
  uprain_rlp(ii) = uprain_rlp(ii) + (rlp_up*mpp/tstep)
  ansnow_rlp(ii) = ansnow_rlp(ii) + ((rip_an+rlp_an)*mpp/tstep)
  rlp_rvv(ii) = rlp_rvv(ii) + (rlp_rv*mpp/tstep)

  if (nl == 1) uprain_1 = uprain_1 + (rlp_up*mpp/tstep)
  if (nl == 2) uprain_2 = uprain_2 + (rlp_up*mpp/tstep)

!------------------------------------------------------------------------
!
!   **************   PRECIPITATE PART II - hm CHANGES ***********************
!
!  Checks for negative values
!  - very unlikely since checked for negative rlp_up/rip_an already.
!  - maybe remove
!
!------------------------------------------------------------------------
    if (uprain_rlp(ii) < 0.) then
      print*,'ii uprain_rlp(ii) = ',ii,uprain_rlp(ii)
      stop 'uprain_rlp negative in precipitate'
    endif
    if (ansnow_rlp(ii) < 0.) then
      print*,'ii ansnow_rlp(ii) = ',ii,ansnow_rlp(ii)
      stop 'ansnow_rlp negative in precipitate'
    endif

!------------------------------------------------------------------------
!
!   **************   PRECIPITATE PART II - hm CHANGES ***********************
!
!  Adding hm of rlp to hm of precip
!  Adjust hmuprain and hmansnow_start: 
!
!  - up to level ii or ii+1?
!  - hmuprain is a half level variable.
!  I think ii, since uprain production at full level ii increase hmuprain
!    at the half level below it.
!
!------------------------------------------------------------------------
  hmuprain_start = hmuprain_start + (dhm_up*mpp/tstep)
  hmansnow_start = hmansnow_start + (dhm_an*mpp/tstep)

  hmuprain_prod(ii) = hmuprain_prod(ii) + (dhm_up*mpp/tstep)
  hmansnow_prod(ii) = hmansnow_prod(ii) + (dhm_an*mpp/tstep)

!------------------------------------------------------------------------
!
!   **************   PRECIPITATE PART II - hm CHANGES ***********************
!
!   Update parcel hmpp/kmpp
!
!   Discussion: 
!    - "precipitate" is called from subroutine updraft and subroutine detrain.
!    - When called from updraft, this hmpp
!     refers to the new hm of the updraft parcel, after taking into account
!     the loss of the excess condensate. 
!    - When precipitate is called from detrain it can be for a partial detrain
!      or a full detrain
!    - For a partial detrain, this hmpp refers to that fraction of the air parcel 
!      that is detraining.
!    - For a full detrain, hmpp refers to the entire air parcel
!    - But in both detrain cases, ALL of the condensate should have been 
!      redistributed (since rlp_max = 0), so hmpp should refer to the hm of the dry 
!      air parcel whose mass is dmd.
!
!   hmpp
!   - New value based on hm losses from precip and net detrainment into background atm
!     due to detrained condensate.
!
!------------------------------------------------------------------------
  hminit = hmpp  ! store old value
  hmpp = hmpp - dhm_up - dhm_rv - dhm_an 
  kmpp = hmpp - (1. + rtpp)*g*zz
! print*,'in precipitate kmpp = ',kmpp

!------------------------------------------------------------------------
!
!   **************   PRECIPITATE PART II - hm CHANGES ***********************
!
!  Call get_rv to determine new parcel properties.
!
!   - Use kmpp in get_rv to determine new T, rlp_in, rvpp
!
!  - In the call to get_rv, rvpp, rlp_in, and tt are purely "out" variables,
!    so can change. However, in practice, it seems rvpp should not change
!    from value given above, since rtpp has gone down. kmpp has also 
!    decreased, however, so may in principle be possible for rvpp to change.
!
!  - The rlp_in below is misnamed. It is purely an output variable
!    for get_rv. The input rlp is implicit in rtpp and rvpp. 
!
!  - Precipitate should  not remove condensate below max_rlp.
!    However, when ice is removed, there is a smaller reduction
!    in parcel enthalpy, so effective heating of the parcel. This can lead to
!    evaporation of condensate, and possibly, removal of additional condensate,
!    so that the final parcel condensate value is below max_rlp.
!
!   - In get_rv: kmpp, rtpp, p(ii): pure in
!             rvpp, rlp_in, tt: pure out
!
!   Note: in the case of full or partial detrainment this call is not needed. All
!     you care about is the hm of the parcel, so that hmflow is correctly
!     calculated in "detrain". You don't need the temperature.
!
!   Calling get_rv from subroutine precipitate
!
!    - in the case of a detraining air parcel, the call to get_rv is for a dry air
!      parcel only, so no change.
!    - In the case of a call from updraft, the rl and ri of the parcel would have
!      been modified from the input values. I think kmpp and rtpp have to be adjusted
!      to remove the ice part, and then again another call to get_t would be needed.
!
!------------------------------------------------------------------------
  allow_sat = 0  ! do not allow saturation: produce rl if needed
  get_rv_call = 4
  call_get_rv_precipitate = call_get_rv_precipitate + 1
  tt_in = tt
! print*,'call_get_rv_precipitate =',call_get_rv_precipitate
  call get_rv(kmpp,  &
              rtpp,  &
             p(ii),  &
              rvpp,  &  !
            rlp_in,  &  ! 
            rip_in,  &  ! updated parcel ri (fixed)
                tt,  &  ! out
       get_rv_call,  &  ! pure in
       allow_sat)

!------------------------------------------------------------------------
!
!  I wanted to find out if this call to get_rv was accomplishing anything.
!  It does seem that rv changes here, i.e. rvpp-rvpp_in = 0.01 typically,
!  and also small changes in T.
!
!------------------------------------------------------------------------

! print*,'Change in temperature call to get_rv in precipitate = ',tt-tt_in
! print*,'Change in rv call to get_rv in precipitate = ',rvpp-rvpp_in
! if (abs(tt-tt_in) > 1.0E-4) then
!   print*,'CHANGE in TT from call get_rv in precipitate'
! endif
! if (abs(rvpp-rvpp_in) > 1.0E-4) then
!   print*,'CHANGE in RVV from call get_rv in precipitate'
! endif

!------------------------------------------------------------------------
!
!   **************   PRECIPITATE PART II - hm CHANGES ***********************
!
!   - Perhaps unneccessary check
!   - Use new tt and rvpp to determine hmpp
!   - Check for inconsistency in parcl hm between incoming value and
!     value with modified rip_in,rvpp
!   hminit: parcel hm prior to condensate removal
!   April/2010: better to define hmpp in terms of absolute changes.
!
!------------------------------------------------------------------------
  hmpp = hmoist(tt,rtpp,rvpp,rip_in,zz)
  hmerr = abs(hminit - hmpp - dhm_up - dhm_rv - dhm_an)/hmpp
  if (hmerr > hmerr_max) then
    print*,'in precipitate: problem with change in hm due to precip'
    print*,'relative hm error hmerr = ',hmerr
    print*,'hminit = ',hminit
    print*,'hmpp = ',hmpp
    print*,'tt = ',tt
    print*,'zz = ',zz*0.001
    print*,'ii = ',ii
    stop 'hm change problem in precipitate'
  endif

!------------------------------------------------------------------------
!
!   Comments on changes in updraft momentum during precipitation (Jan 2010)
!
!     - momentum of updraft condensate that is converted to updraft
!       precipitation is added back to updraft parcel (precipitation
!       is assumed to not carry momentum)
!     - momentum of the parcel = uwindp*(dry mass of parcel)
!     - since the dry mass of the parcel does not change, to conserve
!       momentum the uwindp of the parcel must not change either.
!     - this appears confusing, since the real wind of the parcel must
!       increase when the precipitate returns its momentum to it.
!     - However uwind(real) = uwind/(1 + rtt), so since the rtt of the
!       updraft parcel decreases during precipitation formation, the real
!       wind increases if uwindp of the parcel stays constant
!     - the bottom line: no adjustment of uwindp/vwindp is required for
!       updraft precipitation that carries away no momentum, even though
!       the real u/v winds of the updraft parcel increase.
!
!   Comments on implementing momentum conservation for anvil (condensate)
!
!     - You want the anvil condensate to have the same u/v as the background
!       atmosphere level at which it detrains. This is assumed to be the case
!       in subroutine conserve where the column momentum is calculated.
!     - It turns out that column momentum is conserved when there is no change
!       in uwindp/vwindp when anvil condensate is formed (i.e. rl), and there
!       is no change in uflow/vflow associated with the detrainment of this
!       condensate into the background atmosphere.
!     - Again, it comes down to the fact that the addition or removal of water
!       from a layer, if the water is not assumed to have any momentum, must
!       not changed the uwind/vwind (per unit dry mass of a layer)
!
! The bottom line: ignore momentum considerations in subroutine precipitate
!
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!
!  Update density temperature
!
!------------------------------------------------------------------------
  tdpp = tt*(1. + epsi*rvpp)/(1. + rtpp)

!------------------------------------------------------------------------
!
!  Check for negative rtpp
!
!------------------------------------------------------------------------
  tiny = -1.0E-16
  if ((rtpp < tiny).or.(rip_in < tiny).or.(rlp_in < tiny)) then
    print*,'------------------------------------'
    print*,'rtpp = ',rtpp
    print*,'rip_in = ',rip_in
    print*,'rlp_in = ',rlp_in
    print*,'tt = ',tt
    print*,'zz = ',zz
    print*,'ii = ',ii
    print*,'rvpp = ',rvpp
    print*,'hmpp = ',hmpp
    print*,'mpp = ',mpp
    print*,'nprecip = ',nprecip
    print*,'hmerr = ',hmerr
    stop 'negative water problem in subroutine precipitate'
  endif

!print*,'exiting precipitate'
!------------------------------------------------------------------------
!
!   End
!
!------------------------------------------------------------------------
end subroutine precipitate
!------------------------------------------------------------------------
!**********************************************************************73
!--------------------------------------------------------------------
!
!  subroutine detrain
!
!  Called only from subroutine updrafts
!
!  Adjusts the dry air parcel mass mp the tendencies.
!
!  In the case of entrainment, sigma = dMd/dMd+Mdp
!  In the case of detrainment, the same applies, except dMd < 0 for sigma < 0.
!  Mdp is the initial dry mass of the parcel.
!  The detrained dry mass is dMd = sigma*Mp(n)/(1 - sigma)
!  The detrained vapor mass is dMv = rvp(i)*dMd
!  The detrained condensate mass is dMl = rlp(i)*dMd
!
!  Inputs:
!    mass : the input dry parcel mass. [kg/m2]
!
!  Outputs:
!    mass : the adjusted dry parcel mass. [kg/m2]
!
!  full_detrain = 1: detrain entire air parcel
!  full_detrain = 0: use sigma to calculate dmd
!
!------------------------------------------------------------------------
subroutine detrain(i,ip,it,full_detrain,sigma,nl,tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : hmoist,enthalpy,latent,tkelvin,cl,g,ci,fusion
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In variables
!
!------------------------------------------------------------------------
integer, intent(in) :: i     ! model level at which detrainment is occurring
integer, intent(in) :: ip    ! parcel index (same as i for updrafts)
integer, intent(in) :: it    ! parcel start index 
integer, intent(in) :: full_detrain  !
real(r8), intent(in) :: sigma
integer, intent(in) :: nl
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!  Local Variables
!
!------------------------------------------------------------------------
real(r8) :: kmp,rtp,dmd,small
integer  :: nums,mm,ii,from_detrain
real(r8) :: tp_det, tdp_det, rvp_det, rlp_det, rip_det, hmp_det

!------------------------------------------------------------------------
!
!   Check for sign error in input mass
!
!------------------------------------------------------------------------
  if (mp(ip) < 0.) then
    print*,'mp(ip) = ',mp(ip)
    stop 'mp(ip) negative on entry to detrain'
  endif

!------------------------------------------------------------------------
!
!  ************ In subroutine detrain ********************************
!
!------------------------------------------------------------------------
  if (sigma > 0.) stop 'sigma positive in detrain'

!------------------------------------------------------------------------
!
!  **********  Calculate dmd for partial detrain  ********************
!
!  full_detrain = 0: use sigma to calculate dmd
!  full_detrain = 1: detrain entire air parcel
!
!  dMd = sigma*mp(n)/(1 - sigma)
!  - note that sigma < 0 here since detrain was called, so dmd should
!    be positive
!
!------------------------------------------------------------------------
  small = 0.0000000001
  if (full_detrain == 0) then
!   dmd = - mp(ip)*(sigma/(1.-sigma))
    dmd = - mp(ip)*sigma
    if (dmd < 0.0) then
       print*,'----------------------------------------'
       print*,'in detrain sigma = ',sigma
       print*,'mp(ip) = ',mp(ip)
       print*,'dmd = ',dmd
       stop 'dmd gone negative'
    endif
    mp(ip) = mp(ip) - dmd
    if (dmd < small) then
!      print*,'dmd mp(ip) full_detrain sigma = ',dmd,mp(ip),full_detrain,sigma
!      stop 'dmd less than small in full_detrain = 0'
    endif

!------------------------------------------------------------------------
!
!  **********  Calculate dmd for full detrain  ********************
!
!  full_detrain = 1:
!  Full detrainment applies to updrafts which have reached their LNB, and
!    all downdrafts.
!
!------------------------------------------------------------------------
  elseif (full_detrain == 1) then
!------------------------------------------------------------------------
!
!   Adjust masses.
!
!------------------------------------------------------------------------
    dmd = mp(ip)
    mp(ip) = 0.
!------------------------------------------------------------------------
!
!  full_detrain out of range
!
!------------------------------------------------------------------------
  else
     stop 'full_detrain not defined correctly'
  endif

!------------------------------------------------------------------------
!
!  Initialize detraining parcel properties
!
!  - The mass of the detraining parcel is dmd.
!  - These may be modified in below call to precipitate (e.g. rlp_det should
!    end up being zero always), so not the final values.
!  - These have to be internally consistent; if not flagged on entry to
!    precipitate.
!
!------------------------------------------------------------------------
  tp_det = tp(ip)
  tdp_det = tdp(ip)
  rvp_det = rvp(ip)
  rlp_det = rlp(ip)
  rip_det = rip(ip)
  hmp_det = hmp(ip)

!------------------------------------------------------------------------
!
!  Call precipitate to determine detraining parcel properties.
! 
!  - In the case of a partial detrain, you do not want to change the parcel
!    properties. Changes are applied only to the dmd fraction of the original 
!    air parcel.
!
!  - Coming out of this call to precipitate, the only quantities that matter
!    are those that affect the flow variables below : dmd, rvp_det, and hmp_det.
!    (rlp_det should be zero).
!
!  Call precipitate from detrain  (can be partial or full)
!
!------------------------------------------------------------------------
from_detrain = 1

! print*,'Calling precipitate from detrain full_detrain = ',full_detrain
if (dmd > 1.0E-08) then
  call precipitate(  &
     tp_det,         &  ! inout: modified detraining parcel T
     z(i),           &  ! in:
     i,              &  ! in: height index of background atm
     rvp_det,        &  ! inout: may change
     rlp_det,        &  ! inout: should be zero on exit
     rip_det,        &  ! inout 
     hmp_det,        &  ! inout: modified due to condensate removal
     dmd,            &  ! in: dry mass
     tdp_det,        &  ! inout: updated due to condensate removal
     from_detrain,   &  !
     full_detrain,   &  ! in:
     nl,tstep)          ! in:
endif

!------------------------------------------------------------------------
!
!  Diagnostics: Define detrainment quantities
!
!    - just for diagnostics; not used in temp tendency. 
!    - These are quantities which scale with the mass of detraining parcels.
!      Quantities which scale with the mass of precipitating parcels should
!      be defined in precipitate.
!
!------------------------------------------------------------------------

  tdiff_detrain_up(i) = tdiff_detrain_up(i) + (tp_det - t(i))*dmd
  cape_detrain_up(i) = cape_detrain_up(i) + cape_detrain*dmd
  hm_diff_det(i) = hm_diff_det(i) + (hmp_det-hm(i))*dmd
  mass_detrain_up(i) = mass_detrain_up(i) + dmd

!------------------------------------------------------------------------
!
!  ************ In subroutine detrain ********************************
!
!  Adjust detrainment flows
!  - What is returned from the above call to precipitate are the new parcel
!    properties in which the parcel no longer has condensate. The hm of the parcel
!    should be adjusted in precipitate, and all flows associated with the
!    condensate of detraining air parcels handled there.
!
!  units of rvflow,etc: [kg water/m2*s]
!
!------------------------------------------------------------------------
  dmflow(2,i) = dmflow(2,i) + (dmd/tstep)
  rvflow(2,i) = rvflow(2,i) + (rvp_det*dmd/tstep)
  hmflow(2,i) = hmflow(2,i) + (hmp_det*dmd/tstep)
  uflow(2,i) = uflow(2,i) + (uwindp(ip)*dmd/tstep)
  vflow(2,i) = vflow(2,i) + (vwindp(ip)*dmd/tstep)


!------------------------------------------------------------------------
!
!  ************ In subroutine detrain ********************************
!
!  Update diagnostics
!
!------------------------------------------------------------------------

  if (nl == 1) then
    updet_1(i) = updet_1(i) + (dmd/tstep)
  elseif (nl == 2) then
    updet_2(i) = updet_2(i) + (dmd/tstep)
  endif

!------------------------------------------------------------------------
!
!  Detraining condensate - OBS?
!
!  Jan 2014: All condensate of detraining updraft air parcels is now handled
!    via a call to precipitate, and there assumed to either precipitate or
!    detrain as rv, so rlp/rip here should be zero. All potential
!    sources of rlflow/riflow are now zero, so these arrays should always be
!    zero. The only possible exception is if there is condensate formation in
!    a downdraft.
!
!  - rip should always be zero.
!  - the implied change of phase between water and ice here does not
!    matter. Is included in the temperature via hm conservation. 
!    (Not 100% sure)
!
!------------------------------------------------------------------------
  if (t(i) < tkelvin) then
    riflow(2,i) = riflow(2,i) + ((rlp_det+rip_det)*dmd/tstep) 
  else
    rlflow(2,i) = rlflow(2,i) + ((rlp_det+rip_det)*dmd/tstep)
  endif 

!------------------------------------------------------------------------
!
!   End
!
!------------------------------------------------------------------------
end subroutine detrain
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   function sigmoidal
!
!------------------------------------------------------------------------
real(r8) function sigmoidal(x_in,x_half,x_scale,y_min,y_max)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: x_in
real(r8), intent(in) :: x_half
real(r8), intent(in) :: x_scale
real(r8), intent(in) :: y_min
real(r8), intent(in) :: y_max
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
real(r8) :: tt
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
tt = (x_in - x_half)/x_scale
sigmoidal = y_min + (y_max/(1. + exp(-tt))) 
!------------------------------------------------------------------------
!
!    END     (function sigmoidal)
!
!------------------------------------------------------------------------
end function sigmoidal
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine get_stuff
!  - should be broken up into individual subroutines for clarity
!
!------------------------------------------------------------------------
subroutine get_stuff
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: i,j,ngot_ml,ngot,ngot_strat,ngot_hm_min
real(r8) :: dmdry,rv_thresh,rh_diff,dp
real(r8) :: col_wat_have,col_wat_sat,z_low,z_high
real(r8) :: tot_mass,ff,hm_min,mass_total
!------------------------------------------------------------------------
!
!  Determine i_ml
!  - first level from surface below 0 C (can have i_ml = 1)
!
!------------------------------------------------------------------------
  i_ml = 0
  i_strat = 0

  ngot_strat = 0
  ngot_ml = 0
  ngot_hm_min = 0

  hm_min = 500000.

  do i = 1,nd
!------------------------------------------------------------
!  ml
!------------------------------------------------------------
    if (ngot_ml == 0) then
    if (t(i) < tkelvin) then
      i_ml = i
      ngot_ml = 1
    endif
    endif
!------------------------------------------------------------
!  i_strat
!------------------------------------------------------------
    if (ngot_strat == 0) then
    if (t(i) < t_strat) then
      i_strat = i
      ngot_strat = 1
    endif
    endif
!------------------------------------------------------------------------
!  End loop over i
!------------------------------------------------------------------------
  end do

  if (i_strat == 0) stop 'i_strat not calculated correctly'
  if (i_ml == 0) stop 'i_ml not calculated correctly'

!------------------------------------------------------------------------
!
!  Compute col_rl,colrh
!  Dec 2013: changed definition of colrh
!  62.5 mm col water corresponds to 0.89
!  62.5/0.89 = 70.22
!
!------------------------------------------------------------------------
  colwat = 0.
  col_rl = 0.
  colrh = 0.
  z_low = 2000.
  z_high = 18000.
  col_wat_sat = 0.
  col_wat_have = 0.
  do i = 1,nd
    dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
    col_rl = col_rl + dmdry*rl(i) 
    if (z(i) < z_high) then
      col_wat_sat = col_wat_sat + dmdry*rs(i)
      col_wat_have = col_wat_have + dmdry*rv(i)
    endif
    colwat = colwat + dmdry*rv(i)
  end do
  colrh = col_wat_have/col_wat_sat
!------------------------------------------------------------------------
!
!    END  (subroutine get_stuff)
!
!------------------------------------------------------------------------
end subroutine get_stuff
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
subroutine bl_mix(land_fraction,tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g, hmoist
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: land_fraction
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: i
integer :: i_rem, i_add, nn_print
real(r8) :: dmdry, dmd, lr
real(r8) :: cinf, drv, conv_factor
real(r8) :: rv_mass, rv_rem, mass_rem, hm_rem, ri_rem, cappp

f_1 = 0.
f_2 = 0.
f_3 = 0.
f_4 = 0.

!------------------------------------------------------------------------
!
!  RH mixing
!
!------------------------------------------------------------------------
do i = 1,4
!------------------------------------------------------------------------
!
!  Define conv_factor
!
!------------------------------------------------------------------------
if (i == 1) cappp = cape_1
if (i == 2) cappp = cape_2
if (i == 3) cappp = cape_3
if (i == 4) cappp = cape_4
conv_factor = cappp/cape_rh_scale
if (conv_factor < 0.) conv_factor = 0.
!------------------------------------------------------------------------
!
! If RH exceeds limit move up one layer
!
!------------------------------------------------------------------------
dmd = 0.
dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
if (do_rh_mix == 1) then
if (conv_factor > 0.) then
if (rh(i) > rh_max_bl) then
   dmd = conv_factor*amp_rh_factor*(rh(i) - rh_max_bl)*dmdry
   if (dmd > fraction_max*dmdry) dmd = fraction_max*dmdry
endif
endif
endif
!------------------------------------------------------------------------
!  check for dmd negative
!------------------------------------------------------------------------
   if (dmd < 0.) then
     print*,'--------- dmd negative in if_conv_tend.f90 ----'
     print*,'dmd = ',dmd
     print*,'amp_rh_factor = ',amp_rh_factor
     stop 'dmd gone negative'
   endif
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
   if (i == 1) f_1 = dmd/dmdry
   if (i == 2) f_2 = dmd/dmdry
   if (i == 3) f_3 = dmd/dmdry
   if (i == 4) f_4 = dmd/dmdry
!------------------------------------------------------------------------
!  entrain
!------------------------------------------------------------------------
  dmflow(1,i) = dmflow(1,i) + (dmd/tstep)
  rvflow(1,i) = rvflow(1,i) + (rv(i)*dmd/tstep)
  hmflow(1,i) = hmflow(1,i) + (hm(i)*dmd/tstep)
  uflow(1,i) = uflow(1,i) + (uwind(i)*dmd/tstep)
  vflow(1,i) = vflow(1,i) + (vwind(i)*dmd/tstep)
!------------------------------------------------------------------------
!  detrain
!------------------------------------------------------------------------
  if (i == nd) stop 'i problem in if_conv_tend.f90'
  dmflow(2,i+1) = dmflow(2,i+1) + (dmd/tstep)
  rvflow(2,i+1) = rvflow(2,i+1) + (rv(i)*dmd/tstep)
  hmflow(2,i+1) = hmflow(2,i+1) + (hm(i)*dmd/tstep)
  uflow(2,i+1) = uflow(2,i+1) + (uwind(i)*dmd/tstep)
  vflow(2,i+1) = vflow(2,i+1) + (vwind(i)*dmd/tstep)
!------------------------------------------------------------------------
!
!  end loop over levels
!
!------------------------------------------------------------------------
end do

!------------------------------------------------------------------------
!
!    END  (subroutine bl_mix)
!
!------------------------------------------------------------------------
end subroutine bl_mix
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
subroutine get_km_spectrum
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: it, nl
real(r8) :: km_from_hm

!------------------------------------------------------------------------
!
!  Loop over it from 1 to maxlev
!  - how do I know what km(it) is here? KKKK
!
!------------------------------------------------------------------------
do it = 1,maxlev
do nl = 1,nlaunch
  if (nl_depend == 1) then
    km_inc = km_width/float(nlaunch - 1)
    kmstart(it,nl) = km(it) + float(nl - 1)*km_inc
  else
    kmstart(it,nl) = km(it) + km_width
  endif
  km_from_hm = hm(it) - (1. + rt(it))*z(it)*g
  if (kmstart(it,nl) > km_bad) then
    print*,'======= kmstart too BIG ========='
    print*,'kmstart(it,nl) = ',kmstart(it,nl)
    print*,'km(it) = ',km(it)
    print*,'km_width = ',km_width
    stop 'km too big'
  endif
end do
end do

!------------------------------------------------------------------------
!
!    END  (subroutine get_km_spectrum)
!
!------------------------------------------------------------------------
end subroutine get_km_spectrum
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
subroutine get_cape_spectrum
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : 
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: nl, it,ilnb
real(r8) :: cin, t_start, rv_start

!------------------------------------------------------------------------
!
!  Loop over it from 1 to maxlev
!
!------------------------------------------------------------------------
do it = 1,maxlev
if (t(it) > t_min_cape) then
do nl = 1,nlaunch
!------------------------------------------------------------------------
!
!  Calculate cape of parcel distribution
!
!   - capp(it,nl) = (starting level)
!
!------------------------------------------------------------------------

    t_start = tstart(it,nl)
    rv_start = rvstart(it,nl)

    call cape(t_start,rv_start,it,ilnb,capp(it,nl))

!------------------------------------------------------------------------
! Store lnb
! when ilnb is zero, looks like caped is zero also.
!------------------------------------------------------------------------
!   print*,'ilnb = ',ilnb
    if (ilnb /= 0) then
      lnb(it,nl) = z(ilnb)
    else
!     print*,'ilnb capp(it,nl) = ',ilnb,capp(it,nl)
    endif
!------------------------------------------------------------------------
!
!  End loop over maxlev
!
!------------------------------------------------------------------------
end do
endif
end do

!------------------------------------------------------------------------
!
!    END  (subroutine get_cape_spectrum)
!
!------------------------------------------------------------------------
end subroutine get_cape_spectrum
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine rl_remove
!
!
!------------------------------------------------------------------------
subroutine rl_remove(tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g, ci, cl, latent, fusion, cpd, eps
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!  In/Out Variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: i
real(r8) :: dmdry, drl_uprain, hm_wat, hm_of_wat, drl, cinf, drl_evap
!------------------------------------------------------------------------
!
! 
!------------------------------------------------------------------------
  do i = 1,nd   
!------------------------------------------------------------------------
!  
!     REMOVING RL To Uprain or evap 
!  
!  There is no need to adjust hm here (for evap), since the hm of a grid 
!    level is left unchanged by this internal transformation.
!
!  This is a decision to REMOVE rl via uprain or evap
!  - rl_tar should be known at this point
!
!------------------------------------------------------------------------
   drl = rl(i) + rl_added(i) - rl_tar(i)
   if (drl > f_rl_rem*rl(i)) drl = f_rl_rem*rl(i)
   if ( drl > 0. ) then
!------------------------------------------------------------------------
!  Divide drl into drl_uprain and drl_evap
!------------------------------------------------------------------------
      if (rh(i) < rh_rl_evap_min) then
        drl_uprain = 0.
        drl_evap =  drl
      elseif (rh(i) > rh_rl_evap_max) then
        drl_uprain = drl
        drl_evap =  0.
      else
        cinf = (rh(i)-rh_rl_evap_min)/(rh_rl_evap_max-rh_rl_evap_min)
        drl_uprain = cinf*drl
        drl_evap = (1.-cinf)*drl
      endif
!     print*,'rl_added before rl evap and uprain prod = ',rl_added(i)
      rl_added(i) = rl_added(i) - drl_evap - drl_uprain
!     print*,'drl_evap = ',drl_evap
!     print*,'drl_uprain = ',drl_uprain
!     print*,'rl_added after rl evap and uprain prod = ',rl_added(i)
!------------------------------------------------------------------------
!  Calculate rl_evap
!------------------------------------------------------------------------
      dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
      rl_evap(i) = drl_evap*dmdry/tstep   ! No other source of this term
!------------------------------------------------------------------------
!  Calculate uprain_rl and hmuprain adjustments
!------------------------------------------------------------------------
      uprain_rl(i) = uprain_rl(i) + (drl_uprain*dmdry/tstep)  
      hm_of_wat = cl*t(i) + g*z(i)
      hm_wat = hm_of_wat*drl_uprain*dmdry/tstep  !  removal of hm from layer
      hmuprain_start = hmuprain_start + hm_wat
      hmuprain_prod(i) = hmuprain_prod(i) + hm_wat
!------------------------------------------------------------------------
!  A removal of hm from the layer (represented here as an entrainment).
!  - this is the same process as rv to ansnow. Removing the same amount
!    of hm from the layer, and adding the same amount to hmansnow. That
!    would leave a warming since resulting in a reduction of rv, rather
!    than ri here.
!------------------------------------------------------------------------
     hmflow(1,i) = hmflow(1,i) + hm_wat   ! actually updraft entrain
!------------------------------------------------------------------------
!  end if for rl exceeds rl_tar
!------------------------------------------------------------------------
    endif
!------------------------------------------------------------------------
!
!  End loop over heights
!
!------------------------------------------------------------------------
  end do
!------------------------------------------------------------------------
!
!    END  (subroutine rl_remove)
!
!------------------------------------------------------------------------
end subroutine rl_remove
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine ri_remove
!
!------------------------------------------------------------------------
subroutine ri_remove(tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g, ci, latent, fusion, cpd, eps, cl
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!  In/Out Variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: i
real(r8) :: dmdry, dri, dri_evap, dri_ansnow, hm_ice, hm_of_ice, cinf

!------------------------------------------------------------------------
!
!  This is a decision to REMOVE "excess" condensate to evap or ansnow
!  - ri_tar must be known
!
!------------------------------------------------------------------------
  do i = 1,nd   
    dri = (tstep/tscale_ri_rem(i))*(ri(i) + ri_added(i) - ri_tar(i))
    if (dri > 1.0E-12) then
!------------------------------------------------------------------------
!  To prevent negative ri
!------------------------------------------------------------------------
      if (dri < 0.) dri = 0.
      if (dri > f_ri_remove_max*ri(i)) then 
        dri = f_ri_remove_max*ri(i)
      endif
!------------------------------------------------------------------------
!  Divide dri into dri_ansnow and dri_evap
!------------------------------------------------------------------------
      if (rh(i) < rh_min(i)) then
        dri_ansnow = 0.
        dri_evap =  dri
      elseif (rh(i) > rh_max(i)) then
        dri_ansnow = dri
        dri_evap =  0.
      else
        cinf = (rh(i)-rh_min(i))/(rh_max(i)-rh_min(i))
        dri_ansnow = cinf*dri
        dri_evap = (1.-cinf)*dri
      endif
      ri_added(i) = ri_added(i) - dri_ansnow - dri_evap
!------------------------------------------------------------------------
!  Calculate ri_evap
!  There is no need to adjust hm here, since the hm of a grid level is
!    left unchanged by this internal transformation.
!------------------------------------------------------------------------
      dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
      ri_evap(i) = dri_evap*dmdry/tstep   ! No other source of this term
!------------------------------------------------------------------------
!  Add dri_ansnow to ansnow_ri (to represent ri sink)
!------------------------------------------------------------------------
       dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
       ansnow_ri(i) = ansnow_ri(i) + (dri_ansnow*dmdry/tstep)  
!------------------------------------------------------------------------
!  should use function for hm of ice here since used elsewhere
!------------------------------------------------------------------------
       hm_of_ice = cl*t(i) - fusion(t(i)) + g*z(i)
       hm_ice = hm_of_ice*dri_ansnow*dmdry/tstep  ! removal of hm from a layer

       if (hm_ice < 0.) then
          print*,'dri_ansnow = ',dri_ansnow
          print*,'t(i) = ',t(i)
          print*,'ri(i) = ',ri(i)
          print*,'dri = ',dri
          print*,'hm_ice = ',hm_ice
          print*,'cl*t(i) = ',cl*t(i)
          print*,'fusion(t(i)) = ',fusion(t(i))
          print*,'g*z(i) = ',g*z(i)
          print*,'z(i) = ',z(i)
          print*,'dri_ansnow = ',dri_ansnow
          stop 'hm_ice negative'
       endif
       hmansnow_start = hmansnow_start + hm_ice
       hmansnow_prod(i) = hmansnow_prod(i) + hm_ice
!------------------------------------------------------------------------
!  Removal of hm from the layer (represented here as an entrainment).
!  - this is the same process as rv to ansnow. Removing the same amount
!    of hm from the layer, and adding the same amount to hmansnow. That
!    would leave a warming since resulting in a reduction of rv, rather
!    than ri here.
!------------------------------------------------------------------------
       hmflow(1,i) = hmflow(1,i) + hm_ice   ! actually updraft entrain
    endif
!------------------------------------------------------------------------
!
!  end loop over heights
!
!------------------------------------------------------------------------
   end do
!------------------------------------------------------------------------
!
!    END  (subroutine ri_remove)
!
!------------------------------------------------------------------------
end subroutine ri_remove
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine rv_to_ansnow_or_ri
!
!  - better to invoke after updrafts since could then enforce a 
!    more sophisticated attempt to avoid negative rvnew.
!
!------------------------------------------------------------------------
subroutine rv_to_ansnow_or_ri(tstep,this_lat)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g, ci, latent, fusion, cpd, eps, cl
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!  In/Out Variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: tstep
real(r8), intent(in) :: this_lat
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: i,jj
real(r8) :: dmdry,rh_fraction,drv_ansnow,drv_ansnow_max,dri_prod,hm_ice,xx,pp
real(r8) :: ee,rv_want,dri_prod_max,hm_of_ice
real(r8) :: drv_remove(nd),rh_remove,cinf,ri_need(nd)

!------------------------------------------------------------------------
!
!  Calculate preliminary drv_remove profile
!  - based on local RH exceeding prescribed dvalue
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  Calculate preliminary drv_remove profile
!  - based on local RH exceeding prescribed value
!  - drv_remove is later allocated into parts going to ansnow and another
!    part going to ri
!
!------------------------------------------------------------------------
do i = 1,nd
  if (i >= i_strat) then
!   print*,'----------' 
    if (abs(this_lat) < lat_rh_wat) then
      rv_want = rh_max(i)*rs(i)
!     print*,'using ice' 
    else
      rv_want = rh_max(i)*rs_wat(i)
!     print*,'rs(i) = ',rs(i)
!     print*,'rs_wat(i) = ',rs_wat(i)
!     print*,'using water rv_want = ',rv_want
    endif
!   print*,'this_lat = ',this_lat
!   print*,'i t(i) = ',i,t(i)
!   print*,'rs(i) rs_wat(i) = ',rs(i),rs_wat(i)
!   rv_want = rh_max(i)*rs(i)
    if (rv(i) > rv_want) then
      drv_remove(i) = rv(i) - rv_want
    else
      drv_remove(i) = 0.
    endif
  else
    drv_remove(i) = 0.
  endif
end do

!------------------------------------------------------------------------
!
!  Add to drv_remove an amount needed to bring ri + ri_added to ri_tar.
!  - Added to ri later.
!  - Never remove more than rv_rem_max fraction
!
!------------------------------------------------------------------------
do i = 1,nd
  ri_need(i) = ri_tar(i) - ri(i) - ri_added(i)
  if (ri_need(i) > 0.) then
!   drv_remove(i) = ri_need(i)
    drv_remove(i) = drv_remove(i) + ri_need(i)
!   if (drv_remove(i) < 0.) drv_remove(i) = 0.
  endif
  if (drv_remove(i) > rv_rem_max*rv(i)) drv_remove(i) = rv_rem_max*rv(i)
  dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
end do

!print*,'-----------------'
!print*,'Numbers here in mm per day'
!print*,'ansnow_strat_target = ',ansnow_strat_target*3600.*24.
!print*,'change from prev in mm per day ansnow_strat_diff = ',ansnow_strat_diff*3600.*24.
!print*,'-----------------'

!------------------------------------------------------------------------
!
!  - to avoid negative rvnew
!
!------------------------------------------------------------------------
do i = 1,nd
  if (drv_remove(i) > rv_fraction_remove*rv(i)) then
    drv_remove(i) = rv_fraction_remove*rv(i)
  endif
end do

!------------------------------------------------------------------------
!
!  If drv_remove ends up being positive
!  - partition into ansnow and ri prod
!
!------------------------------------------------------------------------
do i = 1,nd
  if (drv_remove(i) > 0.) then
!------------------------------------------------------------------------
!  Partition drv_remove
!
!  If ri + ri_added exceeds target, allocate all drv to ansnow
!------------------------------------------------------------------------
     if (ri_need(i) < 0.)  then 
       drv_ansnow = drv_remove(i)
       dri_prod = 0.
     else
!------------------------------------------------------------------------
!  Otherwise, allocate drv_remove to dri_prod as long as condensate doesn't
!    exceed ri_tar: however detrainment also trying to keep ri at max value.
!------------------------------------------------------------------------
       if (drv_remove(i) > ri_need(i)) then
         dri_prod = ri_need(i)
         drv_ansnow = drv_remove(i) - ri_need(i)
       else
         dri_prod = drv_remove(i)
         drv_ansnow = 0.
       endif
       if (dri_prod < 0.) dri_prod = 0.
       if (drv_ansnow < 0.) drv_ansnow = 0.
     endif
!------------------------------------------------------------------------
!  Add drv_ansnow to ansnow_rv
!------------------------------------------------------------------------
      dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
      ansnow_rv(i) = ansnow_rv(i) + (drv_ansnow*dmdry/tstep)
      ri_prod(i) = ri_prod(i) + (dri_prod*dmdry/tstep)
      ri_added(i) = ri_added(i) + dri_prod
!------------------------------------------------------------------------
!  should use function for hm of ice here (but maybe not used elsewhere?)
!------------------------------------------------------------------------
      hm_of_ice = cl*t(i) - fusion(t(i)) + g*z(i)
      hm_ice = hm_of_ice*drv_ansnow*dmdry/tstep ! removal of hm from a layer

      if (hm_ice < 0.) then
        print*,'hm_ice = ',hm_ice
        print*,'cl*t(i) = ',cl*t(i)
        print*,'fusion(t(i)) = ',fusion(t(i))
        print*,'g*z(i) = ',g*z(i)
        print*,'t(i) = ',t(i)
        print*,'hm_of_ice = ',hm_of_ice
        print*,'drv_ansnow = ',drv_ansnow
        stop 'hm_ice negative rv_to_ansnow_or_ri'
      endif

      hmansnow_start = hmansnow_start + hm_ice
      hmansnow_prod(i) = hmansnow_prod(i) + hm_ice
!------------------------------------------------------------------------
!  A removal of hm from the layer (represented here as an entrainment).
!  - note that hm_ice < 0 typically, so removal of negative is a source.
!------------------------------------------------------------------------
      hmflow(1,i) = hmflow(1,i) + hm_ice   ! actual updraft entrain
!------------------------------------------------------------------------
!  - ri_prod is just internal redistribution of hm in the layer, so no
!    need to calculate hm change
!------------------------------------------------------------------------
    endif
!------------------------------------------------------------------------
!
!  End loop over heights
!
!------------------------------------------------------------------------
  end do

!------------------------------------------------------------------------
!
!    END  (subroutine rv_to_ansnow_or_ri)
!
!------------------------------------------------------------------------
end subroutine rv_to_ansnow_or_ri
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine rv_to_uprain
!  - drizzle
!
!------------------------------------------------------------------------
subroutine rv_to_uprain(tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g,ci,latent,fusion,cpd,eps,cl,tkelvin
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!  In/Out Variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: i,jj
real(r8) :: dmdry,dhm_up,drv_uprain,hm_wat
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
do i = 1,nd
!------------------------------------------------------------------------
!
!  Calculate drv_uprain from Boundary Layer
!
!------------------------------------------------------------------------
  drv_uprain = 0.
  if (p(i) > press_rv_to_uprain) then
  if (rv(i) > rh_rv_to_uprain*rs(i)) then
    drv_uprain = (tstep/tscale_rv_to_uprain)*(rv(i) - rh_rv_to_uprain*rs(i))
  endif
  endif
!------------------------------------------------------------------------
!  Define dhm_up: hm going to hmuprain from both processes
!------------------------------------------------------------------------
   dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
   uprain_rv(i) = uprain_rv(i) + (drv_uprain*dmdry/tstep)
   hm_wat = cl*t(i) + g*z(i)
   dhm_up = drv_uprain*hm_wat
!------------------------------------------------------------------------
! hm adjustments
! hmflow: is for updraft entrainment, as a way of representing removal
!   of hm from the layer
!------------------------------------------------------------------------
   hmuprain_start = hmuprain_start + (dhm_up*dmdry/tstep)
   hmflow(1,i) = hmflow(1,i) + (dhm_up*dmdry/tstep)

   hmuprain_prod(i) = hmuprain_prod(i) + (dhm_up*dmdry/tstep)
!------------------------------------------------------------------------
!
!  End loop over heights
!
!------------------------------------------------------------------------
end do
!------------------------------------------------------------------------
!
!    END  (subroutine rv_to_uprain)
!
!------------------------------------------------------------------------
end subroutine rv_to_uprain
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine rv_to_rl
!
!------------------------------------------------------------------------
subroutine rv_to_rl(tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!  In/Out Variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
integer :: i
real(r8) :: drl_diff, dmdry
!------------------------------------------------------------------------
!
!    PRODUCING RL FROM RV
!
!  Calculate rl_prod
!  - no need to calculate hmflow here since no change in hm of layer
!
!------------------------------------------------------------------------
do i = 1,nd
  drl_diff = rl_tar(i) - rl(i) - rl_added(i)
  if (drl_diff > 0.) then
    dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
    rl_prod(i) = rl_prod(i) + (drl_diff*dmdry/tstep)
    rl_added(i) = rl_added(i) + drl_diff
  endif
end do
!------------------------------------------------------------------------
!
!    END  (subroutine rv_to_rl)
!
!------------------------------------------------------------------------
end subroutine rv_to_rl
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   subroutine define_uprain_start
!
!------------------------------------------------------------------------
subroutine define_uprain_start
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g,latent,cl
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
integer :: i
real(r8) :: t_uprain
!------------------------------------------------------------------------
!
!   Define uprain
!
!   units of uprain prod: kg water/m2*s
!   units of rain (here): kg water/m2*s
!
!  uprain prod: defined at full levels.
!  uprain: defined at half levels.
!
!  I don't evaporate uprain_rv
!
!------------------------------------------------------------------------
  uprain_surf_rv = 0.
  uprain_start_rlp = 0.
  do i = 1,nd
    uprain_surf_rv = uprain_surf_rv + uprain_rv(i)
    uprain_start_rlp = uprain_start_rlp + uprain_rlp(i) + uprain_rl(i)
  end do
!------------------------------------------------------------------------
!
!  Define have_precip
!  - this number must be small for conservation reasons. uprain is very 
!    unlikely to be ever this small anyway, due to threshold on the mass 
!    fraction
!
!------------------------------------------------------------------------
  if (uprain_surf_rv + uprain_start_rlp > 0.000000000001) have_precip = 1
!------------------------------------------------------------------------
!
!    END  (subroutine define_uprain_start)
!
!------------------------------------------------------------------------
end subroutine define_uprain_start
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------

!
!   subroutine define_ansnow
!
!------------------------------------------------------------------------
subroutine define_ansnow
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : 
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
integer :: i

!------------------------------------------------------------------------
!
!  It is mainly for diagnostic reasons that retain two distinct ansnow
!
!------------------------------------------------------------------------
  ansnow_conv = 0.
  ansnow_strat = 0.
  ansnow_strat_rv = 0.
  ansnow_strat_ri = 0.
  do i = 1,nd
    ansnow_conv = ansnow_conv + ansnow_rlp(i)
    ansnow_strat = ansnow_strat + ansnow_rv(i) + ansnow_ri(i)
    ansnow_strat_rv = ansnow_strat_rv + ansnow_rv(i) 
    ansnow_strat_ri = ansnow_strat_ri + ansnow_ri(i) 
    if ((ansnow_conv < 0.).or.   &
        (ansnow_strat < 0.)) then
      print*,'========= NEGATIVE ansnow ======='
      print*,'ansnow_conv = ',ansnow_conv
      print*,'ansnow_strat = ',ansnow_strat
      print*,'ansnow_rv(i) = ',ansnow_rv(i)
      print*,'ansnow_ri(i) = ',ansnow_ri(i)
      print*,'ansnow_rlp(i) = ',ansnow_rlp(i)
      stop 'negative precip'
    endif
  end do
  ansnow_start = ansnow_conv + ansnow_strat  

!------------------------------------------------------------------------
!
!  Define have_precip
!  - maybe need very small number here or else have problems with dperr in
!    conv_tendencies
!
!------------------------------------------------------------------------
  if (ansnow_conv > 1.0E-14 ) have_precip = 1
  if (ansnow_strat > 1.0E-14 ) have_precip = 1
!------------------------------------------------------------------------
!
!    END  (subroutine define_ansnow)
!
!------------------------------------------------------------------------
end subroutine define_ansnow
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   subroutine evap
!
!------------------------------------------------------------------------
subroutine evap(tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : cl,ci,cpd,fusion,g,latent,epsi,rsat,t_from_km
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!  In/Out Variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local integer variables
!
!------------------------------------------------------------------------
integer :: i, i_top, ii, nprint, nr, ngot
integer :: do_melt
integer :: nn_want, nn
!------------------------------------------------------------------------
!
!   Local real variables
!
!------------------------------------------------------------------------
real(r8) :: uu,cinf
real(r8) :: ansnow, anrain
real(r8) :: hmansnow, hmanrain
real(r8) :: hmuprain
real(r8) :: dmdry, dtt, drv, drvp, drh
real(r8) :: mass_ansnow
real(r8) :: km_wat, km_ice, hm_of_wat, hm_of_ice
real(r8) :: precip_mm_day, cc, small
real(r8) :: tot_prod, rh_factor
real(r8) :: mass_uprain_evap
real(r8) :: mass_anrain_evap
real(r8) :: mass_ansnow_subl
real(r8) :: ansnow_prod(nd)
real(r8) :: error_up, error_anrain, error_ansnow
real(r8) :: slope, pp, p1, p2, uprain_from_rlp, lv
real(r8) :: rh_diff
real(r8) :: rh1, error_small, rh_guess, error
real(r8) :: mass_ratio, f_evap, drv_evap, drv_subl
real(r8) :: drh_uprain
real(r8) :: drh_anrain
real(r8) :: fraction_remove, f_evap_down_max, rain_down_start
real(r8) :: hm_down, km_down, rl_down, ri_down, drh_down, rate_down, f_mass
real(r8) :: td_down, t_down, rain_down, dHm_precip, b_down, rh_down, rs_down
real(r8) :: mass_water_evap, b_lower, km_lower, t_lower, td_lower
real(r8) :: td_low, b_low, km_low, t_low, rs_grid, rh_grid, fff
real(r8) :: b_av, b_prev, water_error, hm_error,drv_tot,rv_error,f_drv
real(r8) :: rv_down, drh_add, mass_down,drv_down,drv_max,mass_down_old
real(r8) :: rv_mass_avail, fraction_evap, hmansnow_old, hmanrain_old
integer :: do_down, do_down_uprain, do_down_anrain, i_put, iter
integer :: have_down, nprint_stuff, i_low, i_down, i_lower, i_lowest
real(r8) :: rv_gain_anrain_down, rv_gain_uprain_down
real(r8) :: rv_gain_anrain_melt, rv_gain_ansnow_subl
real(r8) :: rv_gain_uprain_evap, rv_gain_anrain_evap
real(r8) :: tot_rv_gain,precip_loss


!------------------------------------------------------------------------
!
!  *************  BEFORE ALL LOOPS  *************************
!
!  Diagnostics for water conservation
!
!------------------------------------------------------------------------
  rv_gain_anrain_down = 0.
  rv_gain_uprain_down = 0.
  rv_gain_anrain_melt = 0.
  rv_gain_ansnow_subl = 0.
  rv_gain_uprain_evap = 0.
  rv_gain_anrain_evap = 0.

!------------------------------------------------------------------------
!
!  *************  BEFORE ALL LOOPS  *************************
!
!  Determine whether to do detailed printing
!
!------------------------------------------------------------------------
if (nprint_evap == 1) then 
  nprint = 1
  print*,'ENTERING EVAP: Setting nprint equal 1'
else
  nprint = 0
endif

aaa = ansnow_strat*3600.*24.
if (aaa > 400.0) then
  print*,'ansnow_strat n mm per day = ',aaa
  print*,'Turning on printing'
  nprint = 1
endif

!print*,'========  ENTERING EVAP for new timestep ============'

!------------------------------------------------------------------------
!
!  *************  BEFORE ALL LOOPS  *************************
!
!  Define total uprain/ansnow prod
!  - purely local array
!
!------------------------------------------------------------------------
do i = 1,nd
  ansnow_prod(i) = ansnow_rlp(i) + ansnow_rv(i) + ansnow_ri(i)
end do

!------------------------------------------------------------------------
!
!  *************  BEFORE ALL LOOPS  *************************
!
!   Determine i_top: 
!   - top level at which one of precip sources is non-zero.
!   - may get i_top = -1 here, in which case, don't do loop below.
!
!------------------------------------------------------------------------
i_top = -1
small = 1.0E-14

do i = 1,nd
  tot_prod = uprain_rlp(i) + ansnow_prod(i)
  if ((tot_prod > small).and.(z(i) < z_evap_max)) i_top = i
end do

!------------------------------------------------------------------------
!
!  *************  BEFORE ALL LOOPS  *************************
!
!  Initialize variables:
!  uprain_from_rlp, ansnow, and anrain are continously updated.
!  Current values at that level.
!
!------------------------------------------------------------------------
uprain_from_rlp = 0.
ansnow = 0.
anrain = 0.

hmanrain = 0.
hmansnow = 0.
hmuprain = 0.

if (i_top /= -1) then
  do i = i_top+1,nd
    uprain_from_rlp = uprain_from_rlp + uprain_rlp(i) + uprain_rl(i)
    ansnow = ansnow + ansnow_prod(i)
    hmuprain = hmuprain + hmuprain_prod(i)
    hmansnow = hmansnow + hmansnow_prod(i)
  end do
else
  uprain_from_rlp = uprain_start_rlp
  ansnow = ansnow_conv + ansnow_strat
  anrain = 0.
  hmuprain = hmuprain_start
  hmansnow = hmansnow_start
  hmanrain = 0.
endif

!------------------------------------------------------------------------
!
!  *************  BEFORE ALL LOOPS  *************************
!
!   Diagnostics
!
!------------------------------------------------------------------------
if (nprint == 1) then
  print*,'-- AT START OF EVAP --'
  print*,'i_top = ',i_top
  print*,'starting hmansnow (above i_top) = ',hmansnow
  print*,'starting ansnow (above i_top) = ',ansnow
endif

!------------------------------------------------------------------------
!
!  *************  BEFORE ALL LOOPS  *************************
!
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!
!  *************  START LOOP OVER HEIGHTS *************************
!
!------------------------------------------------------------------------
if (i_top >= 1) then
do i = i_top,1,-1

!print*,'======= AT TOP OF STARTING ANOTHER LEVEL i = ',i

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!
!------------------------------------------------------------------------
mass_uprain_evap = 0.
mass_anrain_evap = 0.
mass_ansnow_subl = 0.

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  Subroutine evap
!  - update precp with local contribution
!
!  Note on precipitation budget:
!  Before all loops, ansnow is initialized as the total production above
!  i_top. At the start of every loop (here below), I add the local production.
!  Within a loop, I subtract the melting and sublimation. So the variable
!  ansnow is continuously updated throughout all loops (and similarly for
!  other rain variables). This is a bit odd. Suppose that I have
!  initiated a downdraft at level i. Then I iterate evaporation, so that the
!  parcel may descend from level i to the surface. At a level below i at which
!  evaporation occurs, the only rain available is the rain that started at
!  level i, minus all evaporation that occurred after that. I do not include
!  the production that occurs between the evaporation level and level i.
!  This rain would of course be made available later as I loop downward in i.
!
!------------------------------------------------------------------------
uprain_from_rlp = uprain_from_rlp + uprain_rlp(i) + uprain_rl(i)
uprain_rlp_rvv = uprain_rlp_rvv + rlp_rvv(i)
ansnow = ansnow + ansnow_prod(i)

!print*,'New ansnow after ansnow_prod added = ',ansnow
!print*,'ansnow_prod(i) = ',ansnow_prod(i)

hmuprain = hmuprain + hmuprain_prod(i)
hmansnow = hmansnow + hmansnow_prod(i)

if (nprint == 1) then
  print*,'adjusted hmansnow from hmansnow_prod = ',hmansnow
  print*,'t(i) = ',t(i)
  print*,'i hmansnow_prod(i) = ',i,hmansnow_prod(i)
endif


!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  ********   Start Melting Downdrafts  **************
!
!  Should likely enter around T equal 260
!  - immediately determine mass_down from ansnow
!  - start evaporating ansnow and testing for negative buoyancy
!  - when T goes below 270 or so start melting
!
!
!  Determine do_melt
!
!  MMM
!
!------------------------------------------------------------------------
do_melt = 1

if (ansnow < 1.0E-12) do_melt = 0
if (switch_melt == 0) do_melt = 0
if (t(i) < t_melt) do_melt = 0

if (nprint == 1) then
  print*,'i do_melt = ',i,do_melt
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  If do_melt = 1:
!
!------------------------------------------------------------------------
if (do_melt == 1) then

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Needed thermo
!
!  I am using the following definitions for the specific enthalpy of
!    water and ice (J/kg):
!
!  km_wat = cl*t(i)                  [J/kg wat]
!  km_ice = km_wat - fusion(t(i))    [J/kg ice]
!
!------------------------------------------------------------------------

km_wat = cl*t(i) 
km_ice = km_wat - fusion(t(i))
hm_of_wat = km_wat + g*z(i)
hm_of_ice = km_ice + g*z(i)

if (km_ice < 0.) then
  print*,'WARNING: km_ice negative'
  print*,'km_ice in do_melt = ',km_ice
  print*,'t(i) = ',t(i)
  print*,'fusion(t(i)) = ',fusion(t(i))
! stop
endif

!print*,'*****  ENTERING do_melt = 1 for starting level i = ',i
!print*,'starting ansnow = ',ansnow
!print*,'starting anrain = ',anrain
!print*,'t(i) = ',t(i)
!print*,'-----'

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Calculate mass_ansnow
!
!------------------------------------------------------------------------
mass_ansnow = tstep*ansnow   ! to get [kg snow/m2]

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Set ansnow equal zero and add to anrain.
!
!------------------------------------------------------------------------
anrain = anrain + (mass_ansnow/tstep)
ansnow = ansnow - (mass_ansnow/tstep)
ansnowmelt(i) = mass_ansnow/tstep

!print*,'DOING MELTING: anrain = ',anrain
!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Adjust hmansnow, hmanrain
!
!  - The net increase in the precipitation MSE (i.e. hmansnow+hmanrain)
!    must be balanced by removal of MSE from the background atm.
!
!   To get hmansnow [J/m2]:
!   hmansnow = mass_ansnow*(km_ice + g*z(i))
!            = [kg ice/m2]*[J/kg ice]
!            = [J/m2]
!
!  - Assume that all of the ansnow is melted to anrain. Assume also that
!    when the ansnow is melted it has a temperature equal to t(i).
!  - This forces you to write:
!    hmanrain = hmanrain_old + (mass_ansnow/tstep)*hm_of_wat
!  - Then the amout of hm subtracted from the downdraft parcel becomes
!    determined by Hm conservation, i.e. equal to the net change in the
!    Hm of anrain + ansnow.
!  - Change in HM of precip
!   tstep*(hmanrain_new + hmansnow_new - hmanrain_old - hmansnow_old)   
!   = tstep*[hmanrain_old + (mass_ansnow/tstep)*hm_of_wat - 
!            hmanrain_old - hmansnow_old]
!   = mass_ansnow*hm_of_wat - hmansnow_old*tstep
!   dHm_precip = mass_ansnow*hm_of_wat - hmansnow_old*tstep
!
!------------------------------------------------------------------------
hmansnow_old = hmansnow
hmanrain_old = hmanrain
hmanrain = hmanrain_old + (mass_ansnow/tstep)*hm_of_wat
hmansnow = 0.

dHm_precip = mass_ansnow*hm_of_wat - hmansnow_old*tstep

if (nprint == 1) then
  print*,'---'
  print*,'adjusted hmansnow = ',hmansnow
  print*,'hm_of_ice = ',hm_of_ice
  print*,'km_ice = ',km_ice
  print*,'km_wat = ',km_wat
  print*,'fusion(t(i)) = ',fusion(t(i))
  print*,'g*z(i) = ',g*z(i)
  print*,'Added amount = ',-(mass_ansnow/tstep)*hm_of_ice
  print*,'mass_ansnow = ',mass_ansnow
  print*,'---'
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Calculate mass_down
!
!  mass_down =  mass_melt_ratio*mass_ansnow
!  Suppose mass_melt_ratio = 100
!  lf/cpd = 3.3E+05/1004 = 330000/1000 = 330 
!  Then cooling would be 3.3 degrees
!
!------------------------------------------------------------------------
mass_down =  mass_melt_ratio*mass_ansnow
dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))
if ((mass_down/dmdry) > 0.2) mass_down = 0.2*dmdry

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Initialize parcel properties before evap loop
!
!  What is hm_down?
!  Change in Hm of parcel:
!  -dHm_precip = mass_down*(hm_down - hm(i))
!  -dHm_precip/mass_down = hm_down - hm(i)
!  hm_down = hm(i) - (dHm_precip/mass_down) 
!
!------------------------------------------------------------------------
hm_down = hm(i) - (dHm_precip/mass_down)
rv_down = rv(i)
km_down = hm_down - (1. + rv_down)*g*z(i)
get_t_call = 1
ri_down = 0.
rl_down = 0.
call t_from_km(km_down,p(i),rv_down,ri_down,rl_down,t_down,get_t_call)
td_down = t_down*(1. + epsi*rv_down)/(1. + rv_down)
b_down = g*(td_down - td(i))/td(i)
rs_down = rsat(t_down,p(i))
rh_down = rv_down/rs_down
i_put = i

!print*,'---------  in MELTING after melting but before evap loop ---------'
!print*,'i z(i) = ',i,z(i)
!print*,'t_down = ',t_down
!print*,'t(i) = ',t(i)
!print*,'b_down = ',b_down
!print*,'rh_down = ',rh_down
!print*,'rh(i) = ',rh(i)
!print*,'-------  finished parcel properties before melting evap loop '

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Start evap loop
!  - loop over iter_evap
!  - get drv_evap/mass_anrain_evap
!  - find buoyancy at next lowest level: b_below
!  - move down if b_av < 0.
!  - end loop
!  - detrain at i_put
!  MMM
!
!  Find mass_anrain_evap and update anrain
!
!------------------------------------------------------------------------
do iter = 1,iter_evap
if (i_put > 1) then

!print*,'STARTING iter = ',iter

rh_diff = rh_melt_target - rh_down
!print*,'rh_diff = ',rh_diff
if (rh_diff > 0.) then 
  if (anrain > 1.0E-10) then
!    print*,'DOING anrain evap'
    drv_evap = rh_melt_target*rs_down - rv_down
    drv_subl = 0.
    mass_anrain_evap = drv_evap*mass_down
    mass_ansnow_subl = 0.
    f_evap = mass_anrain_evap/(tstep*anrain)
    if (f_evap > f_evap_melt_max) then 
      mass_anrain_evap = f_evap_melt_max*tstep*anrain
      drv_evap = mass_anrain_evap/mass_down
    endif
  elseif (ansnow > 1.0E-10) then
!   print*,'DOING ansnow subl'
    drv_evap = 0.
    drv_subl = rh_melt_target*rs_down - rv_down
    mass_anrain_evap = 0.
    mass_ansnow_subl = drv_subl*mass_down
    f_evap = mass_ansnow_subl/(tstep*ansnow)
    if (f_evap > f_evap_melt_max) then 
      mass_ansnow_subl = f_evap_melt_max*tstep*ansnow
      drv_subl = mass_ansnow_subl/mass_down
    endif
  else
!    print*,'NOT DOING evap or subl'
    mass_anrain_evap = 0.
    mass_ansnow_subl = 0.
    drv_evap = 0.
    drv_subl = 0.
  endif
else
  mass_anrain_evap = 0.
  mass_ansnow_subl = 0.
  drv_evap = 0.
  drv_subl = 0.
endif

anrain = anrain - (mass_anrain_evap/tstep)
ansnow = ansnow - (mass_ansnow_subl/tstep)
anrainevap(i_put) = anrainevap(i_put) + (mass_anrain_evap/tstep)
ansnowsubl(i_put) = ansnowsubl(i_put) + (mass_ansnow_subl/tstep)
hmanrain = hmanrain - (mass_anrain_evap/tstep)*hm_of_wat
hmansnow = hmansnow - (mass_ansnow_subl/tstep)*hm_of_ice

!print*,'--- adjusting precip for evap or subl ----'
!print*,'drv_subl = ',drv_subl
!print*,'drv_evap = ',drv_evap
!print*,'new anrain = ',anrain
!print*,'new ansnow = ',ansnow
!print*,'----'


!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Re-determine local parcel properties with evap
!
!------------------------------------------------------------------------
rv_down = rv_down + drv_evap + drv_subl
hm_down = hm_down + drv_evap*hm_of_wat + drv_subl*hm_of_ice
km_down = hm_down - (1. + rv_down)*g*z(i_put)
get_t_call = 1
ri_down = 0.
rl_down = 0.
call t_from_km(km_down,p(i_put),rv_down,ri_down,rl_down,t_down,get_t_call)
td_down = t_down*(1. + epsi*rv_down)/(1. + rv_down)
b_down = g*(td_down - td(i_put))/td(i_put)
rs_down = rsat(t_down,p(i_put))
rh_down = rv_down/rs_down

!print*,'local b_down = ',b_down
!print*,'td_down = ',td_down
!print*,'td(i_put) = ',td(i_put)
!print*,'t_down = ',t_down
!print*,'t(i_put) = ',t(i_put)
!print*,'rv_down = ',rv_down
!print*,'rv(i_put) = ',rv(i_put)
!print*,'i_put = ',i_put
!print*,'i = ',i
!print*,'-----'

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Find parcel properties at lower level
!
!------------------------------------------------------------------------
i_low = i_put - 1
if (i_low == 0) stop 'i_low is zero'
km_low = hm_down - (1. + rv_down)*g*z(i_low)
get_t_call = 1
ri_down = 0.
rl_down = 0.
call t_from_km(km_low,p(i_low),rv_down,ri_down,rl_down,t_low,get_t_call)
td_low = t_low*(1. + epsi*rv_down)/(1. + rv_down)
b_low = g*(td_low - td(i_low))/td(i_low)
b_av = 0.5*(b_low + b_down) + b_precip_melt

!print*,'b_low = ',b_low
!print*,'b_av = ',b_av

if (b_av < 0.) then 
  i_put = i_low
  t_down = t_low
  rs_down = rsat(t_down,p(i_put))
  rh_down = rv_down/rs_down
!  print*,'Moving to lower level i_put = ',i_put
endif


!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE MELTING DOWNDRAFTS ********************
!
!  end loop over iter_evap
!
!------------------------------------------------------------------------
endif
end do

!print*,'---  ENDING iter loop ---'

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE MELTING DOWNDRAFTS ********************
!
!  Define fmass_dn
!  - want to consider fmass_dn to be positive at half level heights between
!    i and i_put.
!  - Suppose i_put = i-1, then the half level height between is i
!
!------------------------------------------------------------------------

if (i_put < i) then
  do ii = i_put+1,i
    fmass_dn(ii) = fmass_dn(ii) + (mass_down/tstep)
  end do
endif

!------------------------------------------------------------------------
!
!  Melting Downdraft Entrainment/Detrainment
!
!------------------------------------------------------------------------

  dmflow(1,i) = dmflow(1,i) + (mass_down/tstep)
  rvflow(1,i) = rvflow(1,i) + (rv(i)*mass_down/tstep)
  hmflow(1,i) = hmflow(1,i) + (hm(i)*mass_down/tstep)
  uflow(1,i) = uflow(1,i) + (uwind(i)*mass_down/tstep)
  vflow(1,i) = vflow(1,i) + (vwind(i)*mass_down/tstep)

!  detrain

  dmflow(2,i_put) = dmflow(2,i_put) + (mass_down/tstep)
  rvflow(2,i_put) = rvflow(2,i_put) + (rv_down*mass_down/tstep)
  hmflow(2,i_put) = hmflow(2,i_put) + (hm_down*mass_down/tstep) 
  uflow(2,i_put) = uflow(2,i_put) + (uwind(i)*mass_down/tstep)
  vflow(2,i_put) = vflow(2,i_put) + (vwind(i)*mass_down/tstep)

! Melting diagnostics
  mtent(i) = mtent(i) + (mass_down/tstep)
  mtdet(i_put) = mtdet(i_put) + (mass_down/tstep)

! This actually includes ansnow subl and anrain evap
  rv_gain_anrain_melt = rv_gain_anrain_melt + (rv_down-rv(i))*mass_down

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE do_melt = 1 *************************
!
!  Stage 1: ansnow melting 
!
!  Diagnostics
!
!------------------------------------------------------------------------

if (nprint == 1) then
  print*,'-- IN do_melt = 1 i = ',i 
  print*,'mass_ansnow = ',mass_ansnow
  print*,'ansnow_start = ',ansnow_start*3600.*24
  print*,'ansnow = ',ansnow*3600.*24
  print*,'anrain = ',anrain*3600.*24
  if (ansnow_start > small) then
    print*,'Fraction of original ansnow still frozen = ',ansnow/ansnow_start
  endif
  print*,'hmanrain = ',hmanrain
  print*,'hmansnow = ',hmansnow
endif

!------------------------------------------------------------------------
! endif for do_melt = 1
!------------------------------------------------------------------------
endif

!------------------------------------------------------------------------
!
!     ---------  END MELTING DOWNDRAFTS ----------
!
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  Stage 2: uprain/ansnow/anrain evaporation
!
!  Calculate mass_evap
!
!  - mass_evap should have units: kg water/m2
!
!------------------------------------------------------------------------
dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rv(i)+rl(i)+ri(i)))

if (rh(i) < rh_uprain ) then

  rh_diff = rh_uprain - rh(i)

!  mass_uprain_evap = rs(i)*rh_diff*dmdry*tstep*uprain_from_rlp*rate_evap_uprain

if (i == 1) then
  mass_uprain_evap = rs(i)*rh_diff*dmdry*tstep*uprain_from_rlp*rate_evap_surf
else
  mass_uprain_evap = rs(i)*rh_diff*dmdry*tstep*uprain_from_rlp*rate_evap_uprain
endif

  rv_mass_avail = uprain_from_rlp*tstep

  if (rv_mass_avail > 0.000001) then
    fraction_evap = mass_uprain_evap/rv_mass_avail
    if (fraction_evap > max_evap) then 
      fraction_evap = max_evap
      mass_uprain_evap = max_evap*rv_mass_avail
    endif
  endif

  nprint_stuff = 0
  if (nprint_stuff == 1) then
! if ((p(i) < 500.*100.).and.(rh(i) < 0.4)) then
    print*,'----------'
    print*,'Evaporating uprain'
    print*,'z(i) = ',z(i)
    print*,'fraction_evap = ',fraction_evap
    print*,'rh(i) = ',rh(i)
    print*,'rs(i) = ',rs(i)
    print*,'rv(i) = ',rv(i)
    print*,'rate_evap_uprain = ',rate_evap_uprain
    print*,'uprain_from_rlp = ',uprain_from_rlp
    print*,'rv mass in layer = ',rv(i)*dmdry
    print*,'rv_mass_avail = ',rv_mass_avail
    print*,'water evaporated mass_uprain_evap = ',mass_uprain_evap
    print*,'fractional increase in water = ',mass_uprain_evap/(rv(i)*dmdry)
    print*,'----------'
! endif
  endif
else
  mass_uprain_evap = 0.
endif

if (nprint == 1) then
 print*,'--- FIRST CALC ---'
 print*,'dmdry = ',dmdry
 print*,'rs(i) = ',rs(i)
 print*,'mass_uprain_evap = ',mass_uprain_evap
 print*,'uprain_from_rlp = ',uprain_from_rlp
 print*,'rate_evap_uprain = ',rate_evap_uprain
 print*,'Expected reduction in uprain_from_rlp = ',mass_uprain_evap/tstep
 print*,'rh(i) = ',rh(i)
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  Stage 2: uprain/ansnow/anrain evaporation
!
!  anrain/ansnow evaporation
!
!
!------------------------------------------------------------------------
if (t(i) > tkelvin) then
  if (rh(i) < rh_ansnow ) then

    rh_diff = rh_ansnow - rh(i)

    if (i == 1) then
      mass_anrain_evap = rs(i)*rh_diff*dmdry*tstep*anrain*rate_evap_surf
    else
      mass_anrain_evap = rs(i)*rh_diff*dmdry*tstep*anrain*rate_evap_anrain
    endif

    rv_mass_avail = anrain*tstep

    if (rv_mass_avail > 0.000001) then
      fraction_evap = mass_anrain_evap/rv_mass_avail
      if (fraction_evap > max_evap) then 
        mass_anrain_evap = max_evap*rv_mass_avail
      endif
    endif
  else
    mass_anrain_evap = 0.
  endif
  mass_ansnow_subl = 0.
else
  if ((rh(i) < rh_ansnow ).and.(ansnow > 1.0E-10)) then

    rh_diff = rh_ansnow - rh(i)
    mass_ansnow_subl = rs(i)*rh_diff*dmdry*tstep*ansnow*rate_evap_ansnow
    rv_mass_avail = ansnow*tstep

    fraction_evap = mass_ansnow_subl/rv_mass_avail
    if (fraction_evap > max_evap) then 
      mass_ansnow_subl = max_evap*rv_mass_avail
    endif
  else
    mass_ansnow_subl = 0.
  endif
  mass_anrain_evap = 0.
endif

if (mass_anrain_evap < -0.0000001) then
  print*,'anrain = ',anrain
  print*,'rh(i) = ',rh(i)
  stop 'stop mass_anrain_evap negative in evap'
endif

if (mass_ansnow_subl < -0.0000001) stop 'stop mass_ansnow_subl negative in evap'
if (mass_uprain_evap < -0.0000001) stop 'stop mass_uprain_evap negative in evap'

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  Stage 2: uprain/ansnow/anrain evaporation
!
!  Determine rain mass evaporated:
!
!  [mass_anrain_evap] = kg water/m2
!  [mass_uprain_evap] = kg water/m2
!
!------------------------------------------------------------------------
uprain_from_rlp = uprain_from_rlp - (mass_uprain_evap/tstep)
anrain = anrain - (mass_anrain_evap/tstep)
ansnow = ansnow - (mass_ansnow_subl/tstep)

if (uprain_from_rlp < -0.0000001) then 
  print*,'---------- NEGATIVE uprain_from_rlp ------'
  print*,'mass_uprain_evap = ',mass_uprain_evap 
  print*,'uprain_from_rlp = ',uprain_from_rlp
  print*,'old value uprain_from_rlp = ',uprain_from_rlp+(mass_uprain_evap/tstep)
  print*,'i z(i) = ',i,z(i)
  stop 'uprain_from_rlp negative in evap'
endif

if (anrain < -0.0000001) stop 'anrain negative in evap'

if (ansnow < -0.0000001) then 
  print*,'-------------'
  print*,'i = ',i
  print*,'t(i) = ',t(i)
  print*,'rh(i) = ',rh(i)
  print*,'ansnow = ',ansnow
  print*,'old ansnow = ',ansnow + (mass_ansnow_subl/tstep)
  print*,'mass_ansnow_subl = ',mass_ansnow_subl
  print*,'rv_mass_avail = ',rv_mass_avail
  print*,'rv_mass_avail/tstep = ',rv_mass_avail/tstep
  print*,'max_evap = ',max_evap
  stop 'ansnow negative in evap'
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  Stage 2: uprain/ansnow/anrain evaporation
!
!  Diagnostics:
!
!------------------------------------------------------------------------
uprainevap(i) = uprainevap(i) + mass_uprain_evap/tstep
anrainevap(i) = anrainevap(i) + mass_anrain_evap/tstep
ansnowsubl(i) = ansnowsubl(i) + mass_ansnow_subl/tstep

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  Stage 2: uprain/ansnow/anrain evaporation
!
!  Adjust hmuprain
!
!------------------------------------------------------------------------
km_wat = cl*t(i)
hm_of_wat = km_wat + g*z(i)

km_ice = cl*t(i) - fusion(t(i))
hm_of_ice = km_ice + g*z(i)

if (km_ice < 0.) then
  print*,'------- STOPPING: km_ice negative -----'
  print*,'km_ice in evap = ',km_ice 
  print*,'t(i) = ',t(i)
  print*,'fusion(t(i)) = ',fusion(t(i))
! stop
endif

hmuprain = hmuprain - (mass_uprain_evap/tstep)*hm_of_wat
hmanrain = hmanrain - (mass_anrain_evap/tstep)*hm_of_wat
hmansnow = hmansnow - (mass_ansnow_subl/tstep)*hm_of_ice

if (nprint == 1) then
  print*,'New hmansnow value after subl = ',hmansnow
  print*,'Amount added = ',-(mass_ansnow_subl/tstep)*hm_of_ice
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  Adjust hm of background atmosphere
!  - source of MSE from evap
!  - note that anrain downdrafts handled above
!
!  hmflow(2,i): is for detrainment.
!  [hmflow] = J/m2/s
!  [hmuprain] = J/m2/s
!
!  Local detrainment of hm from precip to the background atm.
!    Downdraft part not included (detrained
!    into surface layer). 
!
!------------------------------------------------------------------------
hmflow(2,i) = hmflow(2,i) + (mass_uprain_evap/tstep)*hm_of_wat
hmflow(2,i) = hmflow(2,i) + (mass_anrain_evap/tstep)*hm_of_wat
hmflow(2,i) = hmflow(2,i) + (mass_ansnow_subl/tstep)*hm_of_ice

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  Source of rv to background atmosphere (added to detrainment numrv = 2)
!  rvflow(numrv,nd)   F    [kg water/m2/s]
!  - actually no justification for including as separate
!
!------------------------------------------------------------------------
rvflow(2,i) = rvflow(2,i) + (mass_uprain_evap/tstep)   &
                          + (mass_anrain_evap/tstep)   &
                          + (mass_ansnow_subl/tstep) 

  rv_gain_uprain_evap = rv_gain_uprain_evap + mass_uprain_evap
  rv_gain_anrain_evap = rv_gain_anrain_evap + mass_anrain_evap
  rv_gain_ansnow_subl = rv_gain_ansnow_subl + mass_ansnow_subl

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
! Evaporation Diagnostics
!
!------------------------------------------------------------------------
if (nprint == 1) then
  print*,'-------------- IN EVAP level i = ',i
  print*,'rh(i) = ',rh(i)
  print*,'z(i) = ',z(i)
  print*,'t(i) = ',t(i)
  print*,'mass_uprain_evap = ',mass_uprain_evap
  print*,'mass_anrain_evap = ',mass_anrain_evap
  print*,'mass_ansnow_subl = ',mass_ansnow_subl
  print*,'Current value uprain uprain_from_rlp = ',uprain_from_rlp
  print*,'In mm per day Current value uprain uprain_from_rlp = ',uprain_from_rlp*3600.*24.
  print*,'uprain_start_rlp = ',uprain_start_rlp
  if (uprain_start_rlp > small) then
    print*,'Fraction of uprain_start_rlp = ',uprain_from_rlp/uprain_start_rlp
  endif
  print*,'ansnow = ',ansnow
  print*,'anrain = ',anrain
  print*,'ansnow_start = ',ansnow_start
  if (ansnow_start > small) then
    print*,'Fraction of original total ansnow = ',(ansnow+anrain)/ansnow_start
  endif
  print*,'---------------'
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!
!  *************  START PENETRATIVE DOWNDRAFTS ********************
!
!  Stage 3: anrain/uprain downdrafts
!
!  Set to zero if:
!  - T colder than 0 C
!  - downdraft mass exceeds allowable fraction of mass_start_dp
!
!------------------------------------------------------------------------

do_down_uprain = 1
do_down_anrain = 1

if (t(i) < t_down_min) then 
  do_down_uprain = 0
  do_down_anrain = 0
endif

if (i < i_min_down ) then 
  do_down_uprain = 0
  do_down_anrain = 0
endif

if (rh(i) > 1. ) then 
  do_down_uprain = 0
  do_down_anrain = 0
endif

if (rate_down_anrain < 0.00001 ) then 
  do_down_anrain = 0
endif

if (rate_down_uprain < 0.00001 ) then 
  do_down_uprain = 0
endif

if (anrain*3600.*24. < anrain_thresh) do_down_anrain = 0
if (uprain_from_rlp*3600.*24. < anrain_thresh) do_down_uprain = 0

drh_uprain = rh_uprain - rh(i)
drh_anrain = rh_ansnow - rh(i)

if (drh_uprain < 0.) do_down_uprain = 0
if (drh_anrain < 0.) do_down_anrain = 0

if (rh(i) > rh_down_allow) do_down_uprain = 0
if (rh(i) > rh_down_allow) do_down_anrain = 0

if (nprint == 1) then
  print*,'---------   for i = ',i
  print*,'do_down_uprain = ',do_down_uprain
  print*,'do_down_anrain = ',do_down_anrain
endif

!------------------------------------------------------------------------
!  
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!
!   Loop over anrain/uprain downdrafts
!
!------------------------------------------------------------------------
do nn = 1,2

if (nn == 1) then
  do_down = do_down_anrain 
  rate_down = rate_down_anrain
  md_rat = md_rat_anrain
  rain_down = anrain
else
  do_down = do_down_uprain
  rate_down = rate_down_uprain
  md_rat = md_rat_uprain
  rain_down = uprain_from_rlp
endif
rain_down_start = rain_down

!------------------------------------------------------------------------
!
!  Initialize downdraft properties as the same as the backgound atmosphere
!
!------------------------------------------------------------------------
if (do_down == 1) then

rv_down = rv(i)
rs_down = rs(i)
rh_down = rh(i)
hm_down = hm(i)
i_down = i
mass_down = 0.
drv_tot = 0.

!print*,'--------  Starting downdraft loop for nn i = ',nn,i
!print*,'mass_down = ',mass_down

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!
!  Loop over iterations
!  - only enter iteration if have moved down on previous iteration 
!    (Otherwise downdraft rh at max value anyway).
!
!------------------------------------------------------------------------
do iter = 1,iter_down
if (i_down > 1) then

!print*,'--------  Entering downdraft iteration = ',iter
!------------------------------------------------------------------------
!  
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!
!  Define drh_down
!  Here, drh_down should be WITH RESPECT TO THE BACKGROUND ATMOSPHERE at
!    the downdraft parcel level. Otherwise, as a parcel sinks you can have
!    lots of evaporation (and lots of induced upward motion) at height
!    levels where the RH of the atmosphere is very high (i.e. top of the
!    BL) and it gets driven to a higher and higher value.
!
!------------------------------------------------------------------------

if (nn == 1) then
  if (use_bg_rh == 1) then
    drh_down = rh_ansnow - rh(i_down)
  else
    drh_down = rh_ansnow - rh_down
  endif
else
  if (use_bg_rh == 1) then
    drh_down = rh_uprain - rh(i_down)
  else
    drh_down = rh_uprain - rh_down
  endif
endif

if (drh_down < 0.) then 
  drh_down = 0.
  drv_down = 0.
endif

!print*,'rh_down = ',rh_down
!print*,'rv_down = ',rv_down
!print*,'rs_down = ',rs_down
!print*,'drh_down = ',drh_down
!print*,'mass_down = ',mass_down

!------------------------------------------------------------------------
!  
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!
!  Do another iteration if still have rh_down < rh_ansnow/rh_uprain
!
!------------------------------------------------------------------------
if (drh_down > 0.0001) then

! print*,'ENTERED iteration iter i_down = ',iter,i_down

if (use_rh_scale == 1) then
  rh_factor = drh_down_scale*log(1. + (drh_down/drh_down_scale))
else
  rh_factor = drh_down
endif

!------------------------------------------------------------------------
!  
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!  *************  INSIDE drh_down > 0           ********************
!
!  calculate mass_water_evap (kg water/m2)
!  - rain_down is initialized as the anrain or uprain at the starting
!    level, and then progressively reduced as evaporation proceeds. It
!    is not increased due to rain production below the starting level.
!
!------------------------------------------------------------------------

  mass_water_evap = 0.
  f_evap = rs_down*rh_factor*dmdry*rate_down
  f_evap_down_max = 0.0001*f_evap_down_press*(phalf(i)-phalf(i+1))
  if (f_evap > f_evap_down_max) f_evap = f_evap_down_max
  if (f_evap > 0.5) stop 'WARNING f_evap too big'
  mass_water_evap = f_evap*rain_down*tstep
  if (mass_water_evap < -0.000000001) stop 'mass_water negative'

if (nprint == 1) then
  print*,'--------------'
  print*,'CALCULATING mass_water_evap'
  print*,'f_evap_down_max = ',f_evap_down_max
  print*,'f_evap = ',f_evap
  print*,'(phalf(i)-phalf(i+1)) = ',(phalf(i)-phalf(i+1))
  print*,'initial mass_water_evap = ',mass_water_evap
  print*,'Current value of mass_down = ',mass_down
endif

!------------------------------------------------------------------------
!  
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!  *************  INSIDE drh_down > 0           ********************
!
!  Define mass_down on first iteration 
!  - possibly modify if too large a fraction of layer mass
!  - in this case proportionately reduce mass_water_evap also
!
!------------------------------------------------------------------------
  fff = f_mass_max*mass_start_tot
  if (fff < fff_min) fff = fff_min

  if (iter == 1) then 
    mass_down = md_rat*mass_water_evap
    dmdry = phalf(i)-phalf(i+1)/(g*(1.+rv(i)+rl(i)+ri(i)))
    f_mass = mass_down/dmdry
    if (nprint == 1) then
       print*,'--------------'
       print*,'Calculating mass_down on first iteration'
       print*,'md_rat = ',md_rat
       print*,'mass_water_evap = ',mass_water_evap
       print*,'mass_down = ',mass_down
       print*,'f_mass = ',f_mass
    endif
    if (f_mass > fff) then
      mass_down_old = mass_down
      mass_down = fff*dmdry
      mass_water_evap = (mass_down/mass_down_old)*mass_water_evap
      if (nprint == 1) then
         print*,'--------------'
         print*,'New mass_down = ',mass_down
         print*,'mass_start_tot = ',mass_start_tot
         print*,'Redefining mass_down since f_mass = ',f_mass
         print*,'fff = ',fff
      endif
      f_mass = mass_down/dmdry
    endif
    if (nprint == 1) then
      print*,'--------------'
      print*,'mass_down defined on first iteration as ',mass_down
      print*,'fraction of mass removed f_mass = ',f_mass
    endif
  endif

!------------------------------------------------------------------------
!  
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!  *************  INSIDE drh_down > 0           ********************
!
!   Define drv_down
!   - check if drv_max exceeds limit.
!   - If so, adjust mass_water_evap
!   - In principle should calculate new temperature here, and get real RH.
!
!------------------------------------------------------------------------
  if (iter == 1) then
    drv_down = 1./md_rat
  else
    if (mass_down > 1.0E-8) then
      drv_down = mass_water_evap/mass_down
    else
      drv_down = 0.
    endif
  endif
  drv_max = (rh_down_max - rh_down)*rs_down

  f_drv = drv_max/rv_down
  if (f_drv > f_drv_max) then
    drv_max = f_drv_max*rs_down
  endif

if (nprint == 1) then
  print*,'--------------'
  print*,'Preliminary drv_down = ',drv_down
  print*,'drv_max = ',drv_max
endif

  if (drv_down > drv_max) then
    if (nprint == 1) then
      print*,'--------------'
      print*,'Adjusting mass_water_evap by ratio = ',(drv_max/drv_down)
    endif
    mass_water_evap = (drv_max/drv_down)*mass_water_evap
    f_evap = mass_water_evap/(rain_down*tstep)
    drv_down = drv_max
    if (nprint == 1) then
      print*,'New drv_down = ',drv_down
    endif
  endif

if (nprint == 1) then
  print*,'--------------'
  print*,'f_evap = ',f_evap
  print*,'final mass_water_evap = ',mass_water_evap
endif

!------------------------------------------------------------------------
!  
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!  *************  INSIDE drh_down > 0           ********************
!
!   Initialize parcel
!
!  What is hm_down?
!  hm_of_wat = km_wat + g*z(i) = cl*t(i) + g*z(i)
!  mass_water = mass_down*drv_down
!  HM of parcel = hm(i)*mass_down + hm_of_wat*mass_water
!               = hm(i)*mass_down + hm_of_wat*mass_down*drh_parcel*rs(i)
!  hm_down = HM/mass_down
!          = hm(i) + hm_of_wat*drh_parcel*rs(i)
!          = hm(i) + (rv_down - rv(i))*hm_of_wat
!
!------------------------------------------------------------------------

  drv_tot = drv_tot + drv_down  ! update accumulated evap
  rv_down = rv_down + drv_down
  km_wat = cl*t(i_down)
  hm_of_wat = km_wat + g*z(i_down)
  hm_down = hm_down + drv_down*hm_of_wat
  km_down = hm_down - (1. + rv_down)*g*z(i_down)

water_error = abs(drv_down*mass_down - mass_water_evap)

if (nprint == 1) then
  print*,'--------------'
  print*,'New rv_down = ',rv_down
  print*,'water_error = ',water_error
endif

if (water_error > 0.000001) then
  print*,'---------  water error ----------'
  print*,'drv_down*mass_down = ',drv_down*mass_down
  print*,'drv_down = ',drv_down
  print*,'mass_down = ',mass_down
  print*,'mass_water_evap = ',mass_water_evap
  stop 'water error'
endif

hm_error = abs(drv_down*hm_of_wat*mass_down - mass_water_evap*hm_of_wat)
!print*,'hm_error = ',hm_error
if (hm_error > 0.0000001) then
  print*,'--------------'
  print*,'hm_of_wat = ',hm_of_wat
  print*,'drv_down*hm_of_wat*mass_down = ',drv_down*hm_of_wat*mass_down
  print*,'mass_water_evap*hm_of_wat = ',mass_water_evap*hm_of_wat
  stop 'hm_error'
endif


if (nprint == 1) then
  print*,'---------------------------'
  print*,'Entering downdraft loop i nn = ',i,nn
  print*,'fraction_remove = ',fraction_remove
  print*,'New rv_down = ',rv_down 
  print*,'preliminary rs_down before evap cooling = ',rs_down 
  print*,'Increase in parcel RH = ',drv_down/rs_down
  print*,'preliminary new parcel rh before evap cooling = ',rv_down/rs_down
  print*,'rh_down_max = ',rh_down_max
  print*,'drh_down = ',drh_down
  print*,'rh(i) = ',rh(i)
  print*,'mass_down = ',mass_down
  print*,'drv_down = ',drv_down
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!  *************  INSIDE drh_down > 0           ********************
!
!  Update precipitation diagnostics.
!  - need current value of hm_of_wat here
!  - a bit strange here: at the starting level, anrain initialized as the 
!    remaining anrain at that height. anraindown is where evaporation actually
!    occurs ...
!
!------------------------------------------------------------------------

if (nn == 1) then
  anrain = anrain - (mass_water_evap/tstep)
  rain_down = anrain
  anraindown(i_down) = anraindown(i_down) + (mass_water_evap/tstep)
  hmanrain = hmanrain - (mass_water_evap/tstep)*hm_of_wat
  if (anrain < -0.0000000001) then
    print*,'anrain = ',anrain
    stop 'anrain negative'
  endif
  if (nprint == 1) then
    print*,'--------  Adjusting anrain  ----------'
    print*,'New anrain = ',anrain
    print*,'Subtracted amount = ',mass_water_evap/tstep
    print*,'----------------------------------------'
  endif
else
  uprain_from_rlp = uprain_from_rlp - (mass_water_evap/tstep)
  rain_down = uprain_from_rlp
  upraindown(i_down) = upraindown(i_down) + (mass_water_evap/tstep)
  hmuprain = hmuprain - (mass_water_evap/tstep)*hm_of_wat
  if (uprain_from_rlp < 0.) then
    print*,'=========  NEGATIVE uprain_from_rlp after DOWNDRAFTS ========'
    print*,'uprain_from_rlp = ',uprain_from_rlp
    stop  'NEGATIVE uprain'
  endif
  if (nprint == 1) then
    print*,'--------  Adjusting uprain_from_rlp  ----------'
    print*,'New uprain_from_rlp = ',rain_down
    print*,'Subtracted amount = ',mass_water_evap/tstep
    print*,'----------------------------------------'
  endif
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!  *************  INSIDE drh_down > 0           ********************
!
!  Part 2: Calculate Buoyancy at current and lower level
!
! DDD
!
!------------------------------------------------------------------------

!
!  buoyancy and RH at current level
!
  get_t_call = 1
  ri_down = 0.
  rl_down = 0.
  call t_from_km(km_down,p(i_down),rv_down,ri_down,rl_down,t_down,get_t_call)
  td_down = t_down*(1. + epsi*rv_down)/(1. + rv_down)
  b_down = g*(td_down - td(i_down))/td(i_down)
  rs_down = rsat(t_down,p(i_down))
  rh_down = rv_down/rs_down

if (nprint == 1) then
  print*,'-------------------'
  print*,'Buoyancy at current level '
  print*,'km_down = ',km_down
  print*,'t(i_down) = ',t(i_down)
  print*,'t_down = ',t_down
  print*,'Buoyancy at current level = ',b_down
  print*,'----------------'
end if

  if (rh_down > 1.20) then
    print*,'-----  WARNING: IN DOWNDRAFTS RH exceeds MAX VALUE ----'
    print*,'rh_down = ',rh_down
    print*,'drv_down = ',drv_down
    print*,'Fractional increase = ',drv_down/rv_down
    print*,'b_down = ',b_down
    print*,'t_down = ',t_down
    print*,'------------------------------------'
  endif

  rh_dn(i_down) = rh_dn(i_down) + rh_down*mass_down
  mass_rh_dn(i_down) = mass_rh_dn(i_down) + mass_down

!
!  Loop over all lower levels and find buoyancy
!
i_lowest = i_down
 
!do i_lower = 1,i_down-1
do i_lower = i_down-1,i_down-1
  if (i_down == 1) stop 'i_down is one'
!------------------------------------------------------------------
!  Find b_lower
!------------------------------------------------------------------
  km_lower = hm_down - (1. + rv_down)*g*z(i_lower)
  call t_from_km(km_lower,p(i_lower),rv_down,ri_down,rl_down,t_lower,get_t_call)
  td_lower = t_lower*(1. + epsi*rv_down)/(1. + rv_down)
  b_lower = g*(td_lower - td(i_lower))/td(i_lower)
!------------------------------------------------------------------
!  Diagnostics
!------------------------------------------------------------------
if (nprint == 1) then
  print*,'----- in i_lower = ',i_lower
  print*,'km_lower = ',km_lower
  print*,'t(i_lower) = ',t(i_lower)
  print*,'t_lower = ',t_lower
  print*,'Buoyancy at lower level = ',b_lower
endif
!------------------------------------------------------------------
!  Find new i_lowest
!------------------------------------------------------------------
  b_av = 0.5*(b_down + b_lower) + b_precip_down
  dmdry = phalf(i_lower)-phalf(i_lower+1)
  dmdry = dmdry/(g*(1.+rv(i_lower)+rl(i_lower)+ri(i_lower)))
  f_mass = (dndet(i_lower)*tstep + mass_down)/dmdry
  if (nprint == 1) then
     print*,'f_mass = ',f_mass
     print*,'fff = ',fff
  endif
  if ((b_av < 0.).and.(i_lower < i_lowest).and.(f_mass < fff)) then
    i_lowest = i_lower
    rs_down = rsat(t_lower,p(i_lowest))
    rh_down = rv_down/rs_down
    if (nprint == 1) then
      print*,'---------------'
      print*,'Moving down to level i_lower = ',i_lower
      print*,'new rs_down at lower level = ',rs_down
      print*,'new rh_down at lower level = ',rh_down
    endif
  endif
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
end do
i_down = i_lowest

if (nprint == 1) then
  print*,'ENDING an iteration'
endif


!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE LOOP OVER ITERATIONS   ********************
!  *************  INSIDE drh_down > 0           ********************
!
!  endif for drh_down > 0
!  endif for i_down > 1
!  end loop over iterations
!
!------------------------------------------------------------------------
endif
endif
end do

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!
!  Define fmass_dn
!  - purely diagnostic
!  - want to consider fmass_dn to be positive at half level heights between
!    i and i_down.
!  - Suppose i_down = i-1, then the half level height between is i
!
!------------------------------------------------------------------------

if (i_down < i) then
  do ii = i_down+1,i
    fmass_dn(ii) = fmass_dn(ii) + (mass_down/tstep)
  end do
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!
!  Stage 3: anrain/uprain downdrafts
!
!  Entrain and detrain
!
!  Units of hmflow:   [J/m2*s] 
!  Units of hm_of_wat: J/kg water
!
!------------------------------------------------------------------------
 
if (nprint == 1) then
  print*,'------- FINAL DETRAIN OF PARCEL: Updating flow variables for nn = ',nn
  print*,'i = ',i
  print*,'i_down = ',i_down
  print*,'mass_down = ',mass_down
  print*,'rv_down = ',rv_down
  print*,'drv_tot = ',drv_tot
  print*,'rv(i) = ',rv(i)
  print*,'------------'
end if

rv_error = abs(rv(i) + drv_tot - rv_down)
!print*,'rv_error = ',rv_error
if (rv_error > 0.00000001) then
  stop 'rv_error'
endif

  dmflow(1,i) = dmflow(1,i) + (mass_down/tstep)
  rvflow(1,i) = rvflow(1,i) + (rv(i)*mass_down/tstep)
  hmflow(1,i) = hmflow(1,i) + (hm(i)*mass_down/tstep)
  uflow(1,i) = uflow(1,i) + (uwind(i)*mass_down/tstep)
  vflow(1,i) = vflow(1,i) + (vwind(i)*mass_down/tstep)

!  detrain

  dmflow(2,i_down) = dmflow(2,i_down) + (mass_down/tstep)
  rvflow(2,i_down) = rvflow(2,i_down) + (rv_down*mass_down/tstep)
  hmflow(2,i_down) = hmflow(2,i_down) + (hm_down*mass_down/tstep) 
  uflow(2,i_down) = uflow(2,i_down) + (uwind(i_down)*mass_down/tstep)
  vflow(2,i_down) = vflow(2,i_down) + (vwind(i_down)*mass_down/tstep)

if (nn == 1) then
  rv_gain_anrain_down = rv_gain_anrain_down + (rv_down-rv(i))*mass_down
else if (nn == 2) then
  rv_gain_uprain_down = rv_gain_uprain_down + (rv_down-rv(i))*mass_down
endif

if (nprint == 1) then
  print*,'===========  CHECK ==============='
  print*,'rv gained by atmosphere = ',(rv_down-rv(i))*mass_down
  print*,'This should equal lost by rain = ',(rain_down_start-rain_down)*tstep
  print*,'===========  CHECK ==============='
endif

!------------------------------------------------------------------------
!
!  *************  INSIDE LOOP OVER HEIGHTS *************************
!  *************  INSIDE PENETRATIVE DOWNDRAFTS ********************
!
! Downdraft diagnostics
!
!------------------------------------------------------------------------
  dnent(i) = dnent(i) + (mass_down/tstep)
  dndet(i_down) = dndet(i_down) + (mass_down/tstep)

  rs_grid = rsat(t(i_down),p(i_down))
  rh_grid = rv_down/rs_grid
  drv_dndet(i_down) = drv_dndet(i_down) + ((rv_down-rv(i_down))*mass_down/tstep)
  dt_dndet(i_down) = dt_dndet(i_down) + ((t_down-t(i_down))*mass_down/tstep)
  if (mass_down > 0.000000001) then
!   print*,'i_down rh_grid mass_down = ',i_down,rh_grid,mass_down
  endif

  if (rh_grid > 1.20) then
    print*,'WARNING RH_GRID OF DETRAINING DOWNDRAFT PARCEL LARGE -----'
    print*,'RH with respect to grid temperature rh_grid = ',rh_grid
    print*,'i_down = ',i_down
    print*,'p(i_down) = ',p(i_down)
    print*,'t(i_down) = ',t(i_down)
    print*,'rs_grid = ',rs_grid
    print*,'downdraft t_down = ',t_down
    print*,'downdraft rv_down = ',rv_down
    print*,'downdraft rs_down = ',rs_down
    print*,'downdraft rh_down = ',rh_down
    print*,'-----------------------------'
  endif

if (nn == 1) drv_an(i) = drv_tot
if (nn == 2) drv_up(i) = drv_tot

!------------------------------------------------------------------------
!  
!  *********** END OF DOWNDRAFTS   *************************
!
!  endif for do_down = 1  (anrain or uprain)
!  end loop over nn
!
!------------------------------------------------------------------------
endif
end do

if (nprint == 1) then
   print*,'END OF EVAP/DOWN LEVEL'
endif

!------------------------------------------------------------------------
!
!  End loop over levels
!
!------------------------------------------------------------------------
end do
endif

if (nprint == 1) then
  print*,'END LOOP OVER LEVELS'
endif

!------------------------------------------------------------------------
!
!  *************  AFTER LOOP OVER HEIGHTS *************************
!
!  Diagnostics
!  - all have units [kg water/m2/s]
!  - initialized as zero earlier
!
!------------------------------------------------------------------------
do i = 1,nd
  uprain_evap = uprain_evap + uprainevap(i)
  uprain_down = uprain_down + upraindown(i)
  anrain_down = anrain_down + anraindown(i)
  anrain_evap = anrain_evap + anrainevap(i)
  ansnow_subl = ansnow_subl + ansnowsubl(i)
  ansnow_melt = ansnow_melt + ansnowmelt(i)
  if (dndet(i) > 0.00000001) then
    drv_dndet(i) = drv_dndet(i)/dndet(i)
    dt_dndet(i) = dt_dndet(i)/dndet(i)
  else
    drv_dndet(i) = bad
    dt_dndet(i) = bad
  endif
  if (mass_rh_dn(i) > 0.0001) then
    rh_dn(i) = rh_dn(i)/mass_rh_dn(i)
!   print*,'i rh_dn(i) = ',i,rh_dn(i)
  else
    rh_dn(i) = bad
  endif
end do

!------------------------------------------------------------------------
!
!  *************  AFTER LOOP OVER HEIGHTS *************************
!
!  Define surface variables
!  - should all still be initialized as zero
!
!------------------------------------------------------------------------
uprain_surf_rlp = uprain_from_rlp

anrain_surf = anrain
ansnow_surf = ansnow

hmuprain_surf = hmuprain
hmanrain_surf = hmanrain
hmansnow_surf = hmansnow

!------------------------------------------------------------------------
!
!  *************  AFTER LOOP OVER HEIGHTS *************************
!
!  ansnow Conservation check:
!
!------------------------------------------------------------------------
error_ansnow = abs(ansnow_start - ansnow_melt - ansnow_subl - ansnow_surf)

if (error_ansnow > small) then
  print*,'=============== ANSNOW WATER ERROR in EVAP ============='
  print*,'error_ansnow = ansnow_start - ansnow_melt - ansnow_subl - ansnow_surf'
  print*,'i_top = ',i_top
  print*,'error_ansnow = ',error_ansnow
  print*,'ansnow_start (mm/day) = ',ansnow_start*3600.*24.
  print*,'ansnow_surf (mm/day) = ',ansnow_surf*3600.*24.
  print*,'ansnow_melt (mm/day) = ',ansnow_melt*3600.*24.
  print*,'ansnow_subl (mm/day) = ',ansnow_subl*3600.*24.
  print*,'hmansnow_start = ',hmansnow_start
  print*,'hmansnow_surf = ',hmansnow_surf
  stop 'ansnow error in evap'
endif

!------------------------------------------------------------------------
!
!  *************  AFTER LOOP OVER HEIGHTS *************************
!
!  anrain Conservation check:
!
!------------------------------------------------------------------------
error_anrain = abs(ansnow_melt - anrain_down - anrain_evap - anrain_surf)

if (error_anrain > small) then
  print*,'=============== ANRAIN WATER ERROR in EVAP ============='
  print*,'i_top = ',i_top
  print*,'error_anrain = ',error_anrain
  print*,'Fractional error = ',error_anrain/ansnow_melt
  print*,'ansnow_melt = ',ansnow_melt
  print*,'anrain_down = ',anrain_down
  print*,'anrain_evap = ',anrain_evap
  print*,'anrain_surf = ',anrain_surf
  print*,'ansnow_melt (mm/day) = ',ansnow_melt*3600.*24.
  print*,'anrain_down (mm/day) = ',anrain_down*3600.*24.
  print*,'anrain_evap (mm/day) = ',anrain_evap*3600.*24.
  print*,'anrain_surf (mm/day) = ',anrain_surf*3600.*24.
  print*,'hmansnow_start = ',hmansnow_start
  print*,'hmansnow_surf = ',hmansnow_surf
  print*,'hmanrain_surf = ',hmanrain_surf
  stop 'anrain water error in evap'
endif

!------------------------------------------------------------------------
!
!  *************  AFTER LOOP OVER HEIGHTS *************************
!
!  uprain Conservation check:
!
!------------------------------------------------------------------------
error_up = abs(uprain_evap + uprain_down - (uprain_start_rlp-uprain_surf_rlp))

if (error_up > small) then
  print*,'=============== UP WATER ERROR in EVAP ============='
  print*,'error_up = ',error_up
  print*,'uprain_down = ',uprain_down
  print*,'uprain_evap = ',uprain_evap
  print*,'i_top = ',i_top
  print*,'uprain_start_rlp = ',uprain_start_rlp*3600.*24.
  print*,'uprain_surf_rlp = ',uprain_surf_rlp*3600.*24.
  print*,'Decrease in uprain_rlp (kg/m2/s) = ',uprain_start_rlp-uprain_surf_rlp
  print*,'Should equal uprain_evap + uprain_down = ',uprain_evap + uprain_down
  if (uprain_start_rlp > small) then
    print*,'Fraction uprain reaching surface = ',uprain_surf_rlp/uprain_start_rlp
  endif
  print*,'hmuprain_start = ',hmuprain_start
  print*,'hmuprain_surf = ',hmuprain_surf
  if (hmuprain_start > small) then
    print*,'Fraction hmuprain reaching surface = ',hmuprain_surf/hmuprain_start
  endif
  print*,'mass_uprain_evap = ',mass_uprain_evap
  stop 'water error in evap'
endif

!------------------------------------------------------------------------
!
!  *************  AFTER LOOP OVER HEIGHTS *************************
!
! Diagnostics
!------------------------------------------------------------------------
if (nprint == 1) then
  print*,'-------------- LEAVING EVAP  ------------'
  print*,'i_top = ',i_top
  print*,'uprain_start_rlp = ',uprain_start_rlp*3600.*24.
  print*,'uprain_surf_rlp = ',uprain_surf_rlp*3600.*24.
  print*,'Decrease in uprain (kg/m2/s) = ',uprain_start_rlp-uprain_surf_rlp
  print*,'Should equal uprain_evap+anrain_down = ',uprain_evap+anrain_down
  if (uprain_start_rlp > small) then
    print*,'Fraction uprain reaching surface = ',uprain_surf_rlp/uprain_start_rlp
  endif
  print*,'hmuprain_start = ',hmuprain_start
  print*,'hmuprain_surf = ',hmuprain_surf
  if (hmuprain_start > small) then
    print*,'Fraction hmuprain reaching surface = ',hmuprain_surf/hmuprain_start
  endif
  print*,'ansnow_start = ',ansnow_start*3600.*24.
  print*,'ansnow_surf = ',ansnow_surf*3600.*24.
  print*,'anrain_surf = ',anrain_surf*3600.*24.
  print*,'Decrease in ansnow (kg/m2/s) = ',ansnow_start-ansnow_surf-anrain_surf
  print*,'Should equal anrain_evap = ',anrain_evap
  if (ansnow_start > small) then
    print*,'Fraction ansnow at surface = ',(ansnow_surf+anrain_surf)/ansnow_start
  endif
  print*,'hmansnow_start = ',hmansnow_start
  print*,'hmansnow_surf = ',hmansnow_surf
  print*,'hmanrain_surf = ',hmanrain_surf
  print*,'cape_1 = ',cape_1
  print*,'cape_2 = ',cape_2
  print*,'cape_3 = ',cape_3
  print*,'cape_4 = ',cape_4
  print*,'----------------------------------'
endif

!------------------------------------------------------------------------
!
!  *************  AFTER LOOP OVER HEIGHTS *************************
!
! Diagnostics for water conservation
! - does total gain by atm equal precip loss?
!
!------------------------------------------------------------------------
precip_loss = ansnow_start + uprain_start_rlp - anrain_surf -  &
              ansnow_surf - uprain_surf_rlp
precip_loss = precip_loss*tstep
tot_rv_gain = rv_gain_anrain_down + rv_gain_uprain_down +  &
              rv_gain_anrain_melt + rv_gain_ansnow_subl +  &
              rv_gain_uprain_evap + rv_gain_anrain_evap

if (nprint == 1) then
  print*,'------  FINAL WATER CHECKS ---------------'
  print*,'rv_gain_anrain_down = ',rv_gain_anrain_down
  print*,'rv_gain_uprain_down = ',rv_gain_uprain_down
  print*,'rv_gain_anrain_melt = ',rv_gain_anrain_melt
  print*,'rv_gain_ansnow_subl = ',rv_gain_ansnow_subl
  print*,'rv_gain_uprain_evap = ',rv_gain_uprain_evap
  print*,'rv_gain_anrain_evap = ',rv_gain_anrain_evap
  print*,'tot_rv_gain = ',tot_rv_gain
  print*,'Should equal precip_loss = ',precip_loss
  print*,'-------------------------------------------'
endif

!------------------------------------------------------------------------
!
!  Average downward movement 
!  Int(p(det-ent))/Int(ent)
! 
!------------------------------------------------------------------------
xxx = 0.
yyy = 0.
zzz = 0.
do i = 1,nd
  xxx = xxx + dndet(i)*p(i)
  yyy = yyy + dnent(i)*p(i)
  zzz = zzz + dnent(i)
end do
dp_dn = 0.
if (zzz > 0.00000001) then
  dp_dn = 0.01*(xxx - yyy)/zzz
! print*,'dp_dn = ',dp_dn
endif

!------------------------------------------------------------------------
!
!  END  (subroutine evap)
!
!------------------------------------------------------------------------
end subroutine evap
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Local real variables
!
!   subroutine calc_vert
!
!------------------------------------------------------------------------
subroutine calc_vert(tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: rvv,totentrain,totrl_dest,totrv_ansnow,totri_ansnow,totrv_uprain
real(r8) :: drymass,vapmass,condmass,small,totrl_uprain,totriflow
real(r8) :: dp_fraction,ri_move,om_extra
integer  :: i,nums,npp

!------------------------------------------------------------------------
!
!  Initialize changes in surface pressure due to updrafts/downdrafts as xero.
!
!------------------------------------------------------------------------
  dpsurf = 0.
!------------------------------------------------------------------------
! 
!   Brewer Dobson adjustment:
!   (updrafts only)
!
!   Add in Brewer Dobson circulation by adjusting dmflow/hmflow/
!   dmflow(1,i) : updraft entrainment
!   dmflow(2,i) : updraft detrainment
!
!   BD transport:
!   positive change in updraft detrainment in the troposphere (BD circulation
!     is a source of dry mass/rv), and negative change in updraft detrainment
!     in the stratosphere (removal).
!
!   Look at bdflow and make sure changes sign.
!
!   momentum: assume that there is no input or removal of momentum
!     associated with the BD circulation (uflow/vflow unchanged)
!
!------------------------------------------------------------------------
  do i = 1,nd
    dmflow(2,i) = dmflow(2,i) + bdflow(i)
    hmflow(2,i) = hmflow(2,i) + hm(i)*bdflow(i)
    rvflow(2,i) = rvflow(2,i) + rv(i)*bdflow(i)
  end do
!------------------------------------------------------------------------
!
!  Need criteria for doing vert calculations.
!
!  All vert arrays should be previously initialized as zero.
!
!  Should calculate vert arrays here if:
!
!  (i) updraft dry mass entrainment exceeds some threshold
!  (ii) downdraft dry mass entrainment exceeds some threshold
!  (iii) BD circulation on
!  (iv) conversion of rp to ansnow  (results in background vertical motions
!      since consider precipitation to be removal of mass from column)
!
!  Do not need to calculate vert arrays for changes which do not contribute
!    to vertical motions:
!
!  (i) conversion of rl to rp 
!  (ii) aging of rl
!
!  Still need to calculate tendencies from these non-zero flow arrays, however.
!
!------------------------------------------------------------------------
  totriflow = 0.
  totentrain = 0.
  totentrain = 0.
  totrl_dest = 0.
  totrv_ansnow = 0.
  totri_ansnow = 0.
  totrv_uprain = 0.
  totrl_uprain = 0.
  do i = 1,nd
    totriflow = totriflow + abs(riflow(1,i)) 
    totentrain = totentrain + dmflow(1,i) 
    totentrain = totentrain + dmflow(1,i) 
    totrv_ansnow = totrv_ansnow + ansnow_rv(i)
    totri_ansnow = totri_ansnow + ansnow_ri(i)
    totrv_uprain = totrv_uprain + uprain_rv(i)
    totrl_uprain = totrl_uprain + uprain_rl(i)
  end do
! print*,'totentrain = ',totentrain
! print*,'totentrain = ',totentrain
  small = 1.0E-14
!------------------------------------------------------------------------
!
!  return if nothing happening
!
!------------------------------------------------------------------------
  if ((abs(totentrain) < small).and.   &
      (abs(totriflow) < small).and.    &
      (abs(totrv_ansnow) < small).and. &
      (abs(totri_ansnow) < small).and. &
      (abs(totrv_uprain) < small).and. &
      (abs(totrl_uprain) < small).and. &
      (abs(totentrain) < small).and.   &
      (abs(totrl_dest) < small)) return   
! print*,'ENTERING fmass calculation'
!------------------------------------------------------------------------
!
!   Defining Induced Mass Fluxes fmass(1):
!   ------------------------------------
!
!   Retain separate up/down mass fluxes for diagnostic purposes, even though
!     could combine.
!
!   fmass(nd+1): Updraft mass flux: [kg/m2*s]
!
!   negative fmass: downward
!
!   ------------   zh(nd+1)  fmass(nd+1) -----------------  TOP HALF LEVEL
!
!   +++++++++++++ z(nd), dmflow(1,nd), dmflow(2,nd) +++++++++++++++
!
!   -------------   zh(nd)  fmass(nd)   -------------------
!
!   Guiding principles in defining mass flux:
!   -----------------------------------------
!   (1) Keep the top model level phalf(nd+1) locked in place.
!   (2) Set the cloud mass fluxes at the top of the atmosphere = 0.
!   (3) Define the mass fluxes downward from the top in such a way that
!       that the dp of every layer is fixed.
!   (4) In the case of precipitation, there will be more water vapor removed
!       from the column then added. This will give rise to an induced upward
!       mass flux at the surface. This is the reduction in surface pressure
!       due to rainfall.
!
!------------------------------------------------------------------------
  fmass(nd+1) = 0.            ! induced mass flux zero at the top
!------------------------------------------------------------------------
!
!  Loop over levels
!
!------------------------------------------------------------------------
  do i = nd,1,-1
!------------------------------------------------------------------------
!
!   Basic Total Cloud Mass Flux Formula: fmass(1)
!   ---------------------------------------------
!   This loop finds the mass flux
!   at the bottom of every grid box, that keeps the mass in every grid
!   box constant, except for the bottom level, in which case, there is
!   an adjustment in surface pressure phalf(1).
!
!   Uses conservation of mass: for fixed pressure levels, the sum of all horizontal
!    and vertical inflows of mass must equal to zero. With mass flux at the top of
!    a layer assumed known, and the horizontal mass flows known, it solves for
!    the mass flux through the bottom of a layer.
!
!   Units of dmflow/rvflow/rlflow/riflow and fmass should all be kg/m2/sec 
!     - no need to adjust units.
!
!   drymass : net dry mass source (detrainment - entrainment)
!   dmflow(2,i): net detrainment
!   dmflow(1,i): net entrainment
!
!------------------------------------------------------------------------
  drymass = dmflow(2,i) - dmflow(1,i)  ! dry mass (detrain-entrain)
  vapmass = rvflow(2,i) - rvflow(1,i)  ! vapor mass (detrain-entrain)
  condmass = rlflow(2,i) + riflow(2,i) - & !
             rlflow(1,i) - riflow(1,i)     ! rl+ri (detrain-entrain)
!------------------------------------------------------------------------
!
!   "Special" ways updrafts remove mass from a layer
!   - ansnow_rv could be counted as "updraft" or "downdraft"
!   - condmass not a good name, more like "extras"
!
!------------------------------------------------------------------------
    condmass =  condmass              &    ! 
              - ansnow_rv(i)     &    !  removes rv
              - uprain_rv(i)     &    !  removes rv
              - ansnow_ri(i)     &    !  removes ri
              - uprain_rl(i)          !  removes rl
!------------------------------------------------------------------------
!  fmass: net change in mass 
!------------------------------------------------------------------------
  fmass(i) = fmass(i+1) - drymass - vapmass - condmass   ! fmass is positive up
!------------------------------------------------------------------------
! print stuff
!------------------------------------------------------------------------
npp = 0
if ((npp == 1).and.(i < 15)) then  
 print*,'------------- i = ',i
 print*,'fmass(i) = ',fmass(i)
 print*,'dmflow(1,i) = ',dmflow(1,i)
 print*,'dmflow(2,i) = ',dmflow(2,i)
 print*,'rlflow(1,i) = ',rlflow(1,i)
 print*,'rlflow(2,i) = ',rlflow(2,i)
 print*,'rl(i) = ',rl(i)
 print*,'drymass = ',drymass
 print*,'vapmass = ',vapmass
 print*,'condmass = ',condmass
endif
!------------------------------------------------------------------------
!
!  
!
!------------------------------------------------------------------------
  dpsurf = dpsurf + g*(drymass + vapmass + condmass)*tstep

!------------------------------------------------------------------------
! print stuff
!------------------------------------------------------------------------
if ((npp == 1).and.(i < 15)) then  
  print*,'dpsurf = ',dpsurf
endif
!------------------------------------------------------------------------
!
!  End loop over levels.
!
!------------------------------------------------------------------------
  end do

!------------------------------------------------------------------------
!
!   fmass(i) now defined.
!
!   Want to define vertical advections of various species
!
!   Since there is a lot of repetition here, should define a subroutine
!     and call it for each species.
!
!   Loop over nd:
!
!   Define: dmvert, rvvert, rlvert, rivert, hmvert
!
!   All these "vert" arrays should be initialized to zero in init_zero
!   These are the net local changes in a grid cell due to vertical advection
!   fmass has units of kg/m2*s, and is positive upward
!   Could move these initializations to a subroutine
!   hm(i) refers to the total mse of background air at level i
!     (dry air + rv/rl), so hmvert refers to vertical mse transport by
!     both terms and similarly for hmflow.
!
!------------------------------------------------------------------------
!  print*,'start loop over nd to find fmass'
  do i = 1,nd
!----------------------------------------------------------------------73
!
!  Upward flux at top of layer i into layer i+1: 
!  Calculate loss from layer i. 
!  The upward mass flux through the top of the layer is assumed to have the
!    properties of the layer.
!  For (i == nd), this flux should be zero by definition.
!    It might be non-zero due to the BD circulation.
!    Maybe because top pressure layer of BD outflow is
!    above model top. Should have an error message. 
!
!----------------------------------------------------------------------73
    if (fmass(i+1) >= 0.) then
      if (turn_off_cond_adv == 2) then 
        rvv = 1.+rv(i)
        rivert(i) = 0.
        rlvert(i) = 0.
      elseif (turn_off_cond_adv == 1) then 
        rvv = 1.+rv(i)+rl(i)
        rivert(i) = 0.
        rlvert(i) = -rl(i)*fmass(i+1)/rvv
      else
        rvv = 1.+rv(i)+rl(i)+ri(i)
        rivert(i) = -ri(i)*fmass(i+1)/rvv
        rlvert(i) = -rl(i)*fmass(i+1)/rvv
      endif
      dmvert(i) = -fmass(i+1)/rvv
      rvvert(i) = -rv(i)*fmass(i+1)/rvv
      hmvert(i) = -hm(i)*fmass(i+1)/rvv
      uvert(i) = -uwind(i)*fmass(i+1)/rvv
      vvert(i) = -vwind(i)*fmass(i+1)/rvv
!----------------------------------------------------------------------73
!  Downward flux at top of layer i from layer i+1: GAIN at layer i
!  Have to multiply by negative since fmass < 0.
!  (i == nd): stop: Flux should never be downward
!  air has properties of layer above.
!----------------------------------------------------------------------73
    else
      if (i == nd) then
        stop 'downward mass flux at top level not good'
      else
        if (turn_off_cond_adv == 2) then 
          rvv = 1. + rv(i+1) 
          rivert(i) = 0.
          rlvert(i) = 0.
        elseif (turn_off_cond_adv == 1) then 
          rvv = 1. + rv(i+1) + rl(i+1)
          rivert(i) = 0.
          rlvert(i) = -rl(i+1)*fmass(i+1)/rvv
        else
          rvv = 1. + rv(i+1) + rl(i+1) + ri(i+1)
          rivert(i) = -ri(i+1)*fmass(i+1)/rvv
          rlvert(i) = -rl(i+1)*fmass(i+1)/rvv
        endif
        dmvert(i) = -fmass(i+1)/rvv
        rvvert(i) = -rv(i+1)*fmass(i+1)/rvv
        hmvert(i) = -hm(i+1)*fmass(i+1)/rvv
        uvert(i) = -uwind(i+1)*fmass(i+1)/rvv
        vvert(i) = -vwind(i+1)*fmass(i+1)/rvv
      endif
!----------------------------------------------------------------------73
!  endif for fmass(i+1) > or < 0
!----------------------------------------------------------------------73
    endif
!----------------------------------------------------------------------73
!
!  Above endif completes adjustments in the vert quantities due to fluxes
!   at the top of layer i, fmass(i+1).
!
!  Next do adjustments in vert arrays due to fluxes through the bottom 
!   of layer i.
!
!  Upward flux at bottom of layer i from i-1: GAIN at layer i
!
!----------------------------------------------------------------------73
    if (fmass(i) >= 0.) then
!----------------------------------------------------------------------73
! not the surface layer.
!----------------------------------------------------------------------73
      if (i /= 1) then
        if (turn_off_cond_adv == 2) then 
          rvv = 1. + rv(i-1) 
          rlvert(i) = 0.
          rivert(i) = 0.
        elseif (turn_off_cond_adv == 1) then 
          rvv = 1. + rv(i-1) + rl(i-1)
          rivert(i) = 0.
          rlvert(i) = rlvert(i) + rl(i-1)*fmass(i)/rvv
        else
         rvv = 1. + rv(i-1) + rl(i-1) + ri(i-1)
         rivert(i) = rivert(i) + ri(i-1)*fmass(i)/rvv
         rlvert(i) = rlvert(i) + rl(i-1)*fmass(i)/rvv
        endif
        dmvert(i) = dmvert(i) + fmass(i)/rvv
        rvvert(i) = rvvert(i) + rv(i-1)*fmass(i)/rvv
        hmvert(i) = hmvert(i) + hm(i-1)*fmass(i)/rvv
        uvert(i) = uvert(i) + uwind(i-1)*fmass(i)/rvv
        vvert(i) = vvert(i) + vwind(i-1)*fmass(i)/rvv
!----------------------------------------------------------------------73
!  i = 1 surface layer
!  When precipitation occurs, for updrafts, the induced vertical mass flux
!  at the bottom of the surface layer will be upward. In this case, leave dmvert,
!  rvvert,rlvert,rivert,hmvert, etc unchanged. The i-1 quantities that would be needed
!  to define an upward flux are not available.
!----------------------------------------------------------------------73
      else
      endif
!----------------------------------------------------------------------73
!  fmass(i) < 0.
!  Downward flux at bottom of layer i to i-1: Loss from layer i
!  Note that fmass < 0 so the vert arrays are reduced.
!----------------------------------------------------------------------73
    else
      if (i /= 1) then
        if (turn_off_cond_adv == 2) then 
          rvv = 1. + rv(i) 
          rivert(i) = 0.
          rlvert(i) = 0.
        elseif (turn_off_cond_adv == 1) then 
          rvv = 1. + rv(i) + rl(i)
          rivert(i) = 0.
          rlvert(i) = rlvert(i) + rl(i)*fmass(i)/rvv
        else
          rvv = 1. + rv(i) + rl(i) + ri(i)
          rivert(i) = rivert(i) + ri(i)*fmass(i)/rvv
          rlvert(i) = rlvert(i) + rl(i)*fmass(i)/rvv
        endif
        dmvert(i) = dmvert(i) + fmass(i)/rvv
        rvvert(i) = rvvert(i) + rv(i)*fmass(i)/rvv
        hmvert(i) = hmvert(i) + hm(i)*fmass(i)/rvv
        uvert(i) = uvert(i) + uwind(i)*fmass(i)/rvv
        vvert(i) = vvert(i) + vwind(i)*fmass(i)/rvv
!----------------------------------------------------------------------73
!  i = 1:
!  There is a downward mass flux at the bottom of the surface layer. This
!    will occur when there is a net increase in column mass due to the sum
!    of all exchanges of mass at different levels.
!  This should not occur for updrafts, which produce precipitation and remove
!    mass from the column. For downdrafts however, which take mass from precip
!    and add to column, it should be expected. In this case, allow for an increase
!    in surface pressure so that these "exiting" fluxes are retained in the layer.
!    You do not have to add these fluxes as a source. Simply don't register
!    as a loss (i.e. do nothing). Implicitly assumes surface pressure will
!    increase.
!----------------------------------------------------------------------73
      else
!----------------------------------------------------------------------73
! endif for i/= 1
!----------------------------------------------------------------------73
      endif
!----------------------------------------------------------------------73
! endif for fmass(i) > 0 or < 0.
!----------------------------------------------------------------------73
    endif
!------------------------------------------------------------------------
!
!   End loop over nd
!
!------------------------------------------------------------------------
!     print*,'i rivert(i) ri(i) = ',i,rivert(i),ri(i)
  end do

!------------------------------------------------------------------------
!
!  Corrections to rvvert to avoid negative rv
!  - look for rvvert that would cause rv to become negative.
!  - negative rv causes lots of downstream problems
!  - In this case dmvert must also be adjusted to keep total mass
!    fluxes between layers constant.
!  - hmflow should also be adjusted
!  - rvvert is the net local change in rv due to vertical advection (from
!    above and below).
!  -
!  rvvert(nd)                     F  [kg vapor/m2/s]
!
!  ZZZZ
!
!------------------------------------------------------------------------
! do i = 1,nd
! end do

!------------------------------------------------------------------------
!
!  Diagnostics only: defining sums of hmvert,dmvert.
!  These vert arrays should sum to zero. 
!  Vertical advection should just move stuff around, not create it.
!
!------------------------------------------------------------------------
  sum_hmvert = 0.
  sum_dmvert = 0.
  sum_rlvert = 0.
  sum_rivert = 0.
  sum_rvvert = 0.
  do i = 1,nd
    sum_hmvert = sum_hmvert + hmvert(i)
    sum_dmvert = sum_dmvert + dmvert(i)
    sum_rvvert = sum_rvvert + rvvert(i)
    sum_rlvert = sum_rlvert + rlvert(i)
    sum_rivert = sum_rivert + rivert(i)
  end do
!------------------------------------------------------------------------
!
!  END  (calc_vert)
!
!------------------------------------------------------------------------
end subroutine calc_vert
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   subroutine conv_tendencies
!
!   Calculate tendencies from "vert" and "flow" arrays
!
!   Called ONCE for both updrafts+downdrafts
!
!   Basic state not updated here.
!
!  All p(i) levels are treated as fixed, but surface pressure may change.
!  When evaluating changes to basic state, need to know changes in surface
!    pressure from updrafts/downdrafts.
!
!  pnew and phalfnew are being used here only in a diagnostic sense to make sure that
!    the vertical advection is done in such a way as to keep the p(i) fixed.
!    Check for changes (errors).
!
!------------------------------------------------------------------------
subroutine conv_tendencies(land_fraction,tstep)
!------------------------------------------------------------------------
!
!   Use Statements
!   Aug 2010: erased all "if_conv_params" statements since should all be
!     at top of module?
!
!------------------------------------------------------------------------
use if_conv_solvers, only : g,dperr_max,enthalpy,fusion,latent,t_from_km
use if_conv_solvers, only : cl,latent,cpd,cpvmcl,epsi,rdd
use if_conv_params, only: km,z,p
use if_conv_params, only: rvtend,rltend,ritend,utend,vtend,kmtend
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: land_fraction
real(r8), intent(in) :: tstep
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: hmerr,dmold,hmm,tconvert,dz,ddv,ddu,dperr,dmnew,dp,dm
real(r8) :: hmtotnew,hmtotold,sum_hmall
real(r8) :: dhm,drv
real(r8) :: drl,dri,dag
real(r8) :: specadd_up,specadd_dn
real(r8) :: ddt,rtold,dmdry,sec_in_day
real(r8) :: pnew(nd),phalfnew(nd+1)
real(r8) :: dhmflow,dppnew,colhm_old,colhm_new,dcolhm,totevap
integer :: i,ii,nums,nsat
real(r8) :: dtt, km_g, hm_g, dzz, zzz, rho_g, tdpp, rtt_old
real(r8) :: flux, te_col_old, te_col_new, kkk, lff, lvv
real(r8) :: se_col_old, se_col_new, precip_energy, tot_surf_precip
real(r8) :: col_rv_new,col_rl_new,col_ri_new,col_wat_new
real(r8) :: col_rv_old,col_rl_old,col_ri_old,col_wat_old
real(r8) :: col_hm_new, col_hm_old, col_hm_expect, hm_rel_err
real(r8) :: rel_te_error, te_error_watts, te_col_xpd

integer :: niter_max, iter, nt, ngot
!------------------------------------------------------------------------
!
!   Local arrays
!
!------------------------------------------------------------------------
real(r8) :: rvnew(nd)
real(r8) :: rtnew(nd)
real(r8) :: rlnew(nd)
real(r8) :: rinew(nd)
real(r8) :: kmnew(nd)
real(r8) :: dpnew(nd)

!------------------------------------------------------------------------
!
! These are defined as purely local variables, since only tendencies 
!   used in driver. Same names as in driver, but hopefully not causing
!   a problem.
!
!------------------------------------------------------------------------
real(r8) :: uwindnew(nd),vwindnew(nd),hmnew(nd)
real(r8) :: rlnew_small,rvv
!------------------------------------------------------------------------
!
!   All vert quantities now defined for updrafts and downdrafts.
!
!   - dmold: old dry mass of layer
!   - dmnew: new dry mass of layer (all up and down motions + ent/det)
!   - rvnew: new rv of layer.
!
!------------------------------------------------------------------------
  tconvert = (3600.*24./tstep)
  hmtotold = 0.
  hmtotnew = 0.
  sum_hmall = 0.
  do i = 1,nd
!----------------------------------------------------------
!  dmnew calculation: drymass budgets
!----------------------------------------------------------
!    print*,'start of loop over nd i = ',i
    dmold = (phalf(i) - phalf(i+1))/(1.+rv(i)+rl(i)+ri(i))
    dmold = dmold/g
    dm = (dmflow(2,i) - dmflow(1,i) + dmvert(i))*tstep
    dmnew = dmold + dm
!--------------------------------------------------------------------------
!
!  Could have a check here to see if any of these mass flow terms are larger 
!    than dmold. Can get problems when this is the case. But hopefully
!    previous limits on updrafts/downdrafts have lowered this possibility
!
!--------------------------------------------------------------------------
!----------------------------------------------------------
!  rvnew calculation
!----------------------------------------------------------
    drv = rvflow(2,i) - rvflow(1,i)    &  ! rvflow may include rlp det
             + rvvert(i)               &  ! vertical transport
             - ansnow_rv(i)            &
             - uprain_rv(i)            &        
             + rl_evap(i)              & 
             + ri_evap(i)              &
             - rl_prod(i)              &
             - ri_prod(i)      
    drv = drv*tstep
    rvnew(i) = (rv(i)*dmold + drv)/dmnew
!----------------------------------------------------------
!
!  rlnew calculation
!
!----------------------------------------------------------
    drl =   rlflow(2,i) - rlflow(1,i)  &  ! 
            + rlvert(i)                &  ! vertical transport
            - rl_evap(i)               &      
            + rl_prod(i)               &      
            - uprain_rl(i)                      
    drl = drl*tstep
    rlnew(i) = (rl(i)*dmold + drl)/dmnew
!----------------------------------------------------------
!
!  rinew calculation
!  ICE
!
!----------------------------------------------------------
    dri =  riflow(2,i) - riflow(1,i)  &  ! 
         + rivert(i)                  &  ! vertical transport
         - ansnow_ri(i)               &  ! ri to ansnow
         - ri_evap(i)                 &  ! ri evap
         + ri_prod(i)                    ! ri evap
    dri = dri*tstep
    rinew(i) = (ri(i)*dmold + dri)/dmnew
!----------------------------------------------------------
!
!  Check for negative or too large change in rvnew:
!  - This is where a lot of the problems start ...
!    If you have negative rvnew, then you can end up with very high T
!    when tnew is later calculated from km.
!  Possible solutions:
!  (1) Just set rvnew = 0. This causes dperr problems and km conservation
!     problems.
!
!   Basically happens when a mass flux in/out of a layer is significantly
!    larger than the mass of the layer, and the rv of the mass leaving is
!    larger than the rv of the mass entering.
!    (i) when stratiform precipitation is very large, downdraft detrainment
!     into the lowest layer can be very large. (SHOULD BE FIXED)
!    (ii) large values of entrainment of updraft air parcels (mass_ratio)
!
!  Feb 2017: changed threshold to reduce error messages
!
!----------------------------------------------------------
!    if (rvnew(i) < 0.) then
    if (rvnew(i) < -1.0E-04) then
      print*,'==== WARNING: NEGATIVE RVNEW ======' 
      print*,'nstep = ',nstep
      print*,'May be associated with extremely strong vertical gradients'
      print*,'in rv combined with large vertical mass fluxes.'
      print*,'Originating problem may be large vertical motion'
      print*,'leading to large RH in the upper troposphere.'
      print*,'rvnew(i) = ',rvnew(i)
      print*,'rv(i) = ',rv(i)
      if (i < nd) print*,'rv(i+1) = ',rv(i+1)
      if (i /= 1) print*,'rv(i-1) = ',rv(i-1)
      print*,'rh(i) = ',rh(i)
      if (i < nd) print*,'rh(i+1) = ',rh(i+1)
      if (i /= 1) print*,'rh(i-1) = ',rh(i-1)
      print*,'abs(rvnew(i)-rv(i)) = ',abs(rvnew(i)-rv(i))
      print*,'tstep = ',tstep
      print*,'i z(i) = ',i,z(i)
      print*,'t(i) = ',t(i)
      print*,'p(i) = ',p(i)
      if (i < nd) print*,'phalf(i) phalf(i+1) = ',phalf(i),phalf(i+1)
      print*,'rv(i) times dmold/dmnew= ',rv(i)*(dmold/dmnew)
      print*,'Relative change in rv = ',(rvnew(i)-rv(i))/rv(i)
      print*,'dmnew = ',dmnew
      print*,'dmold = ',dmold
      print*,'original kg of water vapor mass in layer rv(i) times dmold = ',rv(i)*dmold
      print*,'************* drv = ',drv
      print*,'drv/dmnew = ',drv/dmnew
      print*,'source: up vapor detrainment rvflow(2,i) = ',rvflow(2,i)
      print*,'sink: up vapor entrainment rvflow(1,i) = ',rvflow(1,i)
      print*,'sink or source vertical advection rvvert(i) = ',rvvert(i)
      print*,'kg of vapor added via up rv detrain rvflow(2,i) = ',rvflow(2,i)*tstep
      print*,'kg of vapor added via up rlflow(2,i) = ',rlflow(2,i)*tstep
      print*,'kg of vapor lost via rv up entrain rvflow(1,i) = ',rvflow(1,i)*tstep
      print*,'kg of vapor lost via cond ansnow_rv(i) = ',ansnow_rv(i)*tstep
      print*,'kg of vapor lost via cond uprain_rv(i) = ',uprain_rv(i)*tstep
      print*,'kg of air added or lost via rvvert(i) = ',rvvert(i)*tstep
      print*,'kg of vapor added via dn rlflow(2,i) = ',rlflow(2,i)*tstep
      print*,'kg of vapor lost via rv dn entrain rvflow(1,i) = ',rvflow(1,i)*tstep
      print*,'********* UPDRAFT MASS FLUXES'
      print*,'negative fmass is downward'
      print*,'fmass(i) = ',fmass(i)
      print*,'rv(i) = ',rv(i)
      if (i < nd) print*,'fmass(i+1) = ',fmass(i+1)
      if (i < nd) print*,'rv(i+1) = ',rv(i+1)
      print*,'************* MASS FLUXES TIMES tstep to get kg/m2 ************'
      print*,'mass fluxes in a tstep larger than dmold could be the problem'
      print*,'dmold = ',dmold
      print*,'updraft fmass(i) times tstep = ',fmass(i)*tstep
      if (i < nd) print*,'updraft fmass(i+1) times tstep = ',fmass(i+1)*tstep
      print*,'*************** MASS DET/ENT TIMES tstep to get kg/m2 **************'
      if ((i == nd).or.(i == (nd-1))) then 
        print*,'dmflow has i+1 and i+2'
        stop 'will get segmentation fault from printing this stuff'
      endif
      print*,'kg of dry air removed via up ent dmflow(1,i) = ',dmflow(1,i)*tstep
      print*,'kg of dry air added via up det dmflow(2,i) = ',dmflow(2,i)*tstep
      print*,'************ RAINFALL in mm/day (mult by 24.3600.) ****************'
      print*,'mass_start_1 = ',mass_start_1
      print*,'mass_start_2 = ',mass_start_2
      print*,'km_width = ',km_width
      print*,'land_fraction = ',land_fraction
      print*,'ASSUME FIXED LATER'
!      stop 'rlnew negative in if_conv_tend'
    endif
!----------------------------------------------------------
!  
!  Check for negative rinew.  
!
!  - like rvnew, can be caused by large entrainment or vertical advection losses
!    associated with large mass fluxes or rain rates.
!  - can be fixed by decreasing mass_ratio_max or sigma_max, or by decreasing the
!    pressure level where downdrafts start to detrain
!  
!----------------------------------------------------------
    rlnew_small = 1.0E-05
    if (rinew(i) < -rlnew_small) then
      print*,'WARNING: == NEGATIVE rinew in if_conv_tend.f90 tendency calculation '
      print*,'nstep = ',nstep
      print*,'i z(i) = ',i,z(i)
      print*,'temperature t(i) = ',t(i)
      if (i < nd) print*,'phalf(i) phalf(i+1) = ',phalf(i),phalf(i+1)
      print*,'tstep = ',tstep
      print*,'col_ri = ',col_ri
      print*,'ri(i) rinew(i) = ',ri(i),rinew(i)
      if (i.ne.nd) print*,'ri(i+1) = ',ri(i+1)
      if (i.ne.1) print*,'ri(i-1) = ',ri(i-1)
      print*,'dmold = ',dmold
      print*,'dmnew = ',dmnew
      print*,'original kg of ri mass in layer: ri(i) times dmold = ',ri(i)*dmold
      print*,'************* ri UPDRAFT STUFF **************'
      print*,'dri = dri*tstep'
      print*,'dri = ',dri
      print*,'dri/dmnew = ',dri/dmnew
      print*,'source or sink: kg of ri from updraft rivert(i) = ',rivert(i)*tstep
      print*,'*************** MASS FLUXES TIMES tstep to get kg/m2 **************'
      print*,'mass fluxes in a tstep larger than dmold could be the problem'
      print*,'dmold = ',dmold
      print*,'updraft fmass(i) times tstep = ',fmass(i)*tstep
      if (i < nd) print*,'updraft fmass(i+1) times tstep = ',fmass(i+1)*tstep
      print*,'*************** MASS DET/ENT TIMES tstep to get kg/m2 **************'
      if ((i == nd).or.(i == (nd-1))) then 
        print*,'dmflow has i+1 and i+2'
        stop 'will get segmentation fault from printing this stuff'
      endif
      print*,'kg of dry air removed via up ent dmflow(1,i) = ',dmflow(1,i)*tstep
      print*,'kg of dry air added via up det dmflow(2,i) = ',dmflow(2,i)*tstep
      print*,'************ RAINFALL in mm/day (mult by 24.3600.) ***********'
      print*,'uprain_start_rlp = ',uprain_start_rlp*24.*3600.
      print*,'ansnow_conv = ',ansnow_conv*24.*3600.
      print*,'ansnow_strat = ',ansnow_strat*24.*3600.
      print*,'---------------------------------------------'
      print*,'mass_start_1 = ',mass_start_1
      print*,'mass_start_2 = ',mass_start_2
      print*,'km_width = ',km_width
      print*,'land_fraction = ',land_fraction
      print*,'ASSUME FIXED LATER'
!      stop 'rinew negative in if_conv_tend'
    endif
!----------------------------------------------------------
!  rtnew calculation
!----------------------------------------------------------
    rtnew(i) = rvnew(i) + rlnew(i) + rinew(i)
!------------------------------------------------------------------------
!
!  dpnew calculation:
!
!  - The change dp in pressure levels should be zero, except for the surface.
!  - Should be true by construction of the induced vertical transports in the 
!    background atm
!  - Check for mass consistency with previously defined pressure levels
!
!------------------------------------------------------------------------
    dpnew(i) = dmnew*(1. + rtnew(i))*g
    dp = phalf(i) - phalf(i+1)
    dperr = abs(dpnew(i) - dp)/dp
    if ((dperr > dperr_max).and.(i /= 1)) then
      print*,'========= entering dperr i = ',i
      print*,'i fractional dperr = ',i,dperr
      print*,'dpnew(i) = ',dpnew(i)
      print*,'dp = ',dp
      print*,'uprain_start_rlp = ',uprain_start_rlp
      print*,'ansnow_conv = ',ansnow_conv
      print*,'ansnow_strat = ',ansnow_strat
      print*,'Absolute pressure error (Pa) = ',dpnew(i) - dp
      print*,'updraft bottom fmass(i) = ',fmass(i)
      if (i < nd) print*,'updraft top fmass(i+1) = ',fmass(i+1)
      print*,'net updraft advection dmvert(i) = ',dmvert(i)
      print*,'entrainment updraft dmflow(1,i) = ',dmflow(1,i)
      print*,'detrainment updraft dmflow(2,i) = ',dmflow(2,i)
      print*,'downdraft bottom fmass(i) = ',fmass(i)
      if (i < nd) print*,'downdraft top fmass(i+1) = ',fmass(i+1)
      print*,'net updft rv advection rvvert(i) = ',rvvert(i)
      print*,'net updft rl advection rlvert(i) = ',rlvert(i)
      print*,'net updft ri advection rivert(i) = ',rivert(i)
      print*,'net updraft rv ent rvflow(1,i) = ',rvflow(1,i)
      print*,'net updraft rv det rvflow(2,i) = ',rvflow(2,i)
      print*,'net updraft rc det rlflow(2,i) = ',rlflow(2,i)
      print*,'i z(i) = ',i,0.001*z(i)
      print*,'rv(i) rvnew(i) = ',rv(i),rvnew(i)
      print*,'rl(i) rlnew(i) = ',rl(i),rlnew(i)
      print*,'ri(i) rinew(i) = ',ri(i),rinew(i)
      if (i < nd) print*,'ri(i+1) = ',ri(i+1)
      if (i > 1) print*,'ri(i-1) = ',ri(i-1)
      print*,'The total mass (dry air+water vapor+condensate) is'
      print*,'inconsistent with prescribed fixed pressure layers'
       print*,'rtnew(i) - rt(i) = ',rtnew(i)-rt(i)
       print*,'rvnew(i) - rv(i) = ',rvnew(i)-rv(i)
       print*,'rlnew(i) - rl(i) = ',rlnew(i)-rl(i)
       print*,'rinew(i) - ri(i) = ',rinew(i)-ri(i)
       print*,'dmvert(i) = ',dmvert(i)
       print*,'rivert(i) = ',rivert(i)
       print*,'riflow(1,i) = ',riflow(1,i)
       print*,'riflow(2,i) = ',riflow(2,i)
      stop 'dperr problem in conv_tendencies'
    endif
!------------------------------------------------------------------------
!
!   Find new momentum values
!
!------------------------------------------------------------------------
  ddu = (uflow(2,i) - uflow(1,i) + uvert(i))*tstep
  uwindnew(i) = (uwind(i)*dmold + ddu)/dmnew
  ddv = (vflow(2,i) - vflow(1,i) + vvert(i))*tstep
  vwindnew(i) = (vwind(i)*dmold + ddv)/dmnew
!------------------------------------------------------------------------
!
!   Calculate hmnew(i)
!    numhm = 1: hm entrainment
!    numhm = 2: hm detrainment
!
!------------------------------------------------------------------------
  dhm = (hmflow(2,i) - hmflow(1,i) + hmvert(i))*tstep
  hmnew(i) = (hm(i)*dmold + dhm)/dmnew
!------------------------------------------------------------------------
!
!  ***********   Diagnostics  ****************
!
!
!  sum_hmall : sum of all hm changes due to hmflow and hmvert from both
!         updrafts and downdrafts. Loss of hm in the column should
!         equal escape of hm from precipitation (hmuprain + hmansnow_surf).
!
!------------------------------------------------------------------------
  hmtotold = hmtotold + hm(i)*dmold
  hmtotnew = hmtotnew + hmnew(i)*dmnew
  sum_hmall = sum_hmall + dhm
!------------------------------------------------------------------------
!
!  end loop over nd
!
!------------------------------------------------------------------------
  end do
!------------------------------------------------------------------------
!
!  Diagnostics:
!
!------------------------------------------------------------------------
!  print*,'hmtotold = ',hmtotold
!  print*,'hmtotnew = ',hmtotnew
!  print*,'hmtotnew - hmtotold = ',hmtotnew-hmtotold
!  print*,'sum_hmall = ',sum_hmall
!  print*,'hmuprain*tstep = ',hmuprain*tstep
!------------------------------------------------------------------------
!
!   Define phalfnew:
!
!   Use dpnew to define new phalf levels 
!   Only used for diagnostics.
!   Keep top phalf(nd+1) fixed, so that precip reduces surface pressure.
!   Loop above should ensure that advection+convection does not change the
!      amount of mass in any grid box, except lowest layer.
!   p(i) and pnew(i) should be identical, subject to numerical error,
!      except that the surface pressures phalfnew(1) and phalf(2) should
!      be different.
!
!  What to do in cam4? Presumably want to keep pressure levels
!    constant. How to handle change in surface pressure?
!
!------------------------------------------------------------------------
! print*,'dpnew(1) = ',dpnew(1)
  phalfnew(nd+1) = phalf(nd+1)
  do i = nd,1,-1
    dp = phalf(i) - phalf(i+1)
    phalfnew(i) = phalfnew(i+1) + dpnew(i)
    pnew(i) = 0.5*(phalfnew(i) + phalfnew(i+1))
  end do
!------------------------------------------------------------------------
!
!  Check to make sure: dpsurf = phalfnew(1) - phalf(1)
!  - errors here could be in the calculation of dpsurf
!
!------------------------------------------------------------------------
   dperr = dpsurf - (phalfnew(1) - phalf(1))
!   print*,'absolute Pa pressure dperr at surface = ',dperr
   if (abs(dperr) > 1.0E-06) then
     print*,'----------------------------------------------------------'
     print*,'absolute pressure error dperr at surface = ',dperr
     print*,'phalf(1) = ',phalf(1)
     print*,'phalfnew(1) = ',phalfnew(1)
     print*,'phalfnew(1) - phalf(1) = ',phalfnew(1)-phalf(1)
     print*,'dpsurf = ',dpsurf
     print*,'uprain_start_rlp = ',uprain_start_rlp
     print*,'ansnow_conv = ',ansnow_conv
     print*,'ansnow_strat = ',ansnow_strat
     do i = 1,nd
     if (i < 25) then
       print*,'------------------------------------ i = ',i
       print*,'phalfnew(i) - phalf(i) = ',phalfnew(i) - phalf(i)
       print*,'rtnew(i) - rt(i) = ',rtnew(i)-rt(i)
       print*,'rvnew(i) - rv(i) = ',rvnew(i)-rv(i)
       print*,'rlnew(i) - rl(i) = ',rlnew(i)-rl(i)
       print*,'rinew(i) - ri(i) = ',rinew(i)-ri(i)
       print*,'dmvert(i) = ',dmvert(i)
     endif
     end do
     stop 'dperr problem in conv_tend.f90'
   endif
!------------------------------------------------------------------------
!
!  Diagnostics:
!
!  Find change in column hm 
!  (at this point equal to change in column enthalpy).
!
!  Check to see if sum of changes in hm consistent with hmflow.
!    The hmnew here are not the GPE adjusted final hm. The sum of
!    all horizontal hm inflows/outflows should equal the change in
!    column hm (individually for updrafts/downdrafts).
!
!  If now, it means that there has been some mistake in calculating
!    the hmnew from the hmflows. The hmflows themselves could still
!    be in error.
!
!------------------------------------------------------------------------
  dhmflow = 0.
  colhm_new = 0.
  colhm_old = 0.
  do i = 1,nd
    dhmflow = dhmflow + (hmflow(2,i) - hmflow(1,i))*tstep
    dmnew = (phalfnew(i) - phalfnew(i+1))/(1.+rvnew(i)+rlnew(i))
    dmnew = abs(dmnew/g)
    colhm_new = colhm_new + dmold*hmnew(i)
    dmold = (phalf(i) - phalf(i+1))/(1.+rv(i)+rl(i)+ri(i))
    dmold = abs(dmold/g)
    colhm_old = colhm_old + dmold*hm(i)
  end do
  dcolhm = colhm_new - colhm_old
!  print*,'change in hm column updrafts and downdraftsdcolhm = ',dcolhm

!------------------------------------------------------------------------
!
!  *********    Calculate kmnew  *****************
!
!  The purpose of knowing kmnew is get tnew and solve for the heat() 
!    tendency in if_conv.F90. This is actually a dse tendency, used to
!    update the dse in physics_types.F90.
!    Then there is a "call geopotential_dse" to get the temperature
!     and new gpht's. Not sure how they use the dse tendency to get the
!     mix of T and z response. Must use hydrostatic balance.
!
!   in cam4: 
!   dse = cpair*t + gravit*z
!   dse_new = cpair*tnew + gravit*znew
!   dse_new - dse = cpair(tnew - t) + gravit*(znew - z)
!
!   for me:
!   hm = km + (1. + rt)*z*g             (the z is the cam4 z not mine)
!   hmnew = kmnew + (1. + rtnew)*znew*g (pretty sure would need new z here)
!
!   The problem is that znew is not something I can calculate. (V IMP!!)
!   However, suppose we define an ieffective "TNEW" 
!   defined in terms of the dse tendency. 
!
!   dse_new - dse = cpair(TNEW - t) = cpair(tnew - t) + gravit*(znew - z)  ***
!
!   We can also define an effective "KMNEW" that is defined in terms of TNEW,
!     rvnew, etc.
!
!   KMNEW = (cpd + rtnew*cl)*TNEW + lv*rvnew - lf*rinew
!
!   If we knew what KMNEW was, we could use get_t to get TNEW, and then
!     use this in the cam4 heating (dse) tendency.
!
!   From above ***
!   TNEW - t = tnew - t + gravit*(znew - z)/cpair
!   TNEW = tnew + gravit*(znew - z)/cpair  (i.e. like a dse)
!
!   Setting into KNEW:
!
!   KMNEW = (cpd + rtnew*cl)*[tnew + gravit*(znew - z)/cpair] + lv*rvnew - lf*rinew
!         = (cpd + rtnew*cl)*tnew + (cpd + rtnew*cl)*gravit*(znew - z)/cpair + lv*rvnew - lf*rinew
!         = knew + (cpd + rtnew*cl)*gravit*(znew - z)/cpair
!         ~ knew + g*(znew-z) + (cl/cpair)*rtnew*g*(znew - z)
!
!   Using:
!   hmnew = kmnew + (1. + rtnew)*znew*g   (though don't know znew)
!   kmnew = hmnew - (1. + rtnew)*znew*g   (though don't know znew)
!
!
!   KMNEW = hmnew - (1. + rtnew)*znew*g + g*(znew-z) + (cl/cpair)*rtnew*g*(znew - z)
!         = hmnew - znew*g - rtnew*znew*g + g*znew - g*z + (cl/cpair)*rtnew*g*(znew - z)
!         = hmnew - rtnew*znew*g - g*z + (cl/cpair)*rtnew*g*(znew - z)
!
!   In this expression, after hmnew, g*z >> rtnew*znew*g and znew ~ z, so can write:
!
!   KMNEW ~ hmnew - rtnew*z*g - g*z + (cl/cpair)*rtnew*g*(znew - z)
!         = hmnew - (1. + rtnew)*z*g + (cl/cpair)*rtnew*g*(znew - z)
!
!   The problem now is don't know znew - z. But it should be small.
!
!   (cl/cpair)*rtnew*g*(znew - z)     (cl/cpair)*rtnew*(znew-z)    4*rtnew*(znew-z)
!   -----------------------------  ~  -------------------------  ~ ----------------
!       (1. + rtnew)*z*g                       z                        z
!
!   For small enough timestep, (znew-z)/z << 1 except near the surface, and 4*rtnew
!   less than 0.1
!
!   Therefore can approximately write:
!
!   *********    KMNEW ~ hmnew - (1. + rtnew)*z*g    *************
!
!   I think this is what I should use.
!
!   Have hmnew(i),rvnew(i),rlnew(i),rinew(i): how to define tnew(i)?
!
!   Note: in July 2015, I tried using the hmnew to self consistently
!    calculate znew and tnew. You can do this either by integrating up
!    from the top of bottom. Both gave wierd results for the convective
!    heating profile. The problem is, I have no way of affecting the 
!    z tendency. I am only calculating an effective dse tedency.
!
!------------------------------------------------------------------------
  do i = 1,nd
!------------------------------------------------------------------------
!
!  define dmnew/dmold  (may be OBS)
!
!------------------------------------------------------------------------
    dmnew = (phalfnew(i)-phalfnew(i+1))/(1.+rvnew(i)+rlnew(i)+rinew(i))
    dmnew = abs(dmnew/g)
    dmold = (phalf(i)-phalf(i+1))/(1.+rv(i)+rl(i)+ri(i))
    dmold = abs(dmold/g)
!------------------------------------------------------------------------
!
!  Define kmnew
!
!------------------------------------------------------------------------
    if (use_new_kmnew == 1) then
      kmnew(i) = hmnew(i) - (1. + rtnew(i))*z(i)*g
    else
      kmnew(i) = hmnew(i) - (1.+rv(i)+rl(i)+ri(i))*g*z(i)*(dmold/dmnew)
    endif
!------------------------------------------------------------------------
!
!  Check for bad values of kmnew
!
!  These diagnostics were an attempt to understand why I was getting large
!  values of knew in the upper troposphere, and therefore, extremly high
!  temperatures. It happens when I have a combination of a sharp increase
!  in hm between layer i and i+1, very large descending mass fluxes, and very
!  large detrainment into layer i. For example, if the subsidence
!  mass fluxes are 5 times bigger or so that the mass of layer i, and there
!  is detrainment into layer i comparable with the mass of i, then subsidence
!  net removal of hm will not offset the input of hm from detrainment, and
!  very high hm can result.
!  - How is this to be avoided. The basic problem occurs when you have very
!  high hm in the boundary layer (e.g. 375 kJ/kg in bottom 3 levels), which
!  with quasi undilute ascent with mean detrainment at layers near the
!  tropopause where the layer have too little mass to buffer the large mass
!  fluxes and detrainment rates. One way is to start forcing detrainment lower,
!  by not allowing additional ascent if it exceeds the amount of atmospheric mass
!  between the current level and (e.g.) 100 hPa.
!
!------------------------------------------------------------------------
    if (kmnew(i) > km_bad) then
      print*,'******  kmnew exceeds km_bad  ***************'
      print*,'kmnew(i) = ',kmnew(i)
      print*,'km_bad = ',km_bad
      print*,'p(i) = ',p(i)
      print*,'z(i) = ',z(i)
      print*,'t(i) = ',t(i)
      if (i < nd) print*,'t(i+1) = ',t(i+1)
      if (i > 1) print*,'t(i-1) = ',t(i-1)
      print*,'rvnew(i) = ',rvnew(i)
      print*,'rv(i) = ',rv(i)
      if (i < nd) print*,'rv(i+1) = ',rv(i+1)
      if (i > 1) print*,'rv(i-1) = ',rv(i-1)
      print*,'rl(i) = ',rl(i)
      print*,'rlnew(i) = ',rlnew(i)
      print*,'ri(i) = ',ri(i)
      print*,'rinew(i) = ',rinew(i)
      print*,'rt(i) = ',rt(i)
      print*,'rtnew(i) = ',rtnew(i)
      print*,'km(i) = ',km(i)
      print*,'kmnew(i) = ',kmnew(i)
      if (i < nd) print*,'hm(i+1) = ',hm(i+1)
      print*,'hm(i) = ',hm(i)
      if (i > 1) print*,'hm(i-1) = ',hm(i-1)
      if (i < nd) print*,'hmnew(i+1) = ',hmnew(i+1)
      print*,'hmnew(i) = ',hmnew(i)
      print*,'hm(i) = ',hm(i)
      if (i > 1) print*,'hmnew(i-1) = ',hmnew(i-1)
      print*,'************ RAINFALL in mm/day (mult by 24.3600.) ***********'
      print*,'uprain_start_rlp = ',uprain_start_rlp*24.*3600.
      print*,'ansnow_conv = ',ansnow_conv*24.*3600.
      print*,'ansnow_strat = ',ansnow_strat*24.*3600.
      print*,'******************************************'
      print*,'entrainment dmflow(1,i) = ',dmflow(1,i)
      print*,'detrainment dmflow(2,i) = ',dmflow(2,i)
      print*,'dhm = (hmflow(2,i) - hmflow(1,i) + hmvert(i))*tstep'
      print*,'hm entrainment hmflow(1,i) = ',hmflow(1,i)
      print*,'hm detrainment hmflow(2,i) = ',hmflow(2,i)
      print*,'hmvert(i) = ',hmvert(i)
      print*,'entrainment hmflow(1,i)*tstep = ',hmflow(1,i)*tstep
      print*,'detrainment hmflow(2,i)*tstep = ',hmflow(2,i)*tstep
      print*,'hmvert(i)*tstep = ',hmvert(i)*tstep
      dmold = (phalf(i) - phalf(i+1))/(1.+rv(i)+rl(i)+ri(i))
      dmold = dmold/g
      dm = (dmflow(2,i) - dmflow(1,i) + dmvert(i))*tstep
      dmnew = dmold + dm
      dhm = (hmflow(2,i) - hmflow(1,i) + hmvert(i))*tstep
      print*,'dmold = ',dmold
      print*,'dm = ',dm
      print*,'dmnew = ',dmnew
      print*,'dhm = ',dhm
      print*,'hm(i)*dmold = ',hm(i)*dmold
      print*,'This is how hmnew is calculated:'
      print*,'hmnew(i) = (hm(i)*dmold + dhm)/dmnew'
      print*,'dhm/dmnew = ',dhm/dmnew
      print*,'hm(i)*dmold/dmnew = ',hm(i)*dmold/dmnew
      print*,'hmp_mean(i) = ',hmp_mean(i)
      print*,'-------- hmvert calculation ---------'
      if (i > 1) print*,'fmass(i-1) = ',fmass(i-1)
      print*,'fmass(i) = ',fmass(i)
      print*,'fmass(i+1) = ',fmass(i+1)
      if (i > 1) print*,'fmass(i-1)*tstep = ',fmass(i-1)*tstep
      print*,'fmass(i)*tstep = ',fmass(i)*tstep
      print*,'fmass(i+1)*tstep = ',fmass(i+1)*tstep
      print*,'Assume fmass(i+1) < 0. (default case)'
      print*,'This is downward transport of hm from upper layer to i'
      if (i < nd) then
        print*,'Then first contribution hmvert(i) = -hm(i+1)*fmass(i+1)/rvv'
        rvv = 1.+rv(i+1)
        print*,'then first contribution = ',-hm(i+1)*fmass(i+1)/rvv
      endif
      print*,'This would be a positive source of hm'
      print*,'Assume also fmass(i) < 0. (default case)'
      print*,'Downward transport of hm of layer i to i-1'
      print*,'Then hmvert(i) = hmvert(i) + hm(i)*fmass(i)/rvv'
      if (i < nd) then
        rvv = 1.+rv(i+1)
        print*,'Second contribution = ',hm(i)*fmass(i)/rvv
      endif
      print*,'This should be a loss of hm and reduce hmvert'
      print*,'Check that these two terms are consistent with hmvert'
      print*,'-------- surface hm  ---------'
      print*,'hm(1) = ',hm(1)
      print*,'hm(2) = ',hm(2)
      print*,'hm(3) = ',hm(3)
      print*,'hm(4) = ',hm(4)
      print*,'Fractions of air being replaced by ent/det:'
      print*,'ent dmflow(1,i)*tstep/dmold = ',dmflow(1,i)*tstep/dmold
      print*,'det dmflow(2,i)*tstep/dmold = ',dmflow(2,i)*tstep/dmold
      print*,'hmflow has units J/m2*s'
      print*,'dmflow has units kg/m2*s'
      print*,'hmflow/dmflow should have units J/kg'
      print*,'entraining hm = ',hmflow(1,i)/dmflow(1,i)
      print*,'Detraining hm = ',hmflow(2,i)/dmflow(2,i)
      print*,'entrain dmflow(1,i) = ',dmflow(1,i)
      print*,'detrain dmflow(2,i) = ',dmflow(2,i)
      print*,'entrain rvflow(1,i) = ',rvflow(1,i)
      print*,'detrain rvflow(2,i) = ',rvflow(2,i)
      print*,'entrain rlflow(1,i) = ',rlflow(1,i)
      print*,'detrain rlflow(2,i) = ',rlflow(2,i)
      print*,'entrain riflow(1,i) = ',riflow(1,i)
      print*,'detrain riflow(2,i) = ',riflow(2,i)
      print*,'------ downdrafts -------'
      print*,'dndet(i) = ',dndet(i)
      print*,'dndet(i)*tstep = ',dndet(i)*tstep
      print*,'dnent(i) = ',dnent(i)
      print*,'dnent(i)*tstep = ',dnent(i)*tstep
      print*,'Compare with dmold = ',dmold
      print*,'-------------------------'
      print*,'mass_start_1 = ',mass_start_1
      print*,'mass_start_2 = ',mass_start_2
      print*,'km_width = ',km_width
      print*,'land_fraction = ',land_fraction
      print*,'conv_energy = ',conv_energy
      print*,'cape_1 = ',cape_1
      print*,'cape_2 = ',cape_2
      print*,'cape_3 = ',cape_3
      print*,'cape_4 = ',cape_4
      stop 'kmnew too big: would generate large tnew'
    endif
!------------------------------------------------------------------------
!
!  Check for large values of temperature
!
!------------------------------------------------------------------------
    get_t_call = 2
    if (kmnew(i) > 400000.) then
      print*,'-------  Large value of kmnew ----'
      print*,'Calling t_from_km for kmnew(i) = ',kmnew(i)
      print*,'conv_energy = ',conv_energy
    endif
    call t_from_km(kmnew(i),p(i),rvnew(i),rinew(i),rlnew(i),xxx,get_t_call)
    if (xxx > t_bad) then
      print*,'******  T from kmnew TOO HOT in if_conv_tend.f90 ***************'
      print*,'New temperature = ',xxx
      print*,'t_bad = ',t_bad
      print*,'p(i) = ',p(i)
      print*,'z(i) = ',z(i)
      print*,'t(i) = ',t(i)
      print*,'rh(i) = ',rh(i)
      print*,'rv(i) = ',rv(i)
      print*,'rvnew(i) = ',rvnew(i)
      print*,'rl(i) = ',rl(i)
      print*,'rlnew(i) = ',rlnew(i)
      print*,'ri(i) = ',ri(i)
      print*,'rinew(i) = ',rinew(i)
      print*,'rt(i) = ',rt(i)
      print*,'rtnew(i) = ',rtnew(i)
      print*,'km(i) = ',km(i)
      print*,'kmnew(i) = ',kmnew(i)
      print*,'hm(i) = ',hm(i)
      print*,'hmnew(i) = ',hmnew(i)
      if (i < nd) print*,'hm(i+1) = ',hm(i+1)
      if (i > 1) print*,'hm(i-1) = ',hm(i-1)
      print*,'************ RAINFALL in mm/day (mult by 24.3600.) ***********'
      print*,'uprain_start_rlp = ',uprain_start_rlp*24.*3600.
      print*,'ansnow_conv = ',ansnow_conv*24.*3600.
      print*,'ansnow_strat = ',ansnow_strat*24.*3600.
      print*,'******************************************'
      print*,'dhm = (hmflow(2,i) - hmflow(1,i) + hmvert(i))*tstep'
      print*,'hm detrainment hmflow(2,i) = ',hmflow(2,i)
      print*,'hm entrainment hmflow(1,i) = ',hmflow(2,i)
      print*,'hmvert(i) = ',hmvert(i)
      print*,'hmflow(2,i)*tstep = ',hmflow(2,i)*tstep
      print*,'hmflow(1,i)*tstep = ',hmflow(2,i)*tstep
      print*,'hmvert(i)*tstep = ',hmvert(i)*tstep
      print*,'mass_start_1 = ',mass_start_1
      print*,'mass_start_2 = ',mass_start_2
      print*,'km_width = ',km_width
      print*,'land_fraction = ',land_fraction
      print*,'conv_energy = ',conv_energy
      print*,'cape_1 = ',cape_1
      print*,'cape_2 = ',cape_2
      print*,'cape_3 = ',cape_3
      print*,'cape_4 = ',cape_4
      print*,'t(1) = ',t(1)
      print*,'t(2) = ',t(2)
      print*,'t(3) = ',t(3)
      print*,'t(4) = ',t(4)
      print*,'t(5) = ',t(5)
      print*,'rh(1) = ',rh(1)
      print*,'rh(2) = ',rh(2)
      print*,'rh(3) = ',rh(3)
      print*,'rh(4) = ',rh(4)
      print*,'rh(5) = ',rh(5)
      stop 'New temperature exceeds t_bad '
    endif
!------------------------------------------------------------------------
!
!   Calculate tendencies
!    - convert to /day
!    - Note that at this point the new temperatures have not been calculated
!
!------------------------------------------------------------------------
    kmtend(i) = (kmnew(i) - km(i))*tconvert
    rvtend(i) = (rvnew(i) - rv(i))*tconvert
    rltend(i) = (rlnew(i) - rl(i))*tconvert
    ritend(i) = (rinew(i) - ri(i))*tconvert
    utend(i) = (uwindnew(i) - uwind(i))*tconvert
    vtend(i) = (vwindnew(i) - vwind(i))*tconvert
!------------------------------------------------------------------------
!
!  Diagnostic tendencies for condensate
!  - not needed for calculational purposes
!  - convert to kg/kg/day
!  - units of rl_evap are [kg water/m2*s]
!  - Of the rlflow/riflow variables, likely only rlflow(2,i)/riflow(2,i)
!    are non-zero, which are due to updraft detrainment.
! 
!------------------------------------------------------------------------
    sec_in_day = 3600.*24.
    uprain_rl_tend(i) = -sec_in_day*uprain_rl(i)/dmnew   
    rl_evap_tend(i) = -sec_in_day*rl_evap(i)/dmnew   
    ri_evap_tend(i) = -sec_in_day*ri_evap(i)/dmnew  
    ri_prod_tend(i) = sec_in_day*ri_prod(i)/dmnew  
    rl_prod_tend(i) = sec_in_day*rl_prod(i)/dmnew  
    rl_det_tend(i) = rlflow(2,i) - rlflow(1,i) 
    rl_det_tend(i) = sec_in_day*rl_det_tend(i)/dmnew
    ri_det_tend(i) = riflow(2,i) - riflow(1,i) 
    ri_det_tend(i) = sec_in_day*ri_det_tend(i)/dmnew
    rl_vert_tend(i) = rlvert(i) + rlvert(i)
    rl_vert_tend(i) = sec_in_day*rl_vert_tend(i)/dmnew
    ri_vert_tend(i) = rivert(i)
    ri_vert_tend(i) = sec_in_day*ri_vert_tend(i)/dmnew
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
  end do

!------------------------------------------------------------------------
!
!  *******  Check for column hm conservation *****
!
!    hm = (cpd + rt*cl)*t + lv*rv - lf*ri + (1. + rtt)*zzz*g
!
!------------------------------------------------------------------------
col_hm_new = 0.
col_hm_old = 0.
do i = 1,nd
  dmnew = (phalfnew(i)-phalfnew(i+1))/(1.+rvnew(i)+rlnew(i)+rinew(i))
  dmnew = dmnew/g
  dmdry = (phalf(i)-phalf(i+1))/(1.+rv(i)+rl(i)+ri(i))
  dmdry = dmdry/g
  col_hm_new = col_hm_new + hmnew(i)*dmnew
  col_hm_old = col_hm_old + hm(i)*dmdry
end do
col_hm_expect = col_hm_old - (hmuprain_surf+hmansnow_surf+hmanrain_surf)*tstep

if (abs(col_hm_old-col_hm_new) > 10000.) then
  hm_rel_err = abs(col_hm_expect-col_hm_new)/(col_hm_old - col_hm_new)
! print*,'relative col hm error hm_rel_err = ',hm_rel_err
  if (abs(hm_rel_err) > 0.2) then
    print*,'---------  col hm error in if_conv_tend.f90 ----'
    print*,'relative col hm error hm_rel_err = ',hm_rel_err
    print*,'This is relative to col_hm_old - col_hm_new'
    print*,'Absolute Error in column hm = ',col_hm_expect-col_hm_new
    print*,'Absolute Error in watts = ',(col_hm_expect-col_hm_new)/tstep
    print*,'col_hm_new = ',col_hm_new
    print*,'col_hm_old = ',col_hm_old
    print*,'col_hm_expect = ',col_hm_expect
    print*,'Change in col_hm = ',col_hm_new-col_hm_old
    print*,'Converted to watts = ',(col_hm_new-col_hm_old)/tstep
    print*,'Loss of hm from precip = ',(hmuprain_surf+hmansnow_surf+hmanrain_surf)*tstep
    print*,'Loss of hm from hmuprain_surf = ',hmuprain_surf*tstep
    print*,'Loss of hm from hmansnow_surf = ',hmansnow_surf*tstep
    print*,'Loss of hm from hmanrain_surf = ',hmanrain_surf*tstep
    print*,'------  Rain  ---------'
    print*,'anrain_surf (mm/day)= ',anrain_surf*3600.*24.
    print*,'uprain_surf_rv (mm/day) = ',uprain_surf_rv*3600.*24.
    print*,'uprain_surf_rlp (mm/day) = ',uprain_surf_rlp*3600.*24.
    print*,'ansnow_surf (mm/day) = ',ansnow_surf*3600.*24.
      print*,'mass_start_1 = ',mass_start_1
      print*,'mass_start_2 = ',mass_start_2
      print*,'km_width = ',km_width
      print*,'land_fraction = ',land_fraction
      print*,'conv_energy = ',conv_energy
!   stop 'hm relative error exceeds threshold'
    print*,'---------------------------'
  endif
endif

!------------------------------------------------------------------------
!
!  *******  Check for energy conservation in cam4 definition *****
!
!  Problem here: tnew and znew not knowable by me!!!
!  - I can't calculate dse. However I could do this in if_conv using my tnew
!    there as a surrogate for a dse tendency
!
!  I enforce conservation of column hm.
!  However should likely agree with the cam4 definition to within 1-2 percent.
!
!  cam4 calculates column energy change from vertical integral of:
!
!  TE = cpd*T + g*z + KE + (lv0+lf0)*wv + lf0*wl
!
!  must be consistent with the fluxes:
!
!  TF = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._r8*latice
!
!  latice here is lf
!  flx_cnd refers to total rain+snow
!  flx_vap refers to evaporation from surafce (is zero during convection)
!  Therefore:
!
!  TF = - (flx_cnd(i) - flx_ice(i))*latice   (DEPENDS ON SURFACE RAIN ONLY)
!
!  How to make sense:
!  Suppose have snow from vapor:
!  In this case flx_cnd = flx_ice and TF equals zero.
!  Decrease in TE from decrease in column wv must be balanced by increase
!    in column cpd*T + g*z. Ignoring KE:
!
!  TE = cpd*T + g*z + (lv0+lf0)*wv + lf0*wl
!
!
!------------------------------------------------------------------------
te_col_new = 0.
se_col_new = 0.
col_rv_new = 0.
col_rl_new = 0.
col_ri_new = 0.
do i = 1,nd
  dmnew = (phalfnew(i)-phalfnew(i+1))/(1.+rvnew(i)+rlnew(i)+rinew(i))
  dmnew = dmnew/g
  lvv = latent(t(i))
  lff = fusion(t(i))
  kkk = cpd*t(i) + g*z(i) + (lvv+lff)*rv(i) + lff*rl(i)
  te_col_new = te_col_new + kkk*dmnew
  col_rv_new = col_rv_new + rvnew(i)*dmnew
  col_rl_new = col_rl_new + rlnew(i)*dmnew
  col_ri_new = col_ri_new + rinew(i)*dmnew
end do
col_wat_new = col_rv_new + col_rl_new + col_ri_new

te_col_old = 0.
do i = 1,nd
  dmdry = (phalf(i)-phalf(i+1))/(1.+rv(i)+rl(i)+ri(i))
  dmdry = dmdry/g
  lvv = latent(t(i))
  lff = fusion(t(i))
  kkk = cpd*t(i) + g*z(i) + (lvv+lff)*rv(i) + lff*rl(i)
  te_col_old = te_col_old + kkk*dmdry
  se_col_old = se_col_old + (cpd*t(i) + g*z(i))*dmdry
  col_rv_old = col_rv_old + rv(i)*dmdry
  col_rl_old = col_rl_old + rl(i)*dmdry
  col_ri_old = col_ri_old + ri(i)*dmdry
end do
col_wat_old = col_rv_old + col_rl_old + col_ri_old

flux = -(anrain_surf+uprain_surf_rv+uprain_surf_rlp)*fusion(t(1))
flux = tstep*flux


!------------------------------------------------------------------------
!
!  Define the expected new column energy:
!  - the old plus (surface rain - surface snow)*lf*tstep
!  - from the check_energy definition, assuming flx_vap = flx_sen = 0.
!  - (flx_cnd(i) - flx_ice(i)) is the net surface rain, since flx_cnd is rain+snow
!
!------------------------------------------------------------------------
te_col_xpd = te_col_old - (anrain_surf+uprain_surf_rv+uprain_surf_rlp)*fusion(t(1))*tstep
te_error_watts = (te_col_xpd - te_col_new)/tstep

precip_energy = (anrain_surf+uprain_surf_rv+uprain_surf_rlp)*latent(t(1)) + &
                ansnow_surf*(latent(t(1)) + fusion(t(1)))

tot_surf_precip = anrain_surf + uprain_surf_rv + uprain_surf_rlp + ansnow_surf 


if (precip_energy > 10.) then
  rel_te_error = te_error_watts/precip_energy   ! Both of these should be in watts
else
  rel_te_error = 0.
endif

!print*,'precip_energy (watts) = ',precip_energy
!print*,'Rate of increase in DSE (watts) = ',(se_col_new-se_col_old)/tstep
!print*,'-------'

if (nprint_energy == 1) then
!if (rel_te_error > 0.02) then
print*,'***********  CAM4 ENERGY DEFINITION   **********'
print*,'Relative te error (normalized by precip_energy) = ',rel_te_error
print*,'te_error_watts (watts) = ',te_error_watts
print*,'precip_energy (watts) = ',precip_energy
print*,'------  Rain  ---------'
print*,'anrain_surf (mm/day)= ',anrain_surf*3600.*24.
print*,'uprain_surf_rv (mm/day) = ',uprain_surf_rv*3600.*24.
print*,'uprain_surf_rlp (mm/day) = ',uprain_surf_rlp*3600.*24.
print*,'ansnow_surf (mm/day) = ',ansnow_surf*3600.*24.
print*,'total surf precip (mm/day) = ',tot_surf_precip*3600.*24.
print*,'------  Column Energy  ---------'
print*,'te_col_new = ',te_col_new
print*,'te_col_old = ',te_col_old
print*,'te_col_xpd = ',te_col_xpd
print*,'Change te_col_new - te_col_old = ',te_col_new-te_col_old
print*,'Converted to watts = ',(te_col_new-te_col_old)/tstep
print*,'te_error in watts = ',te_error_watts
print*,'------ Column Dry Static Energy ----'
print*,'se_col_new = ',se_col_new
print*,'se_col_old = ',se_col_old
print*,'Change se_col_new - se_col_old = ',se_col_new-se_col_old
print*,'Converted to watts = ',(se_col_new-se_col_old)/tstep
print*,'--------  Flux   ---------'
print*,'flux (surface rain times lf) = ',flux
print*,'Converted to watts flux = ',flux/tstep
print*,'Fractional error using cam4 definition = ',(te_col_new-te_col_old)/flux
print*,'precip_energy (surface precip times lv in watts) = ',precip_energy
print*,'---------  Surface Pressure  -----'
print*,'Change in surface pressure phalfnew(1) - phalf(1) = ',phalfnew(1)-phalf(1)
print*,'Change in column mass = ',(phalfnew(1)-phalf(1))/g
print*,'------- Water Budget ---------'
print*,'Change in column vapor col_rv_new-col_rv_old = ',col_rv_new-col_rv_old
print*,'Change in column liq col_rl_new-col_rl_old = ',col_rl_new-col_rl_old
print*,'Change in column ice col_ri_new-col_ri_old = ',col_ri_new-col_ri_old
print*,'Change in column water col_wat_new-col_wat_old = ',col_wat_new-col_wat_old
print*,'Change in column water converted to precip (mm/day) = ',3600.*24.*((-col_wat_new+col_wat_old)/tstep)
print*,'Should be consistent with tot_precip_surf'
print*,'tot_surf_precip (mm/day) = ',tot_surf_precip*3600.*24.
print*,'--------------'
print*,'mass_start_1 = ',mass_start_1
print*,'mass_start_2 = ',mass_start_2
print*,'km_width = ',km_width
print*,'land_fraction = ',land_fraction
print*,'conv_energy = ',conv_energy
!endif
endif

!------------------------------------------------------------------------
!
!  Compute col_rl_new
!  - not clear what it is used for, likely diagnostics.
!
!------------------------------------------------------------------------
  col_rl_new = 0.
  do i = 1,nd
    dmdry = (phalf(i)-phalf(i+1))/(g*(1.+rvnew(i)+rlnew(i)+ri(i)))
    dmdry = dmdry
    if (dmdry < 0.) stop 'dmdry negative in subroutine conv_tendencies'
    col_rl_new = col_rl_new + dmdry*rlnew(i)
  end do
!------------------------------------------------------------------------
!
!  Determine colr_rl_prd: col_rl production
!
!------------------------------------------------------------------------
  col_rl_prod = 3600.*24.*(col_rl_new - col_rl)/tstep
!------------------------------------------------------------------------
!
!   END   (conv_tendencies)
!
!------------------------------------------------------------------------
end subroutine conv_tendencies
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
subroutine calc_eff
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: tot_rv_entrain 
real(r8) :: tot_rc_entrain
real(r8) :: tot_rv_detrain
real(r8) :: tot_rc_detrain
real(r8) :: tot_rv_precip
real(r8) :: tot_uprain 
real(r8) :: tot_ansnow
real(r8) :: tot_entrain
integer :: i
!------------------------------------------------------------------------
!
!  For diagnostics
!
!------------------------------------------------------------------------
tot_rv_entrain = 0.
tot_rc_entrain = 0.
tot_rv_detrain = 0.
tot_rc_detrain = 0.
tot_rv_precip = 0.
tot_uprain = 0.
tot_ansnow = 0.
!------------------------------------------------------------------------
!  Do not include all uprain sources here: only if from updrafts
!------------------------------------------------------------------------
do i = 1,nd
  tot_rv_entrain = tot_rv_entrain + rvflow(1,i)
  tot_rv_detrain = tot_rv_detrain + rvflow(2,i)
  tot_rc_entrain = tot_rc_entrain + rlflow(1,i)
  tot_rc_detrain = tot_rc_detrain + rlflow(2,i)
  tot_uprain = tot_uprain + uprain_rlp(i)  ! only the part of uprain prod by updrafts
  tot_ansnow = tot_ansnow + ansnow_rlp(i)  ! only the part of ansnow prod by updrafts
end do
!------------------------------------------------------------------------
!
! Calculate updraft_eff
!
!------------------------------------------------------------------------
tot_entrain = tot_rv_entrain + tot_rc_entrain
if (tot_entrain > 0.000001) then
  updraft_eff = (tot_uprain + tot_ansnow)/tot_entrain
  updraft_det = (tot_rc_detrain + tot_rv_detrain)/tot_entrain
  updraft_rlp = tot_rv_precip/tot_entrain
! print*,'updraft_eff = ',updraft_eff
! print*,'updraft_det = ',updraft_det
! print*,'updraft_rlp = ',updraft_rlp
else
  updraft_eff = bad
  updraft_det = bad
  updraft_rlp = bad
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
  if ((updraft_eff < -450).and.(updraft_eff > -500)) then
    print*,'updraft_eff in if_conv_tend.f90 = ',updraft_eff
  endif
endif
!------------------------------------------------------------------------
!
!   END   (calc_eff)
!
!------------------------------------------------------------------------
end subroutine calc_eff
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
subroutine calc_org_new(om5)
!------------------------------------------------------------------------
!
!   Use Statements
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
use if_conv_solvers, only : g
implicit none
!------------------------------------------------------------------------
!
!   In/Out variables
!
!------------------------------------------------------------------------
real(r8), intent(in) :: om5
!------------------------------------------------------------------------
!
!   Local Variables
!
!------------------------------------------------------------------------
real(r8) :: xx
real(r8) :: small, diff, yyy
integer :: i
!------------------------------------------------------------------------
!
!  The units of precipprev should be mm/day
!  - units of uprain are mm/sec (kg/m2/s)
!  - try using ansnow_conv only ...
!  - In not using ansnow_strat here, I think my thinking was that I didn't
!    not want the org to be initiated by high upper trop humidity.
!
!------------------------------------------------------------------------
xx = tstep_cam4/tscale_org
yyy = 3600.*24.*(uprain_start_rlp + ansnow_conv)
if (yyy < 0.) yyy = 0.
precip_org_new = (1.-xx)*precip_org + xx*yyy

cm5_new = 0.5*(fmass(5) + fmass(6))*g*3600.*24.*0.01

  if (precip_org > 10.) then
!  print*,'-----------'
!  print*,'om5 = ',om5
!  print*,'cm5_prev = ',cm5_prev
!  print*,'cm5_new = ',cm5_new
!  print*,'amp_om5 = ',amp_om5
!  print*,'precip_org = ',precip_org
!  print*,'amp_rain = ',amp_rain
!  print*,'amp = ',amp
  endif
!------------------------------------------------------------------------
!
!   END   (calc_org_new)
!
!------------------------------------------------------------------------
end subroutine calc_org_new
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
! end module if_conv_tend
!
!------------------------------------------------------------------------
end module if_conv_tend

