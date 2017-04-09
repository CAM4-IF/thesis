!------------------------------------------------------------------------
!
!   Changes:
!   --------
!
!
!   module if_conv_solvers
!
!   Contains parameters and subroutines needed to solve for thermodynamic 
!      quantities, do mixing, etc.
!   These should be local routines that do not require any knowledge of the
!     state variables other than inputs at one level, and no adjustment of
!     tendencies/mass fluxes/detrainment/entrainment/hmflow/vpflow/etc.
!     Basically should be no arrays in these subroutines/functions
!
!------------------------------------------------------------------------
module if_conv_solvers
!------------------------------------------------------------------------
!
!   Use statements
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
!use if_conv_tend, only: esat_ice
implicit none
!------------------------------------------------------------------------
!
!   Allowable errors
!   - these are also used in conv_tend
!   - In double precision, errors can be set at least 100 times smaller
!   - these should all be relative errors
!
!  rherr_max: 
!     - only used in subroutine get_rv. 
!     - Max allowed deviation of rv from
!       rvs when look for saturated solution. Error here does not give rise to
!       errors in total water, so not a conservation issue.
!     - Sep 2010: increased rherr_max to 0.1. Since was generating an error. This
!       is not a real "error", so shouldn't really matter. Just strange that parcels
!       could be so supersaturated.
!
!  hmerr_max:
!   - reduced hmerr_max to 1.0E-07 due to problems in subroutine start. It is likely
!     that these problems are due to residual discontinuity in esat at 0 C, if this
!     was fixed could probably lower the allowed error again.
!   - June 2014: tried hmerr_max = 1.0E-06 to see if increased speed
!
!
!------------------------------------------------------------------------
real(r8), parameter :: hmerr_max = 1.0E-07 ! max hm/kg error for double precision
real(r8), parameter :: rherr_max = 1.0E-01  ! max rh error; only used in subroutine get_rv
real(r8), parameter :: dperr_max = 1.0E-09  ! max dp error in tendency calc
real(r8), parameter :: dt_small = 1.0E-04   ! max dt in subroutine get_rv
!------------------------------------------------------------------------
!
!  Public Parameters (Thermodynamic Constants)
! - for cam constants see: /home/folkins/cesm1_0/models/csm_share/shr/shr_const_mod.F90
! - does it matter for conservation properties how p0 is defined? Doubt it ...
!
!------------------------------------------------------------------------
real(r8), parameter :: rvv = 461.50464    ! should be same as cam
real(r8), parameter :: rdd = 287.04231    ! should be same as cam
real(r8), parameter :: eps = rdd/rvv
real(r8), parameter :: epsi = 1./eps
real(r8), parameter :: cpd = 1004.64      ! same as cam Sept 2010
real(r8), parameter :: cvd = 717.
real(r8), parameter :: cl = 4188.0        ! same as cam4
real(r8), parameter :: cpv = 1810.0       ! same as cam Sept 2010
real(r8), parameter :: ci = 2117.27       ! same as cam4
real(r8), parameter :: cpvmcl = cpv - cl
real(r8), parameter :: clmci = cl - ci
real(r8), parameter :: lv0 = 2.501E+06    ! SHR_CONST_LATVAP in cam4
real(r8), parameter :: lf0 = 3.337E+05    ! SHR_CONST_LATICE in cam4
real(r8), parameter :: tkelvin = 273.15
real(r8), parameter :: g = 9.80616         ! same as cam Sept 2010
real(r8), parameter :: ginv = 1.0/g
real(r8), parameter :: p0 = 100000.   ! cam4 has: SHR_CONST_PSTD    = 101325.0_R8  standard pressure 
real(r8), parameter :: rhow = 1000.   ! density of water
real(r8), parameter :: bad = -999.
real(r8), parameter :: km_badd = 450000.
real(r8), parameter :: t_badd = 400.
real(r8), parameter :: t_min = 3.0
!------------------------------------------------------------------------
!
!  Private Parameters
!    - only used within the module. Any?
!
!------------------------------------------------------------------------
!integer, parameter :: esat_ice = 0
integer, parameter :: esat_ice = 1
!------------------------------------------------------------------------
!
!   Public interfaces
!
!------------------------------------------------------------------------
  public t_from_km ! (subroutine) t from fixed input km,p,rc,rl,ri
  public get_rv    ! (subroutine) t,rv,rl,ri from fixed input km,rt,p
  public start     ! (subroutine) given fixed hm and rh, determines t,rv,rl,ri,rn)
  public esat      ! (function)
  public rsat      ! (function)
  public drssdt    ! (function)
  public get_t     ! (function) Calculate a temperature tm from 
                   !            enthalpy kmm assumming fixed rvm/rlm/rim 
  public enthalpy  ! (subroutine)
  public hmoist    ! (function)
  public latent    ! (function)
  public fusion    ! (function)
  public dksdt     ! (function)
  public dhmsdt    ! (function)
!------------------------------------------------------------------------
!
!   Private interfaces
!    - only used within the module. Any?
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!   contains
!
!------------------------------------------------------------------------
contains
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   subroutine t_from_km
!
!   Calculates grid cell temperature tm from km/rvm/rlm/rim
!   - water partitioning is assumed fixed
!   - tm solution may give rise to supersaturation. It is better to allow
!     this rather than change basic state variables. There is no check for
!     supersaturation.
!   - could change to input total water, and assume balance is fixed liquid 
!     condensate
!   - This program now actually does very little except checing for errors.
!
!------------------------------------------------------------------------
subroutine t_from_km(kmm,pb,rvm,rlm,rim,tm,get_t_call)
!------------------------------------------------------------------------
!
!  use statements
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
!------------------------------------------------------------------------
!
!  Input variables
!
!------------------------------------------------------------------------
implicit none
real(r8), intent(in) :: kmm
real(r8), intent(in) :: pb
real(r8), intent(in) :: rvm
real(r8), intent(in) :: rlm  ! 
real(r8), intent(in) :: rim  ! 
integer, intent(in) :: get_t_call
!------------------------------------------------------------------------
!
!  Output variable
!
!------------------------------------------------------------------------
real(r8), intent(out) :: tm
!------------------------------------------------------------------------
!
!  Local variables
!
!------------------------------------------------------------------------
real(r8) :: rtm
!------------------------------------------------------------------------
!
!  Check for negative kmm.
!
!------------------------------------------------------------------------
  if ((kmm < 0._r8).or.(kmm > km_badd)) then
    print*,'==== FAILURE: neg or large kmm on entry to subroutine t_from_km ====='
    print*,'input pb = ',pb
    print*,'input kmm = ',kmm
    print*,'km_badd = ',km_badd
    print*,'input rvm = ',rvm
    print*,'input rlm = ',rlm
    print*,'input rim = ',rim
    print*,'get_t_call = ',get_t_call
    print*,'in if_conv_solvers.f90'
    stop 'negative or large kmm on entry to subroutine t_from_km (larger than km_badd)'
  endif
!------------------------------------------------------------------------
!
!  Check for negative rtm 
!
!------------------------------------------------------------------------
  rtm = rvm + rlm + rim
  if (rtm < 0._r8) then
    print*,'====== FAILURE: negative rtm on entry to subroutine t_from_km ====='
    print*,'rtm = ',rtm
    print*,'input pb = ',pb
    print*,'input kmm = ',kmm
    print*,'input rvm = ',rvm
    print*,'input rlm = ',rlm
    print*,'input rim = ',rim
    print*,'KEEP GOING: often comes from rvnew calculation at large rain'
    print*,'get_t_call = ',get_t_call
!   stop 'negative rtm on entry to subroutine t_from_km'
  endif
!------------------------------------------------------------------------
!
!  Use fixed km,rtm,rvm,rim to get tm
!
!------------------------------------------------------------------------
  tm = get_t(kmm,rtm,rvm,rim)
!------------------------------------------------------------------------
!
!   Check if tm in reasonable physical range
!
!------------------------------------------------------------------------
  if ((tm < t_min).or.(tm > t_badd)) then
    print*,'-----  tm problem after get_t call -------------------'
    print*,'pb = ',pb
    print*,'tm = ',tm
    print*,'t_badd = ',t_badd
    print*,'t_min = ',t_min
    print*,'kmm = ',kmm
    print*,'rvm = ',rvm
    print*,'rlm = ',rlm
    print*,'rim = ',rim
    print*,'rtm = ',rtm
    print*,'get_t_call = ',get_t_call
    print*,'Keep Going'
!   stop '****** tm out of physical range in subroutine t_from_km ****'
  endif
!------------------------------------------------------------------------
!
!   End
!
!------------------------------------------------------------------------
end subroutine t_from_km
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  subroutine get_rv
!
!  - Calculates the temperature of a parcel, based on the assumption
!    that rv is either the saturated CC value, or that rl = 0.
!  - Used for situations in which there is assumed to be an exchange
!    between the vapor and condensate.
!
!  Inputs:
!  (i) kmm: total enthalpy of parcel
!  (ii) rtm: total water of parcel
!  (iii) pb: background pressure
!
!  Outputs:
!  (i) rvm
!  (ii) rlm
!  (iv) tm (but input used as starting guess for sat solution)
!
!  Three Constraints:
!  (i) total water conservation
!  (ii) moist enthalpy conservation
!  (iii) C-C relation
!
!  The saturation vapor pressure does not depend on the presence of
!  ice, i.e. es = esi for T < 0 C, and es = esw for T > 0 C. This
!  is inconsistent with how ice is treated, i.e. allowing for the existence
!  of supercooled water. However, allowing es to depend on the presence
!  of ice would introduce a discontinuity in es which would likely make
!  the numerical solution much more difficult.
!
!  MMM
!
!------------------------------------------------------------------------
subroutine get_rv(kmm_in,rtm_in,pb,rvm,rlm,rim_in,tm,get_rv_call,allow_sat)
!------------------------------------------------------------------------
!
!  use statements
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
!------------------------------------------------------------------------
!
!  Input variables
!
!------------------------------------------------------------------------
implicit none
real(r8), intent(in) :: kmm_in    ! input moist enthalpy
real(r8), intent(in) :: rtm_in    ! input total water mixing ratio
real(r8), intent(in) :: rim_in    ! input ice mixing ratio (fixed)
real(r8), intent(in) :: pb        ! input pressure in Pa
integer, intent(in) :: get_rv_call   ! for debugging, figure out what made problem
integer, intent(in) :: allow_sat  ! if allow sat solution
!------------------------------------------------------------------------
!
!  Output variables
!  - In order to subtract the ice enthalpy, tm must be as close as possible 
!    to a realistic temperature on input.
!
!------------------------------------------------------------------------
real(r8), intent(inout) :: tm    ! output temperature (but input used as guess)
real(r8), intent(out) :: rvm   ! output vapor mixing ratio
real(r8), intent(out) :: rlm   ! output liquid condensate (may be zero)
!------------------------------------------------------------------------
!
!  Local variables
!
!------------------------------------------------------------------------
real(r8) :: tg,kg,terror,kmerror,aerror,dkgdt,dk,khave,rsg,rvg,rcg,rig
real(r8) :: rherror,rhmax,dt,rsm,rterror,tm_old,tm_unsat,tm_in
real(r8) :: rtt,rterr,small,km_of_ice
real(r8) :: kmm,rtm,rim
integer :: niter_max,iter,nprint
!------------------------------------------------------------------------
!
!  Check for negative kmm_in:
!
!------------------------------------------------------------------------
! print*,'entering get_rv kmm_in = ',kmm_in
  if (kmm_in < 0._r8) then
    print*,'=== FAILURE negative kmm_in on entry to subroutine get_rv in if_conv_solvers.f90 ====='
    print*,'input kmm_in = ',kmm_in
    print*,'input rtm_in = ',rtm_in
    print*,'input pb = ',pb
    print*,'get_rv_call = ',get_rv_call
    stop 'negative kmm_in on entry to subroutine get_rv'
  endif
!------------------------------------------------------------------------
!
!  Check for negative rtm_in on entry:
!
!------------------------------------------------------------------------
  if (rtm_in < 0._r8) then
    print*,'====== FAILURE: negative rtm_in on entry to subroutine get_rv ======='
    print*,'input kmm_in = ',kmm_in
    print*,'input rtm_in = ',rtm_in
    print*,'input pb = ',pb
    print*,'get_rv_call = ',get_rv_call
    stop 'negative rtm_in on entry to subroutine get_rv'
  endif
!------------------------------------------------------------------------
!
!  Define the local "noice" values: kmm and rtm
!  This definition of the km of ice must be consistent with the enthalpy
!    subroutine.
!  Note that this relies on some realistic input value of tm.
!
!------------------------------------------------------------------------
  tm_in = tm
  rtm = rtm_in - rim_in
! print*,'tm being used for ice enthalpy estimate = ',tm
  km_of_ice = cl*tm - fusion(tm)
  kmm = kmm_in - rim_in*km_of_ice
  if (rtm < 0.) stop 'rtm negative in get_rv'
  if (kmm < 0.) stop 'kmm negative in get_rv'
!------------------------------------------------------------------------
!
!  From now on the subroutine works with the non-ice part of the parcel,
!    using kmm and rtm.
!
!  *************   Guess the unsaturated solution  *********************
!
!  - rtm defined on entry
!  - assume rvm = rtm: is solution consistent?
!
!------------------------------------------------------------------------
  rvm = rtm
  rim = 0.   ! since non-ice part of parcel
  tm = get_t(kmm,rtm,rvm,rim)
  tm_unsat = tm  !  for diagnostics
!------------------------------------------------------------------------
!
!  Have unsaturated guess for tm:
!  - from tm -> rsm
!  - check rvm < rs(tm). 
!  - If yes, have self-consistent solution
!  - calculate solution enthalpy khave
!  - Check enthalpy error and exit (enthalpy error should be zero)
!
!------------------------------------------------------------------------
  rsm = rsat(tm,pb)
! print*,'in get_rv: rtm = ',rtm
! print*,'in get_rv: rsm = ',rsm
  if ((rvm <= rsm).or.(allow_sat == 1)) then
    rvm = rtm
    rlm = 0._r8
    rim = 0.
    call enthalpy(khave,tm,rtm,rvm,rim)
    kmerror = abs(khave-kmm)/kmm
    if (kmerror > hmerr_max) then
      print*,'----------- large kmerror in get_rv ----------------'
      print*,'relative error in enthalpy = ',kmerror
      print*,'This is for the unsaturated solution'
      print*,'get_rv_call = ',get_rv_call
      stop 'relative error in enthalpy in get_rv too big'
    endif
!------------------------------------------------------------------------
!
!   Check if tm in reasonable physical range
!
!------------------------------------------------------------------------
  if ((tm < 3.).or.(tm > 380.)) then
    print*,'--------- unsaturated solution in subroutine get_rv ----------'
    print*,'pb = ',pb
    print*,'tm = ',tm
    print*,'kmm = ',kmm
    print*,'khave = ',khave
    print*,'rvm = ',rvm
    print*,'rtm = ',rtm
    print*,'rim = ',rim
    print*,'rlm = ',rlm
    print*,'rsm = ',rsm
    print*,'allow_sat = ',allow_sat
    print*,'get_rv_call = ',get_rv_call
    print*,'********** keep going ***********'
!   stop '*******  tm out of physical range in subroutine get_rv ********'
  endif
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
  else
!------------------------------------------------------------------------
!
!   ****************   Saturated Solution *************************
!
!   Keep rtm and kmm fixed
!   Let tg = tm (unsaturated guess)
!   find rsg
!   find kg
!   loop over iter:
!     define dkgdt (assumes parcel saturated)
!     dk = kmm - kg  
!     dt = dk*(dT/dk)
!     new tg = tg + dt
!     find new rsg
!     find new kg
!     get estimate of t error, relative k error
!    end loop
!
!  There is no point defining rcg during this loop (could however
!  check that rsg < rtm)
!
!  Sept 2014:
!  - added dt way to exist loop
!  - converges very quickly, even though dt is often very big in the 
!    first loop (> 30 K). Suggest that using the unsaturated guess is
!    not a very good first guess. With the new "do while" condition,
!    almost always only execute two loops, and rherror is usually
!    still very small.
!
!------------------------------------------------------------------------
!   niter_max = 6
    niter_max = 8
    iter = 0
    tg = tm_in
    rsg = rsat(tg,pb)
    rig = 0.
    call enthalpy(kg,tg,rtm,rsg,rig)   ! get kg from tg/rtm/rsg
    dt = 10.
!   do iter = 1,niter_max
    do while ((abs(dt) > dt_small).and.(iter <= niter_max))
      iter = iter + 1
      dkgdt = dksdt(tg,pb,rtm,rsg,rig)
      dk = kmm - kg
      dt = dk/dkgdt
      tg = tg + dt
      rsg = rsat(tg,pb)
      call enthalpy(kg,tg,rtm,rsg,rig)
!     print*,'------------ iter = ',iter
!     print*,'pb = ',pb
!     print*,'kg = ',kg
!     print*,'tg = ',tg
!     print*,'rtm = ',rtm
!     print*,'rsg = ',rsg
!     print*,'rig = ',rig
!     print*,'dkgdt = ',dkgdt
!     print*,' in get_rv iter dt = ',iter,dt
    end do
!------------------------------------------------------------------------
!
!  ****************  END OF ITERATIONS ********************
!
!  Define final rlm from rvm saturated solution and known rtm
!  - note that rtm here is not rtm_in: is rtm_in minus the ice part
!
!------------------------------------------------------------------------
  rvm = rsg
  rim = 0.
  rlm = rtm - rvm - rim
  if (rlm < 0.) then
    rlm = 0.
    rvm = rtm
  endif
!------------------------------------------------------------------------
!
!  Check for negative rlm
!
!------------------------------------------------------------------------
  small = 1.0E-08
  if (rlm < -small) then
    print*,'************** rlm negative problem in subroutine get_rv ****'
    print*,'this is after iterating for the saturated solution'
    print*,'Why did unsaturated solution not work?'
    print*,'rlm just defined as rtm - rvm with rvm from iteration'
    print*,'In this case can just set rvm = rtm and rlm = 0'
    print*,'This would effectively be an unsaturated solution'
    print*,'Could be a problem with starting T guess tm_in = ',tm_in
    print*,'rlm = ',rlm
    print*,'rvm = ',rvm
    print*,'rtm = ',rtm
    print*,'rim = ',rim
    print*,'get_rv_call = ',get_rv_call
    stop 'rlm negative in subroutine get_rv'
  endif
!------------------------------------------------------------------------
!
!  Check if relative error in rtp is OK (should be zero)
!
!------------------------------------------------------------------------
  rterror = abs(rvm + rlm - rtm)/rtm
  if (rterror > 1.0E-09) then
    print*,'------------- rterror in get_rv should be zero ----------------'
    print*,'This is for the saturated solution'
    print*,'rterror = ',rterror
    print*,'rherror = ',rherror
    print*,'input rtm = ',rtm
    print*,'solution rvm = ',rvm
    print*,'solution rlm = ',rlm
    print*,'rsm = ',rsm
    print*,'tm = ',tm
    print*,'khave = ',khave
    print*,'pb = ',pb
    print*,'allow_sat = ',allow_sat
    print*,'get_rv_call = ',get_rv_call
    print*,'-----------------------------------------------------------------'
    stop 'relative rt error in get_rv non-zero'
  endif
!------------------------------------------------------------------------
!
!  Want a solution which preserves total water and moist enthalpy.
!
!  If believe tg, will be some error in moist enthalpy
!  Instead, believe rtm partitioning: rvg,rlg,rig,rsg. 
!  This guarantees water conservation.
!  Find tm by enforcing km conservation.
!  Error will be that rs(tm) will not exactly equal rvg.
!  Just have to live with the possibility of small super/under saturations
!   for air parcels that have condensate.
!
!  Use saturated rvm and known (input) kmm to find tm
!
!------------------------------------------------------------------------
  tm = get_t(kmm,rtm,rvm,rim)
! print*,'tm = ',tm
! print*,'rvm = ',rvm
! print*,'rtm = ',rtm
! print*,'calling enthalpy'
  call enthalpy(khave,tm,rtm,rvm,rim)   
  kmerror = abs(khave-kmm)/khave
!------------------------------------------------------------------------
!
!   Check if relative error in moist enthalpy OK (should be zero)
!   Designed to exactly conserve enthalpy, so error should be zero.
!   - probably redundant
!
!------------------------------------------------------------------------
  if (kmerror > hmerr_max) then
    print*,'------- enthalpy error in subroutine get_rv exceeds tolerance -----'
    print*,'This routine should exactly conserve total water and enthalpy by construction'
    print*,'If conservation is violated might be a segmentation fault'
    print*,'The only error should be that the desired input rh does not equal actual rh'
    print*,'This is for the saturated solution'
    print*,'---------- km ISSUES ------------'
    print*,'relative error in enthalpy = ',kmerror
    print*,'max relative error hmerr_max = ',hmerr_max
    print*,'kmm = ',kmm
    print*,'current khave = ',khave
    print*,'---------- rtm ISSUES ------------'
    print*,'rterror = ',rterror
    print*,'current sum : rvm + rlm = ',rvm+rlm
    print*,'actual current rtm = ',rtm
    print*,'---------- rvm ISSUES ------------'
    print*,'rvm = ',rvm
    print*,'rvg from iteration loop should be same as rvm rvg = ',rvg
    print*,'rsm = ',rsm
    print*,'effective RH (need not be 1) rvm over rsm = ',rvm/rsm
    print*,'---------- enthalpy calculation ------------'
    print*,'tm = ',tm
    print*,'cpvmcl = ',cpvmcl
    print*,'cl = ',cl
    print*,'cpvmcl times tkelvin - lv0 = ',cpvmcl*tkelvin - lv0
    print*,'get_rv_call = ',get_rv_call
    print*,'---------------------------------------------------------------'
    stop 'relative enthalpy error in subroutine get_rv too big'
  endif
!------------------------------------------------------------------------
!
!   Check if relative error in rh is small (in general will be non-zero)
!   - Forced to allow some small supersaturation in forcing rt/km conservation.
!   - Program should just convert any supersaturation to rl, so should not cause
!     a problem.
!   - Jan 25/2011: getting rh errors at very cold Temperatures (200 K).
!
!------------------------------------------------------------------------
  rsm = rsat(tm,pb)
  rherror = abs(rsm - rvm)/rsm
! print*,'in get_rv rherror = ',rherror
  if (rherror > rherr_max) then
     print*,'------- WARNING: RH error exceeds threshold in subroutine get_rv --------'
     print*,'This is for the saturated solution assumption'
     print*,'The final solution should be close to saturation'
     print*,'May not be that serious: can still have small km and rt errors'
     print*,'May mean that original assumption of saturated solution not valid'
     print*,'Try to figure out why unsaturated solution was not chosen'
     print*,'Or the final guess tg from the iteration loop is too different'
     print*,'from the tm from using the final saturated rv from the loop and'
     print*,'enforcing enthalpy and water conservation.'
     print*,'-------------------'
     print*,'pressure in hPa = ',pb*0.01
     print*,'rherror defined as abs(rsm - rvm)/rsm'
     print*,'relative error in RH = ',rherror
     print*,'solved RH of parcel = ',rvm/rsm
     print*,'------- enthalpy ----------'
     print*,'input moist enthalpy kmm_in = ',kmm_in
     print*,'Enthalpy with ice part subtracted kmm = ',kmm
     print*,'ice enthalpy subtracted km_of_ice = ',km_of_ice
     print*,'ice enthalpy subtracted rim_in times km_of_ice = ',rim_in*km_of_ice
     print*,'final moist enthalpy khave = ',khave
     print*,'relative error in km kmerror = ',kmerror
     print*,'------- temps ----------'
     print*,'Input temperature first guess for sat solution tm_in = ',tm_in
     print*,'Initial unsat guess temp tm_unsat = ',tm_unsat
     print*,'tm_unsat determined from get_t using kmm and rtm (assumed rv)'
     print*,'current solution for temp before ice added back in tm = ',tm
     print*,'Final best tg in the iteration (determines rvm) = ',tg
     print*,'A big difference between these two temps causes large rherror.'
     print*,'A large T increment may imply need to increase niter_max.'
     print*,'Final temperature increment in loop dt = ',dt
     print*,'Number of iterations done = ',iter
     print*,'------- ice ----------'
     print*,'fixed input ice rim_in = ',rim_in
     print*,'------- rv ----------'
     print*,'There is no input rv'
     print*,'Unsat guess rv (assumed to equal rtm_in - rim_in) = ',rtm
     print*,'rs from tm_unsat = ',rsat(tm_unsat,pb)
     print*,'Is the unsat guess rv larger than rs from tm_unsat?'
     print*,'current solution for rv rvm = ',rvm
     print*,'current solution for rsm (from tm) = ',rsm
     print*,'solved RH of parcel = ',rvm/rsm
     print*,'------- rl ----------'
     print*,'current solution rlm = ',rlm
     print*,'rherr_max = ',rherr_max
     print*,'final relative total water = ',rterror
     print*,'------- other stuff ----------'
     print*,'get_rv_call = ',get_rv_call
     stop 'relative error in RH in get_rv too big'
  endif
!------------------------------------------------------------------------
!
!   Check if tm in reasonable physical range
!
!   Aug 11, 2011: got an output temperature here of 39 K, so reduced min temp
!    to 20 K. Presumably very small temp was due to a very small amount of 
!    convection reaching a very high altitude (here 53 hPa). Though strange
!    that km of parcel was 40 kJ; maybe hm more reasonable, but doubt it; would
!    correspond to a very low PT; could have originated from convection at high
!    altitudes. Oct 2011: had similar problem.
!
!------------------------------------------------------------------------
  if ((tm < t_min).or.(tm > t_badd)) then
    print*,'----- problem for saturated solution in subroutine get_rv ----------'
    print*,'---- tm less than t_min or larger than t_badd ---'
    print*,'input kmm = ',kmm
    print*,'input rtm = ',rtm
    print*,'input pb = ',pb
    print*,'output tm = ',tm
    print*,'t_min = ',t_min
    print*,'output rvm = ',rvm
    print*,'output rlm = ',rlm
    print*,'khave = ',khave
    print*,'rsm = ',rsm
    print*,'kmerror = ',kmerror
    print*,'get_rv_call = ',get_rv_call
    print*,'Keep Going'
!   stop '******* tm out of physical range in subroutine get_rv ********'
  endif
!------------------------------------------------------------------------
!
!  endif for sat/unsat solution
!
!------------------------------------------------------------------------
  endif
!------------------------------------------------------------------------
!
!  **********   Find the final temperature   ******************
!
!  - Add the ice back in
!  - use the fixed input kmm_in, rim_in, rtm_in
!  - kmm_in and rtm_in should include the ice parts.
!  - use the rvm as determined above
!  - Note that since have a new temperature rvm will no longer in general
!    be as close to the saturated value.
!
!------------------------------------------------------------------------
  tm_old = tm
  tm = get_t(kmm_in,rtm_in,rvm,rim_in)
  nprint = 0
  if (nprint == 1) then
  if (rim_in > 1.0E-09) then
    print*,'---------  in get_rv  ------------'
    print*,'In subroutine get_rv: tm - tm_old = ',tm-tm_old
    print*,'rim_in = ',rim_in
    print*,'kmm = ',kmm
    print*,'rtm = ',rtm
    print*,'rvm = ',rvm
    print*,'tm = ',tm
    print*,'tm_old = ',tm_old
    print*,'--- leaving get_rv  ----'
  endif
  endif
!------------------------------------------------------------------------
!
!   End
!
!------------------------------------------------------------------------
end subroutine get_rv
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Subroutine start:
!
!   - Given fixed km/rh/pp, determines t,rv
!   - assumes that rh <= 1, so therefore no condensate
!   - Used to help determine spectrum of starting rv/t
!   - The routine iteratively solves for a best guess t. With each
!       iteration, a new total water is re-calculated from rt = rh*rs(t).
!       rh is fixed at the input value.
!
!   This does not use get_rv. It probably could be adapted to use get_rv, but
!     it's a bit different since it is going (km,rh) -> t,rv,rl rather
!     conserving total water, as most other routines do.
!
!   Inputs: tg is first guess (default is ambient temperature)
!
!   Two strategies:
!   - Loop iterates a new t, progressively trying to minimize km error. At then
!     end of the loop, you have two choices:
!     (1) Accept change in km:
!         Believe the final t, use the fixed rh to get the final rv, and solve
!         for the actual km. This will be different from the input km. Accept
!         this difference.
!     (2) Accept change in rh:
!         Again believe the final t, and use the input km to find rv. 
!         Accept change in final rh from input. May mean have final rh > 1.
!
!    If kmerror in iterative loop is small, should not make any difference which
!    option is picked.
!
!    However, here have chosen option (1) so as to avoid possibly supersaturated
!      air parcels. Does not lean to column conservation error, since actual km
!      of starting air parcel need not exactly equal kmstart spectrum.
!
!   Sept 2010: having problems with temps near 0 C since drssdt is discontinuous.
!   Sept 2014: max number of iterations here seems to be 4. So this is not
!      likely to be very slow.
! 
!   Jan 2015: this program was modified to include a fixed input rl called
!             rll, which is the background rl at that height. Intent here
!             is to be able to entrain rl from BL, and prevent buildup
!             of high rl in BL at high colrh.
!
!------------------------------------------------------------------------
subroutine start(kmwant,km_out,rhwant,rll,rh_out,tg_in,t_out,rv_out,pp)
!------------------------------------------------------------------------
!
!  In/out variables
!
!------------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: kmwant,rhwant,rll,pp,tg_in
real(r8), intent(out) :: km_out,rh_out,t_out,rv_out
!------------------------------------------------------------------------
!
!  Local variables
!
!------------------------------------------------------------------------
real(r8) :: rig,rsg,drsdt,rvg,lvg,dkmdt,kmerror,kmg,tg,dtt,rtt,rtg
integer :: niter_max,iter
real(r8) :: small,rs_k,rv_k,rt_k,ri_k,km_k
!------------------------------------------------------------------------
!
!  Check for out of range input rhwant:
!
!------------------------------------------------------------------------
if ((rhwant < 0.).or.(rhwant > 1.)) then
! print*,'input rhwant to subroutine start = ',rhwant
! stop 'input rhwant is out of range (negative or larger than 1)'
endif
tg = tg_in 
!------------------------------------------------------------------------
!
!   When using the version of es(t) where es = eice below 0 C, drsst becomes
!   discontinuous at 0 C. This can give rise to non-convergence for solutions
!   close to 0 C. km_k is used below to determine when the solution for tg is on
!   the "wrong side" of tkelvin and fix it.
!
!------------------------------------------------------------------------
  small = 1.0E-09
  rs_k = rsat(tkelvin,pp)
  rv_k = rhwant*rs_k
  rt_k = rv_k + rll
  ri_k = 0.
  call enthalpy(km_k,tkelvin,rt_k,rv_k,ri_k)
!------------------------------------------------------------------------
!
!   Calculate dkm/dt at fixed rh (rh < 1) (changed March 2011):
!   ===========================================================
!
!   km = (cpd + rv*cl)*t + lv*rv    (Unsaturated case)
!   lv = lv0 + cpvmcl*(t - tkelvin)
!   dlv/dt = cpvmcl
!   Assume that rv = rh*rs 
!
!   km = [cpd + cl*rh*rs]*t + lv*rh*rs 
!   dkmdt = cpd + cl*rh*d/dt(rs*t) + rh*d/dt(lv*rs) 
!         = cpd + cl*rh*(drsdt*t + rs) + rh*(cpvmcl*rs + lv*drsdt
!         = cpd + rh*drsdt*(cl*t + lv) + cl*rh*rs + rh*cpvmcl*rs
!         = cpd + rh*drsdt*(cl*t + lv) + rh*rs*(cl + cpvmcl)
!
!  Setting niter_max: in most cases kmerror reaches a particular value after
!    4 iterations, and does not improve after that. 
!
!------------------------------------------------------------------------
  rig = 0.
  niter_max = 100
  iter = 0
  kmerror = 1.
  do while ((kmerror > hmerr_max).and.(iter < 4))
! do while ((kmerror > hmerr_max).or.(iter < 4))
    iter = iter + 1
!------------------------------------------------------------------------
!
!  Check for non-convergence of solution
!  - use enthalpy error as test
!
!------------------------------------------------------------------------
  if (iter > (niter_max - 10)) then
    print*,'----- kmerror exceeds hmerr_max in subroutine start ---------------'
    print*,'iter = ',iter
    print*,'current relative kmerror = ',kmerror
    print*,'max allowed relative km error hmerr_max = ',hmerr_max
    print*,'current kmg = ',kmg
    print*,'desired input kmwant = ',kmwant
    print*,'absolute km error kmwant - kmg = ',kmwant-kmg
    print*,'input rhwant = ',rhwant
    print*,'input tg_in = ',tg_in
    print*,'current tg = ',tg
    print*,'current rvg = ',rvg
    print*,'input pp = ',pp
    print*,'dkmdt = ',dkmdt
    print*,'drsdt = ',drsdt
    print*,'lvg = ',lvg
    print*,'cpd = ',cpd
    print*,'rhwant*drsdt*(cl*tg + lvg) = ',rhwant*drsdt*(cl*tg + lvg)
    print*,'rhwant*rsg*(cl + cpvmcl) = ',rhwant*rsg*(cl + cpvmcl)
    print*,'tg - tkelvin = ',tg-tkelvin
    print*,'dtt = ',dtt
    print*,'-----------------------------------------------------------------'
    if (iter == niter_max) stop 'relative km error not converging in start'
  endif
!------------------------------------------------------------------------
!
!  From new starting tg: tg -> rsg 
!  rsg + rh(fixed at input) -> rvg 
!  tg + rvg + pp -> kmg
!  calculate current km error (to see when can exit iteration)
!  use dkmdt to get dtt and new tg
!  repeat and try to converge
!
!   From above (fixed March 2011):
!
!  dkmdt = cpd + rh*drsdt*(cl*t + lv) + rh*rs*(cl + cpvmcl)
!
!------------------------------------------------------------------------
    rsg = rsat(tg,pp)
    rvg = rhwant*rsg
    rtg = rvg + rll
    call enthalpy(kmg,tg,rtg,rvg,rig)
    kmerror = abs(kmwant-kmg)/kmwant
    drsdt = drssdt(tg,pp)
    lvg = latent(tg)
    dkmdt = cpd + rhwant*drsdt*(cl*tg + lvg) + rhwant*rsg*(cl + cpvmcl)
!   dkmdt = cpd + rhwant*drsdt*(cl*tg + lvg - lv0) + rhwant*rsg*(cl + cpvmcl)
    dtt = ((kmwant - kmg)/dkmdt)
    tg = tg + dtt
!------------------------------------------------------------------------
!
!  If tg is on the wrong side of tkelvin, fix it
!  - or else will never converge
!
!------------------------------------------------------------------------
  if ((kmwant > km_k).and.(tg < tkelvin)) then
    print*,'----------- tg should be larger than tkelvin -----------------------'
    print*,'before change tg dtt = ',tg,dtt
    tg = tkelvin + small
!   print*,'setting tg = tkelvin + small'
!   print*,'km_k = ',km_k
!   print*,'kmwant = ',kmwant
!   print*,'kmwant - km_k = ',kmwant - km_k
  elseif ((kmwant < km_k).and.(tg > tkelvin)) then
!  print*,'----------- tg should be smaller than t kelvin ---------------------------'
!  print*,'before change tg dtt = ',tg,dtt
    tg = tkelvin - small
!  print*,'setting tg = tkelvin + small'
!   print*,'km_k = ',km_k
!   print*,'kmwant = ',kmwant
!   print*,'kmwant - km_k = ',kmwant - km_k
  else
  endif 
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
!    print*,'iter kmerror = ',iter,kmerror
!    print*,'kmg kmwant = ',kmg,kmwant
!------------------------------------------------------------------------
!
!  End loop over iterations
!
!------------------------------------------------------------------------
  end do
! print*,'Number of iterations required in subroutine start = ',iter
!------------------------------------------------------------------------
!
!   From final best tg,rvg, define: t_out,rv_out,rh_out,km_out 
!   - rh_out should be the same as rhwant
!
!------------------------------------------------------------------------
  t_out = tg
  rv_out = rvg 
  rtt = rv_out + rll
  rh_out = rv_out/rsat(t_out,pp)
  call enthalpy(km_out,t_out,rtt,rv_out,rig)
!------------------------------------------------------------------------
!
!   End
!
!------------------------------------------------------------------------
end subroutine start
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Function esat
!
!   Input: temperature t (K)
!   Output: saturation vapor pressure in Pa (not mb), with respect to water
!      for T > 0 C, and with respect to ice for T < 0 C.
!   This definition of esat may be inconsistent with the the possible allowance
!     of the existence of supercooled water in the treatment of ice. However,
!     allowing for the possibility of a discontinuity in esat likely would
!     lead to numerical difficulties.
!
!   esat formula taken from voemel web site, WMO, 2008
!   gives 6.112 hPa at 0 C
!
!  original expression for esat gave:
!   esat = exp(23.33086 - (6111.72784/t) + 0.15215*log(t))
!  But this did not exactly give es = 6.112 hPa at 0 C, so changed it a bit.
!  Discontinuity in es at 0 C was causing problems in subroutine start.
!
!  In fortran log means natural log, so is the log t in the ice esat expression
!    really ln?? may have made a mistake, and caused by rh gap due to a discontinuity
!    in esat at 0 C.
!
!  need: log(e0) = AA - (BB/t) + CC*log(t)
!        AA = log(e0) + (BB/t) - CC*log(t) 
!
!  If you use the water saturation expression for all T, the model behaviour
!    is very different ... esat is used in so many places that it would be hard
!    to diagnose the origin of the differences.
!
!------------------------------------------------------------------------
real(r8) function esat(t)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: t
real(r8) :: tc,denom
  tc = t - tkelvin
  denom = 243.12 + tc
  if (tc > 0.) then
    esat = 6.112*exp(17.62*tc/denom)
  else
    if (esat_ice == 1) then
      esat = exp(23.33167801703645 - (6111.72784/t) + 0.15215*log(t))
    else
      esat = 6.112*exp(17.62*tc/denom)
!     print*,'Using water saturation'
    endif
  endif
  esat = 100.*esat
end function esat
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  function rsat
!
!  Inputs:
!    t: temperature (K)
!    p: pressure (Pa)
!
!  Outputs:
!    rsat
!
!  Interesting: at low pressures in the stratosphere, es > p,
!  so that rsat becomes negative. In this case define rsat as eps*es/p
!
!  Don't think rsat should go negative at ordinary range of t/p
!
!------------------------------------------------------------------------
real(r8) function rsat(t,p)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: t,p
real(r8) :: es
  es = esat(t)
  if (p > 1000.) then
    rsat = eps*es/(p - es)   ! this is the correct definition
                             ! maybe should use actual e but then rs not a fn
                             ! of t and p anymore
  else
    rsat = eps*es/p          ! fix for very low pressures
  endif
  if (rsat < 0.) then
     print*,'======== rsat negative in function rsat ==========================='
     print*,'sometimes a result of extremely high temperatures of around 370 K'
     print*,'t = ',t
     print*,'p = ',p
     print*,'rsat = ',rsat
     print*,'eps =',eps
     print*,'es =',es
     print*,'Re-define as rsat = eps es/p'
     rsat = eps*es/p
!    stop 'rsat negative'
  endif
end function rsat
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  Function drssdt
!
!  Inputs:
!    t: temperature (K)
!    p: pressure (Pa)
!
!  Outputs:
!    drssdt
!
!  If lf is not included below 0 C, desdt is in error by about 10%. With
!    lf, the error is about 0.1%.
!
!------------------------------------------------------------------------
real(r8) function drssdt(t,p)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: t,p
real(r8) :: lv, lf, es, tc, desdt
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
  lv = latent(t)
  es = esat(t)
  tc = t - tkelvin
  if (tc > 0.) then
    desdt = (lv*es)/(rvv*t*t)
  else
    if (esat_ice == 1) then
      lf = fusion(t)
      desdt = ((lv+lf)*es)/(rvv*t*t)
    else
      desdt = (lv*es)/(rvv*t*t)
    endif
  endif
  drssdt = (desdt*eps*p)/((p - es)**2)
end function drssdt
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  function get_t
!
!  Calculate a temperature tm from enthalpy kmm assumming fixed rvm/rlm/rim
!
!  Use this enthalpy definition and solve for T:
!
!  km = (cpd + rt*cl)*T + lv(T)*rv - lf(T)*ri
!
!  Then:
!
!  (cpd + rt*cl)*T = km - lv*rv + lf*ri
!
!  (cpd + rt*cl)*T = 
!  km - [lv0 + cpvmcl*(T - tkelvin)]*rv + [lf0 + clmci*(T - tkelvin)]*ri
!
!  (cpd + rt*cl)*T + cpvmcl*T*rv - clmci*T*ri = 
!  km - [lv0 - cpvmcl*tkelvin]*rv + [lf0 - clmci*tkelvin]*ri
!
!  T*(cpd + rt*cl + cpvmcl*rv - clmci*ri) = 
!  km - rv*(lv0 - cpvmcl*tkelvin) + ri*(lf0 - clmci*tkelvin)
!
!  Gives:
!  T = AA/BB, where:
!  AA = km - rv*(lv0 - cpvmcl*tkelvin) + ri*(lf0 - clmci*tkelvin)
!  BB = (cpd + rt*cl + cpvmcl*rv - clmci*ri)
!
!  Use kl0
!  
!------------------------------------------------------------------------
real(r8) function get_t(km_in,rt_in,rv_in,ri_in)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: km_in,rt_in,rv_in,ri_in
real(r8) :: AA,BB
!------------------------------------------------------------------------
!
!  Non-zero ri_in = 0.  (June 2013)
!  ZZZ
!------------------------------------------------------------------------
AA = km_in - rv_in*(lv0 - cpvmcl*tkelvin) +              &
     ri_in*(lf0 - clmci*tkelvin) 
BB = cpd + rt_in*cl + cpvmcl*rv_in - clmci*ri_in
get_t = AA/BB
!
!print*,'---in get_t --'
!print*,'lv0 = ',lv0
!print*,'-----'
!
end function get_t
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  Subroutine enthalpy (changed March 2011):
!  
!  HHH  
!   
!  Km(md,ml,mi,T) = md*kd(To) + mv*kv(T) + ml*kl(T) + mi*ki(T)
!
!  Can write:
!
!     kd(T) = kd0 + cpd*T 
!     kv(T) = kl(T) + lv(T)
!     kl(T) = kl0 + cl*T 
!     ki(T) = kl(T) - lf(T)
!
!   Have freedom to define "zero" for dry air and water.
!   Set kd0 = 0 and kl0 = 0
!
!  Then:
!
!     kd(T) = cpd*T
!     kv(T) = cl*T + lv(T)
!     kl(T) = cl*T   
!     ki(T) = cl*T - lf(T)`
!
!  km = Km(md,ml,mi,T)/md 
!     = cpd*T + rv*kv(T) + rl*kl(T) + ri*ki(T)
!     = cpd*T + rv*[cl*T + lv(T)] + rl*cl*T + ri*[cl*T - lf(T)]
!     = cpd*T + cl*T(rv + rl+ ri) + rv*lv(T) - ri*lf(T)
!     = cpd*T + cl*T*rt + rv*lv(T) - ri*lf(T)
!     = (cpd + cl*rt)*T + rv*lv(T) - ri*lf(T)
!
!  km_of_ice = cl*T - lf(T)
!  km_of_wat = cl*T 
!   
!   Same as Emanuel's expression
!
!  =================
!
!  Inputs: t (K), rt, rv, ri
!
!  ZZZ
!------------------------------------------------------------------------
subroutine enthalpy(kk,t,rt,rv,ri)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(out) :: kk
real(r8), intent(in) :: t,rt,rv,ri
real(r8) :: lv,lf
  lv = latent(t)
  lf = fusion(t)
  kk = (cpd + rt*cl)*t + lv*rv - lf*ri 
end subroutine enthalpy
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!  function hmoist:
!
!  Inputs: t (K), rt, rv, ri, z (m)
!
!  hm = km + (1. + rt)*g*z
! 
!  [hm] = [J/kg dry air]
!
!  The input z should be the mass weighted z of a layer. Simple mean height
!  could introduce errors.
!
!------------------------------------------------------------------------
real(r8) function hmoist(t,rt,rv,ri,z)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: t,rt,rv,ri,z
real(r8) :: kk
  call enthalpy(kk,t,rt,rv,ri)
  hmoist = kk + (1. + rt)*g*z
end function hmoist
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Latent Heat of Evaporation
!
!   Input: t (K)
!
!------------------------------------------------------------------------
real(r8) function latent(t)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: t
real(r8) :: tc
  tc = t - tkelvin
  latent = lv0 + cpvmcl*tc
end function latent
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Latent Heat of Fusion
!
!   Input: t (K)
!
!   cl = 4188.0       
!   ci = 2117.27    
!   clmci = 2070.73
!   Sept 2010: changed to set equal to lf0 for T > tkelvin
!   - increasing function of T
!   May 2014: got rid of this: was causing numerical problems
!
!------------------------------------------------------------------------
real(r8) function fusion(t)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: t
real(r8) :: tc
  tc = t - tkelvin
  fusion = lf0 + clmci*tc
end function fusion
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Function dksdt
!
!   - Evaluates partial derivative of km with respect to T, assuming the parcel
!     is saturated, and that rt is fixed.
!   - Since rt fixed, no change from March 2011 modification of km
!
!  kk = (cpd + rt*cl)*t + lv*rv - lf*ri 
!
!   See notes for the derivation of:
!  
!   drs/dT = des/dT*(eps*p)/(p - es)**2
!   des/dT = Lv*es/(Rv*T*T)  (valid only for water?)
!   dkm/dT = cpd + rt*cl + cpvmcl*rsm + lv*drsmdt - clmci*rim -lf*drimdt
!
!   dk/dT = cpd + rt*cl + cpvmcl*rs + lv*drsdt - clmci*ri -lf*dridt
!
!   Nov 2009: removed dridt term, so now not valid for non-zero ri
!   (May 2014: I think this is OK, since only used for non-ri part of 
!   an air parcel).
!
!------------------------------------------------------------------------
real(r8) function dksdt(t,p,rt,rs,ri) 
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: t,p,rt,rs,ri
real(r8) :: lv,lf
  lv = latent(t)
! lf = fusion(t)    ! not needed
  dksdt = cpd + rt*cl + cpvmcl*rs + lv*drssdt(t,p) - clmci*ri
end function dksdt
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   Function dhmsdt
!
!  Inputs: t (K), p, rt, rv, ri, z (m)
!
!  k =  (cpd + rt*cl)*T + lv*rv - lf*ri
!  hm = (cpd + rt*cl)*T + lv*rv - lf*ri + (1. + rt)*g*z
!
!  Assumes parcel is saturated, rt fixed.
!  set equal to dkdt, since gravitational term has no temp dependence.
!
!------------------------------------------------------------------------
real(r8) function dhmsdt(t,p,rt,rs,ri) 
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none
real(r8), intent(in) :: t,p,rt,rs,ri
  dhmsdt = dksdt(t,p,rt,rs,ri)
end function dhmsdt
!------------------------------------------------------------------------
!**********************************************************************73
!------------------------------------------------------------------------
!
!   end module
!
!------------------------------------------------------------------------
end module if_conv_solvers
