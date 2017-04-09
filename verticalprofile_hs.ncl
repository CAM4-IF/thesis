; ***********************************************
; xy_3.ncl
;
; Concepts illustrated:
;   - Reversing the Y axis
;   - Changing the line dash pattern in an XY plot
;   - Creating your own line dash pattern for an XY plot
;   - Changing the line color and thickness in an XY plot
;   - Creating a vertical profile plot
;
; ***********************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
;************************************************
begin
;************************************************
; read in data
;************************************************
latlons=asciiread("../latlon",(/9,2/),"float")
lat=latlons(3,0)
lon=latlons(3,1)

o=round(((90+lat)/1.89473),3)
p=round((lon/2.5),3)
r=round(((90+lat)/2.5),3)

data1=asciiread("KOR",(/576,10/),"float")
Tdata=new((/48,12/),float)
Pdata=new((/48,12/),float)
RHdata=new((/48,12/),float)
Zdata=new((/48,12/),float)
Udata=new((/48,12/),float)
Vdata=new((/48,12/),float)
i=0
j=0
k=47
do while(i.le.11)
  Tdata(0:47,i)=data1(j:k,2)
  Pdata(0:47,i)=data1(j:k,0)
  RHdata(0:47,i)=data1(j:k,6)
  Zdata(0:47,i)=data1(j:k,9)
  Udata(0:47,i)=data1(j:k,7)
  Vdata(0:47,i)=data1(j:k,8)
  i=i+1
  j=j+48
  k=k+48
end do
TmeanANN=dim_avg_n(Tdata,1)
PmeanANN=dim_avg_n(Pdata,1)
RHmeanANN=dim_avg_n(RHdata,1)*100
ZmeanANN=dim_avg_n(Zdata,1)*1000
UmeanANN=dim_avg_n(Udata,1)
VmeanANN=dim_avg_n(Vdata,1)




ntimes=dimsizes(PmeanANN)
dTdZ=new(ntimes,double)
do i=0,ntimes-2
meanT=(TmeanANN(i)+TmeanANN(i+1))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(TmeanANN(i+1)-TmeanANN(i))
dTdZ(i)=dT/dZ
end do


f=addfile("../AMIP/ta_ANN_AMIPmean.nc","r")
f1=addfile("../AMIP_4K_control_ANN_climo1.nc","r")
f2  = addfile("CAM-IFv5_AMIP_ANN_climo.nc","r")
cmd="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,ta ../AMIP/ta_ANN_AMIPmean.nc" + " tmp.nc"
cmd1="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,T,PS,Z3 ../AMIP_4K_control_ANN_climo1.nc" + " tmp1.nc"
cmd2="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,T,PS,Z3 CAM-IFv5_AMIP_ANN_climo.nc" + " tmp2.nc"
cmd3="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,T,Z3,PS MERGED_T_trop_ANN.nc" + " tmp3.nc"

areamean=systemfunc(cmd)
areamean1=systemfunc(cmd1)
areamean2=systemfunc(cmd2)
areamean3=systemfunc(cmd3)
f3=addfile("tmp.nc","r")
f4=addfile("tmp1.nc","r")
f5=addfile("tmp2.nc","r")
f6=addfile("tmp3.nc","r")

T3=f3->ta(0,:,:,:)
T4=f4->T(0,:,:,:)
T5=f5->T(0,:,:,:)
T6=f6->T(0,:,:,:)
ps3=f4->PS(0,:,:)
ps4=f5->PS(0,:,:)
;ps5=f6->PS(0,:,:)

z=f4->Z3(0,:,0,0)
z1=f5->Z3(0,:,0,0)
z2=f6->Z3(0,::-1,0,0)



T=f->ta(0,:,:,:)
T1=f1->T(0,:,:,:)
T2=f2->T(0,:,:,:)
hyam=f1->hyam
hybm=f1->hybm
ps1=f1->PS(0,:,:)
ps2=f2->PS(0,:,:)
plev=f->plev

T@_FillValue=1e20
u=int2p_n_Wrap(plev*0.01,T,PmeanANN,0,0)
u1=vinth2p(T1,hyam,hybm,PmeanANN,ps1,2,1000,1,False)
u2=vinth2p(T2,hyam,hybm,PmeanANN,ps2,2,1000,1,False)

hsarray=new((/2,26/),float)
hsarray1=new((/1,dimsizes(z2)/),float)
Tprofs=new((/2,26/),float)
Tprofs(0,:)=T4(:,0,0)
Tprofs(1,:)=T5(:,0,0)
i=0
do i=0,2

;;;;;;;This code calculates moist static energy from T,q and P. Calculates HS on sigma levels first, then interpolate later;;;;;
; NCL script
; functions_profiles.ncl
; Mark Stevens
;*****************************************************************


; input arguments
;  p[*] - pressure on hybrid levels in millibars (lev)
;  tk[*] - temperature in Kelvin (lev)
;
; returns saturation specific humidity (dimensionless)
; uses several functions from Dennis's thermo.f code
;-------------------------------------------------------------------------
; from the fortran function ssh
;   this function returns saturation specific humidity ssh (grams of
;   water vapor per kilogram of moist air) given the pressure p
;   (millibars) and the temperature t (celsius). the equation is given
;   in standard meteorological texts. if t is dew point (celsius), then
;   ssh returns the actual specific humidity.
;   compute the dimensionless mixing ratio.


if (i.eq.2) then 
    T9=T6(::-1,0,0) - 273.15 
	pres= T6&lev(::-1)
	delete(z)
	z=z2
else
    T9 = Tprofs(i,:) - 273.15           ; convert Kelvin to Celsius (lev)
    q = T4(:,0,0)                  ; copy coords and atrributes (lev)
    pres= T4&lev
end if

;-------------------------------------------------------------------------
; from fortran function esw
;   this function returns the saturation vapor pressure esw (millibars)
;   over liquid water given the temperature t (celsius). the polynomial
;   approximation below is due to herman wobus, a mathematician who
;   worked at the navy weather research facility, norfolk, virginia,
;   but who is now retired. the coefficients of the polynomial were
;   chosen to fit the values in table 94 on pp. 351-353 of the smith-
;   sonian meteorological tables by roland list (6th edition). the
;   approximation is valid for -50 < t < 100c
;   include 'lib_dev:[gudoc]edfvaxbox.for/list'
;   baker, schlatter  17-may-1982     original version.
;-----------------------------------------------------------------------------
; es0 = saturation vapor ressure over liquid water at 0C

es0 = 6.1078
    
  pol = 0.99999683        + T9*(-0.90826951e-02 +      \
        T9*(0.78736169e-04 + T9*(-0.61117958e-06 +   \
        T9*(0.43884187e-08 + T9*(-0.29883885e-10 +   \ 
        T9*(0.21874425e-12 + T9*(-0.17892321e-14 +   \
        T9*(0.11112018e-16 + T9*(-0.30994571e-19)))))))))               
  
esw = es0/pol^8    ; 3D array (millibars)

;-------------------------------------------------------------------------
; from the fortran function wmr
;   this function approximates the mixing ratio wmr (grams of water
;   vapor per kilogram of dry air) given the pressure p (mb) and the
;   temperature t (celsius). the formula used is given on p. 302 of the
;   smithsonian meteorological tables by roland list (6th edition).
;   the next two lines contain a formula by herman wobus for the
;   correction factor wfw for the departure of the mixture of air
;   and water vapor from the ideal gas law. the formula fits values
;   in table 89, p. 340 of the smithsonian meteorological tables,
;   but only for temperatures and pressures normally encountered in
;   in the atmosphere.
;-------------------------------------------------------------------------
;  eps = ratio of the mean molecular weight of water (18.016 g/mole)
;                               to that of dry air (28.966 g/mole)
eps = 0.62197
x = 0.02*(T9-12.5+7500./pres)
wfw = 1.+ 4.5e-06*pres + 1.4e-03*x^2
fwesw = wfw*esw
rs = eps*fwesw/(pres-fwesw)      ; dimensionless (g/g or kg/kg) 2D array

;--------------------------------------------------------------------------
;  compute the dimensionless saturation specific humidity.


q  = doubletofloat(rs/(1.+rs))
   
q@long_name = "saturation specific humidity"
q@units = "kg/kg"
q@op_derive = "computed using the function ssh"


;************************************************************************

; compute the saturation moist static energy (smse) from the 
; model data on the hybrid levels
; t[*] array of temperatures(K) on the hybrid levels
; z[*] array of geopotential heights(m) of the hybrid levels
; p[*] array of pressures(mb or hPa) on the hybrid levels



; define constants
Cp = 1.00464e3 ; specific heat of dry air at constant pressure (J/(K*kg)) 
L = 2.501e6    ; latent heat of vaporization (J/kg)
g = 9.80616    ; acceleration due to gravity at sea level (m/s^2)

; compute the dry static energy on the hybrid levels (J/kg)
s = Cp*T9 + g*z

; compute the stauration specific humidity on the hybrid levels
qs = q

; compute the staturation moist static energy on the hybrid levels (J/kg)
; and convert to kJ/kg
hs = (s + L*qs)/1000.0
hs@long_name = "stat moist static energy"  
hs@units = "kJ/kg"


if (i.eq.2) then 
   hsarray1=hs
else
   hsarray(i,:)=hs
end if
delete(T9)
delete(pres)
delete(pol)
delete(s)
delete(hs)
delete(q)
delete(x)
delete(esw)
delete(wfw)
delete(fwesw)
delete(rs)
delete(qs)
end do
print(hsarray)
print(hsarray1)

elevs=f6->lev


 
 
wks   = gsn_open_wks ("png","dtdz_tropmean")                  ; send graphics to PNG file

res1                   = True                       ; plot mods desired


 res1@xyLineColor     = "blue"
 res1@trYReverse        = True                       ; reverse Y-axis
 res1@tiMainString    = "10S-10N tropical mean saturation moist static energy (kj/kg)"
 res1@tiMainFontHeightF=0.018
 res1@xyLineThicknesses = 4.0
 res1@tiXAxisString=""
 res1@tiXAxisFontHeightF=0.015
 res1@trYMaxF=1000
 res1@trYMinF=100
 res1@trXMaxF=90
 res1@trXMinF=60
 res1@tiYAxisString="Pressure (hPa)"
 res1@tiYAxisFontHeightF=0.018
 res1@xyYStyle = "Log"
 res1@tmYLMode = "Explicit"
 res1@tmYLValues        =(/100,200,300,400,500,600,700,850,1000/)
 res1@tmYLLabels=(/100,200,300,400,500,600,700,850,1000/)
 res1@tiXAxisString="CAM4-IF: red, CAM4: Blue, ERAI/MERRA: black"

 res1@tiXAxisStringFontHeightF=0.008
 res1@gsnDraw=False
 res1@gsnFrame=False
 res1@xyLineColor    = "black"
 ;plot               = gsn_csm_xy (wks,dTdZ(0:45),PmeanANN(0:45),res1) 
 ;res1@xyLineColor    = "orange"
 ;plot1              = gsn_csm_xy (wks,dTdZ4(0:45),PmeanANN(0:45),res1)
 res1@xyLineColor    = "blue"
 plot2= gsn_csm_xy (wks,hsarray(0,:),T4&lev,res1)
 res1@xyLineColor    = "red"
 plot3= gsn_csm_xy (wks,hsarray(1,:),T5&lev,res1)
 res1@xyLineColor    = "black"
 plot4= gsn_csm_xy (wks,hsarray1(0,:),elevs(::-1),res1)
 
 
 ;overlay(plot1,plot2)
 overlay(plot2,plot3)
 overlay(plot2,plot4)
 
draw(plot2)               
frame(wks) 


;crop=systemfunc("/usr/bin/pdfcrop dtdz_koror.pdf")

;crop=systemfunc("/usr/bin/pdfcrop dtdz_tropmean.png")
cleanup=systemfunc("rm tmp*.nc")
end