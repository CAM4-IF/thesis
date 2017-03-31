; ***********************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
;************************************************
;;;This script imports an ascii file with a radiosonde profile and compares it with output from the CAM model at the same location
begin
;************************************************
; read in data
;************************************************
o=0
p=0

system("rm data/*.eps")
latlons=asciiread("latlon",(/9,2/),"float")
print(latlons)
fils=systemfunc("ls data/*")
print(fils)
names=asciiread("names",(/9,1/),"string")
print(names)
numfils=dimsizes(fils)

sperrors=new((/numfils,5/),float)

m=0
print(m)
do while (m.le.numfils-1)

print(m)
lat=latlons(m,0)
lon=latlons(m,1)

o=round(((90+lat)/1.89473),3)
p=round((lon/2.5),3)
r=round(((90+lat)/2.5),3)



data1=asciiread(fils(m),(/576,10/),"float")
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

 s1   = static_stability(PmeanANN*100,TmeanANN,0,0)
 dp=dpres_plevel(PmeanANN*100,110000,PmeanANN(46),0)
 top=42
 bot=4
 fint=dim_sum_n(s1(top:bot)*dp(top:bot),0)/dim_sum_n(dp(top:bot),0)
 sperrors(m,0)=fint
f  = addfile("new_data/CAM-IFv5_AMIP_ANN_climo.nc","r")
f1  = addfile("AMIP_4K_control_ANN_climo1.nc","r")
f2=addfile("AMIP/ta_ANN_AMIPmean.nc","r")
f3=addfile("CAM_IF_AMIP4K_OCT24_ANN_climo.nc","r")
f5=addfile("CAM-IFv5_F2000_ANN_climo.nc","r")
f6=addfile("ta_Amon_AllModels_amip4K_ymonmean.nc","r")
f8=addfile("CAM-IFv4_AMIP_ANN_climo_det88.nc","r")
f10=addfile("timvar_proc_AMIP4K.nc","r")
f11=addfile("timvar_CAM4_monthly.nc","r")
f12=addfile("timvar_proc_CAM-IF1_AMIP4K.nc","r")
f13=addfile("timvar_CAM-IFv1_monthly.nc","r")
f14=addfile("timvar_proc_CAM-IF_AMIP4K_OCT24.nc","r")
f15=addfile("timvar_ANN_proc_CAM-IF_AMIP_OCT24.nc","r")



n=0
T6=new((/17,20/),float)
do while(n.le.19)
 f7=addfile("ta_ANN_allmodels_AMIP.nc","r")
 test=f7->ta(0,:,r,p)
 T6(:,n)=f7->ta(n,:,r,p)
 n=n+1
end do

plev=f7->plev ;; CMIP5 p levs

T     = f->T(0,:,:,:)
T1=f1->T(0,:,:,:)
T2=f2->ta(0,:,:,:)
T3=f3->T(0,:,:,:)
T4=f5->T(0,:,:,:)
T5=f6->ta(1,:,r,p)
T8=f8->T(0,:,:,:)
T10=f10->T(0,:,:,:)
T11=f11->T(0,:,:,:)
T12=f12->T(0,:,:,:)
T13=f13->T(0,:,:,:)
T14=f14->T(0,:,:,:)
T15=f15->T(0,:,:,:)

RH    = f->Q(0,:,:,:)
RH1=f1->Q(0,:,:,:)
RH2=f1->Q(0,:,:,:)
RH3=f3->Q(0,:,:,:)
RH4=f8->Q(0,:,:,:)

RH10=f10->Q(0,:,:,:)
RH11=f11->Q(0,:,:,:)
RH12=f12->Q(0,:,:,:)
RH13=f13->Q(0,:,:,:)
RH14=f14->Q(0,:,:,:)
RH15=f15->Q(0,:,:,:)

Z     = f->Z3(0,:,:,:)
Z1=f1->Z3(0,:,:,:)
Z2=f1->Z3(0,:,:,:)
Z3=f3->Z3(0,:,:,:)

Z10=f10->Z3(0,:,:,:)
Z11=f11->Z3(0,:,:,:)
Z12=f12->Z3(0,:,:,:)
Z13=f13->Z3(0,:,:,:)
Z14=f14->Z3(0,:,:,:)
Z15=f15->Z3(0,:,:,:)


U     = f->U(0,:,:,:)
U1=f1->U(0,:,:,:)
;U2=f2->U(0,:,:,:)
U3=f3->U(0,:,:,:)


V     = f->V(0,:,:,:)
V1=f1->V(0,:,:,:)
;V2=f2->V(0,:,:,:)
V3=f3->V(0,:,:,:)




;vinth2p requires the lev_p to be expressed in mb [hPa] (using AMWG levels - gives best vertical resolution)
lev           = PmeanANN         
lev!0         = "lev"                  ; variable/dim name 
lev&lev     =  lev                  ; create coordinate variable
lev@long_name = "pressure"               ; attach some attributes
lev@units     = "hPa"
lev@positive  = "down"

hyam            = f->hyam         ; read from file
hybm            = f->hybm
P0mb            = 1000
ps              = f->PS(0,:,:)       ; these are in P
u   = vinth2p(T,hyam,hybm,lev,ps,2,P0mb,1,False)
q=vinth2p(RH,hyam,hybm,lev,ps,2,P0mb,1,False)
z=vinth2p(Z,hyam,hybm,lev,ps,2,P0mb,1,False)
u2=vinth2p(U,hyam,hybm,lev,ps,2,P0mb,1,False)
v=vinth2p(V,hyam,hybm,lev,ps,2,P0mb,1,False)
pm    = lev




hyam1           = f1->hyam         ; read from file
hybm1            = f1->hybm
P0mb            = 1000
ps1              = f1->PS(0,:,:)       ; these are in P
u1   = vinth2p(T1,hyam1,hybm1,lev,ps1,2,P0mb,1,False)
q1=vinth2p(RH1,hyam1,hybm1,lev,ps1,2,P0mb,1,False)
z1=vinth2p(Z1,hyam1,hybm1,lev,ps1,2,P0mb,1,False)
u3=vinth2p(U1,hyam1,hybm1,lev,ps1,2,P0mb,1,False)
v1=vinth2p(V1,hyam1,hybm1,lev,ps1,2,P0mb,1,False)


;hyam2           = f2->hyam        ; read from file
;hybm2            = f2->hybm
;P0mb            = 1000
;ps2             = f2->PS(0,:,:)       ; these are in P
;uu   = vinth2p(T2,hyam2,hybm2,lev,ps2,2,P0mb,1,False)
;q2=vinth2p(RH2,hyam2,hybm2,lev,ps2,2,P0mb,1,False)
;z2=vinth2p(Z2,hyam2,hybm2,lev,ps2,2,P0mb,1,False)
;u4=vinth2p(U2,hyam2,hybm2,lev,ps2,2,P0mb,1,False)
;v2=vinth2p(V2,hyam2,hybm2,lev,ps2,2,P0mb,1,False)

uu=int2p_Wrap(plev*0.01,T2(:,r,p),lev,1)


hyam3            = f3->hyam         ; read from file
hybm3            = f3->hybm
P0mb            = 1000
ps3              = f3->PS(0,:,:)       ; these are in P
uu1   = vinth2p(T3,hyam3,hybm3,lev,ps3,2,P0mb,1,False)
q3=vinth2p(RH3,hyam3,hybm3,lev,ps3,2,P0mb,1,False)
z3=vinth2p(Z3,hyam3,hybm3,lev,ps3,2,P0mb,1,False)
u5=vinth2p(U3,hyam3,hybm3,lev,ps3,2,P0mb,1,False)
v3=vinth2p(V3,hyam3,hybm3,lev,ps3,2,P0mb,1,False)

hyam4            = f8->hyam         ; read from file
hybm4            = f8->hybm
P0mb            = 1000
ps4              = f8->PS(0,:,:)       ; these are in P
uu2  = vinth2p(T8,hyam4,hybm4,lev,ps4,2,P0mb,1,False)
q5=vinth2p(RH4,hyam4,hybm4,lev,ps4,2,P0mb,1,False)
;z3=vinth2p(Z3,hyam4,hybm4,lev,ps4,2,P0mb,1,False)
;u5=vinth2p(U3,hyam4,hybm4,lev,ps4,2,P0mb,1,False)
;v3=vinth2p(V3,hyam4,hybm4,lev,ps4,2,P0mb,1,False)


ps5              = f1->PS(0,:,:)       ; these are in P
u10  = vinth2p(T10,hyam4,hybm4,lev,ps5,2,P0mb,1,False)
q10=vinth2p(RH10,hyam4,hybm4,lev,ps5,2,P0mb,1,False)

ps6              = f1->PS(0,:,:)       ; these are in P
u11  = vinth2p(T11,hyam1,hybm1,lev,ps6,2,P0mb,1,False)
q11=vinth2p(RH11,hyam1,hybm1,lev,ps6,2,P0mb,1,False)

ps7             = f12->PS(0,:,:)       ; these are in P
u12  = vinth2p(T12,hyam4,hybm4,lev,ps7,2,P0mb,1,False)
q12=vinth2p(RH12,hyam4,hybm4,lev,ps7,2,P0mb,1,False)

ps8             = f13->PS(0,:,:)       ; these are in P
u13  = vinth2p(T13,hyam4,hybm4,lev,ps8,2,P0mb,1,False)
q13=vinth2p(RH13,hyam4,hybm4,lev,ps8,2,P0mb,1,False)

ps9            = f3->PS(0,:,:)       ; these are in P
u14  = vinth2p(T13,hyam4,hybm4,lev,ps9,2,P0mb,1,False)
q14=vinth2p(RH13,hyam4,hybm4,lev,ps9,2,P0mb,1,False)
;
ps10             = f->PS(0,:,:)       ; these are in P
u15  = vinth2p(T15,hyam4,hybm4,lev,ps10,2,P0mb,1,False)
q15=vinth2p(RH15,hyam4,hybm4,lev,ps10,2,P0mb,1,False)

 s2   = static_stability(lev*100,u(:,o,p),0,0)
 dp=dpres_plevel(lev*100,110000,PmeanANN(46),0)
 top=42
 bot=4
 fint1=dim_sum_n(s2(top:bot)*dp(top:bot),0)/dim_sum_n(dp(top:bot),0)
 sperrors(m,1)=(fint1-fint)

 s3   = static_stability(lev*100,u1(:,o,p),0,0)
 fint2=dim_sum_n(s3(top:bot)*dp(top:bot),0)/dim_sum_n(dp(top:bot),0)
 sperrors(m,2)=(fint2-fint)
 
 s4   = static_stability(lev*100,uu2(:,o,p),0,0)
 fint3=dim_sum_n(s4(top:bot)*dp(top:bot),0)/dim_sum_n(dp(top:bot),0)
 sperrors(m,3)=(fint3-fint)

;amip=int2p_Wrap(plev*0.01,T4,lev,1)
amip= vinth2p(T4,hyam1,hybm1,lev,ps1,2,P0mb,1,False)

 s5   = static_stability(lev*100,amip(:,o,p),0,0)
 fint4=dim_sum_n(s5(top:bot)*dp(top:bot),0)/dim_sum_n(dp(top:bot),0)
 sperrors(m,4)=(fint4-fint)

Tdiff=u(:,o,p)-TmeanANN
Tdiff1=u1(:,o,p)-TmeanANN
Tdiff2=uu-TmeanANN
Tdiff3=uu1(:,o,p)-TmeanANN
Tdiff4=uu2(:,o,p)-TmeanANN
Tdiff10=u10(:,o,p)-TmeanANN
Tdiff11=u11(:,o,p)-TmeanANN
Tdiff12=u12(:,o,p)-TmeanANN
Tdiff13=u13(:,o,p)-TmeanANN
Tdiff14=u14(:,o,p)-TmeanANN
Tdiff15=u15(:,o,p)-TmeanANN


;Tdiff5=amip-TmeanANN
Tdiff5=amip(:,o,p)-TmeanANN

;RHdiff=rh(:,0,0)-RHmeanANN
;RHdiff1=rh1(:,0,0)-RHmeanANN
;RHdiff2=rh1(:,0,0)-RHmeanANN
;RHdiff3=rh1(:,0,0)-RHmeanANN

;amip4k=int2p_Wrap(plev*0.01,T5,lev,1)
T6@_FillValue=1e20
amipmods=int2p_n_Wrap(plev*0.01,T6,lev,0,0)

amipTdiff=amip(:,o,p)-TmeanANN


tk = u-273.15
tk1=u1-273.15
tk2=uu-273.15
tk3=uu1-273.15
tk4=TmeanANN-273.15




;;Calculate RMSE;;
;RMSE1=dim_rmsd(u(0:45),TmeanANN(0:45))
;RMSE2=dim_rmsd(u1(0:45),TmeanANN(0:45))




;q    = mixhum_ptrh(lev, u(:,o,p), rh(:,o,p), 2)
;q@units = "kg/kg"
;q1    = mixhum_ptrh(lev, u1(:,o,p), rh1(:,o,p), 2)
;q1@units = "kg/kg"
;q2    = mixhum_ptrh(lev, uu(:,o,p), rh2(:,o,p), 2)
;q2@units = "kg/kg"
;q3    = mixhum_ptrh(lev, uu1(:,o,p), rh3(:,o,p), 2)
;q3@units = "kg/kg"
q4   = mixhum_ptrh(lev, TmeanANN, RHmeanANN, 2)
q4@units = "kg/kg"
;q5    = mixhum_ptrh(lev, uu2(:,o,p), rh4(:,o,p), 2)
;q5@units = "kg/kg"



;tdc  = dewtemp_trh(u(:,o,p),rh(:,o,p))-273.15
;tdc@units = "degC"            
;tdc1  = dewtemp_trh(u1(:,o,p),rh1(:,o,p))-273.15
;tdc1@units = "degC"
;tdc2  = dewtemp_trh(uu(:,o,p),rh2(:,o,p))-273.15
;tdc2@units = "degC"
;tdc3  = dewtemp_trh(uu1(:,o,p),rh3(:,o,p))-273.15
;tdc3@units = "degC"
;tdc4  = dewtemp_trh(TmeanANN,RHmeanANN)-273.15
;tdc4@units = "degC"


wspd = sqrt(u2^2 + v^2)
wspd1 = sqrt(u3^2 + v1^2)
wspd2 = sqrt(UmeanANN^2 + VmeanANN^2)

wdir = wind_direction(u2,v,0)
wdir1 = wind_direction(u3,v1,0)
wdir2= wind_direction(UmeanANN,VmeanANN,0)

;;;Calculate equivalent potential temeprature
numlevs=dimsizes(pm)
Lv   = 2.51e6  ; latent heat of vaporization at the triple point [J/kg]; Bluestein, p203
;;cpd  = 1005.7  ; specific heat at constant pressure for air [Bolton (1980): MWR]
cpd  = 1004.64 ; specific heat at constant pressure for air
R    = 287.04  ; specific gas constant for air [J/(Kg-K)]
kap  = R/cpd   ; 0.285

q    = q/(1-q) ;
q@long_name = "mixing ratio"
q1    = q1/(1-q1) ;
q1@long_name = "mixing ratio"

q2   = q1/(1-q1) ;
q2@long_name = "mixing ratio"
q3    = q3/(1-q3) ;
q3@long_name = "mixing ratio"
q4    = q4/(1-q4) ;
q4@long_name = "mixing ratio"
q5=q5/(1-q5)
q5@long_name = "mixing ratio"

t_e  = u(:,o,p)+(Lv/cpd)*q(:,o,p)         
t_e1  = u1(:,o,p)+(Lv/cpd)*q1(:,o,p)  
t_e2  = uu+(Lv/cpd)*q2(:,o,p)  
t_e3  = uu2(:,o,p)+(Lv/cpd)*q5(:,o,p)  
t_e4  = TmeanANN+(Lv/cpd)*q4

print(dimsizes(t_e4))
print(dimsizes(t_e))




P0=doubletofloat(P0mb)
theta_e=new((/numlevs/),float)
theta_e1=new((/numlevs/),float)
theta_e2=new((/numlevs/),float)
theta_e3=new((/numlevs/),float)
theta_e4=new((/numlevs/),float)

theta_e = t_e*(P0/pm)^kap    ; p0 and pm must be same units
theta_e1=t_e1*(P0/pm)^kap
theta_e2=t_e2*(P0/pm)^kap
theta_e3=t_e3*(P0/pm)^kap
theta_e4=t_e4*(P0/pm)^kap


;;Calculate d(theta_e)/dz profiles;;
ntimes=dimsizes(lev)
dTdZ=new(ntimes,double)
do i=0,ntimes-2
meanT=(theta_e(i)+theta_e(i+1))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(theta_e(i+1)-theta_e(i))
dTdZ(i)=dT/dZ
end do
dTdZ1=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(theta_e1(i)+theta_e1(i+1))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(theta_e1(i+1)-theta_e1(i))
dTdZ1(i)=dT/dZ
end do
dTdZ2=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(theta_e2(i)+theta_e2(i+1))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(theta_e2(i+1)-theta_e2(i))
dTdZ2(i)=dT/dZ
end do

dTdZ3=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(theta_e3(i)+theta_e3(i+1))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(theta_e3(i+1)-theta_e3(i))
dTdZ3(i)=dT/dZ
end do

dTdZ4=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(theta_e4(i)+theta_e4(i+1))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(theta_e4(i+1)-theta_e4(i))
dTdZ4(i)=dT/dZ
end do

dTdZdiff=dTdZ-dTdZ4
dTdZdiff1=dTdZ1-dTdZ4
dTdZdiff2=dTdZ3-dTdZ4



qdiff=(q(:,o,p)-q4)
qdiff1=(q1(:,o,p)-q4)
qdiff2=(q5(:,o,p)-q4)


;;;;Make series of plots of models vs. observed quantities (T,dT/dZ,Q,Te,dTe/dZ)

;;;T plot (CAM-IFv2,CAM-IFv1,CAM4,CMIP5/AMIP mean)

;;;Calculate RMSE/Bias first:

rmsearray=new(23,float)

tsum=decimalPlaces((sum(Tdiff(1:46))/45),3,True)
tsum1=decimalPlaces((sum(Tdiff1(1:46))/45),3,True)
tsum2=decimalPlaces((sum(Tdiff4(1:46))/45),3, True)



wks   = gsn_open_wks ("eps",  fils(m) + "Tobs")                  ; send graphics to eps file


res                   = True                       ; plot mods desired
res@gsnDraw=False
res@gsnFrame=False
res@trYMaxF=1000
res@trYMinF=80
res@trXMaxF=4
res@trXMinF=-4
res@xyLineColor     = "black"
res@xyYStyle = "Log"
res@gsnXRefLine=0
res@tmYLMode = "Explicit"
res@tmYLValues        =(/100,200,300,400,500,600,700,850,1000/)
res@tmYLLabels=(/100,200,300,400,500,600,700,850,1000/)
res@xyLineColor    = "grey58"

amipmoddiff=amipmods(:,0)-TmeanANN
print(amipmoddiff)
plot4              = gsn_csm_xy(wks,amipmoddiff(1:46),u&lev(1:46),res)
amipmoddiff1=amipmods(:,1)-TmeanANN
plot5              = gsn_csm_xy(wks,amipmoddiff1(1:46),u&lev(1:46),res)
amipmoddiff2=amipmods(:,2)-TmeanANN
plot6              = gsn_csm_xy(wks,amipmoddiff2(1:46),u&lev(1:46),res)
amipmoddiff3=amipmods(:,3)-TmeanANN
plot7              = gsn_csm_xy(wks,amipmoddiff3(1:46),u&lev(1:46),res)
amipmoddiff4=amipmods(:,4)-TmeanANN
plot8              = gsn_csm_xy(wks,amipmoddiff4(1:46),u&lev(1:46),res)
amipmoddiff5=amipmods(:,5)-TmeanANN
plot9             = gsn_csm_xy(wks,amipmoddiff5(1:46),u&lev(1:46),res)
amipmoddiff6=amipmods(:,6)-TmeanANN
plot10              = gsn_csm_xy(wks,amipmoddiff6(1:46),u&lev(1:46),res)
amipmoddiff7=amipmods(:,7)-TmeanANN
plot11              = gsn_csm_xy(wks,amipmoddiff7(1:46),u&lev(1:46),res)
amipmoddiff8=amipmods(:,8)-TmeanANN
plot12              = gsn_csm_xy(wks,amipmoddiff8(1:46),u&lev(1:46),res)
amipmoddiff9=amipmods(:,9)-TmeanANN
plot13              =gsn_csm_xy(wks,amipmoddiff9(1:46),u&lev(1:46),res)
amipmoddiff10=amipmods(:,10)-TmeanANN
plot14              = gsn_csm_xy(wks,amipmoddiff10(1:46),u&lev(1:46),res)
amipmoddiff11=amipmods(:,11)-TmeanANN
plot15              = gsn_csm_xy(wks,amipmoddiff11(1:46),u&lev(1:46),res)
amipmoddiff12=amipmods(:,12)-TmeanANN
plot16              =gsn_csm_xy(wks,amipmoddiff12(1:46),u&lev(1:46),res)
amipmoddiff13=amipmods(:,13)-TmeanANN
plot17              =gsn_csm_xy(wks,amipmoddiff13(1:46),u&lev(1:46),res)
amipmoddiff14=amipmods(:,14)-TmeanANN
plot18              = gsn_csm_xy(wks,amipmoddiff14(1:46),u&lev(1:46),res)
amipmoddiff15=amipmods(:,15)-TmeanANN
plot19              = gsn_csm_xy(wks,amipmoddiff15(1:46),u&lev(1:46),res)
amipmoddiff16=amipmods(:,16)-TmeanANN
plot20              = gsn_csm_xy(wks,amipmoddiff16(1:46),u&lev(1:46),res)
amipmoddiff17=amipmods(:,17)-TmeanANN
plot21              = gsn_csm_xy(wks,amipmoddiff17(1:46),u&lev(1:46),res)
amipmoddiff18=amipmods(:,18)-TmeanANN
plot22              = gsn_csm_xy(wks,amipmoddiff18(1:46),u&lev(1:46),res)
amipmoddiff19=amipmods(:,19)-TmeanANN
plot23              = gsn_csm_xy(wks,amipmoddiff19(1:46),u&lev(1:46),res)
sqe=Tdiff^2
sqe1=Tdiff1^2
sqe2=Tdiff4^2
sqe3=amipmoddiff^2
sqe4=amipmoddiff1^2
sqe5=amipmoddiff2^2
sqe6=amipmoddiff3^2
sqe7=amipmoddiff4^2
sqe8=amipmoddiff5^2
sqe9=amipmoddiff6^2
sqe10=amipmoddiff7^2
sqe11=amipmoddiff8^2
sqe12=amipmoddiff9^2
sqe13=amipmoddiff10^2
sqe14=amipmoddiff11^2
sqe15=amipmoddiff12^2
sqe16=amipmoddiff13^2
sqe17=amipmoddiff14^2
sqe18=amipmoddiff15^2
sqe19=amipmoddiff16^2
sqe20=amipmoddiff17^2
sqe21=amipmoddiff18^2
sqe22=amipmoddiff19^2


rmsearray(0)=decimalPlaces(sqrt(avg(sqe(1:46))),3,True)
rmsearray(1)=decimalPlaces(sqrt(avg(sqe1(1:46))),3,True)
rmsearray(2)=decimalPlaces(sqrt(avg(sqe2(1:46))),3,True)
rmsearray(3)=decimalPlaces(sqrt(avg(sqe3(1:46))),3,True)
rmsearray(4)=decimalPlaces(sqrt(avg(sqe4(1:46))),3,True)
rmsearray(5)=decimalPlaces(sqrt(avg(sqe5(1:46))),3,True)
rmsearray(6)=decimalPlaces(sqrt(avg(sqe6(1:46))),3,True)
rmsearray(7)=decimalPlaces(sqrt(avg(sqe7(1:46))),3,True)
rmsearray(8)=decimalPlaces(sqrt(avg(sqe8(1:46))),3,True)
rmsearray(9)=decimalPlaces(sqrt(avg(sqe9(1:46))),3,True)
rmsearray(10)=decimalPlaces(sqrt(avg(sqe10(1:46))),3,True)
rmsearray(11)=decimalPlaces(sqrt(avg(sqe11(1:46))),3,True)
rmsearray(12)=decimalPlaces(sqrt(avg(sqe12(1:46))),3,True)
rmsearray(13)=decimalPlaces(sqrt(avg(sqe13(1:46))),3,True)
rmsearray(14)=decimalPlaces(sqrt(avg(sqe14(1:46))),3,True)
rmsearray(15)=decimalPlaces(sqrt(avg(sqe15(1:46))),3,True)
rmsearray(16)=decimalPlaces(sqrt(avg(sqe16(1:46))),3,True)
rmsearray(17)=decimalPlaces(sqrt(avg(sqe17(1:46))),3,True)
rmsearray(18)=decimalPlaces(sqrt(avg(sqe18(1:46))),3,True)
rmsearray(19)=decimalPlaces(sqrt(avg(sqe19(1:46))),3,True)
rmsearray(20)=decimalPlaces(sqrt(avg(sqe20(1:46))),3,True)
rmsearray(21)=decimalPlaces(sqrt(avg(sqe21(1:46))),3,True)
rmsearray(22)=decimalPlaces(sqrt(avg(sqe22(1:46))),3,True)

quarts=dim_pqsort_n(rmsearray,2,0)
print(quarts)
index=(ind(quarts.eq.0))+1
print(index)

res@xyLineColor     = "blue"
res@trYReverse        = True                     ; reverse Y-axis
res@tiMainString    =  "ANN mean temperature (model-observed) (K):"
res@gsnRightString  ="~C~CAM4 RMSE: " + rmsearray(1) + "~C~CAM4-IF RMSE: " + rmsearray(0) + " rank: " + index + "~C~CMIP5/AMIP RMSE: " + rmsearray(2)
res@gsnRightStringFontHeightF=0.008
res@gsnLeftString  =names(m,0) + ": (" + lon + "E, " + lat + "N)"
res@gsnLeftStringFontHeightF=0.018
res@tiMainFontHeightF=0.018
res@xyLineThicknesses = 8.0
res@tiXAxisString="CAM4-IF: red, CAM4: Blue, CMIP5 AMIP: orange, CMIP5 AMIP members: black"
res@tiXAxisFontHeightF=0.012
res@tiYAxisString="Pressure (hPa)"
res@tiYAxisFontHeightF=0.018

res@xyLineThicknesses = 8.0
res@xyLineColor    = "red"
plot               = gsn_csm_xy(wks,Tdiff(1:46),u&lev(1:46),res) 
res@xyLineColor    = "blue"
plot1              = gsn_csm_xy(wks,Tdiff1(1:46),u&lev(1:46),res)
res@xyLineColor    = "green"
;plot2              =gsn_csm_xy(wks,Tdiff4(1:46),u&lev(1:46),res)
res@xyLineColor    = "orange"
plot3              = gsn_csm_xy(wks,Tdiff2(1:46),u&lev(1:46),res)
res@xyLineColor    = "grey"
res@xyLineThicknesses = 1.0

mres               = True                     ; marker mods desired
mres@gsMarkerIndex =7                      ; polymarker style
mres@gsMarkerSizeF = 1.                      ; polymarker size
mres@gsMarkerColor = "red"                    ; polymarker color

;dum=new((ntimes-2),typeof(plot))
;error_bar=new((ntimes-2),typeof(plot))
;do i=1,ntimes-3
;dum(i) = gsn_add_polymarker(wks,plot,Tdiff(i),lev(i),mres)
;mres@gsLineColor = "red"                    ; polymarker color
;error_bar(i) = gsn_add_polyline(wks,plot,(/Tdiff(i)+2.056*(sqrt(u15(i,o,p))/sqrt(27)),Tdiff(i)-2.056*(sqrt(u15(i,o,p))/sqrt(27))/),(/lev(i),lev(i)/),mres)
;end do

;dum1=new((ntimes-2),typeof(plot1))
;error_bar1=new((ntimes-2),typeof(plot))
;do i=1,ntimes-3
;dum1(i) = gsn_add_polymarker(wks,plot1,Tdiff1(i),lev(i),mres)
;mres@gsLineColor = "blue"                    ; polymarker color
;error_bar1(i) = gsn_add_polyline(wks,plot1,(/Tdiff1(i)+2.056*(sqrt(u11(i,o,p))/sqrt(27)),Tdiff1(i)-2.056*(sqrt(u11(i,o,p))/sqrt(27))/),(/lev(i),lev(i)/),mres)
;end do



overlay(plot,plot1)
;overlay(plot,plot2)
overlay(plot,plot3)
overlay(plot,plot4)
overlay(plot,plot5)
overlay(plot,plot6)
overlay(plot,plot7)
overlay(plot,plot8)
overlay(plot,plot9)
overlay(plot,plot10)
overlay(plot,plot11)
overlay(plot,plot12)
overlay(plot,plot13)
overlay(plot,plot14)
overlay(plot,plot15)
overlay(plot,plot16)
overlay(plot,plot17)
overlay(plot,plot18)
overlay(plot,plot19)
overlay(plot,plot20)
overlay(plot,plot21)
overlay(plot,plot22)
overlay(plot,plot23)
draw(plot)
;gsn_polyline(wks,plot,(/0,0/),(/0,1000/),False)                  
frame(wks)



;;;dTdZ plot (CAM-IFv2,CAM-IFv1,CAM4,CMIP5/AMIP mean) (u,u1,uu2,amip,TmeanANN)

;;Calculate dT/dZ profiles;;
ntimes=dimsizes(lev)
dTdZ=new(ntimes,double)
do i=0,ntimes-2
meanT=(u(i,o,p)+u(i+1,o,p))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(u(i+1,o,p)-u(i,o,p))
dTdZ(i)=dT/dZ
end do

dTdZ1=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(u1(i,o,p)+u1(i+1,o,p))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(u1(i+1,o,p)-u1(i,o,p))
dTdZ1(i)=dT/dZ
end do

dTdZ2=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(uu(i)+uu(i+1))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(uu(i+1)-uu(i))
dTdZ2(i)=dT/dZ
end do

;dTdZ3=new(dimsizes(lev),double)
;do i=0,ntimes-2
;meanT=(amip(i)+amip(i+1))/2
;dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
;dT=(amip(i+1)-amip(i))
;dTdZ3(i)=dT/dZ
;end do
dTdZ3=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(amip(i,o,p)+amip(i+1,o,p))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(amip(i+1,o,p)-amip(i,o,p))
dTdZ3(i)=dT/dZ
end do

dTdZ4=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(TmeanANN(i)+TmeanANN(i+1))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dT=(TmeanANN(i+1)-TmeanANN(i))
dTdZ4(i)=dT/dZ
end do

dTdZdiff=dTdZ-dTdZ4
dTdZdiff1=dTdZ1-dTdZ4
dTdZdiff2=dTdZ2-dTdZ4
dTdZdiff3=dTdZ3-dTdZ4

;;;Calculate RMSE/Bias first:

sqe24=dTdZdiff^2
sqe25=dTdZdiff1^2
sqe26=dTdZdiff2^2
sqe27=dTdZdiff3^2

rmse4=decimalPlaces(sqrt(avg(sqe24(1:46))),3,True)
rmse5=decimalPlaces(sqrt(avg(sqe25(1:46))),3,True)
rmse6=decimalPlaces(sqrt(avg(sqe26(1:46))),3,True)


tsum4=decimalPlaces((sum(dTdZdiff(1:46))/45),3,True)
tsum5=decimalPlaces((sum(dTdZdiff1(1:46))/45),3,True)
tsum6=decimalPlaces((sum(dTdZdiff2(1:46))/45),3, True)

wks   = gsn_open_wks ("eps",fils(m) + "dTdZobs")                  ; send graphics to eps file

res1                   = True                       ; plot mods desired
res1@xyLineColor     = "blue"
res1@trYReverse        = True                     ; reverse Y-axis
res1@tiMainString    = names(m,0) + ": (" + lon + "E, " + lat + "N) ANN mean"
res1@gsnRightString  ="CAM4-IF RMSE: " + rmse4 +"~C~CAM4 RMSE: " + rmse5 + "~C~CMIP5 AMIP mean RMSE: " + rmse6
res1@gsnRightStringFontHeightF=0.008
res1@gsnLeftString  ="CAM4-IF Bias: " + tsum4 +"~C~CAM4 Bias: " + tsum5 + "~C~CMIP5 AMIP mean Bias: " + tsum6
res1@gsnLeftStringFontHeightF=0.008
res1@tiMainFontHeightF=0.018
res1@xyLineThicknesses = 4.0
res1@tiXAxisString="Temperature lapse rate (model-observed) (K/km): CAM4-IF: red, CAM4: Blue, CMIP5 AMIP: orange"
res1@tiXAxisFontHeightF=0.01
res1@tiYAxisString="Pressure (hPa)"
res1@tiYAxisFontHeightF=0.018

res1@gsnDraw=False
res1@gsnFrame=False
res1@trYMaxF=1000
res1@trYMinF=50
res1@trXMaxF=4
res1@trXMinF=-4
res1@xyLineColor     = "red"
plot               = gsn_csm_xy (wks,dTdZdiff(1:46),u&lev(1:46),res1) 
res1@xyLineColor    = "blue"
plot1              = gsn_csm_xy (wks,dTdZdiff1(1:46),u&lev(1:46),res1)
res1@xyLineColor    = "orange"
plot2              = gsn_csm_xy (wks,dTdZdiff2(1:46),u&lev(1:46),res1)
res1@xyLineColor    = "green"
;plot3              = gsn_csm_xy (wks,dTdZdiff3(1:46),u&lev(1:46),res1)
res1@xyLineColor    = "cyan"
overlay(plot,plot1)
overlay(plot,plot2)
;overlay(plot,plot3)
;overlay(plot,plot4)
;overlay(plot,plot5)
draw(plot)
gsn_polyline(wks,plot,(/0,0/),(/0,1000/),False)                  
frame(wks)






;;;Q plot (CAM-IFv2,CAM-IFv1,CAM4) (qdiff,qdiff1,qdiff2)
;;;Calculate RMSE/Bias first:

sqe8=qdiff^2
sqe9=qdiff1^2
sqe10=qdiff2^2

rmse8=sqrt(avg(sqe8(1:46)))
rmse9=sqrt(avg(sqe9(1:46)))
rmse10=sqrt(avg(sqe10(1:46)))


tsum8=sum(qdiff(1:46))/45
tsum9=sum(qdiff1(1:46))/45
tsum10=sum(qdiff2(1:46))/45

wks   = gsn_open_wks ("eps",fils(m) + "Qobs")                  ; send graphics to eps file

res1@gsnRightString  ="CAM-IFv2 RMSE: " + rmse8 +"~C~CAM-IFv1 RMSE: " + rmse10 +"~C~CAM4 RMSE: " + rmse9 ;+ "~C~CAM-IFv3 RMSE: " + rmse7
res1@gsnRightStringFontHeightF=0.008
res1@gsnLeftString  ="CAM-IFv2 Bias: " + tsum8 +"~C~CAM-IFv1 Bias: " + tsum10 + "~C~CAM4 Bias: " + tsum9 ;+ "~C~CAM-IFv3 Bias: " + tsum7
res1@gsnLeftStringFontHeightF=0.008
res1@tiMainFontHeightF=0.018
res1@xyLineThicknesses = 4.0
res1@tiXAxisString="Q (model-observed) (kg/kg*1000): CAM-IFv2: red, CAM-IFv1: green, CAM4: blue"
res1@tiXAxisFontHeightF=0.012
res1@tiYAxisString="Pressure (hPa)"
res1@tiYAxisFontHeightF=0.018

res1@xyLineColor     = "red"
plot               = gsn_csm_xy (wks,qdiff(1:46)*1000,u&lev(1:46),res1) 
res@xyLineColor    = "blue"
plot1              = gsn_csm_xy (wks,qdiff1(1:46)*1000,u&lev(1:46),res1)
res@xyLineColor    = "green"
plot2              = gsn_csm_xy (wks,qdiff2(1:46)*1000,u&lev(1:46),res1)



dum2=new((ntimes-2),typeof(plot))
error_bar2=new((ntimes-2),typeof(plot))
do i=1,ntimes-3
dum2(i) = gsn_add_polymarker(wks,plot,qdiff(i)*1000,lev(i),mres)
mres@gsLineColor = "red"                    ; polymarker color
error_bar2(i) = gsn_add_polyline(wks,plot,(/qdiff(i)*1000+(2056*(sqrt(q15(i,o,p))/sqrt(27))),qdiff(i)*1000-2056*(sqrt(q15(i,o,p))/sqrt(27))/),(/lev(i),lev(i)/),mres)
end do

dum3=new((ntimes-2),typeof(plot1))
error_bar3=new((ntimes-2),typeof(plot))
do i=1,ntimes-3
dum3(i) = gsn_add_polymarker(wks,plot1,qdiff1(i)*1000,lev(i),mres)
mres@gsLineColor = "blue"                    ; polymarker color
;error_bar3(i) = gsn_add_polyline(wks,plot1,(/qdiff1(i)*1000+2056*(sqrt(q11(i,o,p))/sqrt(27)),qdiff1(i)*1000-2056*(sqrt(q11(i,o,p))/sqrt(27))/),(/lev(i),lev(i)/),mres)
end do


overlay(plot,plot1)
overlay(plot,plot2)
draw(plot)
gsn_polyline(wks,plot,(/0,0/),(/0,1000/),False)                  
frame(wks)

;;;Te plot (CAM-IFv2,CAM-IFv1,CAM4) (theta_e,theta_e1,theta_e3)
tediff=theta_e-theta_e4
tediff1=theta_e1-theta_e4
tediff2=theta_e3-theta_e4

sqe11=tediff^2
sqe12=tediff1^2
sqe13=tediff2^2

rmse11=sqrt(avg(sqe11(1:46)))
rmse12=sqrt(avg(sqe12(1:46)))
rmse13=sqrt(avg(sqe13(1:46)))


tsum11=sum(tediff(1:46))/45
tsum12=sum(tediff1(1:46))/45
tsum13=sum(tediff2(1:46))/45

wks   = gsn_open_wks ("eps",fils(m) + "Teobs")                  ; send graphics to eps file

res@gsnRightString  ="CAM-IFv2 RMSE: " + rmse11 +"~C~CAM-IFv1 RMSE: " + rmse13 +"~C~CAM4 RMSE: " + rmse12 ;+ "~C~CAM-IFv3 RMSE: " + rmse7
res@gsnRightStringFontHeightF=0.008
res@gsnLeftString  ="CAM-IFv2 Bias: " + tsum11 +"~C~CAM-IFv1 Bias: " + tsum13 + "~C~CAM4 Bias: " + tsum12 ;+ "~C~CAm-IFv3 Bias: " + tsum7
res@gsnLeftStringFontHeightF=0.008
res@tiMainFontHeightF=0.018
res@xyLineThicknesses = 4.0
res@tiXAxisString="Theta_e (model-observed) (K): CAM-IFv2: red, CAM-IFv1: green, CAM4: blue"
res@tiXAxisFontHeightF=0.012
res@tiYAxisString="Pressure (hPa)"
res@tiYAxisFontHeightF=0.018

res@trXMaxF=8
res@trXMinF=-8

res@xyLineColor     = "red"
plot               = gsn_csm_xy (wks,tediff(1:46),u&lev(1:46),res) 
res@xyLineColor    = "blue"
plot1              = gsn_csm_xy (wks,tediff1(1:46),u&lev(1:46),res)
res@xyLineColor    = "green"
plot2              = gsn_csm_xy (wks,tediff2(1:46),u&lev(1:46),res)







overlay(plot,plot1)
overlay(plot,plot2)
draw(plot)
gsn_polyline(wks,plot,(/0,0/),(/0,1000/),False)                  
frame(wks)




;;;dTedZ plot (CAM-IFv2,CAM-IFv1,CAM4) (theta_e,theta_e1,theta_e3)


;;Calculate dTe/dZ profiles;;
ntimes=dimsizes(lev)
dTedZ=new(ntimes,double)
do i=0,ntimes-2
meanT=(u(i,o,p)+u(i+1,o,p))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dTe=(theta_e(i+1)-theta_e(i))
dTedZ(i)=dTe/dZ
end do

dTedZ1=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(u1(i,o,p)+u1(i+1,o,p))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dTe=(theta_e1(i+1)-theta_e1(i))
dTedZ1(i)=dTe/dZ
end do

dTedZ2=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(uu2(i,o,p)+uu2(i+1,o,p))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dTe=(theta_e3(i+1)-theta_e3(i))
dTedZ2(i)=dTe/dZ
end do

dTedZ4=new(dimsizes(lev),double)
do i=0,ntimes-2
meanT=(TmeanANN(i)+TmeanANN(i+1))/2
dZ=((287*meanT/9.81)*log(lev(i)/lev(i+1)))/1000
dTe=(theta_e4(i+1)-theta_e4(i))
dTedZ4(i)=dTe/dZ
end do

dtedzdiff=dTedZ-dTedZ4
dtedzdiff1=dTedZ1-dTedZ4
dtedzdiff2=dTedZ2-dTedZ4

sqe34=dtedzdiff^2
sqe35=dtedzdiff1^2
sqe36=dtedzdiff2^2

rmse14=sqrt(avg(sqe34(1:46)))
rmse15=sqrt(avg(sqe35(1:46)))
rmse16=sqrt(avg(sqe36(1:46)))


tsum14=sum(dtedzdiff(1:46))/45
tsum15=sum(dtedzdiff1(1:46))/45
tsum16=sum(dtedzdiff2(1:46))/45


wks   = gsn_open_wks ("eps",fils(m) + "dTedZobs")                  ; send graphics to eps file

res@gsnRightString  ="CAM-IFv2 RMSE: " + rmse14 +"~C~CAM-IFv1 RMSE: " + rmse16 +"~C~CAM4 RMSE: " + rmse15 ;+ "~C~CAM-IFv3 RMSE: " + rmse7
res@gsnRightStringFontHeightF=0.008
res@gsnLeftString  ="CAM-IFv2 Bias: " + tsum14 +"~C~CAM-IFv1 Bias: " + tsum16 + "~C~CAM4 Bias: " + tsum15 ;+ "~C~CAM-IFv3 Bias: " + tsum7
res@gsnLeftStringFontHeightF=0.008
res@tiMainFontHeightF=0.018
res@xyLineThicknesses = 4.0
res@tiXAxisString="d(theta_e)/dZ (model-observed) (K/km): CAM-IFv2: red, CAM-IFv1: green, CAM4: blue"
res@tiXAxisFontHeightF=0.012
res@tiYAxisString="Pressure (hPa)"
res@tiYAxisFontHeightF=0.018

res@trXMaxF=8
res@trXMinF=-8

res@xyLineColor     = "red"
plot               = gsn_csm_xy (wks,dtedzdiff(1:46),u&lev(1:46),res) 
res@xyLineColor    = "blue"
plot1              = gsn_csm_xy (wks,dtedzdiff1(1:46),u&lev(1:46),res)
res@xyLineColor    = "green"
plot2              = gsn_csm_xy (wks,dtedzdiff2(1:46),u&lev(1:46),res)
overlay(plot,plot1)
overlay(plot,plot2)
draw(plot)
gsn_polyline(wks,plot,(/0,0/),(/0,1000/),False)                  
frame(wks)







;;;draw Skew-Ts;;;
;skewtOpts                 = True
;skewtOpts@DrawColAreaFill = True    ; default is False
;skewtOpts@DrawFahrenheit=False
;skewtOpts@DrawColAreaFill=True
;skewtOpts@DrawHeightScale=True
   
;dataOpts   = True
;dataOpts@DrawFahrenheit= False
;dataOpts@PrintZ = True
;dataOpts@DrawColAreaFill=True
;dataOpts@ThermoInfo=True

;skewtOpts@tiMainString = "Koror"


;wks  = gsn_open_wks ("eps", "CAM-IF")
;skewt_bkgd = skewT_BackGround(wks, skewtOpts)
;skewt_data = skewT_PlotData(wks,skewt_bkgd,lev,tk(:,o,p),tdc,z(:,o,p),wspd(:,o,p),wdir(:,o,p),dataOpts)


;draw(skewt_bkgd)
;frame(wks)
;wks  = gsn_open_wks ("eps", "CAM-IF4K")
;skewt_bkgd1 = skewT_BackGround(wks, skewtOpts)
;skewt_data1 = skewT_PlotData(wks,skewt_bkgd1,lev,tk3(:,o,p),tdc3,z3(:,o,p),wspd(:,o,p),wdir(:,o,p),dataOpts)
;draw(skewt_bkgd)
;frame(wks)

m=m+1
delete(lat)
delete(lon)
delete(lev)
delete(res)
;delete(amip4k)
delete(amipmods)
end do
print(sperrors)
system("mv data/*.eps plots/")
system("gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=plots/profiles.pdf plots/*.eps")
end