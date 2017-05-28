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
latlons=asciiread("../../latlon",(/9,2/),"float")
lat=latlons(3,0)
lon=latlons(3,1)

o=round(((90+lat)/1.89473),3)
p=round((lon/2.5),3)
r=round(((90+lat)/2.5),3)

data1=asciiread("../KOR",(/576,10/),"float")
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


f=addfile("../../AMIP/ta_ANN_AMIPmean.nc","r")
f1=addfile("../../AMIP_4K_control_ANN_climo1.nc","r")
f2  = addfile("../CAM-IFv5_AMIP_ANN_climo.nc","r")
f8  = addfile("../CAM-IFv5_AMIP_br_ANN_climo.nc","r")
;f7=addfile("RAOBS.nc","r")
cmd="cdo -fldmean -sellonlatbox,0,360,-20,20 -selname,ta ../../AMIP/ta_ANN_AMIPmean.nc" + " tmp.nc"
cmd1="cdo -fldmean -sellonlatbox,0,360,-20,20 -selname,T,PS,Z3 ../../AMIP_4K_control_ANN_climo1.nc" + " tmp1.nc"
cmd2="cdo -fldmean -sellonlatbox,0,360,-20,20 -selname,T,PS,Z3 ../CAM-IFv5_AMIP_ANN_climo.nc" + " tmp2.nc"
cmd3="cdo -fldmean -sellonlatbox,0,360,-20,20 -selname,T,Z3 ../MERGED_T_trop_ANN.nc" + " tmp3.nc"
cmd4="cdo -fldmean -sellonlatbox,0,360,-20,20 -selname,T,PS,Z3 ../CAM-IFv5_AMIP_br_ANN_climo.nc" + " tmp4.nc"

areamean=systemfunc(cmd)
areamean1=systemfunc(cmd1)
areamean2=systemfunc(cmd2)
areamean3=systemfunc(cmd3)
areamean4=systemfunc(cmd4)
f3=addfile("tmp.nc","r")
f4=addfile("tmp1.nc","r")
f5=addfile("tmp2.nc","r")
f6=addfile("tmp3.nc","r")
f7=addfile("tmp4.nc","r")

T3=f3->ta(0,:,:,:)
T4=f4->T(0,:,:,:)
T5=f5->T(0,:,:,:)
T6=f6->T(0,:,:,:)
T7=f7->T(0,:,:,:)
Zmerra=f6->Z3(0,:,:,:)
Zcam4=f4->Z3(0,:,:,:)
Zcamif=f5->Z3(0,:,:,:)
Zcamifr=f7->Z3(0,:,:,:)
ps3=f4->PS(0,:,:)
ps4=f5->PS(0,:,:)
ps5=f7->PS(0,:,:)
;ps5=f6->PS(0,:,:)

;test=f7->T(:,:,:)
;print(test)

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


dTdZ1=new(dimsizes(PmeanANN),double)
do i=0,ntimes-2
meanT=(u(i,r,p)+u(i+1,r,p))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(u(i+1,r,p)-u(i,r,p))
dTdZ1(i)=dT/dZ
end do

dTdZ2=new(dimsizes(PmeanANN),double)
do i=0,ntimes-2
meanT=(u1(i,o,p)+u1(i+1,o,p))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(u1(i+1,o,p)-u1(i,o,p))
dTdZ2(i)=dT/dZ
end do

dTdZ3=new(dimsizes(PmeanANN),double)
do i=0,ntimes-2
meanT=(u2(i,o,p)+u2(i+1,o,p))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(u2(i+1,o,p)-u2(i,o,p))
dTdZ3(i)=dT/dZ
end do

dtdz1=new(dimsizes(T(:,30,30)),double)
do i=0,dimsizes(T(:,30,30))-2
meanT=(T(i,r,p)+T(i+1,r,p))/2
dz=((287*meanT/9.81)*log(T&plev(i)/T&plev(i+1)))/1000
dt=(T(i+1,r,p)-T(i,r,p))
dtdz1(i)=dt/dz
end do

delete(dz)
dtdz2=new(dimsizes(T1(:,0,0)),double)
do i=0,dimsizes(T1(:,0,0))-2
meanT=(T1(i,o,p)+T1(i+1,o,p))/2
dz=((287*meanT/9.81)*log(T1&lev(i)/T1&lev(i+1)))/1000
dt=(T1(i+1,o,p)-T1(i,o,p))
dtdz2(i)=dt/dz
end do

dtdz3=new(dimsizes(T2(:,0,0)),double)
do i=0,dimsizes(T2(:,0,0))-2
meanT=(T2(i,o,p)+T2(i+1,o,p))/2
dz=((287*meanT/9.81)*log(T2&lev(i)/T2&lev(i+1)))/1000
dt=(T2(i+1,o,p)-T2(i,o,p))
dtdz3(i)=dt/dz
end do




T3@_FillValue=1e20
u3=int2p_n_Wrap(plev*0.01,T3,PmeanANN,0,0)
u4=vinth2p(T4,hyam,hybm,PmeanANN,ps3,2,1000,1,False)
u5=vinth2p(T5,hyam,hybm,PmeanANN,ps4,2,1000,1,False)
u7=vinth2p(T7,hyam,hybm,PmeanANN,ps5,2,1000,1,False)
elevs=f6->lev
u6=int2p_n_Wrap(elevs,T6,PmeanANN,0,0)



dTdZ4=new(dimsizes(PmeanANN),double)
do i=0,ntimes-2
meanT=(u3(i,0,0)+u3(i+1,0,0))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(u3(i+1,0,0)-u3(i,0,0))
dTdZ4(i)=dT/dZ
end do

dTdZ5=new(dimsizes(PmeanANN),double)
do i=0,ntimes-2
meanT=(u4(i,0,0)+u4(i+1,0,0))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(u4(i+1,0,0)-u4(i,0,0))
dTdZ5(i)=dT/dZ
end do

dTdZ6=new(dimsizes(PmeanANN),double)
do i=0,ntimes-2
meanT=(u5(i,0,0)+u5(i+1,0,0))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(u5(i+1,0,0)-u5(i,0,0))
dTdZ6(i)=dT/dZ
end do
delete(dZ)

dTdZ8=new(dimsizes(PmeanANN),double)
do i=0,ntimes-2
meanT=(u7(i,0,0)+u7(i+1,0,0))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(u7(i+1,0,0)-u7(i,0,0))
dTdZ8(i)=dT/dZ
end do
delete(dZ)

dTdZ7=new(dimsizes(PmeanANN),double)
do i=0,ntimes-2
meanT=(u6(i,0,0)+u6(i+1,0,0))/2
dZ=((287*meanT/9.81)*log(PmeanANN(i)/PmeanANN(i+1)))/1000
dT=(u6(i+1,0,0)-u6(i,0,0))
dTdZ7(i)=dT/dZ
end do

dtdz7=new(dimsizes(T6(:,0,0)),double)
do i=0,dimsizes(T6(:,0,0))-2
dz=(Zmerra(i+1,0,0)-Zmerra(i,0,0))/1000
dt=(T6(i+1,0,0)-T6(i,0,0))
dtdz7(i)=dt/dz
end do

delete(dz)
dtdz4=new(dimsizes(T3(:,0,0)),double)
do i=0,dimsizes(T3(:,0,0))-2
meanT=(T3(i,0,0)+T3(i+1,0,0))/2
dz=((287*meanT/9.81)*log(plev(i)/plev(i+1)))/1000
dt=(T3(i+1,0,0)-T3(i,0,0))
dtdz4(i)=dt/dz
end do


dtdz5=new(dimsizes(T4(:,0,0)),double)
do i=0,dimsizes(T4(:,0,0))-2
dz=(Zcam4(i+1,0,0)-Zcam4(i,0,0))/1000
dt=(T4(i+1,0,0)-T4(i,0,0))
dtdz5(i)=dt/dz
end do

dtdz6=new(dimsizes(T5(:,0,0)),double)
do i=0,dimsizes(T5(:,0,0))-2
dz=(Zcamif(i+1,0,0)-Zcamif(i,0,0))/1000
dt=(T5(i+1,0,0)-T5(i,0,0))
dtdz6(i)=dt/dz
end do

dtdz8=new(dimsizes(T7(:,0,0)),double)
do i=0,dimsizes(T7(:,0,0))-2
dz=(Zcamifr(i+1,0,0)-Zcamifr(i,0,0))/1000
dt=(T7(i+1,0,0)-T7(i,0,0))
dtdz8(i)=dt/dz
end do



rmse=decimalPlaces(sqrt(avg((dTdZ2(1:43)-dTdZ(1:43))^2)),3,True)
rmse1=decimalPlaces(sqrt(avg((dTdZ1(1:43)-dTdZ(1:43))^2)),3,True)
rmse2=decimalPlaces(sqrt(avg((dTdZ3(1:43)-dTdZ(1:43))^2)),3,True)

;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("pdf","dtdz_koror")                  ; send graphics to PNG file

 res                   = True                       ; plot mods desired


 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainString    =  ""
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 4.0
 res@tiXAxisFontHeightF=0.018
 res@trYMaxF=1000
 res@trYMinF=100
 res@trXMaxF=0
 res@trXMinF=-9
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.018
 res@xyYStyle = "Log"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/100,200,300,400,500,600,700,850,1000/)
 res@tiXAxisString="dT/dz (K/km)"

 ;res@gsnRightString  ="CAM4-ZM RMSE: " + rmse + " K/km" + "~C~AMIP mean RMSE: " + rmse1 + " K/km" + "~C~CAM4-IF-t RMSE: " + rmse2 + " K/km"
 res@gsnLeftString  ="b) KOR"  + ": (" + lon + "E, " + lat + "N)"
 res@gsnLeftStringFontHeightF=0.018
  res@gsnRightStringFontHeightF=0.008
 res@gsnDraw=False
 res@gsnFrame=False
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 8.0
 plot               = gsn_csm_xy (wks,dTdZ(0:45),PmeanANN(0:45),res) 
 res@xyLineColor    = "orange"
 res@xyLineThicknesses = 4.0
 plot1              = gsn_csm_xy (wks,dtdz1(0:15),T&plev(0:15)*0.01,res)
 res@xyLineColor    = "blue"
 res@xyLineThicknesses = 4.0
 plot2= gsn_csm_xy (wks,dtdz2(0:24),T1&lev(1:25),res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 4.0
 plot3= gsn_csm_xy (wks,dtdz3(0:24),T2&lev(1:25),res)
 
 overlay(plot,plot1)
 overlay(plot,plot2)
 overlay(plot,plot3)
 
gres = True
gres@YPosPercent = 86.    ; expressed as %, 0->100, sets position of top border of legend 
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
gres@XPosPercent = 57    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
gres@ItemSpacePercent = 10.  
lineres = True
lineres@lgLineColors = (/"orange","red","blue","black"/) ; line colors
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres = True
textres@lgLabels = (/"CMIP5/AMIP mean","CAM4-IF-t","CAM4 default","RAOBS"/)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015,0.015,0.015/)   
plot = simple_legend(wks,plot,gres,lineres,textres)

txid = new(1,graphic)                  ; Nine text strings
txres               = True     
txres@txFontHeightF = 0.014
txid(0) = gsn_add_text (wks,plot,"CAM4-ZM RMSE: " + sprintf("%5.3f",rmse) + " K/km" + "~C~CAM4-IF-t RMSE: " + rmse2  + " K/km" + "~C~CMIP5/AMIP RMSE: " + sprintf("%5.3f",rmse1) +  " K/km", -2.2, 300 ,txres)
draw(plot)               
frame(wks) 
 
 
delete(txid)

rmse3=decimalPlaces(sqrt(avg((dTdZ4(1:43)-dTdZ7(1:43))^2)),3,True)
rmse4=decimalPlaces(sqrt(avg((dTdZ5(1:41)-dTdZ7(1:41))^2)),3,True)
rmse5=decimalPlaces(sqrt(avg((dTdZ6(1:43)-dTdZ7(1:43))^2)),3,True)
rmse6=decimalPlaces(sqrt(avg((dTdZ8(1:43)-dTdZ7(1:43))^2)),3,True)


wks   = gsn_open_wks ("pdf","dtdz_tropmean")                  ; send graphics to PNG file

res1                   = True                       ; plot mods desired


 res1@xyLineColor     = "blue"
 res1@trYReverse        = True                       ; reverse Y-axis
 res1@gsnLeftString    = "c) 20S-20N mean lapse rate"
 res1@tiMainFontHeightF=0.022
 res1@xyLineThicknesses = 4.0
 res1@tiXAxisString=""
 res1@tiXAxisFontHeightF=0.015
 res1@trYMaxF=1000
 res1@trYMinF=100
 res1@trXMaxF=0
 res1@trXMinF=-9
 res1@tiYAxisString="Pressure (hPa)"
 res1@xyYStyle = "Linear"
 res1@tmYLMode = "Explicit"
 res1@tmYLValues        =(/100,200,300,400,500,600,700,850,1000/)
 res1@tmYLLabels=(/100,200,300,400,500,600,700,850,1000/)
 res1@tiXAxisString="dT/dz (K/km)"
 res1@tiYAxisFontHeightF=0.018
 res1@tiXAxisFontHeightF=0.018
 res1@gsnDraw=False
 res1@gsnFrame=False
 res1@gsnLeftStringFontHeightF=0.018
 res1@xyLineColor    = "black"
 res1@xyLineThicknesses = 8.0
 plot4= gsn_csm_xy (wks,dtdz7(0:22),T6&lev(0:22),res1)
 res1@xyLineThicknesses = 4.0
 res1@xyLineColor    = "orange"
 plot1              = gsn_csm_xy (wks,dtdz4(0:15),T3&plev(0:15)*0.01,res1)
 res1@xyLineColor    = "blue"
 plot2= gsn_csm_xy (wks,dtdz5(0:24),T4&lev(1:25),res1)
 res1@xyLineColor    = "red"
 plot3= gsn_csm_xy (wks,dtdz6(0:24),T5&lev(1:25),res1)
 res1@xyLineColor    = "green"
 plot5= gsn_csm_xy (wks,dtdz8(0:24),T7&lev(1:25),res1)

 
 
 overlay(plot1,plot4)
 overlay(plot1,plot3)
 overlay(plot1,plot2)
 overlay(plot1,plot5)
 gres = True
gres@YPosPercent = 86.    ; expressed as %, 0->100, sets position of top border of legend 
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
gres@XPosPercent = 62    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
gres@ItemSpacePercent = 10.  
delete(lineres)
lineres = True
lineres@lgLineColors = (/"orange","red","green","blue","black"/) ; line colors
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
delete(textres)
textres = True
textres@lgLabels = (/"CMIP5/AMIP mean","CAM4-IF-t","CAM4-IF-r","CAM4 default","MERRA/ERAI"/)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.012,0.012,0.012,0.012/)   
plot1 = simple_legend(wks,plot1,gres,lineres,textres)
txid = new(1,graphic)                  ; Nine text strings
txres               = True     
txres@txFontHeightF = 0.014
txid(0) = gsn_add_text (wks,plot1,"CAM4 default RMSE: " + rmse4 + " K/km" + "~C~CAM4-IF-r RMSE: " + sprintf("%5.3f",rmse6) + " K/km" + "~C~CAM4-IF-t RMSE: " + sprintf("%5.3f",rmse5) + " K/km" + "~C~CMIP5/AMIP RMSE: " + rmse3 +  " K/km", -2.2, 300 ,txres)
draw(plot1)               
frame(wks) 
crop=systemfunc("pdfcrop dtdz_koror.pdf")
crop=systemfunc("pdfcrop dtdz_tropmean.pdf")
cleanup=systemfunc("rm tmp*.nc")
end
