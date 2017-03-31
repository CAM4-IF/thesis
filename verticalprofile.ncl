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
f2  = addfile("testforian_ANN.nc","r")
f7=addfile("RAOBS.nc","r")
cmd="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,ta ../AMIP/ta_ANN_AMIPmean.nc" + " tmp.nc"
cmd1="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,T,PS ../AMIP_4K_control_ANN_climo1.nc" + " tmp1.nc"
cmd2="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,T,PS testforian_ANN.nc" + " tmp2.nc"
cmd3="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,T,PS ERAI_ANN_climo.nc" + " tmp3.nc"

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
ps5=f6->PS(0,:,:)

test=f7->T(:,:,:)
print(test)

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



T3@_FillValue=1e20
u3=int2p_n_Wrap(plev*0.01,T3,PmeanANN,0,0)
u4=vinth2p(T4,hyam,hybm,PmeanANN,ps3,2,1000,1,False)
u5=vinth2p(T5,hyam,hybm,PmeanANN,ps4,2,1000,1,False)
;u6=vinth2p(T6,hyam,hybm,PmeanANN,ps5,2,1000,1,False)
elevs=f6->lev



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


dTdZ7=new(dimsizes(elevs),double)
do i=0,dimsizes(elevs)-2
meanT=(T6(i,0,0)+T6(i+1,0,0))/2
dZ=((287*meanT/9.81)*log(elevs(i)/elevs(i+1)))/1000
dT=(T6(i+1,0,0)-T6(i,0,0))
dTdZ7(i)=dT/dZ
end do




;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("png","dtdz_koror")                  ; send graphics to PNG file

 res                   = True                       ; plot mods desired


 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainString    = "KOR"  + ": (" + lon + "E, " + lat + "N) ANN mean"
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 4.0
 res@tiXAxisString=""
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
 res@tiXAxisString="dT/dz: CAM4-IF: red, CAM4: Blue, CMIP5 AMIP mean: orange, OBS: black"

 res@tiXAxisStringFontHeightF=0.008
 res@gsnDraw=False
 res@gsnFrame=False
 res@xyLineColor    = "black"
 plot               = gsn_csm_xy (wks,dTdZ(0:45),PmeanANN(0:45),res) 
 res@xyLineColor    = "orange"
 plot1              = gsn_csm_xy (wks,dTdZ1(0:45),PmeanANN(0:45),res)
 res@xyLineColor    = "blue"
 plot2= gsn_csm_xy (wks,dTdZ2(0:45),PmeanANN(0:45),res)
 res@xyLineColor    = "red"
 plot3= gsn_csm_xy (wks,dTdZ3(0:45),PmeanANN(0:45),res)
 
 
 overlay(plot,plot1)
 overlay(plot,plot2)
 overlay(plot,plot3)
 
draw(plot)               
frame(wks) 
 
 
wks   = gsn_open_wks ("png","dtdz_tropmean")                  ; send graphics to PNG file

res1                   = True                       ; plot mods desired


 res1@xyLineColor     = "blue"
 res1@trYReverse        = True                       ; reverse Y-axis
 res1@tiMainString    = "10S-10N tropical mean lapse rate"
 res1@tiMainFontHeightF=0.018
 res1@xyLineThicknesses = 4.0
 res1@tiXAxisString=""
 res1@tiXAxisFontHeightF=0.018
 res1@trYMaxF=1000
 res1@trYMinF=100
 res1@trXMaxF=0
 res1@trXMinF=-9
 res1@tiYAxisString="Pressure (hPa)"
 res1@tiYAxisFontHeightF=0.018
 res1@xyYStyle = "Log"
 res1@tmYLMode = "Explicit"
 res1@tmYLValues        =(/100,200,300,400,500,600,700,850,1000/)
 res1@tmYLLabels=(/100,200,300,400,500,600,700,850,1000/)
 res1@tiXAxisString="dT/dz: CAM4-IF: red, CAM4: Blue, CMIP5 AMIP mean: orange, ERAI: black"

 res1@tiXAxisStringFontHeightF=0.008
 res1@gsnDraw=False
 res1@gsnFrame=False
 res1@xyLineColor    = "black"
 ;plot               = gsn_csm_xy (wks,dTdZ(0:45),PmeanANN(0:45),res1) 
 res1@xyLineColor    = "orange"
 plot1              = gsn_csm_xy (wks,dTdZ4(0:45),PmeanANN(0:45),res1)
 res1@xyLineColor    = "blue"
 plot2= gsn_csm_xy (wks,dTdZ5(0:45),PmeanANN(0:45),res1)
 res1@xyLineColor    = "red"
 plot3= gsn_csm_xy (wks,dTdZ6(0:45),PmeanANN(0:45),res1)
 res1@xyLineColor    = "black"
 plot4= gsn_csm_xy (wks,dTdZ7,elevs,res1)
 
 
 overlay(plot1,plot2)
 overlay(plot1,plot3)
 overlay(plot1,plot4)
 
draw(plot1)               
frame(wks) 


;crop=systemfunc("/usr/bin/pdfcrop dtdz_koror.pdf")

;crop=systemfunc("/usr/bin/pdfcrop dtdz_tropmean.png")
cleanup=systemfunc("rm tmp*.nc")
end
