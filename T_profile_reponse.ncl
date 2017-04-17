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
;************************************************FM
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc2")
print(fils)
print(fils1)
plotfile="T_profile_response" ;;;PDF file name!
datalog=new((/3,26/),"double")
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(2, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
lets=(/"a)","b)","c)","d)"/)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM4_4K_ANN_climo.nc2", "r")
   hyam = h->hyam 
   hybm = h->hybm
   cmd="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,T,PS " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,T,PS " + filen1 + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->T(0,:,:,:)
   T&lev@units="hPa"
   ;T2=f->OMEGA(0,:,:,:)
  ; T2&lev@units="hPa"
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   levs=conform_dims(dimsizes(T),T&lev,0)
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, T, PS/100, 0, opt)/900
   wva!0="lev"
   wva!1="lon"
   wva&lon=T&lon
   opt=(/0,100,1000/)
   ;wva2 = wgt_vert_avg_beta(pm/100, T2, PS/100, 1, opt)/900
   ;wva2!0="lev"
   ;wva2!1="lon"
   ;wva2&lon=T&lon
   fmassavg=avg(wva)
   ;omegavg=avg(wva2)
   
   PS1   = f1->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->T(0,:,:,:)
   T1&lev@units="hPa"
   ;T3=f1->OMEGA(0,:,:,:)
   ;T3&lev@units="hPa"
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt1=(/0,100,1000/)
   wva1 = wgt_vert_avg_beta(levs, T1, PS1/100, 0, opt1)/900
   wva1!0="lev"
   wva1!1="lon"
   wva1&lon=T&lon
   opt1=(/0,100,1000/)
   ;wva3 = wgt_vert_avg_beta(pm1/100, T3, PS/100, 0, opt)/900
   ;wva3!0="lev"
   ;wva3!1="lon"
   ;wva3&lon=T&lon
   fmassavg1=avg(wva1)
   ;omegavg1=avg(wva3)
   response=(fmassavg1-fmassavg)/fmassavg
   ;response1=(omegavg1-omegavg)/omegavg
   print(response)
   ;print(response1)

 res                   = True                       ; plot mods desired
 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 1.0
 res@tiXAxisString=""
 res@tiXAxisFontHeightF=0.02
 res@trYMaxF=1000
 res@trYMinF=30
 ;res@trXMinF=min(T1(:,0,0))-0.0001
 res@trXMinF=-1
 ;res@trXMaxF=0.002
 res@trXMaxF=9
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.02
 res@xyYStyle = "Linear"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmXBMode = "Explicit"
 res@tmXBValues        =(/-2,-1,0,1,2,3,4,5,6,7,8,9,10/)
 res@tmXBLabels=(/-2,-1,0,1,2,3,4,5,6,7,8,9,10/)
 res@tmXBLabelFontHeightF=0.015
 
 ;res@tiXAxisString="Total M~B~c~N~ (kg/m~S~2~N~/s)"
 res@tiXAxisString="~F33~D~F~T (K)"
 res@vpHeightF 	 = 0.7               ; change aspect ratio of plot
 res@vpWidthF 	         = 0.55
 res@gsnXRefLine = 0.0 

 res@tiXAxisStringFontHeightF=0.0015
 res@gsnCenterString  ="a) Temperature response"
 res@gsnLeftString  =""
 res@gsnCenterStringFontHeightF=0.03
 res@gsnRightStringFontHeightF=0.025
 res@gsnDraw=False
 res@gsnFrame=False
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
 T4=T1(:,0,0)-T(:,0,0)
 
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
 datalog(n,:)=T4
 

n=n+1
end do      

res@xyDashPattern  = 0
res@xyLineColor    = "blue"
 res@xyLineThicknesses = 1.0
 plot(0)= gsn_csm_xy  (wks,datalog(0,:),T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 plot2= gsn_csm_xy  (wks,datalog(1,:),T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 plot3= gsn_csm_xy  (wks,datalog(2,:),T1&lev,res)
 overlay(plot(0),plot2)
 overlay(plot(0),plot3)
 gres = True
gres@YPosPercent = 25. 
gres@XPosPercent = 20 
gres@ItemSpacePercent = 5.  
lineres = True
lineres@lgLineColors = (/"black","red","blue"/) ; line colors
lineres@lgDashIndexes=(/0,0,0/)
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres = True
textres@lgLabels = (/"CAM4-IF best T","CAM4-IF best rain","Default CAM4"/)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015/)   
plot3 = simple_legend(wks,plot3,gres,lineres,textres)
 delta=datalog(1,:)-datalog(0,:)
 delta=datalog(2,:)-datalog(0,:)
 delta1=datalog(2,:)-datalog(0,:)
 delta2=datalog(1,:)-datalog(2,:)
 
 
 
 res@trXMinF=-2
 res@trXMaxF=2
 res@xyDashPattern  = 0
 res@xyLineColor    = "blue"
 res@xyLineThicknesses = 1.0
 plot(1)= gsn_csm_xy  (wks,delta,T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 plot4= gsn_csm_xy  (wks,delta1,T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 plot5= gsn_csm_xy  (wks,delta2,T1&lev,res)
 overlay(plot(1),plot4)
 overlay(plot(1),plot5)
 
res4 = True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = ""
res4@lbTitlePosition ="Bottom"
res4@lbTitleFontHeightF=.008
res4@txFontHeightF    = 0.03
res4@gsnMaximize=True


gsn_panel(wks,plot,(/1,2/),res4)
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
cleanup=systemfunc("rm tmp*.nc")
end
