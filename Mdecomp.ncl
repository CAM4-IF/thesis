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
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*FMASS*_CAM*ANN*.nc11")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*FMASS*_CAM*4K_ANN*.nc21")
fils2=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/FMASS_DN*CAM*ANN.nc")
fils3=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/FMASS_DN*CAM*4K*ANN.nc2")

plotfile="M_profile_updown" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods*2, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=1
i=0
lets=(/"c)","a)","b)","c)","d)"/)
names=(/"CAM4 default","CAM4-IF-r","CAM4-IF-t"/)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   filen2=tostring(fils2(i))
   filen3=tostring(fils3(i))
   print(filen2)
   print(filen3)
   cmd="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS " + filen1 + " tmp1.nc"
   cmd2="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,PS ../mass_flux/CAM-IFv5_F2000_ANN_climo.nc tmp2.nc"
   cmd3="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS_DN " + filen2 + " tmp3.nc"
   cmd4="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS_DN " + filen3 + " tmp4.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   areamean2=systemfunc(cmd2)
   areamean2=systemfunc(cmd3)
   areamean2=systemfunc(cmd4)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
   f2=addfile("tmp2.nc","r")
   f3=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/CAM-IFv5_F2000_ANN_climo.nc","r")
   f4=addfile("tmp3.nc","r")
   f5=addfile("tmp4.nc","r")
   hyam = f3->hyam 
   hybm = f3->hybm
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f2->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->FMASS(0,:,:,:)
   T&ilev@units="hPa"
   T2=f4->FMASS_DN(0,:,:,:)
   T2&ilev@units="hPa"
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   levs=conform_dims(dimsizes(T),T&ilev,0)
   opt=(/0,100,1000/)
   updraft=T-T2
   downdraft=T2
   wva = wgt_vert_avg_beta(levs, updraft, PS/100, 0, opt)/900
   wva!0="lev"
   wva!1="lon"
   wva&lon=T&lon
   fmassavg=avg(wva)
   
   PS1   = f2->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->FMASS(0,:,:,:)
   T1&ilev@units="hPa"
   T3=f5->FMASS_DN(0,:,:,:)
   T3&ilev@units="hPa"
   updraft1=T1-T3
   downdraft1=T3
   updraftdiff=updraft1-updraft
   downdraftdiff=downdraft1-downdraft
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt1=(/0,100,1000/)
   wva1 = wgt_vert_avg_beta(levs, T1, PS1/100, 0, opt1)/900
   wva1!0="lev"
   wva1!1="lon"
   wva1&lon=T&lon

 res                   = True                       ; plot mods desired
 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 1.0
 res@tiXAxisString=""
 res@tiXAxisFontHeightF=0.02
 res@trYMaxF=1000
 res@trYMinF=100
 res@trXMinF=min(updraft1(:,0,0))-0.002
 ;res@trXMinF=-0.009
 res@trXMaxF=0.002
 ;res@trXMaxF=max(T1(:,0,0))+0.0001
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.02
 res@xyYStyle = "Linear"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/100,200,300,400,500,600,700,850,1000/)
 res@tmXBMode = "Explicit"
 delete(res@tmXBValues)
 delete(res@tmXBLabels)
 res@tmXBValues        =(/-0.01,-0.008,-0.006,-0.004,-0.002,0,0.002/)
 res@tmXBLabels=(/-0.01,-0.008,-0.006,-0.004,-0.002,0,0.002/)
 res@tmXBLabelFontHeightF=0.015
 
 res@tiXAxisString="M~B~u~N~ (kg m~S~-2~N~ s~S~-1~N~)"
 res@vpHeightF 	 = 0.7               ; change aspect ratio of plot
 res@vpWidthF 	         = 0.5
 res@gsnXRefLine = 0.0 

 res@tiXAxisStringFontHeightF=0.0015
 res@gsnCenterString  =lets(n) + " " + names(n)
 res@gsnLeftString  =""
 res@gsnCenterStringFontHeightF=0.03
 res@gsnRightStringFontHeightF=0.025
 res@gsnDraw=False
 res@gsnFrame=False
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 0
 plot(i)               = gsn_csm_xy (wks,updraft(:,0,0),T1&ilev,res) 
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 plot1              = gsn_csm_xy (wks,updraft1(:,0,0),T1&ilev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
 plot3= gsn_csm_xy  (wks,updraftdiff(:,0,0),T1&ilev,res)
 overlay(plot(i),plot1)
 overlay(plot(i),plot3)
 
 gres = True
gres@YPosPercent = 13.    ; expressed as %, 0->100, sets position of top border of legend
;gres@XPosPercent = 10    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

gres@XPosPercent = 20
gres@ItemSpacePercent = 5.  
lineres = True
lineres@lgLineColors = (/"black","red","black"/) ; line colors
lineres@lgDashIndexes=(/1,0,0/)
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres = True
textres@lgLabels = (/"Response","4K warming","Control"/)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015/)   
plot(i) = simple_legend(wks,plot(i),gres,lineres,textres)

 res@trXMaxF=max(downdraft1(:,0,0))+0.0005
 res@trXMinF=-0.001
 res@tmXBMode = "Explicit"
 delete(res@tmXBValues)
 delete(res@tmXBLabels)
 res@tmXBValues        =(/-0.002,-0.001,0,0.001,0.002,0.003,0.004,0.005/)
 res@tmXBLabels=(/-0.002,-0.001,0,0.001,0.002,0.003,0.004,0.005/)
 res@gsnCenterString  =lets(n+2) + " " + names(n)
res@tiXAxisString="M~B~b~N~ (kg m~S~-2~N~ s~S~-1~N~)"
res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 0
plot(i+2)=gsn_csm_xy (wks,downdraft(:,0,0),T1&ilev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
plot2=gsn_csm_xy (wks,downdraft1(:,0,0),T1&ilev,res) 
res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
plot3=gsn_csm_xy (wks,downdraftdiff(:,0,0),T1&ilev,res) 
overlay(plot(i+2),plot2)
overlay(plot(i+2),plot3)

n=n+1
i=i+1
end do      
resP                            = True
 
 resP@cnInfoLabelOn       = False
 resP@lbLabelBarOn = False
 resP@txString = ""
 resP@txFontHeightF    = 0.03
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 gsn_panel(wks,plot,(/2,2/),resP)  
 delete(wks) 
 print(fils2)
 print(fils3)
;crop=systemfunc("mogrify -trim dtdz_koror.png")
;crop=systemfunc("mogrify -trim dtdz_tropmean.png")
crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
cleanup=systemfunc("rm tmp*.nc")
end
