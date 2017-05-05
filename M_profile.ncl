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
;fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/FMASS_CAM*ANN.nc1")
;fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/FMASS_CAM*4K*ANN.nc2")
print(fils)
print(fils1)
plotfile="M_profile" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods*2, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
lets=(/"a)","b)","c)"/)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   cmd="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS " + filen1 + " tmp1.nc"
   cmd2="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,PS ../mass_flux/CAM-IFv5_F2000_ANN_climo.nc tmp2.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   areamean2=systemfunc(cmd2)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
   f2=addfile("tmp2.nc","r")
   f3=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/CAM-IFv5_F2000_ANN_climo.nc","r")
   hyam = f3->hyam 
   hybm = f3->hybm
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f2->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->FMASS(0,:,:,:)
   T&ilev@units="hPa"
   ;T2=f->OMEGA(0,:,:,:)
  ; T2&lev@units="hPa"
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   levs=conform_dims(dimsizes(T),T&ilev,0)
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
   
   PS1   = f2->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->FMASS(0,:,:,:)
   T1&ilev@units="hPa"
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
 res@trYMinF=100
 res@trXMinF=min(T1(:,0,0))-0.002
 ;res@trXMinF=-0.009
 res@trXMaxF=0.002
 ;res@trXMaxF=max(T1(:,0,0))+0.0001
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.02
 res@xyYStyle = "Linear"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/100,200,300,400,500,600,700,850,1000/)
 ;res@tmXBMode = "Explicit"
 ;res@tmXBValues        =(/-0.01,-0.008,-0.006,-0.004,-0.002,0,0.002/)
 ;res@tmXBLabels=(/-0.01,-0.008,-0.006,-0.004,-0.002,0,0.002/)
 res@tmXBLabelFontHeightF=0.015
 
 ;res@tiXAxisString="Total M~B~c~N~ (kg m~S~-2~N~ s~S~-1~N~)"
 res@tiXAxisString="Updraft M~B~c~N~ (kg m~S~-2~N~ s~S~-1~N~)"
 res@vpHeightF 	 = 0.7               ; change aspect ratio of plot
 res@vpWidthF 	         = 0.55
 res@gsnXRefLine = 0.0 

 res@tiXAxisStringFontHeightF=0.0015
 res@gsnCenterString  =lets(n) + " " + f@case
 res@gsnLeftString  =""
 res@gsnCenterStringFontHeightF=0.03
 res@gsnRightStringFontHeightF=0.025
 res@gsnDraw=False
 res@gsnFrame=False
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 0
 T4=T1(:,0,0)-T(:,0,0)
 plot(n)               = gsn_csm_xy (wks,T(:,0,0),T1&ilev,res) 
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 plot1              = gsn_csm_xy (wks,T1(:,0,0),T1&ilev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
 plot3= gsn_csm_xy  (wks,T4,T1&ilev,res)
 overlay(plot(n),plot1)
 overlay(plot(n),plot3)
 
 gres = True
gres@YPosPercent = 13.    ; expressed as %, 0->100, sets position of top border of legend 
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
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
plot(n) = simple_legend(wks,plot(n),gres,lineres,textres)
 
n=n+1
end do      
resP                            = True
 
 resP@cnInfoLabelOn       = False
 resP@lbLabelBarOn = False
 resP@txString = ""
 resP@txFontHeightF    = 0.03
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 gsn_panel(wks,plot,(/1,3/),resP)  
 delete(wks) 
 print(fils)
 print(fils1)
;crop=systemfunc("mogrify -trim dtdz_koror.png")
;crop=systemfunc("mogrify -trim dtdz_tropmean.png")
crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
cleanup=systemfunc("rm tmp*.nc")
end
