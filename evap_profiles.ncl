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
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/evap/*.nc1")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/evap/*_4K.nc")
;fils2=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/evap/CAM4*ANN*.nc*")
print(fils)
print(fils1)
plotfile="test" ;;;PDF file name!
datalog=new((/3,26/),"double")
datalog1=new((/3,26/),"double")
datalog2=new((/3,26/),"double")
datalog3=new((/3,26/),"double")
datalog4=new((/3,26/),"double")
datalog5=new((/3,26/),"double")
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods*2, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
lets=(/"a)","b)","c)"/)
names=(/"CAM4-IF-r","CAM4-IF-t"/)
do while(n.le.dimsizes(fils)-1)
  
   if (n.eq.0) then
    filen="../CAM/files/CAM4_ANN_climo.nc4"
    filen1="../CAM/files/CAM4_4K_ANN_climo.nc2"
	cmd="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,DCQ " + filen + " tmp.nc"
    cmd1="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,DCQ " + filen1 + " tmp1.nc"
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
   else
   
    filen=tostring(fils(n))
    filen1=tostring(fils1(n))
   cmd="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,UPRAINDOWN,ANRAINDOWN,ANSNOWSUBL,RLP_RVV " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,UPRAINDOWN,ANRAINDOWN,ANSNOWSUBL,RLP_RVV " + filen1 + " tmp1.nc"
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
   end if
  
  if (n.eq.0) then
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f2->PS(0,:,:)                         ; Pa  [punits=1]
   tot=f->DCQ(0,:,:,:)
   tot&lev@units="hPa"
   
   tot1=f1->DCQ(0,:,:,:)
   tot1&lev@units="hPa"
  else
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f2->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->UPRAINDOWN(0,:,:,:)
   T&lev@units="hPa"
   U=f->ANRAINDOWN(0,:,:,:)
   U&lev@units="hPa"
   V=f->ANSNOWSUBL(0,:,:,:)
   V&lev@units="hPa"
   X=f->RLP_RVV(0,:,:,:)
   X&lev@units="hPa"
   tot=T+U+V+X
   tot!0="lev"
   tot&lev=T&lev
   tot&lev@units="hPa"
   
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   levs=conform_dims(dimsizes(T),T&lev,0)
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, T, PS/100, 0, opt)/900
   wva!0="lev"
   wva!1="lon"
   wva&lon=T&lon
   opt=(/0,100,1000/)
   UPRAINDOWNavg=avg(wva)
 
   
   PS1   = f2->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->UPRAINDOWN(0,:,:,:)
   T1&lev@units="hPa"
   U1=f1->ANRAINDOWN(0,:,:,:)
   U1&lev@units="hPa"
   V1=f1->ANSNOWSUBL(0,:,:,:)
   V1&lev@units="hPa"
   X1=f1->RLP_RVV(0,:,:,:)
   X1&lev@units="hPa"
   tot1=T1+U1+V1+X1
   tot1!0="lev"
   tot1&lev=T&lev
   tot1&lev@units="hPa"
   end if

 res                   = True                       ; plot mods desired
 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 1.0
 res@tiXAxisString=""
 res@tiXAxisFontHeightF=0.02
 res@trYMaxF=1000
 res@trYMinF=100
 res@trXMinF=min(tot1(:,0,0))
 ;res@trXMinF=-0.009
 ;res@trXMaxF=0.002
 res@trXMaxF=max(tot1(:,0,0))+max(tot1(:,0,0))*0.1
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.02
 res@tmXBFormat="0*+^e"
 res@xyYStyle = "Linear"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/100,200,300,400,500,600,700,850,1000/)
 res@tmXBLabelFontHeightF=0.015
 
 res@tiXAxisString="kg m~S~-2~N~ s~S~-1~N~"
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
 res@xyLineThicknesses = 3.0
 res@xyDashPattern  = 0
 if (n.eq.0) then
  res@xyLineColor    = "red"
   plot(n)= gsn_csm_xy (wks,tot(:,0,0),tot&lev,res)
  res@xyLineThicknesses = 1.0
  plot1= gsn_csm_xy (wks,tot1(:,0,0),tot1&lev,res)
   overlay(plot(n),plot1)
 else
 T4=T1(:,0,0)-T(:,0,0)
 T5=U1(:,0,0)-U(:,0,0)
 T6=V1(:,0,0)-V(:,0,0)
 T7=X1(:,0,0)-X(:,0,0)
 plot(n)               = gsn_csm_xy (wks,T(:,0,0),T1&lev,res) 
 res@xyLineColor    = "blue"
 plot4              = gsn_csm_xy (wks,U(:,0,0),T1&lev,res)
 res@xyLineColor    = "orange"
 plot5              = gsn_csm_xy (wks,V(:,0,0),T1&lev,res)
 res@xyLineColor    = "cyan"
 plot6              = gsn_csm_xy (wks,X(:,0,0),T1&lev,res)
 res@xyLineColor    = "red"
 plot10              = gsn_csm_xy (wks,tot(:,0,0),T1&lev,res)
 res@xyLineThicknesses = 1.0
 res@xyLineColor    = "black"
 plot1              = gsn_csm_xy (wks,T1(:,0,0),T1&lev,res)
 res@xyLineColor    = "blue"
 plot7              = gsn_csm_xy (wks,U1(:,0,0),T1&lev,res)
 res@xyLineColor    = "orange"
 plot8              = gsn_csm_xy (wks,V1(:,0,0),T1&lev,res)
 res@xyLineColor    = "cyan"
 plot9              = gsn_csm_xy (wks,X1(:,0,0),T1&lev,res)
 res@xyLineColor    = "red"
 plot11              = gsn_csm_xy (wks,tot1(:,0,0),T1&lev,res)
 
 overlay(plot(n),plot1)
 ;overlay(plot(n),plot3)
 overlay(plot(n),plot4)
 overlay(plot(n),plot5)
 overlay(plot(n),plot6)
 overlay(plot(n),plot7)
 overlay(plot(n),plot8)
 overlay(plot(n),plot9)
 overlay(plot(n),plot10)
 overlay(plot(n),plot11)
 ;overlay(plot(n),plot12)
 
 gres = True
gres@YPosPercent = 25    ; expressed as %, 0->100, sets position of top border of legend 
gres@XPosPercent = 40
gres@ItemSpacePercent = 5.  
lineres = True
lineres@lgLineColors = (/"red","black","blue","orange","cyan"/) ; line colors
lineres@lgDashIndexes=(/0,0,0,0,0/)
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres = True
textres@lgLabels = (/"NET","UPRAINDOWN","ANRAINDOWN","ANSNOWSUBL","RLP_RVV"/)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015,0.015,0.015/)   
plot(n) = simple_legend(wks,plot(n),gres,lineres,textres)
 end if
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