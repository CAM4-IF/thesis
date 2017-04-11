begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*_ANN*.nc11")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*4K_ANN*.nc21")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="FMASS_ANN" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods*2, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
lets=(/"a)","b)","c)","d)","e)","f)"/)
i=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile(filen, "r")
   hyam = h->hyam 
   hybm = h->hybm
   cmd="cdo -mermean -sellonlatbox,0,360,-30,30 -selname,FMASS,PS,OMEGA " + filen + " tmp.nc"
   cmd1="cdo -mermean -sellonlatbox,0,360,-30,30 -selname,FMASS,PS,OMEGA " + filen1 + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->FMASS(0,:,:,:)
   T&ilev@units="hPa"
   T2=f->OMEGA(0,:,:,:)
   T2&lev@units="hPa"
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   ;print(pm(4,:,:))
   ;print(dimsizes(T))
   levs=conform_dims(dimsizes(T),T&ilev,0)
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, T, PS/100, 0, opt)/900
   wva!0="lev"
   wva!1="lon"
   wva&lon=T&lon
   opt=(/0,100,1000/)
   wva2 = wgt_vert_avg_beta(pm/100, T2, PS/100, 1, opt)/900
   wva2!0="lev"
   wva2!1="lon"
   wva2&lon=T&lon
   fmassavg=avg(wva)
   omegavg=avg(wva2)
   
   PS1   = f1->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->FMASS(0,:,:,:)
   T1&ilev@units="hPa"
   T3=f1->OMEGA(0,:,:,:)
   T3&lev@units="hPa"
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt1=(/0,100,1000/)
   wva1 = wgt_vert_avg_beta(levs, T1, PS1/100, 0, opt1)/900
   wva1!0="lev"
   wva1!1="lon"
   wva1&lon=T&lon
   opt1=(/0,100,1000/)
   wva3 = wgt_vert_avg_beta(pm1/100, T3, PS/100, 0, opt)/900
   wva3!0="lev"
   wva3!1="lon"
   wva3&lon=T&lon
   fmassavg1=avg(wva1)
   omegavg1=avg(wva3)
   response=(fmassavg1-fmassavg)/fmassavg
   response1=(omegavg1-omegavg)/omegavg
   print(response)
   print(response1)
   
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
 
  Qdiff=T1(:,0,:)-T(:,0,:)
  Qdiff!0="lev"
  Qdiff&lev=T&ilev
  Qdiff!1="lon"
  Qdiff&lon=T&lon
  Qdiff&lev@units="hPa"
  
  Qdiff1=T3(:,0,:)-T2(:,0,:)
  Qdiff1!0="lev"
  Qdiff1&lev=T2&lev
  Qdiff1!1="lon"
  Qdiff1&lon=T&lon
  Qdiff1&lev@units="hPa"
  
;;;;calculate change in omega_up;;;;
  omegup=avg(wva2 < 0)
  omegup1=avg(wva3<0)
  omegresponse=(omegup1-omegup)/omegup
  
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels =(/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)*2e-4
  res@cnFillColors =(/10,15,20,25,30,35,40,45,50,55,60,65,0,0,180,185,190,195,200,205,210,215,220,225,230,235/)
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@lbLabelBarOn = False
  if(i.eq.2) then
   res@lbLabelBarOn = True
   res@lbOrientation =   "horizontal"  
   res@lbLabelAutoStride    = True                 ; optimal labels
   res@lbTitleString    = "kg/m^2/s"
   res@lbTitlePosition ="Bottom"
   res@pmLabelBarOrthogonalPosF=1.15
   res@lbTitleFontHeightF=.015
  else
  end if
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                ; turn off contour info label
  res@tiXAxisString    = ""
  res@tiYAxisSide = "Left"	
  res@tiYAxisString    = "Pressure (hPa)"
  res@trYMaxF = 1000.0
  res@trYMinF = 75.0      ; troposphere only
  res@tiXAxisFontHeightF =0.012
  res@tmYROn = False
  res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7

   
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels =(/-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)*1e-3
 cres@cnLineLabelsOn       = True                 ; turn on line labels
 cres@cnFillOn             = False                 ; turn on color fill
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False                ; turn off contour info label
 cres@gsnRightString = ""
 cres@gsnLeftString = ""
 cres@tiYAxisString    = "Pressure (hPa)"
 cres@tiXAxisString    = "Longitude"
 cres@tmYROn = False
 cres@tmYRLabelsOn = False
 cres@gsnContourZeroLineThicknessF  = 0.0
 cres@cnLineThicknessF =0.5
 ;cres@gsnContourNegLineDashPattern =-1  
 cres@gsnLeftString =""
 cres@gsnContourZeroLineThicknessF = 2
 cres@gsnRightString    = lets(i) + " " + f@case
 cres@trYMinF            = 75.   ; Limits for Y axis.  The limits
 cres@trYMaxF            = 1000.
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelFontHeightF =0.01
 cres@cnLineLabelPlacementMode = "Constant"
 cres@cnLineLabelInterval=1
 cres@gsnLeftString  = "response (%): " + response*100
 cres@gsnCenterString  = "FMASS"
 cres@gsnLeftStringFontHeightF=0.012
 cres@gsnRightStringFontHeightF=0.015
 cres@gsnContourPosLineDashPattern =1
gsn_define_colormap(wks,"BlueWhiteOrangeRed")
 plot(i)  = gsn_csm_pres_hgt(wks, Qdiff(::-1,:), res)  
 plot2  = gsn_csm_pres_hgt(wks, T(::-1,0,:), cres)
 
 
 cres@cnLevels =(/-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)*4e-3
 res@cnLevels =(/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)*1e-3
 res@cnFillColors = (/10,15,20,25,30,35,40,45,50,55,60,65,0,0,180,185,190,195,200,205,210,215,220,225,230,235/)
 cres@gsnContourPosLineDashPattern =1 
 cres@gsnLeftString  = "response (upward only) (%): ~C~" + omegresponse*100
 cres@gsnRightString  = lets(i+1) + " " + f@case
 cres@gsnCenterString  = "OMEGA"
 if(i.eq.2) then
  res@lbTitleString    = "Pa/s"
 else
 end if
 plot(i+1)  = gsn_csm_pres_hgt(wks, Qdiff1(::-1,:), res)  
 plot3  = gsn_csm_pres_hgt(wks, T2(::-1,0,:), cres)

 
 overlay(plot(i),plot2)
 overlay(plot(i+1),plot3)
   
   cleanup=systemfunc("rm tmp*.nc")
   delete(T)
   delete(T1)
   delete(Qdiff)
   delete(T2)
   delete(T3)
   delete(Qdiff1)
  ; delete(pm)
  ; delete(pm1)
  print(n)
n=n+1
i=i+2
end do
T=f->FMASS(0,:,:,:)
 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbLabelBarOn = False
 resP@txString = "ANN 30S-30N mean FMASS response (shaded) and climo (contours)"
 resP@txFontHeightF    = 0.01
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 gsn_panel(wks,plot,(/3,2/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop " + plotfile)  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(fils1)
end