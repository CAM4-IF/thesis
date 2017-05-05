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
  res@cnLevels =(/-14,-12,-10,-8,-6,-4,-2,2,4,6,8,10,12,14/)*1e-4
  res@cnFillColors =(/40,50,60,70,80,90,100,0,140,150,160,170,180,190,200/)
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@lbLabelBarOn = False
  if(i.eq.2) then
   res@lbLabelBarOn = True
   res@lbOrientation =   "horizontal"  
   ;res@lbLabelAutoStride    = True                ; optimal labels
   res@lbTitleString    = "kg m~S~-2~N~ s~S~-1~N~"
   res@lbTitlePosition ="Bottom"
   res@pmLabelBarOrthogonalPosF=1.22
   res@lbTitleFontHeightF=.015
   res@lbBoxEndCapStyle="TriangleBothEnds"
   res@lbAutoManage=True
   ;res@lbLabelStrings =  (/"-2.4e-3","-2e-3","-1.6e-3","-1.2e-3","-8e-4","-4e-4","4e-4","8e-4","1.2e-3","1.6e-3","2e-3","2.4e-3"/)
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
  res@cnFillOpacityF=1
  
res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
  res@vpWidthF 	         = 0.7


   
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
 cres@tiYAxisString    = "Pressure (hPa)"
 cres@tiXAxisString    = "Longitude"
 cres@tmYROn = False
 cres@tmYRLabelsOn = False
 cres@gsnContourZeroLineThicknessF  = 0.0
 cres@cnLineThicknessF =0.5
 ;cres@gsnContourNegLineDashPattern =-1  
 cres@gsnLeftString =""
 cres@gsnContourZeroLineThicknessF = 2
 cres@gsnLeftString    = lets(i) + " " + f@case
 cres@trYMinF            = 75.   ; Limits for Y axis.  The limits
 cres@trYMaxF            = 1000.
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelFontHeightF =0.01
 cres@cnLineLabelPlacementMode = "Constant"
 cres@cnLineLabelInterval=1
 cres@gsnRightString  = "response (%): " + decimalPlaces(response*100,2,True)
 cres@gsnCenterString  = "M~B~c~N~"
 cres@gsnLeftStringFontHeightF=0.03
 cres@gsnRightStringFontHeightF=0.03
 cres@gsnCenterStringFontHeightF=0.035
 cres@gsnContourPosLineDashPattern =1
gsn_define_colormap(wks,"ncl_default")
 plot(i)  = gsn_csm_pres_hgt(wks, Qdiff(::-1,:), res)  
 plot2  = gsn_csm_pres_hgt(wks, T(::-1,0,:), cres)
 
 delete(cres@cnLevels)
 delete(res@cnLevels)
 cres@cnLevels =(/-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)*4e-3
 res@cnLevels =(/-14,-12,-10,-8,-6,-4,-2,2,4,6,8,10,12,14/)*0.5e-3
 res@cnFillColors = (/40,50,60,70,80,90,100,0,140,150,160,170,180,190,200/)
 cres@gsnContourPosLineDashPattern =1 
 ;cres@gsnRightString  = "~F33~w~F~~S~~F34~-~F~~N~ response (%): " + decimalPlaces(omegresponse*100,2,True)
 cres@gsnRightString  = ""
 cres@gsnLeftString  = lets(i+1) + " " + f@case
 cres@gsnCenterString  = "~F33~w~F~"
 if(i.eq.2) then
  res@lbTitleString    = "Pa s~S~-1~N~"
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
 resP@txString = ""
 resP@txFontHeightF    = 0.03
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 gsn_panel(wks,plot,(/3,2/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop " + plotfile)  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(fils1)
end
