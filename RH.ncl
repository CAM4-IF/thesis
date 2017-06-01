begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*4K*ANN*.nc2")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="test_RH" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
lets=(/"a)","b)","c)"/)
names=(/"CAM4 default","CAM4-IF-r","CAM4-IF-t"/)
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   filen2="/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM-IFv4_F2000_4K_ANN_climo.nc1"
   g=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM-IFv4_F2000_4K_ANN_climo.nc1","r")
   hyam = g->hyam 
   hybm = g->hybm
   cmd="cdo -mermean -sellonlatbox,0,360,-30,30 -selname,RELHUM,PS " + filen + " tmp.nc"
   cmd1="cdo -mermean -sellonlatbox,0,360,-30,30 -selname,RELHUM,PS " + filen1 + " tmp1.nc"
   cmd2="cdo -mermean -sellonlatbox,0,360,-30,30 -selname,PS " + filen2 + " tmp2.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   areamean2=systemfunc(cmd2)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
   f2=addfile("tmp2.nc","r")
   T4=f->RELHUM(0,:,:,:)
   T4&lev@units="hPa"
  
   P0=100000.
   ;PS=new((/1,1/),float)
   ;PS(0,:)=1000.
   ;PS(:,0)=1000.
   ;PS1=new((/1,1/),float)
   ;PS1(0,:)=1000.
   ;PS1(:,0)=1000.
   PS=f->PS(0,:,:)
   PS1=f1->PS(0,:,:)
   
   T9=f1->RELHUM(0,:,:,:)
   
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa]
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa]    
   opt=(/0,850,1000/)
   print(pm)
   wva = wgt_vert_avg_beta(pm/100, T9, PS, 0, opt)/900

   
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
 
  Qdiff=(T9)-(T4)
  Qdiff!0="lev"
  Qdiff&lev=T9&lev
  Qdiff!2="lon"
  Qdiff&lon=T9&lon
  Qdiff&lev@units="hPa"
  wva1 = wgt_vert_avg_beta(pm1/100, Qdiff, PS1, 0, opt)/900
  response=100*avg(wva1)/avg(wva)
  
  res                      = True                 ; plot mods desired

  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15/)
  res@cnFillColors = (/12,18,24,40,60,90,106,116,128,130,142,148,156,168,180,208,218,230,250,255/)
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@lbLabelBarOn = False
  res@lbLabelAutoStride    = True                 ; optimal labels
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                ; turn off contour info label
  res@tiXAxisString    = ""
  res@tiYAxisSide = "Left"	
  res@tiYAxisString    = "Pressure (hPa)"
  res@trYMaxF = 1000.0
  res@trYMinF = 100.0      ; troposphere only
  res@tiXAxisFontHeightF =0.012
  res@tmYROn = False
  res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7
  
   ;Qzon=dim_rmvmean_n_Wrap(wva4(0,:),0)
   ;Qzon!0="lon"
   ;Qzon&lon=T&lon
   ;Qzon1=dim_rmvmean_n_Wrap(wva9(0,:),0)
   ;Qzon1!0="lon"
   ;Qzon1&lon=T&lon
   ;Qzonmean=avg(Qzon({120:180}))
   ;Qzonmean1=avg(Qzon1({120:180}))
   ;Qzonmean3=avg(wva4)
   ;Qzonmean4=avg(wva9)
   
   ;response=Qzonmean1-Qzonmean
   ;response1=Qzonmean4-Qzonmean3
   ;percentresponse=(response/(Qzonmean))*100.
   ;percentresponse1=(response1/(Qzonmean3))*100.
 ;  datalog(n)=percentresponse
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels = (/10,20,30,40,50,60,70,80,90/)
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
 cres@gsnContourNegLineDashPattern =1  
 cres@gsnLeftString =""
 cres@gsnContourZeroLineThicknessF = 2
 cres@gsnCenterString    =lets(n) + " " + names(n)
 cres@trYMinF            = 100.   ; Limits for Y axis.  The limits
 cres@trYMaxF            = 1000.
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelPlacementMode = "Constant"
 cres@cnLineLabelInterval=1
 cres@gsnLeftString  =""
 cres@gsnLeftStringFontHeightF=0.022
 cres@gsnRightStringFontHeightF=0.022
 cres@gsnCenterStringFontHeightF=0.03
 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 plot(n)  = gsn_csm_pres_hgt(wks, Qdiff(:,0,:), res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_pres_hgt(wks, T4(:,0,:), cres)
 overlay(plot(n),plot2)
 
   
   cleanup=systemfunc("rm tmp*.nc")
   delete(T4)
   delete(T9)
   delete(Qdiff)
   delete(pm)
   delete(pm1)
n=n+1
end do
 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = ""
 resP@lbTitleString    = "_~H-35~~V-29~~F33~d~F~~F33~H~F~/~F33~H~F~~H-20~~V+29~~N~_" + "~H+10~~V-29~~N~(%)"
 resP@lbTitleFontHeightF=.012
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= True
 resP@txFontHeightF    = 0.0085
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 resP@pmLabelBarOrthogonalPosF=-0.02
 gsn_panel(wks,plot,(/4,3/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop test_RH.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(fils1)
;print(datalog)
end
