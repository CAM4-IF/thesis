begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc1")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="all_radQ_JJA_panel" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   filen2="/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM-IFv4_F2000_4K_ANN_climo.nc1"
   g=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM-IFv4_F2000_4K_ANN_climo.nc1","r")
   hyam = g->hyam 
   hybm = g->hybm
   cmd="cdo -mermean -sellonlatbox,0,360,-15,15 -selname,DTCOND,QRL,QRS,DTV,PS " + filen + " tmp.nc"
   cmd1="cdo -mermean -sellonlatbox,0,360,-15,15 -selname,DTCOND,QRL,QRS,DTV,PS " + filen1 + " tmp1.nc"
   cmd2="cdo -mermean -sellonlatbox,0,360,-15,15 -selname,PS " + filen2 + " tmp2.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   areamean2=systemfunc(cmd2)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
   f2=addfile("tmp2.nc","r")
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)  
   PS1   = f1->PS(0,:,:)  

   ; Pa  [punits=1]
   T=f->DTCOND(0,:,:,:)
   T1=f->QRL(0,:,:,:)
   T2=f->QRS(0,:,:,:)
   T3=f->DTV(0,:,:,:)
   T4=T2(:,0,:)+T1(:,0,:)
   ;T2!0="lev"
   ;T2&lev=T&lev
   ;T2!1="lon"
   ;T2&lon=T&lon
   T1&lev@units="hPa"
   Ttot=T2+T1
   T4!0="lev"
   T4&lev=T&lev
   T4!1="lon"
   T4&lon=T&lon
   T4&lev@units="hPa"
   T4=dim_rmvmean_n_Wrap(T4,1)
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt=(/0,10000,20000/)
   diff=10000
   wva = wgt_vert_avg_beta(pm, T, PS, 1, opt)/diff ;vertically pressure thickness weighted sum
   wva1 = wgt_vert_avg_beta(pm, T1, PS, 1, opt)/diff
   wva2 = wgt_vert_avg_beta(pm, T2, PS, 1, opt)/diff
   wva3 = wgt_vert_avg_beta(pm, T3, PS, 1, opt)/diff
   wva4 = wgt_vert_avg_beta(pm, Ttot, PS, 1, opt)/diff
   
   wva4!0="lev"
   wva4!1="lon"
   wva4&lon=T&lon
   
   
   T5=f1->DTCOND(0,:,:,:)
   T6=f1->QRL(0,:,:,:)
   T7=f1->QRS(0,:,:,:)
   T8=f1->DTV(0,:,:,:)
   T9=T7(:,0,:)+T6(:,0,:)
   Ttot1=T7+T6
   T9!0="lev"
   T9&lev=T&lev
   T9!1="lon"
   T9&lon=T&lon
   T9&lev@units="hPa"
   T9=dim_rmvmean_n_Wrap(T9,1)
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt1=(/0,10000,20000/)
   wva5 = wgt_vert_avg_beta(pm1, T5, PS1, 1, opt1)/diff
   wva6 = wgt_vert_avg_beta(pm1, T6, PS1, 1, opt1)/diff
   wva7 = wgt_vert_avg_beta(pm1, T7, PS1, 1, opt1)/diff
   wva8 = wgt_vert_avg_beta(pm1, T8, PS1, 1, opt1)/diff
   wva9 = wgt_vert_avg_beta(pm1, Ttot1, PS1, 1, opt1)/diff
   wva9!0="lev"
   wva9!1="lon"
   wva9&lon=T&lon
   
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
 
  Qdiff=(T9)-(T4)
  Qdiff!0="lev"
  Qdiff&lev=T&lev
  Qdiff!1="lon"
  Qdiff&lon=T&lon
  Qdiff&lev@units="hPa"
  res                      = True                 ; plot mods desired

  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-30,-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20,30/)*2.0e-6
  res@cnFillColors = (/2,6,12,18,24,40,60,90,106,116,128,130,142,148,156,168,180,208,218,230,250,255/)
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
  res@trYMinF = 50.0      ; troposphere only
  res@tiXAxisFontHeightF =0.012
  res@tmYROn = False
  res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7
  
   Qzon=dim_rmvmean_n_Wrap(wva4(0,:),0)
   Qzon!0="lon"
   Qzon&lon=T&lon
   Qzon1=dim_rmvmean_n_Wrap(wva9(0,:),0)
   Qzon1!0="lon"
   Qzon1&lon=T&lon
   Qzonmean=avg(Qzon({120:180}))
   Qzonmean1=avg(Qzon1({120:180}))
   Qzonmean3=avg(wva4)
   Qzonmean4=avg(wva9)
   
   response=Qzonmean1-Qzonmean
   response1=Qzonmean4-Qzonmean3
   percentresponse=(response/(Qzonmean))*100.
   percentresponse1=(response1/(Qzonmean3))*100.
 ;  datalog(n)=percentresponse
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels = (/-260,-240,-220,-200,-180,160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240,260/)*2.0e-7
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
 cres@gsnRightString    = f@case
 cres@trYMinF            = 50.   ; Limits for Y axis.  The limits
 cres@trYMaxF            = 1000.
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelPlacementMode = "Constant"
 cres@cnLineLabelInterval=1
 cres@gsnLeftString  = "Climatological Q (100-200 hPa):" + Qzonmean3 + "~C~Q response:" + response1 + "~C~Climatological Q* asc. reg. (100-200 hPa): " + Qzonmean + "~C~Q* response:" + response
 cres@gsnLeftStringFontHeightF=0.012
 cres@gsnRightStringFontHeightF=0.015
 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 plot(n)  = gsn_csm_pres_hgt(wks, Qdiff(::-1,:), res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_pres_hgt(wks, T4(::-1,:), cres)
 overlay(plot(n),plot2)
 
   
   cleanup=systemfunc("rm tmp*.nc")
   delete(T)
   delete(T1)
   delete(T2)
   delete(T3)
   delete(T4)
   delete(T5)
   delete(T6)
   delete(T7)
   delete(T8)
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
 resP@txString = "JJA 0-10N radiative Q* (contours) and response (shading)"
 resP@lbTitleString    = "K/s"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= True
 resP@txFontHeightF    = 0.0085
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/4,3/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop" + plotfile +".pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(fils1)
;print(datalog)
end
