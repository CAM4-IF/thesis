begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc1")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="omega_ANN" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile(filen, "r")
   hyam = h->hyam 
   hybm = h->hybm
   cmd="cdo -mermean -sellonlatbox,0,360,-10,10 -selname,OMEGA,PS " + filen + " tmp.nc"
   cmd1="cdo -mermean -sellonlatbox,0,360,-10,10 -selname,OMEGA,PS " + filen1 + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->OMEGA(0,:,:,:)
   T&lev@units="hPa"
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt=(/0,10000,100000/)
   wva = wgt_vert_avg_beta(pm, T, PS, 1, opt)/90000
   wva!0="lev"
   wva!1="lon"
   wva&lon=T&lon
   
   PS1   = f1->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->OMEGA(0,:,:,:)
   T1&lev@units="hPa"
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt1=(/0,10000,100000/)
   wva1 = wgt_vert_avg_beta(pm1, T1, PS1, 1, opt1)/90000
   wva1!0="lev"
   wva1!1="lon"
   wva1&lon=T&lon
   
   
   
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
 
  Qdiff=T1(:,0,:)-T(:,0,:)
  Qdiff!0="lev"
  Qdiff&lev=T&lev
  Qdiff!1="lon"
  Qdiff&lon=T&lon
  Qdiff&lev@units="hPa"
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-30,-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20,30/)*5.0e-4
  res@cnFillColors = (/2,6,12,18,24,40,60,90,106,116,0,0,142,148,156,168,180,208,218,230,250,255/)
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

 
;;;;Calculate 500 hPa upward omega mean;;;; 
   omegstar=T(:,0,:)
   omegstar!1="lon"
   omegstar&lon=T&lon
   omegstar1=T1(:,0,:)
   omegstar1!1="lon"
   omegstar1&lon=T&lon
   
   omegstarwc=dim_rmvmean_n_Wrap(T(:,0,:),1)
   omegstarwc1=dim_rmvmean_n_Wrap(T1(:,0,:),1)
   omegstarwc!1="lon"
   omegstarwc&lon=T&lon
   omegstarwc1!1="lon"
   omegstarwc1&lon=T&lon
   
   omegstarwc500=omegstarwc(18,{60:180}) < 0
   omegstarwc500mean=avg(omegstarwc500)
   omegstarwc5001=omegstarwc1(18,{60:180}) < 0
   omegstarwc500mean1=avg(omegstarwc5001)
   
   omegstarwc500diff=omegstarwc500mean1-omegstarwc500mean
   responsewc=100*(omegstarwc500diff/avg(omegstarwc500mean))
   
   
   omegstar500=omegstar(18,:)
   omegstar500=omegstar500 < 0
   omegstar500mean=avg(omegstar500)
   
   omegstar5001=omegstar1(18,:) < 0
   omegstar500mean1=avg(omegstar5001)
   
   omegstar500diff=omegstar500mean1-omegstar500mean
   
   response=100*(omegstar500diff/avg(omegstar500mean))
   
   print(response)
   
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels = (/-260,-240,-220,-200,-180,160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240,260/)*5.0e-4
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
 cres@gsnLeftString  = "w_500_up response (%):" + response + "~C~w_500_up mean:" + omegstar500mean + "~C~w_500_up* WC response (%):" + responsewc + "~C~w_500_up* WC mean: " + omegstarwc500mean
 cres@gsnLeftStringFontHeightF=0.012
 cres@gsnRightStringFontHeightF=0.015
 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 plot(n)  = gsn_csm_pres_hgt(wks, Qdiff(::-1,:), res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_pres_hgt(wks, T(::-1,0,:), cres)
 overlay(plot(n),plot2)
 
   
   cleanup=systemfunc("rm tmp*.nc")
   delete(T)
   delete(T1)
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
 resP@txString = "CAM4/CAM-IF ANN 0-10N omega response (shaded) and control climatology (contours)"
 resP@lbTitleString    = "K/s"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.0085
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/3,3/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop" + plotfile)  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(fils1)
end
