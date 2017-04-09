begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*.nc2")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="t_ANN_rei" ;;;PDF file name!
wks = gsn_open_wks("png", plotfile) ;; open a workstation file!
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+1, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   
   f= addfile(filen, "r")
   hyam = f->hyam 
   hybm = f->hybm
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)
   Q=f->DTCOND(0,:,:,:)
   Q1=f->QRL(0,:,:,:)
   Q2=f->QRS(0,:,:,:)
   Q3=f->DTV(0,:,:,:)
   prc=f->PRECC(0,:,:)
   prl=f->PRECL(0,:,:)
   precip=prc+prl
   precip!0="lat"
   precip&lat=Q&lat
   precip!1="lon"
   precip&lon=Q&lon
   Qtot=Q+Q1+Q2+Q3
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt=(/0,10000,92500/)
   wva = wgt_vert_avg_beta(pm, Q, PS, 0, opt)/(77500)
   wva!0="lat"
   wva&lat=Q&lat
   wva!1="lon"
   wva&lon=Q&lon
   tropmean=avg(wva({-10:10},{90:210}))
   pmean=avg(precip({-10:10},{90:210}))
   
   
   res                      = True                
   res@tmYRMode = "Automatic"
   res@cnLevelSelectionMode = "ExplicitLevels"   
   res@cnLevels = (/0,5,10,15,20,25,30,35,40/)*3e-6
   ;res@cnLevels = (/-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20/)*3e-5
   ;res@cnLevels =(/-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18/)
   ;res@cnFillColors = (/2,6,12,24,40,60,90,100,0,0,0,0,156,170,190,210,220,230,240,250/)
   res@cnFillColors = (/0,166,176,186,196,206,216,226,236,246/)
   res@cnLinesOn       = False               
   res@cnLineLabelsOn       = False                 
  res@cnFillOn             = True                
  res@lbLabelBarOn = True
  res@lbLabelAutoStride    = True                 
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                
  res@cnFillOpacityF=0.7
  res@cnLineThicknessF = 1.                 
  res@cnLineLabelsOn   = False                
  res@cnInfoLabelOn = False                 
  res@gsnSpreadColors     = False
  res@gsnMaximize      = True
  res@mpDataBaseVersion    = "MediumRes"  
  res@mpMaxLatF =35; specify the plot domain
  res@mpMinLatF = -35
  res@mpCenterLonF         = -200
  res@tmXTOn="False"
  res@gsnZonalMean    = True
  res@gsnZonalMeanXMinF=0e-6
  res@gsnZonalMeanXMaxF=20e-6
  res@gsnZonalMeanYRefLine=0
  res@gsnRightString    = f@case
  res@lbTitleString    = "K/s"
  res@lbTitleFontHeightF=.007   
  res@gsnLeftString  = "10S-0N/90E-210E mean: " + tropmean + "~C~" + pmean
  
    ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels = (/-3.6,-3.4,-3.2,-3.,-2.8,-2.6,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6/)*2e-5
 ;cres@cnLevels=(/-100,-50,-25,-10,0,10,25,50,100/)
 cres@cnLineLabelsOn       = True                 ; turn on line labels
 cres@cnFillOn             = False                 ; turn on color fill
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False                ; turn off contour info label
 cres@cnLineThicknessF =0.5
 cres@gsnRightString    = f@case
 cres@gsnRightStringFontHeightF    = 0.05
 cres@cnLevelSpacingF=200.0
 cres@gsnContourZeroLineThicknessF  = 2.0
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelFontHeightF      = 0.004
 cres@cnLineLabelPlacementMode = "Constant"
 cres@gsnContourNegLineDashPattern =2
 cres@cnLineLabelInterval=1
 ;cres@gsnLeftString  = "10S-10N response pattern correlation: " + pcor + "~C~10S-10N tropical mean Q:" + qstarmean
 ;cres@gsnLeftString  = "10S-10N response pattern correlation: " + pcor
 cres@gsnLeftStringFontHeightF=0.025
 cres@tmXTOn="False"
 cres@gsnZonalMean    = True
 
 plot(n)  = gsn_csm_contour_map(wks, wva, res)  
   n=n+1
end do
   files=tostring(fils(0))
   f1= addfile(files, "r")
   PS1   = f1->PS(0,:,:)
   Q4=f1->DTCOND(0,:,:,:)
   Q5=f1->QRL(0,:,:,:)
   Q6=f1->QRS(0,:,:,:)
   Q7=f1->DTV(0,:,:,:)
   Qtot1=Q4+Q5+Q6
   pm1   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   wva1 = wgt_vert_avg_beta(pm1, Q4, PS1, 0, opt)/(77500)
   wva1!0="lat"
   wva1&lat=Q4&lat
   wva1!1="lon"
   wva1&lon=Q4&lon
   
   diff=wva-wva1
   diff!0="lat"
   diff&lat=Q4&lat
   diff!1="lon"
   diff&lon=Q4&lon
   res1=res
   res1@cnLevels =(/-20,-15,-10,-5,0,5,10,15,20/)*1e-6
   res1@cnFillColors = (/10,20,30,40,0,0,170,180,190,200/)
   res1@gsnZonalMeanXMinF=-5e-6
   res1@gsnZonalMeanXMaxF=5e-6
   res1@gsnZonalMeanYRefLine=0
   res1@gsnRightString    = "CAM-IFv4-CAM4"
   res1@gsnLeftString  =""
   
  
   
   plot(n)  = gsn_csm_contour_map(wks, diff, res1)

 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "ANN mean 925-150 hPa pressure weighted sum of condensational heating"
 resP@lbTitleString    = "K/s"
 resP@lbTitleFontHeightF=.007
  resP@lbLabelBarOn=False
 ;;resP@lbTitlePosition ="Top"
 ;resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.01
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 ;resP@lbLabelStride=1
 gsn_panel(wks,plot,(/3,1/),resP)  
 print(fils)
 ;delete(wks) 
;crop=systemfunc("/usr/bin/pdfcrop test.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
end
