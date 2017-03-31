begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_AMIP/AMIP/psl/*ANN*.nc")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_AMIP/AMIP4K/psl/*ANN*.nc")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="CAM_ANN_WC" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   cmd="cdo -fldmean -sellonlatbox,200,280,-5,5 -selname,psl " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,200,280,-5,5 -selname,psl " + filen1 + " tmp1.nc"
   cmd2="cdo -fldmean -sellonlatbox,100,180,-5,5 -selname,psl " + filen + " tmp2.nc"
   cmd3="cdo -fldmean -sellonlatbox,100,180,-5,5 -selname,psl " + filen1 + " tmp3.nc"
   tahiti=systemfunc(cmd)
   tahiti1=systemfunc(cmd1)
   darwin=systemfunc(cmd2)
   darwin1=systemfunc(cmd3)
   f2=addfile("tmp.nc", "r")
   f3=addfile("tmp1.nc","r")
   f4=addfile("tmp2.nc","r")
   f5=addfile("tmp3.nc", "r")
   print(fils)
   print(fils1)
   
 
   
   psl=f2->psl(0,:,:) 
   psl1=f3->psl(0,:,:) 
   psl2=f4->psl(0,:,:) 
   psl3=f5->psl(0,:,:) 
   
   
   td=psl2-psl
   td1=psl3-psl1
   
   response=100*((td1-td)/td)
   datalog(n)=response(0,0)
   print(response)
   
   f    = addfile(filen, "r")
   f1=addfile(filen1,"r")
   T=f->psl(0,:,:) 
   T1=f1->psl(0,:,:)
   
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
 
  Qdiff=T1-T
  Qdiff!1="lon"
  Qdiff&lon=T&lon
  Qdiff!0="lat"
  Qdiff&lat=T&lat
  
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20/)*10
  res@cnFillColors = (/2,6,12,24,40,60,90,100,0,0,0,0,156,170,190,210,220,230,240,250/)
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@lbLabelBarOn = False
  res@lbLabelAutoStride    = True                 ; optimal labels
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                ; turn off contour info label
  ;res@tiXAxisString    = ""
  ;res@tiYAxisSide = "Left"	
  ;res@tiYAxisString    = "Pressure (hPa)"
  ;res@trYMaxF = 1000.0
  ;res@trYMinF = 100.0      ; troposphere only
  ;res@tiXAxisFontHeightF =0.012
  ;res@tmYROn = False
  ;res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7
  res@cnLineThicknessF = 1.                   ; thicker lines
  res@cnLineLabelsOn   = False                
  res@cnInfoLabelOn = False                   ; turn on contour info label
  res@gsnSpreadColors     = False
  res@gsnMaximize      = True
  res@mpDataBaseVersion    = "MediumRes"  
  res@mpMaxLatF = 30                    ; specify the plot domain
  res@mpMinLatF = -30
  res@mpCenterLonF         = -200
  res@tmXTOn="False"

 
   ;control=dim_avg_n_Wrap(wva4,1)
   ;case=dim_avg_n_Wrap(wva9,1)
   ;response=case-control
   ;percentresponse=(response/control)*100.
   ;datalog(n)=percentresponse
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@cnLevelSelectionMode = "AutomaticLevels"
 cres@cnLineLabelsOn       = True                 ; turn on line labels
 cres@cnFillOn             = False                 ; turn on color fill
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False                ; turn off contour info label
 cres@cnLineThicknessF =1
 cres@gsnRightString    = f@model_id
 cres@gsnRightStringFontHeightF    = 0.05
 cres@cnLevelSpacingF=200.0
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelFontHeightF      = 0.004
 cres@cnLineLabelPlacementMode = "Constant"
 cres@cnLineLabelInterval=1
 cres@gsnLeftString  = "SOI response (%):" + response(0,0) +"~C~SOI mean:" + td(0,0)
 cres@gsnLeftStringFontHeightF=0.017
 cres@tmXTOn="False"
 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 
 
 plot(n)  = gsn_csm_contour_map(wks, Qdiff, res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_contour(wks, T, cres)
 ;plot3=gsn_csm_vector(wks,udiff,vdiff,vcres)
 overlay(plot(n),plot2)
 ;overlay(plot(n),plot3)
 
   
   cleanup=systemfunc("rm tmp*.nc")
   delete(T)
   delete(T1)
   delete(Qdiff)
   
   
n=n+1
end do
 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "ANN AMIP4K-AMIP_control SLP response (shaded) and control climatology (contours)"
 resP@lbTitleString    = "Pa"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.0085
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/5,2/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop CAM_ANN_WC.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(fils1)
delete(fils)
delete(fils1)
print(datalog)
end
