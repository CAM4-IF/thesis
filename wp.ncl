begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="IWP_diff" ;;;PDF file name!
wks = gsn_open_wks("png", plotfile) ;; open a workstation file!
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+1, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   
   f= addfile(filen, "r")
   f1=addfile("CERES_ANN.nc","r")
   hyam = f->hyam 
   hybm = f->hybm
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)
   Q=f->TGCLDIWP(0,:,:)
   Q1=f1->iwp_total_mon(0,:,:)
   Qdiff=Q-Q1
   Qdiff!0="lat"
   Qdiff&lat=Q&lat
   Qdiff!1="lon"
   Qdiff&lon=Q&lon
   
   
   res                      = True                
   res@tmYRMode = "Automatic"
   res@cnLevelSelectionMode = "ExplicitLevels"   
   ;res@cnLevels = (/5,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,220,240/)*1.5
  ; res@cnFillColors = (/0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100/)   
  res@cnLevels = (/-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100/)*2
 res@cnFillColors = (/2,6,12,24,40,60,70,90,100,0,0,0,0,156,170,180,190,210,220,230,240,250/)
   res@cnLineLabelsOn       = False                 
  res@cnFillOn             = True                
  res@lbLabelBarOn = False
  res@lbLabelAutoStride    = True                 
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False  
  res@cnLineThicknessF =0.5                 
  res@cnLineLabelsOn   = False                
  res@cnInfoLabelOn = False                 
  res@gsnSpreadColors     = False
  res@gsnMaximize      = True
  res@mpDataBaseVersion    = "MediumRes"  
  res@mpMaxLatF =35; specify the plot domain
  res@mpMinLatF = -35
  res@mpCenterLonF         = -200
  res@tmXTOn="False"
  res@gsnRightString    = f@case
  res@gsnRightStringFontHeightF=0.01
  res@lbTitleString    = "K/s"
  res@lbTitleFontHeightF=.007   
  res@gsnLeftString  = ""
 ; res@cnFillPalette="NCV_jaisnd"
  
  
  plot(n)  = gsn_csm_contour_map(wks, Qdiff, res)  
   n=n+1
end do

 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "ANN IWP (model minus CERES)"
 resP@lbLabelBarOn = True
 resP@lbTitleString    = "g/m^2"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 ;resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.015
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 ;resP@lbLabelStride=1
 gsn_panel(wks,plot,(/5,1/),resP)  
 print(fils)
 ;delete(wks) 
;crop=systemfunc("/usr/bin/pdfcrop test.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
end
