begin
;fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc")
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc2")
plotfile="ANN_pr_response"
;wtype="pdf"
;wtype@wkHeight =1400
wks = gsn_open_wks("epsi", plotfile)                      ;Open a workstation
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+2, graphic)    ;;;Create a panel plot with panels = number of models!
lets=(/"a)","b)","c)","d)"/)
n=0
do while(n.le.dimsizes(fils)-1)
filen=tostring(fils(n))
filen1=tostring(fils1(n))
pr1= addfile(filen, "r")
pr2=addfile(filen1,"r")
if (n.eq.3) then
  p=pr1->pr(0,{-30:30},:)
  p3=p*86400
  p7=pr2->pr(0,{-30:30},:)
  p6=p7*86400
  test=p6-p3
  print(test)
else
  p1 = pr1->PRECC(0,{-30:30},:)
  p2= pr1->PRECL(0,{-30:30},:)
  p4=pr2->PRECC(0,{-30:30},:)
  p5=pr2->PRECL(0,{-30:30},:)
  p3=(p1+p2)*86400000
  p6=(p4+p5)*86400000
end if
p3!0="lat"
p3&lat=pr1->lat({-30:30})
print(dimsizes(pr1))
print(dimsizes(p3))
p3!1="lon"
p3&lon=pr1->lon




p6!0="lat"
p6&lat=pr2->lat({-30:30})
p6!1="lon"
p6&lon=pr1->lon

diff=p6-p3
diff!0="lat"
diff&lat=pr1->lat({-30:30})
diff!1="lon"
diff&lon=pr1->lon


res3                           = True
res3@gsnDraw =  False
res3@gsnFrame = False 
res3@cnLinesOn                 = False
res3@mpProjection = "Mercator"
res3@mpCenterLonF         = -200              ; center plot at 180
res3@mpLimitMode="LatLon"
res3@mpMaxLatF = 30                      ; specify the plot domain
res3@mpMinLatF = -30
res3@mpGridAndLimbOn      = False          ; turn on limb and grid lines
res3@mpGridLatSpacingF = 15
res3@mpGridLonSpacingF = 45
res3@pmTickMarkDisplayMode  = "Always" 
res3@mpOutlineOn = True
res3@mpFillOn = False
res3@cnLevelSelectionMode = "ExplicitLevels" ; use explicit levels
res3@cnLevels = (/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12/)*0.5
res3@cnFillColors = (/240,230,220,210,200,190,180,170,160,150,140,130,0,120,110,105,95,90,85,80,70,60,50,40,30/)
res3@cnFillOn = True
res3@cnFillMode= "RasterFill"
res3@cnRasterSmoothingOn = "True"
res3@cnLineThicknessF = 0.25                   ; thicker lines
res3@cnLineLabelsOn   = False            
res3@cnInfoLabelOn = False                   ; turn on contour info label
res3@gsnSpreadColors     = False
res3@gsnMaximize      = True
res3@mpDataBaseVersion    = "MediumRes"    
res3@lbLabelBarOn        = False  
res3@gsnCenterStringFontHeightF= 0.012
res3@tmXTOn="False"
res3@tmYROn="False"
res3@tmYLOn="True"
res3@gsnCenterString= pr1@case
;res3@gsnRightString= "RMSE: " + decimalPlaces(rmse,2,True) + " Ocean: " + decimalPlaces(rmseocean,2,True) + " Land: " + decimalPlaces(rmseland,2,True) + "~C~Pat. Cor.: " + decimalPlaces(pcor,2,True) + " Ocean: " + decimalPlaces(pcorocean,2,True) + " Land: " + decimalPlaces(pcorland,2,True)
res3@gsnLeftString= lets(n)
res3@gsnRightStringFontHeightF= 0.0065
res3@gsnLeftStringFontHeightF= 0.012
res3@mpGridLineThicknessF=0.5

cres=True
cres@cnLevelSelectionMode = "ExplicitLevels"
cres@cnLevels = (/0,4,6,8,10,12/)
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False   
 cres@cnLineLabelsOn   = False 
 cres@cnLineThicknessF = 0.5

plot(n)= gsn_csm_contour_map(wks,diff,res3)
plot1=gsn_csm_contour(wks,p3,cres)

;;;;attach zonal means;;;
  zres=True
  zres@trXMinF              = -2.      ; Could also use gsnZonalMeanXMinF
  zres@trXMaxF              = 2.      ; Could also use gsnZonalMeanXMaxF
  zres@trYMinF              = -30
  zres@trYMaxF              =  30

  zres@tmXBMode             = "Explicit"
  zres@tmXBValues           = (/-2,-1,0,1,2/)
  zres@tmXBLabels           = zres@tmXBValues + ""   ; Make it a string
  zres@tmXBLabelFontHeightF = 0.01                   ; Make slightly smaller.
  zres@vpWidthF             = 0.1
  zonal_id = gsn_csm_attach_zonal_means(wks,plot(n),diff,zres)
overlay(plot(n),plot1)
n=n+1
 delete(p3)
 delete(p6)
end do


;;;;Add GPCP rainfall as 3rd panel;;;


res4 = True
res4@gsnPanelLabelBar=True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = "Precipitation rate (mm/day)"
res4@lbTitlePosition ="Bottom"
res4@lbTitleFontHeightF=.008
res4@txFontHeightF    = 0.03
;res4@pmLabelBarWidthF=0.3
;res4@txString = "ANN"
;res4@pmLabelBarOrthogonalPosF=0.95
res4@gsnMaximize=True
res4@gsnPaperOrientation="landscape"


gsn_panel(wks,plot,(/4,1/),res4)
delete(wks)
string1="eps2eps " + plotfile +".epsi" + " " + plotfile +".eps"
print(string1)
system(string1) 
system("epstopdf " + plotfile +".eps") 
system("rm " + plotfile +".eps") 
system("rm " + plotfile +".epsi")
;crop=systemfunc("/usr/bin/pdfcrop ANN_pr.pdf")
end
