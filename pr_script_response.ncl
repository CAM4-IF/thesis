begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc2")
plotfile="ANN_pr_response"
wks = gsn_open_wks("epsi", plotfile)                      ;Open a workstation
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+2, graphic)    ;;;Create a panel plot with panels = number of models!
lets=(/"a)","b)","c)","d)"/)
names=(/"CAM4 default","CAM4-IF-r","CAM4-IF-t","AMIP4K"/)
n=0
do while(n.le.dimsizes(fils)-1)
filen=tostring(fils(n))
filen1=tostring(fils1(n))
print(filen1)
pr1= addfile(filen, "r")
pr2=addfile(filen1,"r")
pr3=addfile("stip.nc","r")
pr4=addfile("stip1.nc","r")
pr5=addfile("amip4k_pr_maxresponse.nc","r")
pr6=addfile("amip4k_pr_minresponse.nc","r")
stip=pr3->pr(0,{-30:30},:)
stip1=pr4->pr(0,{-30:30},:)
if (n.eq.3) then
  p=pr1->pr(0,{-30:30},:)
  p3=p*86400
  p7=pr2->pr(0,{-30:30},:)
  p6=p7*86400
  diff=p6
else
  p1 = pr1->PRECC(0,{-30:30},:)
  p2= pr1->PRECL(0,{-30:30},:)
  p4=pr2->PRECC(0,{-30:30},:)
  p5=pr2->PRECL(0,{-30:30},:)
  p3=(p1+p2)*86400000
  p6=(p4+p5)*86400000
  maximum=pr5->pr(0,{-30:30},:)
  m=maximum*86400
  minimum=pr6->pr(0,{-30:30},:)
  m1=minimum*86400
end if
p3!0="lat"
p3&lat=pr1->lat({-30:30})
p3!1="lon"
p3&lon=pr1->lon



if (n.eq.3) then
else
p6!0="lat"
p6&lat=pr2->lat({-30:30})
p6!1="lon"
p6&lon=pr1->lon
diff=p6-p3
diffmax=diff-m
diffmin=diff-m1
diffmax!0="lat"
diffmax&lat=pr2->lat({-30:30})
diffmax!1="lon"
diffmax&lon=pr2->lon
diffmin!0="lat"
diffmin&lat=pr2->lat({-30:30})
diffmin!1="lon"
diffmin&lon=pr2->lon
end if
diff!0="lat"
diff&lat=pr2->lat({-30:30})
diff!1="lon"
diff&lon=pr2->lon



;;;;Calculate tropical means (grid-area weighted);;;;;

pi=4.*atan(1.0) 
rad=(pi/180.) 
coswgt=cos(rad*p3&lat) 
coswgt!0 = "lat" 
coswgt&lat= p3&lat 
llats = -30. 
llatn = 30. 
llonw = 0. 
llone = 360. 
aa=wgt_areaave(p3({llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)
delete(coswgt)
coswgt=cos(rad*diff&lat) 
coswgt!0 = "lat" 
coswgt&lat= diff&lat 
aa1=wgt_areaave(diff({llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)
response=100*(aa1/aa)





;;;;Create the plot resources;;;;


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
res3@cnFillColors = (/245,240,230,220,210,200,190,180,170,160,150,140,0,110,105,95,90,85,80,70,60,55,50,45,40/)
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
res3@gsnCenterString= names(n)
res3@gsnRightString= "Response: " + sprintf("%5.2f",response) + "%"
res3@gsnLeftString= lets(n)
res3@gsnRightStringFontHeightF= 0.012
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

if (n.eq.3) then
else
;;;Add stipling where response is larger than max response of AMIP4K members;;;
gt=diffmax>0
gt!0="lat"
gt&lat=pr2->lat({-30:30})
gt!1="lon"
gt&lon=pr2->lon

lt=diffmin<0
lt!0="lat"
lt&lat=pr2->lat({-30:30})
lt!1="lon"
lt&lon=pr2->lon
print(lt)
res30 = True                            
  res30@gsnDraw=False
  res30@gsnFrame=False
  res30@cnFillOn = True                   ; turn on color/pattern fill
  res30@cnMonoFillPattern = False   ; allow different patterns
  res30@cnMonoFillColor = True       ; only use one color (black)
  res30@cnLevelSelectionMode = "ExplicitLevels" ; set explicit cnlev
  res30@cnLevels   = (/0,0.1/)    ; 
  res30@cnFillPatterns = (/-1,-1,17/) ;
  res30@gsnAddCyclic = True   ; may or may not be needed
  res30@cnInfoLabelOn       = False       ; turn off info label
  res30@cnLinesOn           = False       ; do not draw contour lines
  res30@cnLineLabelsOn      = False    ; do not draw contour labels
  res30@cnFillScaleF        = 0.75         ; add extra density
  res30@gsnRightString   = ""    
  res30@gsnLeftString    = ""
  res30@gsnCenterString  = ""
  res30@lbLabelBarOn=False
  
  plot2   = gsn_csm_contour(wks,gt, res30)
  
  
  res31=res30
  res31@cnLevels   = (/-0.1,0/)    ; 
  res31@cnFillPatterns = (/6,-1,-1/) ; 
  plot3   = gsn_csm_contour(wks,lt, res31)
 overlay(plot(n),plot2)
 overlay(plot(n),plot3)
  delete(m1)
 delete(diffmin)
end if


n=n+1
 delete(p3)
 delete(p6)
 delete(coswgt)
 delete(diff)

end do

;;
;;
;Add stipling where 10 or more models agree on sign
delete(res3)
res3 = True                            
  res3@gsnDraw=False
  res3@gsnFrame=False
  res3@cnFillOn = True                   ; turn on color/pattern fill
  res3@cnMonoFillPattern = False   ; allow different patterns
  res3@cnMonoFillColor = True       ; only use one color (black)
  res3@cnLevelSelectionMode = "ExplicitLevels" ; set explicit cnlev
  res3@cnLevels   = (/0,1,2/)    ; only have 1 contour level
  res3@cnFillPatterns = (/-1,17,17,17/) ; don't fill <0.95, stipple >=0.95
  res3@gsnAddCyclic = True   ; may or may not be needed
  res3@cnInfoLabelOn       = False       ; turn off info label
  res3@cnLinesOn           = False       ; do not draw contour lines
  res3@cnLineLabelsOn      = False    ; do not draw contour labels
  res3@cnFillScaleF        = 0.75         ; add extra density
  res3@gsnRightString   = ""    
  res3@gsnLeftString    = ""
  res3@gsnCenterString  = ""
  res3@lbLabelBarOn=False
  
  plot2   = gsn_csm_contour(wks,stip, res3)
 
  res5 = res3                           
  res5@cnFillPatterns = (/-1,6,6,6/) 
  res5@cnLevels   = (/10,11,12/)
  
  plot3   = gsn_csm_contour(wks,stip1, res5)
  
overlay(plot(3),plot2)
overlay(plot(3),plot3)

;;;;Add GPCP rainfall as 3rd panel;;;


res4 = True
res4@gsnPanelLabelBar=True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = "P (mm day~S~-1~N~)"
res4@lbTitlePosition ="Bottom"
res4@lbTitleFontHeightF=.012
res4@txFontHeightF    = 0.03
res4@gsnMaximize=True
res4@gsnPaperOrientation="landscape"
res4@pmLabelBarOrthogonalPosF=-0.02


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
