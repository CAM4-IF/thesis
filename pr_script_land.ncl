begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc")
plotfile="ANN_pr"
wks = gsn_open_wks("pdf", plotfile)                      ;Open a workstation
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+2, graphic)    ;;;Create a panel plot with panels = number of models!
n=0
do while(n.le.dimsizes(fils)-1)
filen=tostring(fils(n))
pr1= addfile(filen, "r")
GPCP=addfile("PRECL_ANN_remap_climo.nc", "r")
land=pr1->LANDFRAC(0,{-30:30},:)
;pr2= addfile("CAM/files/", "r")
;pr3=addfile("CAM/files/","r")
;------------------------------------------------------------------------------------
;surface temp (ts)
;------------------------------------------------------------------------------------
p1 = pr1->PRECC(0,{-30:30},:)
p2= pr1->PRECL(0,{-30:30},:)
p3=(p1+p2)*86400000
pr2= GPCP->PRECT_LAND(0,{-30:30},:)
pr2!0="lat"
pr2&lat=GPCP->lat({-30:30})
pr2!1="lon"
pr2&lon=GPCP->lon
p3!0="lat"
p3&lat=pr1->lat({-30:30})
p3!1="lon"
p3&lon=pr1->lon

diff=p3-pr2
diff!0="lat"
diff&lat=pr1->lat({-30:30})
diff!1="lon"
diff&lon=pr1->lon

ocean_only=mask(p3,pr2,-999)
ocean_only!0="lat"
ocean_only&lat=pr1->lat({-30:30})
ocean_only!1="lon"
ocean_only&lon=pr1->lon

land_only=mask(p3,(pr2.ne.-999),0.9999999)

land_only!0="lat"
land_only&lat=pr1->lat({-30:30})
land_only!1="lon"
land_only&lon=pr1->lon


GPCPmax=max(pr2({-30:30},:))
GPCPmean=avg(pr2({-30:30},:))
GPCPvar=variance(pr2({-30:30},:))

rmse=(wgt_arearmse(p3({-30:30},:),pr2({-30:30},:),1.0,1.0,0))
print(rmse)
pcor=pattern_cor(p3({-30:30},:),pr2({-30:30},:),1.0,0)

;rmseocean=(wgt_arearmse(ocean_only({-30:30},:),pr2({-30:30},:),1.0,1.0,0))
;pcorocean=pattern_cor(ocean_only({-30:30},:),pr2({-30:30},:),1.0,0)

rmseland=(wgt_arearmse(land_only({-30:30},:),pr2({-30:30},:),1.0,1.0,0))
pcorland=pattern_cor(land_only({-30:30},:),pr2({-30:30},:),1.0,0)

;print(rmseocean)
;print(pcorocean)
rmseocean=0
pcorocean=0

meanpr=avg(p3({-30:30},:))
maxpr=max(p3({-30:30},:))
variancepr=variance(p3({-30:30},:))


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
res3@cnLevels = (/3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12/)
res3@cnFillColors = (/0,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240/)
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
res3@gsnCenterString=  pr1@case
res3@gsnRightString= "RMSE : " + decimalPlaces(rmse,2,True) + " Ocean: " + decimalPlaces(rmseocean,2,True) + " Land: " + decimalPlaces(rmseland,2,True) + "~C~Pat. Cor. : " + decimalPlaces(pcor,2,True) + " Ocean: " + decimalPlaces(pcorocean,2,True) + " Land: " + decimalPlaces(pcorland,2,True)
res3@gsnLeftString= "Mean : " + decimalPlaces(meanpr,2,True) + " mm/day" + "~C~Max : " + decimalPlaces(maxpr,2,True) + " mm/day"  + "~C~Variance : " + decimalPlaces(variancepr,2,True) + " mm/day"
res3@gsnRightStringFontHeightF= 0.0065
res3@gsnLeftStringFontHeightF= 0.0065
res3@mpGridLineThicknessF=0.5

;res3@gsnZonalMean    = True
;res3@gsnZonalMeanXMinF=0
;res3@gsnZonalMeanXMaxF=8
;res3@gsnZonalMeanYRefLine=0
;res3@tmXBLabelFontHeightF=0.006
;res3@mpGridLineColor="black"





gsn_define_colormap(wks,"WhiteGreen")


plot(n)= gsn_csm_contour_map(wks,p3,res3)
;;;;attach zonal means;;;
  zres=True
  zres@trXMinF              = 0.      ; Could also use gsnZonalMeanXMinF
  zres@trXMaxF              = 7.      ; Could also use gsnZonalMeanXMaxF
  zres@trYMinF              = -30
  zres@trYMaxF              =  30

  zres@tmXBMode             = "Explicit"
  zres@tmXBValues           = (/0,3,6/)
  zres@tmXBLabels           = zres@tmXBValues + ""   ; Make it a string
  zres@tmXBLabelFontHeightF = 0.01                   ; Make slightly smaller.
  zres@vpWidthF             = 0.1
  zonal_id = gsn_csm_attach_zonal_means(wks,plot(n),p3,zres)

n=n+1
end do


;;;;Add GPCP rainfall as 3rd panel;;;



res30                           = res3

res30@gsnCenterString="PREC/L (1948 - 2000)"
res30@gsnRightString= "Mean : " + GPCPmean + " mm/day"
res30@gsnLeftString=  "Max : " + GPCPmax + " mm/day" + "~C~Variance : " + GPCPvar + " mm/day"


plot(n)=gsn_csm_contour_map(wks,pr2,res30)

  zonal_id = gsn_csm_attach_zonal_means(wks,plot(n),pr2,zres)



res4 = True
res4@gsnPanelLabelBar=True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = "Precipitation rate (mm/day)"
res4@lbTitlePosition ="Top"
res4@lbTitleFontHeightF=.004 
res4@txFontHeightF    = 0.015
;res4@pmLabelBarWidthF=0.3
res4@txString = "ANN - Land Only"
;res4@pmLabelBarOrthogonalPosF=0.95


gsn_panel(wks,plot,(/4,2/),res4)
;draw(wks)
;frame(wks)
delete(wks)
crop=systemfunc("/usr/bin/pdfcrop ANN_pr.pdf")
end