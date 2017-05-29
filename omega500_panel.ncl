begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc4")
plotfile="ANN_omega500"
wks = gsn_open_wks("pdf", plotfile)                      ;Open a workstation
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+2, graphic)    ;;;Create a panel plot with panels = number of models!
n=0
do while(n.le.dimsizes(fils)-1)
filen=tostring(fils(n))
pr1= addfile(filen, "r")
MERRA=addfile("MERRA_ANN_omega_climo_remap.nc", "r")
land=pr1->LANDFRAC(0,{-30:30},:)
;pr2= addfile("CAM/files/", "r")
;pr3=addfile("CAM/files/","r")
;------------------------------------------------------------------------------------
;surface temp (ts)
;------------------------------------------------------------------------------------
p1 = pr1->OMEGA(0,18,{-30:30},:)
pr2= MERRA->OMEGA500(0,{-30:30},:)
;pr2!0="lat"
;pr2&lat=MERRA->lat({-30:30})
;pr2!1="lon"
;pr2&lon=MERRA->lon
;p3!0="lat"
;p3&lat=pr1->lat({-30:30})
;p3!1="lon"
;p3&lon=pr1->lon



print(dimsizes(pr2))
print(dimsizes(p1))
diff=p1-pr2
diff!0="lat"
diff&lat=pr1->lat({-30:30})
diff!1="lon"
diff&lon=pr1->lon

ocean_only=mask(p1,land,0)
ocean_only!0="lat"
ocean_only&lat=pr1->lat({-30:30})
ocean_only!1="lon"
ocean_only&lon=pr1->lon

land_only=mask(p1,(land.ne.0),0.9999999)

land_only!0="lat"
land_only&lat=pr1->lat({-30:30})
land_only!1="lon"
land_only&lon=pr1->lon


MERRAmax=max(pr2({-15:15},:))
MERRAmean=avg(pr2({-15:15},:))
MERRAvar=variance(pr2({-15:15},:))

rmse=(wgt_arearmse(p1({-15:15},:),pr2({-15:15},:),1.0,1.0,0))
print(rmse)
pcor=pattern_cor(p1({-15:15},:),pr2({-15:15},:),1.0,0)

rmseocean=(wgt_arearmse(ocean_only({-15:15},:),pr2({-15:15},:),1.0,1.0,0))
pcorocean=pattern_cor(ocean_only({-15:15},:),pr2({-15:15},:),1.0,0)

rmseland=(wgt_arearmse(land_only({-15:15},:),pr2({-15:15},:),1.0,1.0,0))
pcorland=pattern_cor(land_only({-15:15},:),pr2({-15:15},:),1.0,0)

print(rmseocean)
print(pcorocean)


meanpr=avg(p1({-15:15},:))
maxpr=max(p1({-15:15},:))
variancepr=variance(p1({-15:15},:))


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
res3@gsnStringFontHeightF=0.005
res3@mpOutlineOn = True
res3@mpFillOn = False
res3@cnLevelSelectionMode = "ExplicitLevels" ; use explicit levels
res3@cnLevels = (/-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12/)*1e-2
res3@cnFillColors = (/4,12,20,28,40,60,0,0,150,170,180,190,200,210,240/)
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
res3@tmXTOn="False"
res3@tmYROn="False"
res3@tmYLOn="True"
res3@gsnCenterString=  pr1@case
res3@gsnCenterStringFontHeightF=0.008
res3@gsnRightString= "RMSE (MERRA, 15S-5N): " + decimalPlaces(rmse,3,True) + " Ocean: " + decimalPlaces(rmseocean,3,True) + " Land: " + decimalPlaces(rmseland,3,True) + "~C~Pat. Cor. (MERRA, 15S-5N): " + decimalPlaces(pcor,3,True) + " Ocean: " + decimalPlaces(pcorocean,3,True) + " Land: " + decimalPlaces(pcorland,3,True)
res3@gsnLeftString= "Mean (15S-5N): " + decimalPlaces(meanpr,3,True) + " Pa/s" + "~C~Max (15S-5N): " + decimalPlaces(maxpr,3,True) + " Pa/s" +  "~C~Variance (15S-5N): " + decimalPlaces(variancepr,3,True) + " Pa/s" 
res3@gsnRightStringFontHeightF= 0.006
res3@gsnLeftStringFontHeightF= 0.006
res3@mpGridLineThicknessF=0.5
res3@gsnCenterStringFontHeightF=0.012

res3@gsnZonalMean    = True
res3@gsnZonalMeanXMinF=-0.05
res3@gsnZonalMeanXMaxF=0.05
res3@gsnZonalMeanYRefLine=0
res3@tmXBLabelFontHeightF=0.006
res3@mpGridLineColor="black"

res3@gsnRightStringFontHeightF= 0.006
res3@gsnLeftStringFontHeightF= 0.006
res3@tmXBLabelFontHeightF=0.006

gsn_define_colormap(wks,"NCV_jaisnd")


plot(n)= gsn_csm_contour_map(wks,p1,res3)

n=n+1
end do


;;;;Add MERRA rainfall as 3rd panel;;;



res30                           = res3

res30@gsnCenterString="MERRA (1979-2013)"
res30@gsnRightString=  "Variance (15S-5N):" + decimalPlaces(MERRAvar,3,True) + " Pa/s"
res30@gsnLeftString= "Mean (15S-5N): " + decimalPlaces(MERRAmean,3,True) + " Pa/s" + "~C~Max (15S-5N): " + decimalPlaces(MERRAmax,3,True) + " Pa/s"


plot(n)=gsn_csm_contour_map(wks,pr2,res30)


res4 = True
res4@gsnPanelLabelBar=True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = "500 hPa Omega (Pa/s)"
res4@lbTitlePosition ="Top"
res4@lbTitleFontHeightF=.004 
res4@txFontHeightF    = 0.01
;res4@pmLabelBarWidthF=0.3
res4@txString = "ANN"
;res4@pmLabelBarOrthogonalPosF=0.95


gsn_panel(wks,plot,(/4,2/),res4)
;draw(wks)
;frame(wks)
delete(wks)
crop=systemfunc("/usr/bin/pdfcrop ANN_omega500.pdf")
end
