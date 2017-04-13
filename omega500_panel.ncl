begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc4")
plotfile="ANN_omega500"
wks = gsn_open_wks("epsi", plotfile)                      ;Open a workstation
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+2, graphic)    ;;;Create a panel plot with panels = number of models!
lets=(/"a)","b)","c)","d)"/)
n=0
do while(n.le.dimsizes(fils)-1)
filen=tostring(fils(n))
pr1= addfile(filen, "r")
MERRA=addfile("../merged_omega500_climo_remap.nc", "r")
land=pr1->LANDFRAC(0,{-30:30},:)
;pr2= addfile("../CAM/files/", "r")
;pr3=addfile("../CAM/files/","r")
;------------------------------------------------------------------------------------
;surface temp (ts)
;------------------------------------------------------------------------------------
p1 = pr1->OMEGA(0,18,{-30:30},:)
pr2= MERRA->OMEGA(0,0,{-30:30},:)
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


MERRAmax=max(pr2({-30:30},:))
MERRAmean=avg(pr2({-30:30},:))
MERRAvar=variance(pr2({-30:30},:))

rmse=(wgt_arearmse(p1({-30:30},:),pr2({-30:30},:),1.0,1.0,0))
print(rmse)
pcor=pattern_cor(p1({-30:30},:),pr2({-30:30},:),1.0,0)

rmseocean=(wgt_arearmse(ocean_only({-30:30},:),pr2({-30:30},:),1.0,1.0,0))
pcorocean=pattern_cor(ocean_only({-30:30},:),pr2({-30:30},:),1.0,0)

rmseland=(wgt_arearmse(land_only({-30:30},:),pr2({-30:30},:),1.0,1.0,0))
pcorland=pattern_cor(land_only({-30:30},:),pr2({-30:30},:),1.0,0)

print(rmseocean)
print(pcorocean)


meanpr=avg(p1({-30:30},:))
maxpr=max(p1({-30:30},:))
variancepr=variance(p1({-30:30},:))


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
res3@cnLevels = (/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)*1e-2
res3@cnFillColors = (/10,20,30,40,50,60,70,80,90,100,110,120,0,0,140,150,160,170,180,190,200,210,220,230,240,250/)
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
res3@gsnCenterString=  lets(n) + " " +  pr1@case
res3@gsnCenterStringFontHeightF=0.001
res3@gsnRightString= "RMSE: " + decimalPlaces(rmse,4,True) + " Pat. Cor.: " + decimalPlaces(pcor,3,True)
res3@gsnLeftString= "Variance: " + decimalPlaces(variancepr,5,True) + " Pa/s" 
res3@gsnRightStringFontHeightF= 0.01
res3@gsnLeftStringFontHeightF= 0.01
res3@mpGridLineThicknessF=0.5
res3@gsnCenterStringFontHeightF=0.012

;res3@gsnZonalMean    = True
;res3@gsnZonalMeanXMinF=-0.05
;res3@gsnZonalMeanXMaxF=0.05
;res3@gsnZonalMeanYRefLine=0
res3@tmXBLabelFontHeightF=0.006
res3@mpGridLineColor="black"

res3@tmXBLabelFontHeightF=0.006

gsn_define_colormap(wks,"nice_gfdl")


plot(n)= gsn_csm_contour_map(wks,p1,res3)
;;;;attach zonal means;;;
  zres=True
  zres@trXMinF              = -0.05      ; Could also use gsnZonalMeanXMinF
  zres@trXMaxF              = 0.05      ; Could also use gsnZonalMeanXMaxF
  zres@trYMinF              = -30
  zres@trYMaxF              =  30

  zres@tmXBMode             = "Explicit"
  zres@tmXBValues           = (/-0.05,0,0.05/)
  zres@tmXBLabels           = zres@tmXBValues + ""   ; Make it a string
  zres@tmXBLabelFontHeightF = 0.01                   ; Make slightly smaller.
  zres@vpWidthF             = 0.1
  zonal_id = gsn_csm_attach_zonal_means(wks,plot(n),p1,zres)

n=n+1
end do


;;;;Add MERRA rainfall as 3rd panel;;;



res30                           = res3

res30@gsnCenterString=lets(n) + " MERRA/ERAI"
res30@gsnRightString=  ""
res30@gsnLeftString= "Variance: " + decimalPlaces(MERRAvar,5,True) + " Pa/s"


plot(n)=gsn_csm_contour_map(wks,pr2,res30)
zonal_id = gsn_csm_attach_zonal_means(wks,plot(n),pr2,zres)

res4 = True
res4@gsnPanelLabelBar=True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = "Omega (Pa/s)"
res4@lbTitlePosition ="Bottom"
res4@lbTitleFontHeightF=.008
res4@txFontHeightF    = 0.03
res4@txString = "ANN"
res4@gsnMaximize=True
res4@gsnPaperOrientation="landscape"


gsn_panel(wks,plot,(/4,1/),res4)
;draw(wks)
;frame(wks)
delete(wks)
string1="eps2eps " + plotfile +".epsi" + " " + plotfile +".eps"
print(string1)
system(string1) 
system("epstopdf " + plotfile +".eps") 
system("rm " + plotfile +".eps") 
system("rm " + plotfile +".epsi")
end
