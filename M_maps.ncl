begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*_ANN*.nc11")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*4K_ANN*.nc21")
plotfile="ANN_M_omega"
wks = gsn_open_wks("epsi", plotfile)                      ;Open a workstation
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+2, graphic)    ;;;Create a panel plot with panels = number of models!
lets=(/"a)","b)","c)","d)"/)
n=0
do while(n.le.dimsizes(fils)-1)
filen=tostring(fils(n))
filen1=tostring(fils1(n))
f= addfile(filen, "r")
f1=addfile(filen1,"r")
T=f->FMASS(0,:,:,:)
T1=f1->FMASS(0,:,:,:)

O=f->OMEGA(0,:,:,:)
O1=f1->OMEGA(0,:,:,:)

PS   = f->PS(0,:,:)
PS1   = f1->PS(0,:,:)
hyam = f->hyam 
hybm = f->hybm
P0   = 100000. 
levs=conform_dims(dimsizes(T),T&ilev,0)
opt=(/0,100,1000/)
wva = wgt_vert_avg_beta(levs, T, PS/100, 0, opt)/900
wva1 = wgt_vert_avg_beta(levs, T1, PS1/100, 0, opt)/900
MDIFF=wva1-wva
MDIFF!0="lat"
MDIFF&lat=f->lat
MDIFF!1="lon"
MDIFF&lon=f->lon
wva!0="lat"
wva&lat=f->lat
wva!1="lon"
wva&lon=f->lon
pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)
pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm) 
wva2 = wgt_vert_avg_beta(pm/100, O, PS/100, 0, opt)/900
wva3 = wgt_vert_avg_beta(pm/100, O1, PS1/100, 0, opt)/900
odiff=wva3-wva2
odiff!0="lat"
odiff&lat=f->lat
odiff!1="lon"
odiff&lon=f->lon
wva2!0="lat"
wva2&lat=f->lat
wva2!1="lon"
wva2&lon=f->lon

fmassavg=avg(wva)
omegavg=avg(odiff)
avgomeg=avg(omegavg)/avg(wva2)
avgm=avg(fmassavg)/avg(wva)
print(avgomeg)
print(avgm)
print(avg(wva))
pcor=pattern_cor(wva({-30:30},:),wva2({-30:30},:),1.0,0)

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
;res3@cnLevels = (/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12/)*0.5e-3
res3@cnLevels = (/-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15/)*1e-3
res3@cnFillColors = (/170,160,150,140,130,120,110,100,90,80,70,60,50,40,20,2,0/)
;res3@cnFillColors = (/240,230,220,210,200,190,180,170,160,150,140,130,0,120,110,105,95,90,85,80,70,60,50,40,30/)
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
res3@gsnCenterString= f@case
res3@gsnLeftString= lets(n)
res3@gsnRightStringFontHeightF= 0.012
res3@gsnLeftStringFontHeightF= 0.012
res3@mpGridLineThicknessF=0.5
res3@gsnRightString="r~B~pat~N~: " + sprintf("%5.2f",pcor)

cres=True
cres@cnLevelSelectionMode = "ExplicitLevels"
cres@cnLevels = (/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12/)*8e-3
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False   
 cres@cnLineLabelsOn   = False 
 cres@cnLineThicknessF = 2
 cres@gsnLeftString=""
 cres@gsnRightString=""
 cres@gsnContourZeroLineThicknessF  = 4.0
  cres@gsnContourNegLineDashPattern =2

gsn_define_colormap(wks,"nice_gfdl")
plot(n)= gsn_csm_contour_map(wks,wva,res3)
;plot1=gsn_csm_contour(wks,wva2,cres)
;;;;attach zonal means;;;
  zres=True
  zres@trXMinF              = -2*6e-3      ; Could also use gsnZonalMeanXMinF
  zres@trXMaxF              = 2*6e-3    ; Could also use gsnZonalMeanXMaxF
  zres@trYMinF              = -30
  zres@trYMaxF              =  30

  zres@tmXBMode             = "Explicit"

  zres@tmXBValues           = (/-2,0,2/)*6e-3
  zres@tmXBLabels           = zres@tmXBValues + ""   ; Make it a string
  zres@tmXBLabelFontHeightF = 0.01                   ; Make slightly smaller.
  zres@vpWidthF             = 0.1
 zonal_id = gsn_csm_attach_zonal_means(wks,plot(n),wva,zres)
;overlay(plot(n),plot1)
n=n+1
 delete(T)
 delete(T1)
 delete(MDIFF)
 delete(wva)
 delete(wva1)
end do


res4 = True
res4@gsnPanelLabelBar=True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = "M~B~int~N~ (kg m~S~-2~N~s~S~-1~N~)"
res4@lbTitlePosition ="Bottom"
res4@lbTitleFontHeightF=.015
res4@txFontHeightF    = 0.03
res4@pmLabelBarOrthogonalPosF=-0.02
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
end
