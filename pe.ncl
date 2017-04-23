begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*_ANN*.nc11")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*4K_ANN*.nc21")
fils2=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc4")
fils3=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc2")
plotfile="ANN_pe"
wks = gsn_open_wks("epsi", plotfile)                      ;Open a workstation
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+2, graphic)    ;;;Create a panel plot with panels = number of models!
lets=(/"a)","b)","c)","d)"/)
n=0
do while(n.le.dimsizes(fils)-1)
filen=tostring(fils(n))
filen1=tostring(fils1(n))
filen2=tostring(fils2(n))
filen3=tostring(fils3(n))

pr1= addfile(filen2, "r")
pr2=addfile(filen3,"r")
 p1 = pr1->PRECC(0,{-30:30},:)
  p2= pr1->PRECL(0,{-30:30},:)
  p4=pr2->PRECC(0,{-30:30},:)
  p5=pr2->PRECL(0,{-30:30},:)
  p3=(p1+p2)*1000
  p6=(p4+p5)*1000
p6!0="lat"
p6&lat=pr2->lat({-30:30})
p6!1="lon"
p6&lon=pr2->lon
p3!0="lat"
p3&lat=pr1->lat({-30:30})
p3!1="lon"
p3&lon=pr1->lon
q1=pr1->Q(0,23:25,{-30:30},:)
q2=pr2->Q(0,23:25,{-30:30},:)
q1avg=dim_avg_n_Wrap(q1,0)
q2avg=dim_avg_n_Wrap(q2,0)
mprime=(p3)/(q1avg)
mprime!0="lat"
mprime&lat=pr2->lat({-30:30})
mprime!1="lon"
mprime&lon=pr2->lon

mprime1=(p6)/(q2avg)

mprime1!0="lat"
mprime1&lat=pr2->lat({-30:30})
mprime1!1="lon"
mprime1&lon=pr2->lon


f= addfile(filen, "r")
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
wva!0="lat"
wva&lat=pr1->lat
wva!1="lon"
wva&lon=pr1->lon
wva1 = wgt_vert_avg_beta(levs, T1, PS1/100, 0, opt)/900

wva1!0="lat"
wva1&lat=pr1->lat
wva1!1="lon"
wva1&lon=pr1->lon


pemap=wva({-30:30},:)/(mprime)
pemap1=wva1({-30:30},:)/(mprime1)
pemap!0="lat"
pemap&lat=pr1->lat({-30:30})
pemap!1="lon"
pemap&lon=pr1->lon

pemap1!0="lat"
pemap1&lat=pr1->lat({-30:30})
pemap1!1="lon"
pemap1&lon=pr1->lon


preciponly=mask(pemap,(p3*86400.gt.2),True)
preciponly!0="lat"
preciponly&lat=pr1->lat({-30:30})
preciponly!1="lon"
preciponly&lon=pr1->lon



pe=avg(preciponly)
pe1=avg(preciponly)
print(avg(pe))
print(avg(pe1))
;;;;compute 



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
res3@cnLevels = (/1,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-20,-25,-30/)*1e-1
res3@cnFillColors = (/0,200,190,180,170,160,150,140,130,120,110,100,90,80,70,60,50,40,20,2,0/)
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
res3@gsnRightString="r~B~pat~N~: " + decimalPlaces(pcor,2,True)
res3@cnMissingValFillColor = "white"

cres=True
cres@cnLevelSelectionMode = "ExplicitLevels"
cres@cnLevels = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)/86400
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
plot(n)= gsn_csm_contour_map(wks,preciponly,res3)
plot1=gsn_csm_contour(wks,p3,cres)
;;;;attach zonal means;;;
  zres=True
  zres@trXMinF              = -3      ; Could also use gsnZonalMeanXMinF
  zres@trXMaxF              = 1   ; Could also use gsnZonalMeanXMaxF
  zres@trYMinF              = -30
  zres@trYMaxF              =  30

  zres@tmXBMode             = "Explicit"

  zres@tmXBValues           = (/-3,-2,-1,0,1/)
  zres@tmXBLabels           = zres@tmXBValues + ""   ; Make it a string
  zres@tmXBLabelFontHeightF = 0.01                   ; Make slightly smaller.
  zres@vpWidthF             = 0.1
 zonal_id = gsn_csm_attach_zonal_means(wks,plot(n),preciponly,zres)
overlay(plot(n),plot1)
n=n+1
 delete(T)
 delete(T1)
 delete(MDIFF)
 delete(wva)
 delete(wva1)
 delete(q1avg)
 delete(q2avg)
end do


res4 = True
res4@gsnPanelLabelBar=True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = "M~B~int~N~/M~S~'~N~"
res4@lbTitlePosition ="Bottom"
res4@lbTitleFontHeightF=.015
res4@txFontHeightF    = 0.03
res4@pmLabelBarOrthogonalPosF=-0.02
res4@gsnMaximize=True
res4@gsnPaperOrientation="landscape"

gsn_panel(wks,plot,(/4,1/),res4)
delete(wks)

string1="eps2eps " + plotfile +".epsi" + " " + plotfile +".eps"
system(string1) 
system("epstopdf " + plotfile +".eps") 
system("rm " + plotfile +".eps") 
system("rm " + plotfile +".epsi")
end
