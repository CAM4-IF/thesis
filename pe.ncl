begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*_ANN*.nc11")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*4K_ANN*.nc21")
fils2=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc4")
fils3=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc2")
plotfile="ANN_pe"
wks = gsn_open_wks("pdf", plotfile)                      ;Open a workstation
datalog=new((/3,4608/),"double")
datalog1=new((/3,4608/),"double")
datalog2=new((/3,4608/),"double")
datalog3=new((/3,4608/),"double")
datalog4=new((/3,4608/),"double")
datalog5=new((/3,4608/),"double")
modelnames=new(3,"string")
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods+2, graphic)    ;;;Create a panel plot with panels = number of models!
plots=new(9, graphic)
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

pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)
opt=(/0,10000,100000/)
wva2 = wgt_vert_avg_beta(pm, O, PS, 0, opt)/90000
wva3 = wgt_vert_avg_beta(pm, O1, PS, 0, opt)/90000
wva2!0="lat"
wva2&lat=pr1->lat
wva2!1="lon"
wva2&lon=pr1->lon
wva3!0="lat"
wva3&lat=pr1->lat
wva3!1="lon"
wva3&lon=pr1->lon


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


preciponly=-1/(mask(pemap,(p3*86400.gt.2),True))
preciponly!0="lat"
preciponly&lat=pr1->lat({-30:30})
preciponly!1="lon"
preciponly&lon=pr1->lon



preciponly1=-1/(mask(pemap1,(p6*86400.gt.2),True))
preciponly1!0="lat"
preciponly1&lat=pr1->lat({-30:30})
preciponly1!1="lon"
preciponly1&lon=pr1->lon

peresponse=(preciponly1-preciponly)*-1
peresponse!0="lat"
peresponse&lat=pr1->lat({-30:30})
peresponse!1="lon"
peresponse&lon=pr1->lon


pe=avg(preciponly)
pe1=avg(preciponly1)


mprimeponly=(mask(mprime,(p3*86400.gt.4),True))
mprimeponly!0="lat"
mprimeponly&lat=pr1->lat({-30:30})
mprimeponly!1="lon"
mprimeponly&lon=pr1->lon

mprimeponly1=(mask(mprime1,(p6*86400.gt.4),True))
mprimeponly1!0="lat"
mprimeponly1&lat=pr1->lat({-30:30})
mprimeponly1!1="lon"
mprimeponly1&lon=pr1->lon


mintponly=(mask(wva({-30:30},:),(p3*86400.gt.4),True))
mintponly!0="lat"
mintponly&lat=pr1->lat({-30:30})
mintponly!1="lon"
mintponly&lon=pr1->lon

mintponly1=(mask(wva1({-30:30},:),(p6*86400.gt.4),True))
mintponly1!0="lat"
mintponly1&lat=pr1->lat({-30:30})
mintponly1!1="lon"
mintponly1&lon=pr1->lon


oponly=(mask(wva2({-30:30},:),(p3*86400.gt.4),True))
oponly!0="lat"
oponly&lat=pr1->lat({-30:30})
oponly!1="lon"
oponly&lon=pr1->lon


oponly1=(mask(wva3({-30:30},:),(p3*86400.gt.4),True))
oponly1!0="lat"
oponly1&lat=pr1->lat({-30:30})
oponly1!1="lon"
oponly1&lon=pr1->lon





mprimeresponse=mprime1-mprime
mintresponse=wva1({-30:30},:)-wva({-30:30},:)

x=ndtooned(mprimeponly)
y=ndtooned(mintponly)*-1
x1=ndtooned(mprimeponly1)
y1=ndtooned(mintponly1)*-1
x2=ndtooned(oponly)*-1
x3=ndtooned(oponly1)*-1

datalog(n,:)=x
datalog1(n,:)=y
datalog2(n,:)=x1
datalog3(n,:)=y1
datalog4(n,:)=x2
datalog5(n,:)=x3

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
wva2 = wgt_vert_avg_beta(pm, O, PS, 1, opt)/900
wva3 = wgt_vert_avg_beta(pm, O1, PS1, 1, opt)/900
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

modelnames=(/"CAM4 default","CAM4-IF-r","CAM4-IF-t"/)

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
res3@cnLevels = (/1,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-20,-25,-30/)*-1e-1
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
res3@gsnRightString=""
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
  zres@trXMinF              = -1      ; Could also use gsnZonalMeanXMinF
  zres@trXMaxF              = 3  ; Could also use gsnZonalMeanXMaxF
  zres@trYMinF              = -30
  zres@trYMaxF              =  30

  zres@tmXBMode             = "Explicit"

  zres@tmXBValues           = (/-1,0,1,2,3/)
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
res4@lbTitleString    = "_~H-20~~V-28~M~B~int~N~/M~S~'~N~"
res4@lbTitlePosition ="Bottom"
res4@lbTitleFontHeightF=.015
res4@txFontHeightF    = 0.03
res4@pmLabelBarOrthogonalPosF=-0.02
res4@gsnMaximize=True
res4@gsnPaperOrientation="landscape"

gsn_panel(wks,plot,(/4,1/),res4)
delete(wks)
n=0
delete(lets)
lets=(/"a)","b)","c)"/)
wks = gsn_open_wks("pdf", "pe_scat") 
do while(n.le.dimsizes(fils)-1)
 lets=(/"a)","b)","c)"/)
 names=(/"CAM4 default","CAM4-IF-r","CAM4-IF-t"/)
 res                   = True                     ; plot mods desired
  res@gsnMaximize       = True                     ; maximize plot
  res@tiMainString      = ""           ; add title
  res@xyMarkLineMode    = "Markers"                ; choose to use markers
  res@xyMarkers         =  16                      ; choose type of marker  
  res@xyMarkerColor     = "Black"               ; Marker color
  res@xyMarkerSizeF     = 0.0005                     ; Marker size (default 0.01)
  res@gsnDraw=False
  res@gsnFrame=False
  res@trXMinF=0.001
  res@trXMaxF=0.02
  res@trYMinF=0.001
  res@gsnXRefLine = -3. 
  res@gsnXRefLine = 0.0 
  res@gsnYRefLine = 0.0  
  ;res@trYMaxF=0.022
  res@gsnLeftString  =lets(n) + " " + modelnames(n)
  res@gsnLeftStringFontHeightF=0.03
  res@tiYAxisString = "_~H-30~~V-28~-M~B~int~N~ (kg m~S~-2~N~ s~S~-1~N~)"
  res@tiXAxisString = "_~H-20~~V-28~M~S~'~N~ (kg m~S~-2~N~ s~S~-1~N~)"
  res@tmXBFormat="0@*+^sg"
  rc=regline(datalog(n,:)<0.008,datalog1(n,:)<0.02)
  rc1=regline(datalog2(n,:)<0.008,datalog3(n,:)<0.02)
  yreg=rc*datalog(n,:)+rc@yintercept
  yreg1=rc1*datalog2(n,:)+rc1@yintercept
  cor=escorc(datalog(n,:),datalog1(n,:))
  cor1=escorc(datalog2(n,:),datalog3(n,:))
  plots(n)  = gsn_csm_xy (wks,datalog(n,:),datalog1(n,:),res)                    ; create plot
  
  res@xyMarkerColor     = "Red" 
  plot1 = gsn_csm_xy (wks,datalog2(n,:),datalog3(n,:),res)                    ; create plot
  res@xyMarkLineMode="Lines"
  res@xyLineColor="Grey"
  plot2=gsn_csm_xy(wks,datalog(n,:),yreg,res)
  res@xyLineColor="Orange"
  plot3=gsn_csm_xy(wks,datalog2(n,:),yreg1,res)
  overlay(plots(n),plot1)
  overlay(plots(n),plot2)
  overlay(plots(n),plot3)
  txid1 = new(3,graphic)                  ; Nine text strings
 txres               = True     
 txres@txFontHeightF = 0.022
 txid1(0) = gsn_add_text (wks,plots(n),"Control r: " + sprintf("%5.2f",cor) + "~C~" + "+4K r: " + sprintf("%5.2f",cor1), 0.005,max(datalog1(n,:))-0.002, txres)
 gres=True
  gres@YPosPercent = 50. 
gres@XPosPercent = 75
gres@ItemSpacePercent = 5.  
lineres1 = True
lineres1@lgLineColors = (/"Grey","Orange"/) ; line colors
lineres1@lgDashIndexes=(/0,0,0,0,0,0/)
lineres1@lgLineThicknesses = 2.5                        ; line thicknesses
lineres1@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres1 = True
textres1@lgLabels = (/"Control","+4K" /)  ; legend labels (required)
textres1@lgLabelFontHeights = (/0.018,0.018,0.018/)   
plots(n) = simple_legend(wks,plots(n),gres,lineres1,textres1)

delete(lets)
lets=(/"d)","e)","f)"/)
names=(/"CAM4 default","CAM4-IF-r","CAM4-IF-t"/)
res@xyMarkLineMode    = "Markers" 
res@tiXAxisString = "_~H-27~~V-22~-~F33~w~F~~B~int~N~ (Pa s~S~-1~N~)"
res@trXMinF=-0.1
res@trXMaxF=0.15
res@trYMinF=-0.001
res@gsnLeftString  =lets(n) + " " + modelnames(n)
;res@trYMaxF=0.022
;res@tmXBFormat="0*+^e"
res@xyMarkerColor     = "Black" 
rc2=regline((datalog4(n,:)>0),datalog1(n,:))
rc3=regline((datalog5(n,:)>0),datalog3(n,:))
r    = escorc(datalog4(n,:),datalog1(n,:)) 
r1    = escorc(datalog5(n,:),datalog3(n,:)) 
yreg2=rc2*datalog4(n,:)+rc2@yintercept
yreg3=rc3*datalog5(n,:)+rc3@yintercept
plots(n+3)  = gsn_csm_xy (wks,datalog4(n,:),datalog1(n,:),res)                    ; create plot
res@xyMarkerColor     = "Red" 
plot4= gsn_csm_xy (wks,datalog5(n,:),datalog3(n,:),res)                    ; create plot
overlay(plots(n+3),plot4)
res@xyMarkLineMode="Lines"
res@xyLineColor="Grey"
res@gsnXRefLine = 0.0 
res@gsnYRefLine = 0.0  
plot5=gsn_csm_xy(wks,datalog4(n,:),yreg2,res)
res@xyLineColor="Orange"
plot6=gsn_csm_xy(wks,datalog5(n,:),yreg3,res)
overlay(plots(n+3),plot5)
overlay(plots(n+3),plot6)
cmeanmint=avg(datalog1(n,:) > 0.005)
cmeanomeg=avg(datalog4(n,:) > 0.01)
fmeanmint=avg(datalog3(n,:) > 0.005)
fmeanomeg=avg(datalog5(n,:) > 0.01)

mintresponse1=(datalog3(n,:)-datalog1(n,:))/datalog1(n,:)
omegresponse=(datalog5(n,:)-datalog4(n,:))/datalog4(n,:)


 ;txres               = True     
 ;txres@txFontHeightF = 0.018
 ;txid1(0) = gsn_add_text (wks,plots(n+3),"Mean M~B~int~N~ response (%): " + sprintf("%5.1f",mintresponse1) + "~C~~C~" + "Mean (~F33~w~F~~B~int~N~)~S~~F34~-~F~~N~ response (%): " + sprintf("%5.1f",omegresponse), 0,max(datalog1(n,:))-0.001, txres)
 gres=True
  gres@YPosPercent = 50. 
gres@XPosPercent = 75
gres@ItemSpacePercent = 5.  
lineres1 = True
lineres1@lgLineColors = (/"Grey","Orange"/) ; line colors
lineres1@lgDashIndexes=(/0,0,0,0,0,0/)
lineres1@lgLineThicknesses = 2.5                        ; line thicknesses
lineres1@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres1 = True
textres1@lgLabels = (/"Control","+4K" /)  ; legend labels (required)
textres1@lgLabelFontHeights = (/0.018,0.018,0.018/)   
plots(n+3) = simple_legend(wks,plots(n+3),gres,lineres1,textres1)
  
  
  
lets=(/"g)","h)","i)"/)
names=(/"CAM4 default","CAM4-IF-r","CAM4-IF-t"/)
res1=True
res1@tiMainString      = ""           ; add title
res1@xyMarkers         =  16                      ; choose type of marker   
res1@xyMarkerColor     = "Black"               ; Marker color
res1@xyMarkerSizeF     = 0.0005                     ; Marker size (default 0.01)
res1@xyMarkLineMode    = "Markers" 
res1@tiYAxisString = "_~H-34~~V-28~~F33~d~F~M~B~int~N~~N~/M~B~int~N~~H-38~~V+28~~N~_"
res1@tiXAxisString = "_~H-34~~V-22~~F33~d~F~~F33~w~F~~B~int~N~/~F33~w~F~~B~int~N~~H-38~~V+22~~N~_"
res1@gsnDraw=False
res1@gsnFrame=False
res1@trXMinF=-1
res1@trXMaxF=1
res1@trYMinF=-1
res1@trYMaxF=1
res1@gsnXRefLine = 0.0
res1@gsnYRefLine = 0.0  
res1@gsnMaximize=True
res1@gsnLeftString  =lets(n) + " " + modelnames(n)
 res1@gsnLeftStringFontHeightF=0.03
;res@tmXBFormat="0*+^e"
res1@xyMarkerColor     = "Black" 
x10=omegresponse<1
x20=x10>-1
y10=mintresponse1<1
y20=y10>-1
rc4=regline(x20,y20)
r1    = escorc(x20,y20) 
print(r1)
yreg4=rc4*x20+rc4@yintercept
plots(n+6)  = gsn_csm_xy (wks,x20,y20,res1)                    ; create plot
res1@xyMarkLineMode="Lines"
res1@xyLineColor="Grey"
plot7=gsn_csm_xy(wks,x20,yreg4,res1)
overlay(plots(n+6),plot7)
delete(mintresponse1)
delete(omegresponse)
txres               = True     
txres@txFontHeightF = 0.025
txid1(0) = gsn_add_text (wks,plots(n+6),"r: " + sprintf("%5.2f",r1), -0.5,0.5, txres)
  
n=n+1
delete(lets)
end do

resP                            = True
resP@txString = ""
resP@txFontHeightF    = 0.03
resP@gsnMaximize=True
resP@gsnPaperOrientation = "portrait"
gsn_panel(wks,plots,(/3,3/),resP)  
 delete(wks) 

;string1="eps2eps " + plotfile +".epsi" + " " + plotfile +".eps"
;system(string1) 
;system("epstopdf " + plotfile +".eps") 
;system("rm " + plotfile +".eps") 
;system("rm " + plotfile +".epsi")
system("pdfcrop pe_scat.pdf")
end
