load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc1")
datalog=new(dimsizes(fils),double) ; create list to store output
plotfile="200hpachistar_ANN" ;;;PDF file name!
wks = gsn_open_wks("epsi", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods+1, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
lets=(/"a)","b)","c)","d)"/)
n=0
do while(n.le.dimsizes(fils)-1)
  f=addfile(fils(n),"r") ; add file
  f1=addfile(fils1(n),"r") ; add file
  u     = f->T(:,:,:,:)
  u@units="K"
  q     = f->Q(:,:,:,:)
  u1     = f1->T(:,:,:,:)
  w=f->OMEGA(:,:,:,:)
  w1=f1->OMEGA(:,:,:,:)
  u1@units="K"
  q1     = f1->Q(:,:,:,:)
  ;q    = q/(1-q) ;
  q@long_name = "mixing ratio"
  ;q1    = q1/(1-q1) ;
  q1@long_name = "mixing ratio"
  U=f->U
  V=f->V
  U1=f1->U
  V1=f1->V
  Q=f->Q
  Q1=f1->Q
  Z=f->Z3
  lons=f1->lon
  
  
  
  hyam = f1->hyam ; read from a file the mid-layer coef
  hybm = f1->hybm ; read from a file
  psfc = f->PS(:,:,:)   ; surface pressure [Pa]
  psfc1=f->PS(:,:,:)
  P0   = 100000. ; since ps is in Pa or [ f->P0]
  pm  = pres_hybrid_ccm(psfc,P0,hyam,hybm)/100.
  pm!1="lev"
  pm!2="lat"
  pm&lat=u&lat
  pm!3="lon"
  pm&lon=u&lon
  
  
  
  hyam1 = f1->hyam ; read from a file the mid-layer coef
  hybm1 = f1->hybm ; read from a file
 ; psfc1 = f1->PS(:,{-30:30},:)   ; surface pressure [Pa]
  P01   = 100000. ; since ps is in Pa or [ f->P0]
  pm1  = pres_hybrid_ccm(psfc1,P01,hyam1,hybm1)/1000.
  
  
  lev_p           = (/  10, 20, 30, 50, 70,100,150,200,250,300,350,400,500,600,700,800,875,925,950,975,1000 /)                
  lev_p!0         = "lev_p"                  ; variable/dim name 
  lev_p&lev_p     =  lev_p                   ; create coordinate variable
  lev_p@long_name = "pressure"               ; attach some attributes
  lev_p@units     = "hPa"
  lev_p@positive  = "down"

; convert to pressure levels
 vp   = vinth2p(V(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
 UU=vinth2p(U(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
 vp1   = vinth2p(V1(0,:,:,:),hyam,hybm,lev_p,psfc1(0,:,:),2,1000.,1,False)
 UU1=vinth2p(U1(0,:,:,:),hyam,hybm,lev_p,psfc1(0,:,:),2,1000.,1,False)
 up=vinth2p(u(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
 ZZ=vinth2p(Z(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
 qq=vinth2p(q(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
 ww=vinth2p(w(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
 
 wwstar=dim_rmvmean_n_Wrap(ww, 1)
 
;;;calculate 200 hPa Chi;;;;
nlat=dimsizes(vp&lat)
mlon=dimsizes(up&lat)
uvmsg=1e+36
vpot=new((/96,144/),float,uvmsg)
sf=new((/96,144/),float,uvmsg)
vpot1=new((/96,144/),float,uvmsg)
sf1=new((/96,144/),float,uvmsg)

uv2sfvpf(UU({200},:,:),vp({200},:,:),sf,vpot)
uv2sfvpf(UU1({200},:,:),vp1({200},:,:),sf,vpot1)

chistar=dim_rmvmean_n_Wrap(vpot,0)
chistar1=dim_rmvmean_n_Wrap(vpot1,0)

chistar!0="lat"
chistar&lat=u&lat
chistar!1="lon"
chistar&lon=u&lon

chistar1!0="lat"
chistar1&lat=u&lat
chistar1!1="lon"
chistar1&lon=u&lon
chidiff=chistar1-chistar

chidiff!0="lat"
chidiff&lat=u&lat
chidiff!1="lon"
chidiff&lon=u&lon

chistarmin=min(chistar({-30:30},:))
chistarmin1=min(chistar1({-30:30},:))
response=100*((chistarmin1-chistarmin)/(chistarmin))
  
te=new((/dimsizes(lev_p),dimsizes(u&lat),dimsizes(u&lon)/), "double")
dTdP=new((/dimsizes(lev_p),dimsizes(u&lat),dimsizes(u&lon)/), "double")

i=0
do while(i.le.dimsizes(u&lon)-1)
 k=0
  do while(k.le.dimsizes(u&lat)-1)
   j=0
   do while(j.le.dimsizes(lev_p)-1)
    te(j,k,i)=1004.*up(j,k,i)+9.81*ZZ(j,k,i)+2.5e6*qq(j,k,i)
    j=j+1 
   end do
   j=0
   do while(j.le.dimsizes(lev_p)-2)
    meanT=(up(j,k,i)+up(j+1,k,i))/2
    dP=lev_p(j)-lev_p(j+1)
    dT=te(j,k,i)-te(j+1,k,i)
    dTdP(j,k,i)=dT/dP
    j=j+1
   end do
  k=k+1
 end do
 i=i+1
end do


dTdP!0="lev"
dTdP&lev=lev_p
dTdP&lev@units="hPa"
dTdP!1="lat"
dTdP&lat=u&lat
dTdP!2="lon"
dTdP&lon=u&lon

te!0="lev"
te&lev=lev_p
te&lev@units="hPa"
te!1="lat"
te&lat=u&lat
te!2="lon"
te&lon=u&lon

  
  f2=addfile("../MERGED_UV_ANN.nc","r")

U2=f2->U
V2=f2->V

uvmsg=1e+36
vpot=new((/96,144/),float,uvmsg)
sf=new((/96,144/),float,uvmsg)

uv2sfvpf(U2(0,{200},:,:),V2(0,{200},:,:),sf,vpot)

chistar2=dim_rmvmean_n_Wrap(vpot,0)
chistar2!0="lat"
chistar2&lat=U2&lat
chistar2!1="lon"
chistar2&lon=U2&lon
chistarmin2=min(chistar2({-30:30},:))

  
  pcor=pattern_cor(chistar({-30:30},:),chistar2({-30:30},:),1.0,0)
  rmse=(wgt_arearmse(chistar({-30:30},:),chistar2({-30:30},:),1.0,1.0,0))
  
  
   res                      = True        
   res@cnLevelSelectionMode = "ExplicitLevels"   
   res@cnLevels = (/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)*1e6
   res@cnFillColors = (/40,50,60,70,80,85,90,95,100,105,110,120,0,0,130,140,150,160,170,180,190,200,210,220,230,240/)
   res@cnLinesOn       = False               
   res@cnLineLabelsOn       = False                 
  res@cnFillOn             = True                
  res@lbLabelBarOn = False                 
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                
  res@cnLineLabelsOn   = False                
  res@cnInfoLabelOn = False                 
  res@gsnSpreadColors     = False
  res@gsnMaximize      = True
  res@mpDataBaseVersion    = "MediumRes"  
  res@mpMaxLatF =30; specify the plot domain
  res@mpMinLatF = -30
  res@mpCenterLonF         = -200
  res@tmXTOn="False"
  res@gsnZonalMeanYRefLine=0
  res@mpGridLatSpacingF = 15
  res@mpGridLonSpacingF = 45
  res@pmTickMarkDisplayMode  = "Always" 
  res@gsnLeftString    = "" 
  res@gsnCenterString    =lets(n) + " " + f@case
  res@gsnCenterStringFontHeightF=0.012
  ;res@gsnLeftString  = "chi* min: " + decimalPlaces(chistarmin/1e6,2,True) + "chi* response (%): " + decimalPlaces(response,2,True) + "~C~GMS: " + GMS2
  res@gsnLeftString  = "chi* min: " + decimalPlaces(chistarmin/1e6,2,True) + "x1e6 m^2/s" + "~C~Pattern correlation: " + decimalPlaces(pcor,2,True)
  res@gsnRightString="RMSE: " + decimalPlaces(rmse/1e6,2,True) + "x1e6 m^2/s"
  
res@tmXBLabelFontHeightF=0.006
res@gsnRightStringFontHeightF= 0.01
res@gsnLeftStringFontHeightF= 0.01
res@mpGridLineThicknessF=0.5
  
    ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels = (/-3.6,-3.4,-3.2,-3.,-2.8,-2.6,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6/)*5e6
 ;cres@cnLevels=(/-100,-50,-25,-10,0,10,25,50,100/)
 cres@cnLineLabelsOn       = True                 ; turn on line labels
 cres@cnFillOn             = False                 ; turn on color fill
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False                ; turn off contour info label
 cres@cnLineThicknessF =0.5
 cres@gsnRightString    = ""
 cres@gsnRightStringFontHeightF    = 0.05
 cres@cnLevelSpacingF=200.0
 cres@gsnContourZeroLineThicknessF  = 2.0
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelFontHeightF      = 0.004
 cres@cnLineLabelPlacementMode = "Constant"
 cres@gsnContourNegLineDashPattern =2
 cres@cnLineLabelInterval=1
 cres@gsnLeftStringFontHeightF=0.025
 cres@gsnLeftString    = ""
 cres@tmXTOn="False"
; cres@gsnZonalMean    = True


 
 plot(n)  = gsn_csm_contour_map(wks, chistar, res) 
 plot1=gsn_csm_contour(wks,chistar2,cres)
 overlay(plot(n),plot1)
  n=n+1
  delete(u)
 ; delete(dp)
  delete(u1)
  delete(lons)
  delete(te)
  delete(dTdP)
  delete(chistar)
  delete(chistar1)
  delete(chidiff)
  
end do




print(n)

;res@gsnLeftString  = "chi* min: " + decimalPlaces(chistarmin2/1e6,2,True) + "x1e6 m^2/s"
;res@gsnCenterString=lets(n-1) + " MERRA/ERAI"
;res@gsnRightString =""
 
;plot(n)=gsn_csm_contour_map(wks, chistar2, res) 


resP                            = True
resP@gsnPanelLabelBar=True
resP@cnInfoLabelOn       = False
resP@lbOrientation =   "horizontal"          ; 
resP@lbLabelFontHeightF    = 0.006
resP@txString = "ANN"
resP@lbTitleString    = "m^2/s"
resP@gsnMaximize=True
resP@gsnPaperOrientation = "landscape"
resP@lbTitlePosition ="Bottom"
resP@lbTitleFontHeightF=.008
resP@txFontHeightF    = 0.03
gsn_panel(wks,plot,(/4,1/),resP)  
delete(wks) 
string1="eps2eps " + plotfile +".epsi" + " " + plotfile +".eps"
print(string1)
system(string1) 
system("epstopdf " + plotfile +".eps") 
system("rm " + plotfile +".eps") 
system("rm " + plotfile +".epsi")
end
