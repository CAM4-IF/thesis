load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*4K*ANN*.nc2")
datalog=new(dimsizes(fils),double) ; create list to store output
plotfile="200hpachistar_response_ANN" ;;;PDF file name!
wks = gsn_open_wks("epsi", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods+1, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
lets=(/"a)","b)","c)","d)"/)
f4=addfile("chistar200_maxresponse_ANN.nc","r")
f5=addfile("chistar200_minresponse_ANN.nc","r")
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
  maxchi=f4->velopot
  minchi=f5->velopot
  
  
  
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
maxdiff=chidiff-(maxchi(0,0,:,:))
mindiff=chidiff-(minchi(0,0,:,:))

maxdiff!0="lat"
maxdiff&lat=u&lat
maxdiff!1="lon"
maxdiff&lon=u&lon

mindiff!0="lat"
mindiff&lat=u&lat
mindiff!1="lon"
mindiff&lon=u&lon



chistarmin=min(chistar({-30:30},:))
chistarmin1=min(chistar1({-30:30},:))
response=((chistarmin1-chistarmin)/(chistarmin))
  
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
   res@cnLevels = (/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12/)*0.5e6
   res@cnFillColors = (/40,50,60,70,80,85,90,95,100,105,110,120,0,130,140,150,160,170,180,190,200,210,220,230,240/)
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
  res@gsnLeftString    = lets(n)
  res@gsnCenterString    = f@case
  res@gsnCenterStringFontHeightF=0.018
  res@gsnLeftStringFontHeightF=0.018
  ;res@gsnRightString    = "~F33~D~F~(~F33~c~F~~S~*~N~) response: " + sprintf("%5.2f",response)
res@tmXBLabelFontHeightF=0.006
res@gsnRightStringFontHeightF= 0.018
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
 cres@cnLineThicknessF =2

 cres@cnLevelSpacingF=200.0
 cres@gsnContourZeroLineThicknessF  = 4.0
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelFontHeightF      = 0.004
 cres@cnLineLabelPlacementMode = "Constant"
 cres@gsnContourNegLineDashPattern =2
 cres@cnLineLabelInterval=1
 cres@gsnLeftString    = ""
 cres@tmXTOn="False"

 
 plot(n)  = gsn_csm_contour_map(wks, chidiff, res) 
 plot1=gsn_csm_contour(wks,chistar,cres)
 overlay(plot(n),plot1)
 
 ;;;Add stipling where response is larger than max response of AMIP4K members;;;
gt=maxdiff>0
gt!0="lat"
gt&lat=u&lat
gt!1="lon"
gt&lon=u&lon

lt=mindiff<0
lt!0="lat"
lt&lat=u&lat
lt!1="lon"
lt&lon=u&lon
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
  res30@cnFillScaleF        = 1        ; add extra density
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
delete(U)
delete(V)

delete(f)
delete(f1)

;;;;Add AMIP4K chi* plot;;;;
f=addfile("chi200_response_allmodels_ANN.nc","r") ; add file
f1=addfile("../CMIP5_AMIP/ANN_AMIP_chi200_climo.nc","r") ; add file
f2=addfile("chistarstip.nc","r") 
f3=addfile("chistarstip1.nc","r")

U=f->velopot
V=f1->velopot
stip=f2->velopot(0,0,:,:)
stip1=f3->velopot(0,0,:,:)
res@gsnLeftString    = lets(n)
res@gsnRightString    =""
cres@gsnRightString    = ""
chistar=dim_rmvmean_n_Wrap(U(:,0,:,:),2)
chistar!1="lat"
chistar&lat=U&lat
chistar!2="lon"
chistar&lon=U&lon
chistarmean=dim_avg_n(chistar,0)
chistarmean!0="lat"
chistarmean&lat=U&lat
chistarmean!1="lon"
chistarmean&lon=U&lon


chistar1=dim_rmvmean_n_Wrap(V(0,0,:,:),0)
chistar1!0="lat"
chistar1&lat=U&lat
chistar1!1="lon"
chistar1&lon=U&lon
 res@gsnCenterString    = "AMIP4K"
plot(n)  = gsn_csm_contour_map(wks,chistarmean(:,:), res) 
plot1=gsn_csm_contour(wks,chistar1,cres)
overlay(plot(n),plot1)

res30 = True                            
  res30@gsnDraw=False
  res30@gsnFrame=False
  res30@cnFillOn = True                   ; turn on color/pattern fill
  res30@cnMonoFillPattern = False   ; allow different patterns
  res30@cnMonoFillColor = True       ; only use one color (black)
  res30@cnLevelSelectionMode = "ExplicitLevels" ; set explicit cnlev
  res30@cnLevels   = (/7,8/)    ; 
  res30@cnFillPatterns = (/-1,17,17/) ;
  res30@gsnAddCyclic = True   ; may or may not be needed
  res30@cnInfoLabelOn       = False       ; turn off info label
  res30@cnLinesOn           = False       ; do not draw contour lines
  res30@cnLineLabelsOn      = False    ; do not draw contour labels
  res30@cnFillScaleF        = 1         ; add extra density
  res30@gsnRightString   = ""    
  res30@gsnLeftString    = ""
  res30@gsnCenterString  = ""
  res30@lbLabelBarOn=False
  
  plot2   = gsn_csm_contour(wks,stip, res30)
  res30@cnFillPatterns = (/-1,6,6/) ;
  plot3   = gsn_csm_contour(wks,stip1, res30)
 overlay(plot(n),plot2)
 overlay(plot(n),plot3)

print(n)

;res@gsnLeftString  = "chi* min: " + decimalPlaces(chistarmin2/1e6,2,True) + "x1e6 m^2/s"
;res@gsnCenterString=lets(n-1) + " MERRA/ERAI"
;res@gsnRightString =""
 
;plot(n)=gsn_csm_contour_map(wks, chistar2, res) 
;

resP                            = True
resP@gsnPanelLabelBar=True
resP@cnInfoLabelOn       = False
resP@lbOrientation =   "horizontal"          ; 
resP@lbLabelFontHeightF    = 0.006
resP@txString = ""
resP@lbTitleString    = "~F33~D~F~(~F33~c~F~~S~*~N~) m~S~2~N~/s"
resP@gsnMaximize=True
resP@gsnPaperOrientation = "landscape"
resP@lbTitlePosition ="Bottom"
resP@pmLabelBarOrthogonalPosF=-0.02
resP@lbTitleFontHeightF=.015
resP@txFontHeightF    = 0.03
gsn_panel(wks,plot,(/4,1/),resP)  
delete(wks) 
string1="eps2eps " + plotfile +".epsi" + " " + plotfile +".eps"
print(string1)
system(string1) 
system("epstopdf " + plotfile +".eps") 
system("rm " + plotfile +".eps") 
system("rm " + plotfile +".epsi")
print(fils)
print(fils1)
end
