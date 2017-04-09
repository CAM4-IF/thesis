load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc1")
datalog=new(dimsizes(fils),double) ; create list to store output
plotfile="GMS_CAM4_ANN" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1+1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
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
 ; psfc1 = f1->PS(:,{-10:10},:)   ; surface pressure [Pa]
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

chistarmin=min(chistar({-20:20},:))
chistarmin1=min(chistar1({-20:20},:))
print(chistarmin)
print(chistarmin1)
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



;;;Calculate meridional wind * moist entropy;;;
;mev=V(0,:,:,:)*te
;mev!0="lev"
;mev&lev=U&lev
;mev!1="lat"
;mev&lat=U&lat
;mev!2="lon"
;mev&lon=U&lon

 
;  vp1=vinth2p(V(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
  
;mzmpsi=zonal_mpsi(vp,mev&lat,lev_p*100,psfc(0,:,:))
;mzmpsi!0="lev"
;mzmpsi&lev=lev_p
;mzmpsi!1="lat"
;mzmpsi&lat=U&lat

;zmpsi=zonal_mpsi(vp1,mev&lat,lev_p*100,psfc(0,:,:))
;zmpsi!0="lev"
;zmpsi&lev=lev_p
;zmpsi!1="lat"
;zmpsi&lat=U&lat

;print(max(mzmpsi(:,{-20:20})))
;print(max(zmpsi(:,{-20:20})))
;GMSHC=(max(mzmpsi(:,{-20:20})))/(max(zmpsi(:,{-20:20})))
;print(GMSHC)


;;;;calculate moisture flux convergence;;;
;;advection term (neglect in deep tropics);;;
;advQ=advect_var(U(0,:,:,:),V(0,:,:,:),te,1,"advection of Q","1/m",0)
;advQ!0="lev"
;advQ!1="lat"
;advQ&lat=U&lat
;advQ!2="lon"
;advQ&lon=U&lon
;advQ_trop=advQ(:,{-10:10},:)
;;divergence term;;
;;first set up empty arrays for divergence scalar and vector components;;
;dv=new((dimsizes(U)),float,1e+36)
;ud=new((dimsizes(U)),float,1e+36)
;vd=new((dimsizes(U)),float,1e+36)
;uv2dvf(U,V,dv) ;;calculates divergence scalar
;dv2uvf(dv,ud,vd) ;;calculates divergent wind components
;dv!1="lev"
;dv&lev=lev_p
;dv!2="lat"
;dv&lat=U&lat
;dv!3="lon"
;dv&lon=U&lon
;mdiv=Q(0,:,{-10:10},:)*dv(0,:,{-10:10},:) ;;scale divergence by Q (just deep tropics only)
;;;total Mdiv;;
;mfc=mdiv
;mfc!0="lev"
;mfc&lev=lev_p
;mfc!1="lat"
;troplat=f->lat({-10:10})
;mfc&lat=troplat
;mfc!2="lon"
;mfc&lon=U&lon
;;;horizontal divergence of MSE (ignore in deep tropics);;;
;msediv=te(:,{-10:10},:)*dv(0,:,{-10:10},:)
;;;do vertical integral of numerator and denominator, then do spatial averaging;;;
;;vertically integrate then average denominator (925-100 hPa);;;
opt=(/0,90,1000/)
;wvam = wgt_vert_avg_beta(pm(0,:,{-10:10},{:}),mfc(:,{-10:10},{:}),psfc(0,{-10:10},{:}),0,opt)/910
;mfcmean=avg(wvam)
;;;vertically integrate numerator (vertical advection of MSE);;;
  wva = wgt_vert_avg_beta(lev_p,(dTdP*ww),psfc(0,:,:),0,opt)/910
  wva!0="lat"
  wva&lat=u&lat
  wva!1="lon"
  wva&lon=u&lon
  mask1=f->LANDFRAC(0,:,:)
  wvamask=mask(wva,mask1,0)
  wvamask!0="lat"
  wvamask&lat=u&lat
  wvamask!1="lon"
  wvamask&lon=u&lon
  GMS1=avg(wva({-10:10},{:})) ;;numerator (MSE transport)
;;;;Calculate mass flux in this region (use omega_500 as an estimate);;;;
  opt1=(/0,150,850/)
  mflux=wgt_vert_avg_beta(pm(0,:,:,:),w(0,:,:,:),psfc(0,:,:),0,opt1)/700
  mflux!0="lat"
  mflux&lat=u&lat
  mflux!1="lon"
  mflux&lon=u&lon
 ; mflux=w(0,18,:,:)
  mfluxmask=mask(mflux,mask1,0)
  mfluxmask!0="lat"
  mfluxmask&lat=u&lat
  mfluxmask!1="lon"
  mfluxmask&lon=u&lon
  mfluxavg=avg(mflux({-10:10},{:}))
  print(mfluxavg)
  
  
  GMS2=9.81*GMS1/mfluxavg
  print(GMS1)
  print(GMS2)
  print(fils(n))
  print(fils1(n))
  
  
   res                      = True                
   res@tmYRMode = "Automatic"
   res@cnLevelSelectionMode = "ExplicitLevels"   
   ;res@cnLevels = (/0,5,10,15,20,25,30,35,40/)*3e-4
   res@cnLevels = (/-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20/)*1e6
   ;res@cnLevels =(/-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18/)*3e-5
   res@cnFillColors = (/2,6,12,24,40,60,90,100,0,0,0,0,156,170,190,210,220,230,240,250/)
   ;res@cnFillColors = (/0,166,176,186,196,206,216,226,236,246/)
   res@cnLinesOn       = False               
   res@cnLineLabelsOn       = False                 
  res@cnFillOn             = True                
  res@lbLabelBarOn = False
  res@lbLabelAutoStride    = True                 
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                
  res@cnFillOpacityF=0.7
  res@cnLineThicknessF = 1.                 
  res@cnLineLabelsOn   = False                
  res@cnInfoLabelOn = False                 
  res@gsnSpreadColors     = False
  res@gsnMaximize      = True
  res@mpDataBaseVersion    = "MediumRes"  
  res@mpMaxLatF =35; specify the plot domain
  res@mpMinLatF = -35
  res@mpCenterLonF         = -200
  res@tmXTOn="False"
  res@gsnZonalMeanYRefLine=0
  res@gsnRightString    = f@case
  res@gsnLeftString    = ""
  res@lbTitleFontHeightF=.007   
  res@gsnRightString    = f@case
  ;res@gsnLeftString  = "chi* min: " + chistarmin + "chi* response (%): " + response + "~C~GMS: " + GMS2
  res@gsnLeftString  = "chi* min: " + chistarmin
  
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
 ;cres@gsnLeftString  = "10S-10N response pattern correlation: " + pcor + "~C~10S-10N tropical mean Q:" + qstarmean
 ;cres@gsnLeftString  = "10S-10N response pattern correlation: " + pcor
 cres@gsnLeftStringFontHeightF=0.025
 cres@gsnLeftString    = ""
 cres@tmXTOn="False"
; cres@gsnZonalMean    = True
 
 plot(n)  = gsn_csm_contour_map(wks, chistar, res) 
 ;plot1=gsn_csm_contour(wks,chistar,cres)
 ;overlay(plot(n),plot1)
  n=n+1
  delete(u)
 ; delete(dp)
  delete(u1)
  delete(lons)
  delete(te)
  delete(dTdP)
  delete(wva)
  delete(mflux)
  delete(chistar)
  delete(chistar1)
  delete(chidiff)
  
end do



f2=addfile("MERRA_ANN_UV200_climo_remap.nc","r")

U2=f1->U
V2=f1->V

uvmsg=1e+36
vpot=new((/96,144/),float,uvmsg)
sf=new((/96,144/),float,uvmsg)

uv2sfvpf(U2(0,{200},:,:),V2(0,{200},:,:),sf,vpot)

chistar2=dim_rmvmean_n_Wrap(vpot,0)

chistar2!0="lat"
chistar2&lat=U2&lat
chistar2!1="lon"
chistar2&lon=U2&lon
chistarmin2=min(chistar2({-20:20},:))

print(n)
print(dimsizes(plot))
;plot(n)=gsn_csm_contour_map(wks, chistar2, res) 


resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "ANN 200 hPa chi* (contours) and response (shading)"
 resP@lbTitleString    = "m^2/s"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.0085
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/4,2/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop " + plotfile)  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
end
