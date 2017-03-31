begin

fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*DJF*.nc1")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*DJF*.nc2")
datalog=new(dimsizes(fils),double) ; create list to store output
plotfile="v4_spdotomegbar" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
  f=addfile(fils(n),"r") ; add file
  f1=addfile(fils1(n),"r") ; add file
  u     = f->T(:,:,{-15:5},:)
  u@units="K"
  omeg= f->OMEGA(:,:,{-15:5},:)
  q     = f->Q(:,:,{-15:5},:)
  u1     = f1->T(:,:,{-15:5},:)
  u1@units="K"
  omeg1=f1->OMEGA(:,:,{-15:5},:)
  q1     = f1->Q(:,:,{-15:5},:)
  q    = q/(1-q) ;
  q@long_name = "mixing ratio"
  q1    = q1/(1-q1) ;
  q1@long_name = "mixing ratio"
  lons=f1->lon
  
  hyam = f->hyam ; read from a file the mid-layer coef
  hybm = f->hybm ; read from a file
  psfc = f->PS(:,{-15:5},:)   ; surface pressure [Pa]
  P0   = 100000. ; since ps is in Pa or [ f->P0]
  pm  = pres_hybrid_ccm(psfc,P0,hyam,hybm)
  
  
  
  hyam1 = f1->hyam ; read from a file the mid-layer coef
  hybm1 = f1->hybm ; read from a file
  psfc1 = f1->PS(:,{-15:5},:)   ; surface pressure [Pa]
  P01   = 100000. ; since ps is in Pa or [ f->P0]
  pm1  = pres_hybrid_ccm(psfc1,P01,hyam1,hybm1)
  

  omeg&lev@units="hPa"
  omeg1&lev@units="hPa"
  
  s1=static_stability(pm, u, 1, 0)
  s2=static_stability(pm1,u1, 1, 0)
  
  spo=s1*omeg
  spo1=s2*omeg1
  
  omegdiff=omeg1-omeg
  sdiff=spo1(0,:,:,:)-spo(0,:,:,:)
  sdiff1=dim_avg_n(sdiff,1)
  sdiff1!0="lev"
  sdiff1&lev=u&lev
  sdiff1!1="lon"
  sdiff1&lon=u&lon
  sdiff1&lev@units="hPa"
  
  sbaromeg=s1(0,:,:,:)*omegdiff(0,:,:,:)
  sbaromegavg=dim_avg_n(sbaromeg,1)
  sbaromegavg!0="lev"
  sbaromegavg&lev=u&lev
  sbaromegavg!1="lon"
  sbaromegavg&lon=u&lon
  sbaromegavg&lev@units="hPa"
  
  sdotomegbar=sdiff*omeg(0,:,:,:)
  sdotomegbaravg=dim_avg_n(sdotomegbar,1)
  sdotomegbaravg!0="lev"
  sdotomegbaravg&lev=u&lev
  sdotomegbaravg!1="lon"
  sdotomegbaravg&lon=u&lon
  sdotomegbaravg&lev@units="hPa"
  
  
  opt=(/1,8000,92500/)
 wva5 = wgt_vert_avg_beta(pm(0,:,:,:),sdotomegbar, psfc(0,:,:), 1, opt)
 wva5avg=avg(wva5)
 
  
  
  szmean=dim_avg_n(spo(0,:,:,:),1)
   szmean!0="lev"
   szmean&lev=u&lev
   szmean!1="lon"
   szmean&lon=u&lon
   szmean&lev@units="hPa"
   
 opt=(/1,8000,92500/)
 wva = wgt_vert_avg_beta(pm, spo, psfc, 1, opt)
 wva1 = wgt_vert_avg_beta(pm, spo1, psfc1, 1, opt)
 wvaavg=avg(wva)
 wvaavg1=avg(wva1)
   response=100*(wvaavg1-wvaavg)/(wvaavg)
   print(response)
   
  
  
 
 filen=tostring(fils(n))
 filen1=tostring(fils1(n))
 cmd="cdo -fldmean -sellonlatbox,0,360,-15,5 -selname,T,Q,PS,OMEGA " + filen + " tmp.nc"
 areamean=systemfunc(cmd)
 cmd1="cdo -fldmean -sellonlatbox,0,360,-15,5 -selname,T,Q,PS,OMEGA " + filen1 + " tmp1.nc"
 areamean1=systemfunc(cmd1)
 f2=addfile("tmp.nc","r") ; add file
 T=f2->T(0,:,:,:) ; extract T data from file
 Q=f2->Q(0,:,:,:) ; extract T data from file
 Q    = Q/(1-Q) ;
 Q@long_name = "mixing ratio"

 f3=addfile("tmp1.nc","r") ; add file
 T1=f3->T(0,:,:,:) ; extract T data from file
 Q1=f3->Q(0,:,:,:) ; extract T data from file
 Q1    = Q1/(1-Q1) ;
 Q1@long_name = "mixing ratio"
 omega=f2->OMEGA(0,:,:,:)
 omega1=f3->OMEGA(0,:,:,:)
 PS=f2->PS(0,:,:)
 PS1=f3->PS(0,:,:)
 pm2   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa]
 s3   = static_stability(pm2, T, 0, 0)
 s3!0="lev"
 s3&lev=T&lev
 s3lev=T&lev
 
 pm3   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
 s4   = static_stability(pm3, T1, 0, 0)
 s4!0="lev"
 s4&lev=T&lev
 s4lev=T&lev
 
 spomeg=s3*omega
 spomeg1=s4*omega1
 
 
 print(pm2)
 opt=(/1,8000,92500/)
 wva2 = wgt_vert_avg_beta(pm2, spomeg, PS, 1, opt)
 fint=avg(wva2)
 
 opt=(/1,8000,92500/)
 wva3 = wgt_vert_avg_beta(pm3, spomeg1, PS1, 1, opt)
 fint1=avg(wva3)
 
 print(fint)
 print(fint1)
 
 datalog(n)=fint
 percentresponse=((fint1-fint)/fint)*100
  
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-30,-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20,30/)*2e-7
  res@cnFillColors = (/2,6,12,18,24,40,60,90,106,116,128,130,142,148,156,168,180,208,218,230,250,255/) 
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@lbLabelBarOn = False
  res@lbLabelAutoStride    = True                 ; optimal labels
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                ; turn off contour info label
  res@tiXAxisString    = ""
  res@tiYAxisSide = "Left"	
  res@tiYAxisString    = "Pressure (hPa)"
  res@trYMaxF = 1000.0
  res@trYMinF = 50.0      ; troposphere only
  res@tiXAxisFontHeightF =0.012
  res@tmYROn = False
  res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7

 
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels = (/-260,-240,-220,-200,-180,160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240,260/)*5e-7
 cres@cnLineLabelsOn       = True                 ; turn on line labels
 cres@cnFillOn             = False                 ; turn on color fill
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False                ; turn off contour info label
 cres@gsnRightString = ""
 cres@gsnLeftString = ""
 cres@tiYAxisString    = "Pressure (hPa)"
 cres@tiXAxisString    = "Longitude"
 cres@tmYROn = False
 cres@tmYRLabelsOn = False
 cres@gsnContourZeroLineThicknessF  = 0.0
 cres@cnLineThicknessF =0.5
 cres@gsnContourNegLineDashPattern =1  
 cres@gsnLeftString =""
 cres@gsnContourZeroLineThicknessF = 2
 cres@gsnRightString    = f@case
 cres@trYMinF            = 50.   ; Limits for Y axis.  The limits
 cres@trYMaxF            = 1000.
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelPlacementMode = "Constant"
 cres@cnLineLabelInterval=1
 cres@gsnLeftString  = "925-80 hPa mean:" + decimalPlaces(wva5avg,8,True) 
 ;+ "~C~% response:" + decimalPlaces(response,4,True)
 ;cres@gsnLeftString  = tediff
 cres@gsnLeftStringFontHeightF=0.017

 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 plot(n)  = gsn_csm_pres_hgt(wks, sdotomegbaravg(:,:), res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_pres_hgt(wks, szmean(:,:), cres)
 overlay(plot(n),plot2)
 
  n=n+1
  delete(u)
  delete(s1)
 ; delete(dp)
  delete(u1)
  delete(szmean)
  delete(lons)
  delete(sdiff)
  delete(s2)
  delete(T)
end do
resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "CAM4/CAM-IF DJF (1999-2002) 15S-5N AMIP4K sp_dot*omega_bar (shaded) and sp*omega climatology (contours)"
 resP@lbTitleString    = "K/Pa"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride=True
 resP@txFontHeightF    = 0.0085
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/3,3/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop" + plotfile)  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(datalog)
end
