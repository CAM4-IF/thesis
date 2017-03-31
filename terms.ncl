begin

fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc2")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*ANN*.nc1")
datalog=new(dimsizes(fils),double) ; create list to store output
plotfile="test_sdob" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
  f=addfile(fils(n),"r") ; add file
  f1=addfile(fils1(n),"r") ; add file
  u     = f->T(:,:,{-10:10},:)
  u@units="K"
  omeg= f->OMEGA(:,:,{-10:10},:)
  q     = f->Q(:,:,{-10:10},:)
  u1     = f1->T(:,:,{-10:10},:)
  u1@units="K"
  omeg1=f1->OMEGA(:,:,{-10:10},:)
  q1     = f1->Q(:,:,{-10:10},:)
  q    = q/(1-q) ;
  q@long_name = "mixing ratio"
  q1    = q1/(1-q1) ;
  q1@long_name = "mixing ratio"
  lons=f1->lon
  
  hyam = f->hyam ; read from a file the mid-layer coef
  hybm = f->hybm ; read from a file
  psfc = f->PS(:,{-10:10},:)   ; surface pressure [Pa]
  P0   = 100000. ; since ps is in Pa or [ f->P0]
  pm  = pres_hybrid_ccm(psfc,P0,hyam,hybm)
  
  
  
  hyam1 = f1->hyam ; read from a file the mid-layer coef
  hybm1 = f1->hybm ; read from a file
  psfc1 = f1->PS(:,{-10:10},:)   ; surface pressure [Pa]
  P01   = 100000. ; since ps is in Pa or [ f->P0]
  pm1  = pres_hybrid_ccm(psfc1,P01,hyam1,hybm1)
  

  opt=(/0,15000,92500/)

  
  s1   = static_stability(pm, u, 1, 0)
  s2   = static_stability(pm1, u1, 1, 0)
  s1!1="lev"
  s1&lev=pm(0,:,0,0)
  s1&lat=u&lat
  s1lev=pm(0,:,0,0)
  s2!1="lev"
  s2&lev=pm1(0,:,0,0)
  s2&lat=u1&lat
  s2lev=pm1(0,:,0,0)
  s1s2=s2-s1
  
  spdotomegbar=-1*(s2-s1)*omeg
  sdob=dim_avg_n(spdotomegbar,2)
  sdob!1="lev"
  sdob&lev=u&lev
  sdob!2="lon"
  sdob&lon=lons
  sdob&lev@units="hPa"
  omegdotspbar=-1*(omeg1-omeg)*s1
  obsd=dim_avg_n(omegdotspbar,2)
  obsd!1="lev"
  obsd&lev=u&lev
  obsd!2="lon"
  obsd&lon=lons
  obsd&lev@units="hPa"



  
  top=(ind(s1lev.eq.13877.97))
  bot=(ind(s1lev.eq.92493.44))
  


  szmean=dim_avg_n(s1,2)
  szmean!1="lev"
  szmean&lev=u&lev
  szmean!2="lon"
  szmean&lon=lons
  szmean1=dim_avg_n(s2,2)
  szmean1!1="lev"
  szmean1&lev=u1&lev
  szmean1!2="lon"
  szmean1&lon=lons
  szmean&lev@units="hPa"
  sdiff=szmean1(0,:,:)-szmean(0,:,:)
  sdiff!0="lev"
  sdiff&lev=u&lev
  sdiff&lev@units="hPa"
  sdiff!1="lon"
  sdiff&lon=lons
  printVarSummary(sdiff)
  
  
 filen=tostring(fils(n))
 filen1=tostring(fils1(n))
 cmd="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,T,Q " + filen + " tmp.nc"
 areamean=systemfunc(cmd)
 cmd1="cdo -fldmean -sellonlatbox,0,360,-10,10 -selname,T,Q " + filen1 + " tmp1.nc"
 areamean1=systemfunc(cmd1)
 f2=addfile("tmp.nc","r") ; add file
 T=f2->T(0,:,0,0) ; extract T data from file
 Q=f2->Q(0,:,0,0) ; extract T data from file
 Q    = Q/(1-Q) ;
 Q@long_name = "mixing ratio"

 f3=addfile("tmp1.nc","r") ; add file
 T1=f3->T(0,:,0,0) ; extract T data from file
 Q1=f3->Q(0,:,0,0) ; extract T data from file
 Q1    = Q1/(1-Q1) ;
 Q1@long_name = "mixing ratio"
 
 
 s3   = static_stability(T&lev, T, 0, 0)
 s3!0="lev"
 s3&lev=T&lev
 s3lev=T&lev


 s4   = static_stability(T1&lev, T1, 0, 0)
 s4!0="lev"
 s4&lev=T&lev
 s4lev=T&lev
 dp   = dpres_plevel(T&lev, 1013., 30., 0) 
 fint = dim_sum_n(s3(13:20)*dp(13:20), 0)/dim_sum_n(dp(13:20),0)
 fint1 = dim_sum_n(s4(13:20)*dp(13:20), 0)/dim_sum_n(dp(13:20),0)

 sdob1=dim_avg_n_Wrap(sdob,2)
 obsd1=dim_avg_n_Wrap(obsd,2)
 term1 = -1*dim_sum_n(sdob1(0,13:20)*dp(13:20), 0)/dim_sum_n(dp(13:20),0)
 term2 = -1*dim_sum_n(obsd1(0,13:20)*dp(13:20), 0)/dim_sum_n(dp(13:20),0)
 print(term1)
 printVarSummary(fint)
 printVarSummary(fint1)
 datalog(n)=fint
 percentresponse=((fint1-fint)/fint)*100
  
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-30,-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20,30/)*2e-6
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
 cres@cnLevels = (/-260,-240,-220,-200,-180,160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240,260/)*2e-7
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
 cres@gsnLeftString  = "925-150 hPa mean term 1: " + decimalPlaces(term1,8,True) + "~C~925-150 hPa mean term 2: " + decimalPlaces(term2,8,True)
 ;cres@gsnLeftString  = tediff
 cres@gsnLeftStringFontHeightF=0.012

 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 plot(n)  = gsn_csm_pres_hgt(wks, sdob(0,:,:), res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_pres_hgt(wks, obsd(0,:,:), cres)
 overlay(plot(n),plot2)
 
  n=n+1
  delete(u)
  delete(s1)
 ; delete(dp)
  delete(top)
  delete(bot)
  delete(s1lev)
  delete(u1)
  delete(szmean)
  delete(lons)
  delete(sdiff)
  delete(s2)
  delete(szmean1)
  delete(T)
end do
resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "CAM4/CAM-IF ANN (1998-2005) 10S-10N AMIP4K energy balance term 1 (shaded) and term 2 (contours)"
 resP@lbTitleString    = "K/s"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.0085
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/3,3/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop test_sdob.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(datalog)
end
