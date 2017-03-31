begin

fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_AMIP/AMIP/ta/ANN*")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_AMIP/AMIP4K/ta/ANN*")
datalog=new(dimsizes(fils),double) ; create list to store output
plotfile="AMIP4K_Sp_ANN" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
 f1=addfile(fils(n),"r")
 T1=f1->ta(0,:,{-10:10},:)
 lons=f1->lon
 s   = static_stability(T1&plev, T1, 0, 0)
 s!0="lev"
 s&lev=T1&plev
 szmean=dim_avg_n(s,1)
 szmean!0="lev"
 szmean&lev=T1&plev
 szmean!1="lon"
 szmean&lon=lons
 
 
 f2=addfile(fils1(n),"r")
 T2=f2->ta(0,:,{-10:10},:)
 lons1=f2->lon
 s2   = static_stability(T2&plev, T2, 0, 0)
 s2!0="lev"
 s2&lev=T2&plev
 szmean1=dim_avg_n(s2,1)
 szmean1!0="lev"
 szmean1&lev=T2&plev
 szmean1!1="lon"
 szmean1&lon=lons1
 
 sdiff=szmean1-szmean
 
 sdiff!0="lev"
 sdiff&lev=T2&plev
 sdiff!1="lon"
 sdiff&lon=lons1
 
 filen=tostring(fils1(n))
 cmd="cdo -fldmean -sellonlatbox,60,200,-10,10 " + filen + " tmp.nc"
 areamean=systemfunc(cmd)
 f=addfile("tmp.nc","r") ; add file
 T=f->ta(0,:,0,0) ; extract T data from file
 ;;;calculate sp and set dimension attributes;;;
 s1   = static_stability(T&plev, T, 0, 0)
 s1!0="lev"
 s1&lev=T&plev
 
 s1lev=T&plev
 top=(ind(s1lev.eq.15000))
 bot=(ind(s1lev.eq.92500))
 print(top)
 dp   = dpres_plevel(T&plev, 101300., 3000., 0) 

 fint = dim_sum_n(s1(top:bot)*dp(top:bot), 0)/dim_sum_n(dp(top:bot),0)
 printVarSummary(fint)
 datalog(n)=fint
 
 ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-30,-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20,30/)*1e-5
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
  res@trYMinF = 100.0      ; troposphere only
  res@tiXAxisFontHeightF =0.012
  res@tmYROn = False
  res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7

 
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels = (/-260,-240,-220,-200,-180,160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240,260/)*1e-5
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
 cres@gsnRightString    = f@model_id
 cres@gsnLeftString=datalog(n)
 cres@trYMinF            = 100.   ; Limits for Y axis.  The limits
 cres@trYMaxF            = 1000.
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelPlacementMode = "Constant"
 cres@cnLineLabelInterval=1


 
 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 plot(n)  = gsn_csm_pres_hgt(wks, sdiff, res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_pres_hgt(wks, szmean, cres)
 overlay(plot(n),plot2)
 
 
 n=n+1
 delete(T)
 delete(s1)
 delete(dp)
 delete(top)
 delete(bot)
 delete(s1lev)
 delete(areamean)
 delete(s)
 delete(T1)
 delete(szmean)
 delete(lons)
 delete(lons1)
 delete(sdiff)
 delete(s2)
 delete(szmean1)
 delete(T2)
 cleanup=systemfunc("rm tmp.nc")
end do
;;;;build the panel plot!;;;;
 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "AMIP4K 10S-10NN mean ANN Sp response (shaded) with control run climatology (contours) (1980-2005)"
 resP@lbTitleString    = "K/Pa"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.01
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/5,5/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop AMIP4K_Sp_ANN.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(datalog)
print(fils)
end
