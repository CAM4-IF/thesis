begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_AMIP/AMIP/ua/*ANN*.nc")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_AMIP/AMIP/va/*ANN*.nc")
fils2=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_AMIP/AMIP4K/ua/*ANN*.nc")
fils3=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_AMIP/AMIP4K/va/*ANN*.nc")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="CAM_ANN_WCchi" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   filen2=tostring(fils2(n))
   filen3=tostring(fils3(n))
   c="cdo merge " + filen + " " + filen1 + " uvtmp.nc"
   c1="cdo merge " + filen2 + " " + filen3 + " uvtmp1.nc"
   print(c)
   cmd=systemfunc(c)
   cmd1=systemfunc(c1)
   c2="cdo selname,velopot,stream -sp2gp -dv2ps -uv2dv -remapbil,n45 -sellevel,20000 -chname,ua,U -chname,va,V uvtmp.nc chi.nc"
   c3="cdo selname,velopot,stream -sp2gp -dv2ps -uv2dv -remapbil,n45 -sellevel,20000 -chname,ua,U -chname,va,V uvtmp1.nc chi1.nc"
   cmd2=systemfunc(c2)
   cmd3=systemfunc(c3)
   c4="cdo sub -selname,velopot chi.nc -enlarge,chi.nc -zonmean -selname,velopot chi.nc chi2.nc" ; zonally anomalous chi
   c5="cdo sub -selname,velopot chi1.nc -enlarge,chi1.nc -zonmean -selname,velopot chi1.nc chi3.nc" ; zonally anomalous chi
   cmd4=systemfunc(c4)
   cmd5=systemfunc(c5)

   f=addfile("chi2.nc", "r")
   f1=addfile("chi3.nc","r")
   f2=addfile(filen,"r")
   
   wc=f->velopot(0,0,:,:)
   wcstrength=min(wc({-20:20},:))
   wc1=f1->velopot(0,0,:,:)
   wcstrength1=min(wc1({-20:20},:))
   
   wcresponse=wc1-wc
   wcrp=100*((wcstrength1-wcstrength)/wcstrength)
   
   wcresponse!0="lat"
   wcresponse&lat=wc&lat
   wcresponse!1="lon"
   wcresponse&lon=wc&lon
   print("hi!")
   
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
 
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20/)*0.5e6
  res@cnFillColors = (/2,6,12,24,40,60,90,100,110,0,0,150,156,170,190,210,220,230,240,250/)
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@lbLabelBarOn = False
  res@lbLabelAutoStride    = True                 ; optimal labels
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                ; turn off contour info label
  ;res@tiXAxisString    = ""
  ;res@tiYAxisSide = "Left"	
  ;res@tiYAxisString    = "Pressure (hPa)"
  ;res@trYMaxF = 1000.0
  ;res@trYMinF = 100.0      ; troposphere only
  ;res@tiXAxisFontHeightF =0.012
  ;res@tmYROn = False
  ;res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7
  res@cnLineThicknessF = 1.                   ; thicker lines
  res@cnLineLabelsOn   = False                
  res@cnInfoLabelOn = False                   ; turn on contour info label
  res@gsnSpreadColors     = False
  res@gsnMaximize      = True
  res@mpDataBaseVersion    = "MediumRes"  
  res@mpMaxLatF = 30                    ; specify the plot domain
  res@mpMinLatF = -30
  res@mpCenterLonF         = -200
  res@tmXTOn="False"

 
   ;control=dim_avg_n_Wrap(wva4,1)
   ;case=dim_avg_n_Wrap(wva9,1)
   ;response=case-control
   ;percentresponse=(response/control)*100.
   ;datalog(n)=percentresponse
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
  cres                      = True                 ; plot mods desired
  cres@cnLevelSelectionMode = "AutomaticLevels"
  cres@cnLineLabelsOn       = True                 ; turn on line labels
  cres@cnFillOn             = False                 ; turn on color fill
  cres@gsnDraw = False
  cres@gsnFrame = False
  cres@cnInfoLabelOn = False                ; turn off contour info label
  cres@cnLineThicknessF =1
  cres@gsnRightString    = f2@model_id
  cres@gsnRightStringFontHeightF    = 0.05
  cres@cnLevelSpacingF=200.0
  cres@cnLineLabelDensityF      = 3.0
  cres@cnLineLabelFontHeightF      = 0.004
  cres@cnLineLabelPlacementMode = "Constant"
  cres@cnLineLabelInterval=1
  cres@gsnLeftString  ="chi* response (%): " + wcrp + "~C~control chi*: " + wcstrength
  cres@gsnLeftStringFontHeightF=0.04
  cres@tmXTOn="False"
  cres1=cres
  cres1@cnLineColor="red"
  cres1@gsnContourPosLineDashPattern=1
 
 
  plot(n)  = gsn_csm_contour_map(wks, wcresponse, res)  
  ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
  plot2  = gsn_csm_contour(wks, wc, cres)
  ;plot3=gsn_csm_vector(wks,udiff,vdiff,vcres)
  overlay(plot(n),plot2)
  ;overlay(plot(n),plot3)
 
   
  cleanup=systemfunc("rm *tmp*.nc *chi*.nc")
  
  datalog(n)=wcstrength
   
n=n+1
end do
 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "ANN AMIP4K-AMIP_control 200 hPa Chi* response (shaded) and control climatology (contours)"
 resP@lbTitleString    = "m^2 s^-1"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.0085
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/5,2/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop CAM_ANN_WCchi.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(datalog)
print(fils)
print(fils1)
end
