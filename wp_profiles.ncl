begin
;************************************************
; read in data
;************************************************FM
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc2")
print(fils)
print(fils1)
plotfile="wp_profiles_response" ;;;PDF file name!
datalog=new((/3,26/),"double")
datalog1=new((/3,26/),"double")
datalog2=new((/3,26/),"double")
datalog3=new((/3,26/),"double")
datalog4=new((/3,26/),"double")
datalog5=new((/3,26/),"double")
vertmean=new((/3/),"double")
vertmean1=new((/3/),"double")
vertmean2=new((/3/),"double")
vertmean3=new((/3/),"double")
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(2, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
lets=(/"a)","b)","c)"/)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM4_4K_ANN_climo.nc2", "r")
   hyam = h->hyam 
   hybm = h->hybm
   cmd="cdo -fldmean -sellonlatbox,0,360,-15,15 -selname,T,PS,CLDICE,CLDLIQ,QRL,QRS " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,0,360,-15,15 -selname,T,PS,CLDICE,CLDLIQ,QRL,QRS " + filen1 + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->T(0,:,:,:)
   T&lev@units="hPa"
   lw=f->CLDICE(0,:,:,:)
   sw=f->CLDLIQ(0,:,:,:)
   l=f->QRL(0,:,:,:)
   s=f->QRS(0,:,:,:)
   
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   ;levs=conform_dims(dimsizes(T),T&lev,0)
  
   PS1   = f1->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->T(0,:,:,:)
   T1&lev@units="hPa"
   lw1=f1->CLDICE(0,:,:,:)
   sw1=f1->CLDLIQ(0,:,:,:)
   l1=f1->QRL(0,:,:,:)
   s1=f1->QRS(0,:,:,:)
   
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   
   Qdiff=lw1-lw
   Qdiff1=sw1-sw
   Qraddiff=(l1+s1)-(l+s)
   datalog(n,:)=Qdiff(:,0,0)*1000
   datalog1(n,:)=sw(:,0,0)*1000
   datalog2(n,:)=lw(:,0,0)*1000
   datalog3(n,:)=Qdiff1(:,0,0)*1000
   datalog4(n,:)=Qraddiff(:,0,0)
   datalog5(n,:)=l(:,0,0)+s(:,0,0)
   
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(pm/100, Qdiff, PS/100, 0, opt)/900
   wva1 = wgt_vert_avg_beta(pm1/100, Qdiff1, PS1/100, 0, opt)/900
   wva2 = wgt_vert_avg_beta(pm/100, lw, PS/100, 0, opt)/900
   wva3 = wgt_vert_avg_beta(pm1/100,sw, PS1/100, 0, opt)/900
   wva4 = wgt_vert_avg_beta(pm/100, l+s, PS/100, 0, opt)/900
   wva5 = wgt_vert_avg_beta(pm1/100,Qraddiff, PS1/100, 0, opt)/900
   print(wva)
   iceres=wva(0,0)
   watres=wva1(0,0)
   ice=wva2(0,0)
   wat=wva3(0,0)
   vertmean(n)=(iceres/ice)*100
   vertmean1(n)=(watres/wat)*100
   vertmean2(n)=(wva5/wva4)*100
   
   
   print(vertmean(n))
   print(vertmean1(n))
   print(vertmean2(n))
   print(fils(n))
   
   

 res                   = True                       ; plot mods desired
 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 1.0
 res@tiXAxisString=""
 res@tiXAxisFontHeightF=0.02
 res@trYMaxF=1000
 res@trYMinF=30
 res@trXMinF=-10*1e-3
 res@trXMaxF=50*1e-3
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.02
 res@xyYStyle = "Linear"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmXBMode = "Explicit"
 res@tmXBValues        =(/-2,0,2,4,6,8,10,12/)*4e-3
 res@tmXBLabels=(/-2,0,2,4,6,8,10,12/)*4e-3
 res@tmXBLabelFontHeightF=0.015
 
 res@tiXAxisString="q~B~liq/ice~N~ (g/kg)"
 res@vpHeightF 	 = 0.7               ; change aspect ratio of plot
 res@vpWidthF 	         = 0.55
 res@gsnXRefLine = 0.0 
 res@gsnLeftString  ="a) Cloud ice/water mixing ratio"
 res@gsnLeftStringFontHeightF=0.03
 res@gsnRightStringFontHeightF=0.025
 res@gsnDraw=False
 res@gsnFrame=False
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
 
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
 

n=n+1
end do      
lets=(/"a)","b)","c)"/)
res@xyDashPattern  = 1
res@xyLineColor    = "blue"
 res@xyLineThicknesses = 1.0
 plot(0)= gsn_csm_xy  (wks,datalog(0,:),T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 plot2= gsn_csm_xy  (wks,datalog(1,:),T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 plot3= gsn_csm_xy  (wks,datalog(2,:),T1&lev,res)
 res@xyDashPattern  = 0
 res@xyLineColor    = "blue"
 res@xyLineThicknesses = 3.0
 plot4= gsn_csm_xy  (wks,datalog2(0,:),T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 3.0
 plot5= gsn_csm_xy  (wks,datalog2(1,:),T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 3.0
 plot6= gsn_csm_xy  (wks,datalog2(2,:),T1&lev,res)
 
 
 res@xyDashPattern  = 0
 res@xyLineColor    = "cyan"
 res@xyLineThicknesses = 3.0
 plot7= gsn_csm_xy  (wks,datalog1(0,:),T1&lev,res)
 res@xyLineColor    = "purple"
 res@xyLineThicknesses = 3.0
 plot8= gsn_csm_xy  (wks,datalog1(1,:),T1&lev,res)
 res@xyLineColor    = "orange"
 res@xyLineThicknesses = 3.0
 plot9= gsn_csm_xy  (wks,datalog1(2,:),T1&lev,res)
 
 
 res@xyDashPattern  = 1
 res@xyLineColor    = "cyan"
 res@xyLineThicknesses = 1.0
 plot10= gsn_csm_xy  (wks,datalog3(0,:),T1&lev,res)
 res@xyLineColor    = "purple"
 res@xyLineThicknesses = 1.0
 plot11= gsn_csm_xy  (wks,datalog3(1,:),T1&lev,res)
 res@xyLineColor    = "orange"
 res@xyLineThicknesses = 1.0
 plot12= gsn_csm_xy  (wks,datalog3(2,:),T1&lev,res)
 
 
 overlay(plot(0),plot2)
 overlay(plot(0),plot3)
 overlay(plot(0),plot4)
 overlay(plot(0),plot5)
 overlay(plot(0),plot6)
 overlay(plot(0),plot7)
 overlay(plot(0),plot8)
 overlay(plot(0),plot9)
 overlay(plot(0),plot10)
 overlay(plot(0),plot11)
 overlay(plot(0),plot12)
 gres = True
gres@YPosPercent = 85. 
gres@XPosPercent = 47
gres@ItemSpacePercent = 5.  
lineres = True
lineres@lgLineColors = (/"black","red","blue","orange","purple","cyan"/) ; line colors
lineres@lgDashIndexes=(/0,0,0,0,0,0/)
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres = True
textres@lgLabels = (/"CAM4-IF best T q~B~ice","CAM4-IF best rain q~B~ice","Default CAM4 q~B~ice","CAM4-IF best T q~B~liq","CAM4-IF best rain q~B~liq","Default CAM4 q~B~liq" /)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015/)   
plot3 = simple_legend(wks,plot3,gres,lineres,textres)




 
;;;;;Radiative heating;;;;;
 res1=res
 res1@xyLineThicknesses = 1.0
 res1@tiXAxisString="Q~B~rad~N~ (K/s)"
 res1@trXMinF=-20e-6
 res1@trXMaxF=16e-6
 res1@gsnLeftString  ="b) Radiative heating"
 res1@tmXBLabelFontHeightF=0.012
 delete(res1@tmXBValues)
 delete(res1@tmXBLabels)
 res1@tmXBValues        =(/-24,-20,-16,-12,-8,-4,0,4,8,12,16/)*1e-6
 res1@tmXBLabels=(/-24,-20,-16,-12,-8,-4,0,4,8,12,16/)*1e-6
 
 
 res1@xyDashPattern  = 0
  res1@xyLineColor    = "blue"
   res1@xyLineThicknesses = 3.0
 plot(1)= gsn_csm_xy  (wks,datalog5(0,:),T1&lev,res1)
 res1@xyLineColor    = "red"
 res1@xyLineThicknesses = 3.0
 plot16= gsn_csm_xy  (wks,datalog5(1,:),T1&lev,res1)
 res1@xyLineColor    = "black"
 res1@xyLineThicknesses = 3.0
 plot17= gsn_csm_xy  (wks,datalog5(2,:),T1&lev,res1)
 res1@xyDashPattern  = 1
 res1@xyLineThicknesses = 1.0
 res1@xyLineColor    = "blue"
 plot13= gsn_csm_xy  (wks,datalog4(0,:),T1&lev,res1)
 res1@xyLineColor    = "red"
 res1@xyLineThicknesses = 1.0
 plot14= gsn_csm_xy  (wks,datalog4(1,:),T1&lev,res1)
 res1@xyLineColor    = "black"
 res1@xyLineThicknesses = 1.0
 plot15= gsn_csm_xy  (wks,datalog4(2,:),T1&lev,res1)
 overlay(plot(1),plot13)
 overlay(plot(1),plot14)
 overlay(plot(1),plot15)
 overlay(plot(1),plot16)
 overlay(plot(1),plot17)
 txid1 = new(3,graphic)                  ; Nine text strings
 txres               = True     
 txres@txFontHeightF = 0.015
 txid1(0) = gsn_add_text (wks,plot(1),"CAM4: " + decimalPlaces(vertmean2(0),1,True) + " %", 4.45e-6, 400,txres)
 txid1(1) = gsn_add_text (wks,plot(1),"CAM4-IF best rain: " + decimalPlaces(vertmean2(1),1,True) + " %", 8e-6, 430,txres)
 txid1(2) = gsn_add_text (wks,plot(1),"CAM4-IF best T: " + decimalPlaces(vertmean2(2),1,True) + " %", 7.3e-6, 460,txres)
 gres@YPosPercent = 85. 
gres@XPosPercent = 57
gres@ItemSpacePercent = 5.  
lineres1 = True
lineres1@lgLineColors = (/"black","red","blue"/) ; line colors
lineres1@lgDashIndexes=(/0,0,0,0,0,0/)
lineres1@lgLineThicknesses = 2.5                        ; line thicknesses
lineres1@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres1 = True
textres1@lgLabels = (/"CAM4-IF best T","CAM4-IF best rain","Default CAM4" /)  ; legend labels (required)
textres1@lgLabelFontHeights = (/0.015,0.015,0.015/)   
plot(1) = simple_legend(wks,plot(1),gres,lineres1,textres1)

  
 
res4 = True
res4@cnInfoLabelOn       = False
res4@lbOrientation =   "horizontal"          ; 
res4@gsnPanelMainString    = " "
res4@lbLabelFontHeightF    = 0.006
res4@lbTitleString    = ""
res4@lbTitlePosition ="Bottom"
res4@lbTitleFontHeightF=.008
res4@txFontHeightF    = 0.03
res4@gsnMaximize=True


gsn_panel(wks,plot,(/1,2/),res4)
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
cleanup=systemfunc("rm tmp*.nc")
end
