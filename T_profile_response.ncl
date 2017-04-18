begin
;************************************************
; read in data
;************************************************FM
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc2")
print(fils)
print(fils1)
plotfile="T_profile_response" ;;;PDF file name!
datalog=new((/3,26/),"double")
datalog1=new((/3,26/),"double")
datalog2=new((/3,26/),"double")
datalog3=new((/3,26/),"double")
vertmean=new((/3/),"double")
vertmean1=new((/3/),"double")
vertmean2=new((/3/),"double")
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(4, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
lets=(/"a)","b)","c)","d)"/)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM4_4K_ANN_climo.nc2", "r")
   hyam = h->hyam 
   hybm = h->hybm
   cmd="cdo -fldmean -sellonlatbox,0,360,-15,15 -selname,T,PS,DTCOND,QRL,QRS,DTV,OMEGA " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,0,360,-15,15 -selname,T,PS,DTCOND,QRL,QRS,DTV,OMEGA " + filen1 + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->T(0,:,:,:)
   T&lev@units="hPa"
   cond=f->DTCOND(0,:,:,:)
   lw=f->QRL(0,:,:,:)
   sw=f->QRS(0,:,:,:)
   diff=f->DTV(0,:,:,:)
   omeg=f->OMEGA(0,:,:,:) 
   
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   levs=conform_dims(dimsizes(T),T&lev,0)
  
   
   PS1   = f1->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->T(0,:,:,:)
   T1&lev@units="hPa"
   cond1=f1->DTCOND(0,:,:,:)
   lw1=f1->QRL(0,:,:,:)
   sw1=f1->QRS(0,:,:,:)
   diff1=f1->DTV(0,:,:,:)
   omeg1=f1->OMEGA(0,:,:,:)
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   
   Q=cond+lw+sw+diff
   Q1=cond1+lw1+sw1+diff1
   Qdiff=Q1-Q
   datalog1(n,:)=Qdiff(:,0,0)
   
  
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, Q, PS/100, 0, opt)/950
   Qtot=wva(0,0)
   
   opt1=(/0,100,1000/)
   wva3 = wgt_vert_avg_beta(pm1/100, Qdiff, PS/100, 0, opt)/950
   print(wva3)
   Qtot1=wva3(0,0)
   vertmean(n)=(Qtot1/Qtot)*100
  
   s=static_stability(pm,T,0,0)
   s1=static_stability(pm1,T1,0,0)
   sdiff=s1-s
   opt=(/0,150,1000/)
   wva3 = wgt_vert_avg_beta(levs, s, PS/100, 0, opt)/850
   wva4 = wgt_vert_avg_beta(levs, sdiff, PS/100, 0, opt)/850
   datalog2(n,:)=sdiff(:,0,0)
   vertmean1(n)=(wva4/wva3)*100
   
   omegdiff=omeg1-omeg
   datalog3(n,:)=omegdiff(:,0,0)
   opt=(/0,100,1000/)
   print(omegdiff)
   wva5 = wgt_vert_avg_beta(levs, omeg, PS/100, 0, opt)/900
   wva6 = wgt_vert_avg_beta(levs, omegdiff, PS/100, 0, opt)/900
   vertmean2(n)=(wva6/wva5)*100
   
   

 res                   = True                       ; plot mods desired
 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 1.0
 res@tiXAxisString=""
 res@tiXAxisFontHeightF=0.02
 res@trYMaxF=1000
 res@trYMinF=30
 ;res@trXMinF=min(T1(:,0,0))-0.0001
 res@trXMinF=-1
 ;res@trXMaxF=0.002
 res@trXMaxF=10
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.02
 res@xyYStyle = "Linear"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmXBMode = "Explicit"
 res@tmXBValues        =(/-2,-1,0,1,2,3,4,5,6,7,8,9,10/)
 res@tmXBLabels=(/-2,-1,0,1,2,3,4,5,6,7,8,9,10/)
 res@tmXBLabelFontHeightF=0.015
 
 ;res@tiXAxisString="Total M~B~c~N~ (kg/m~S~2~N~/s)"
 res@tiXAxisString="~F33~D~F~T (K)"
 res@vpHeightF 	 = 0.6               ; change aspect ratio of plot
 res@vpWidthF 	         = 0.55
 res@gsnXRefLine = 0.0 

 res@tiXAxisStringFontHeightF=0.0015
 res@gsnLeftString  ="a) Temperature response"
 res@gsnLeftStringFontHeightF=0.03
 res@gsnRightStringFontHeightF=0.025
 res@gsnDraw=False
 res@gsnFrame=False
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
 T4=T1(:,0,0)-T(:,0,0)
 
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 res@xyDashPattern  = 1
 datalog(n,:)=T4
 

n=n+1
end do      

res@xyDashPattern  = 0
res@xyLineColor    = "blue"
 res@xyLineThicknesses = 1.0
 plot(0)= gsn_csm_xy  (wks,datalog(0,:),T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 plot2= gsn_csm_xy  (wks,datalog(1,:),T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 plot3= gsn_csm_xy  (wks,datalog(2,:),T1&lev,res)
 overlay(plot(0),plot2)
 overlay(plot(0),plot3)
 gres = True
gres@YPosPercent = 25. 
gres@XPosPercent = 20 
gres@ItemSpacePercent = 5.  
lineres = True
lineres@lgLineColors = (/"black","red","blue"/) ; line colors
lineres@lgDashIndexes=(/0,0,0/)
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres = True
textres@lgLabels = (/"CAM4-IF best T","CAM4-IF best rain","Default CAM4"/)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015/)   
plot3 = simple_legend(wks,plot3,gres,lineres,textres)
 delta=datalog(1,:)-datalog(0,:)
 delta=datalog(2,:)-datalog(0,:)
 delta1=datalog(2,:)-datalog(0,:)
 delta2=datalog(1,:)-datalog(2,:)
 
 
 res@gsnLeftString  ="b) Diabatic heating response"
 delete(res@trXMinF)
 delete(res@trXMaxF)
 res@trXMinF=-8e-6
 res@trXMaxF=8e-6
 delete(res@tmXBValues)
 delete(res@tmXBLabels)
 delete(res@tiXAxisString)
 res@tiXAxisString="~F33~D~F~Q~B~tot~N~ (K/s)"
 res@tmXBValues        =(/-8,-6,-4,-2,0,2,4,6,8/)*1e-6
 res@tmXBLabels=(/-8,-6,-4,-2,0,2,4,6,8/)*1e-6
 res@tmXBLabelFontHeightF=0.015
 res@xyDashPattern  = 0
 res@xyLineColor    = "blue"
 
 res@xyLineThicknesses = 1.0
 plot(1)= gsn_csm_xy  (wks,datalog1(0,:),T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 1.0
 plot4= gsn_csm_xy  (wks,datalog1(1,:),T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 1.0
 plot5= gsn_csm_xy  (wks,datalog1(2,:),T1&lev,res)
 overlay(plot(1),plot4)
 overlay(plot(1),plot5)
 txid = new(3,graphic)                  ; Nine text strings
 txres               = True     
 txres@txFontHeightF = 0.015
 txid(0) = gsn_add_text (wks,plot(1),"CAM4: " + decimalPlaces(vertmean(0),1,True) + " %", -5.5e-6, 225,txres)
 txid(1) = gsn_add_text (wks,plot(1),"CAM4-IF best rain: " + decimalPlaces(vertmean(1),1,True) + " %",-4.05e-6, 255,txres)
 txid(2) = gsn_add_text (wks,plot(1),"CAM4-IF best T: " + decimalPlaces(vertmean(2),1,True) + " %", -4.35e-6, 285,txres)
 
 
;;;;;Static stability;;;;;
res1=res
 res1@xyLineThicknesses = 1.0
 res1@tiXAxisString="~F33~D~F~~F33~s~F~ (K/Pa)"
 res1@trXMinF=-8e-4
 res1@trXMaxF=8e-4
 res1@gsnLeftString  ="c) Static stability response"
 res1@tmXBLabelFontHeightF=0.012
 delete(res1@tmXBValues)
 delete(res1@tmXBLabels)
 res1@tmXBValues        =(/-10,-8,-6,-4,-2,0,2,4,6,8,10/)*1e-4
 res1@tmXBLabels=(/-10,-8,-6,-4,-2,0,2,4,6,8,10/)*1e-4
 res1@xyLineColor    = "blue"
 plot(2)= gsn_csm_xy  (wks,datalog2(0,:),T1&lev,res1)
 res1@xyLineColor    = "red"
 res1@xyLineThicknesses = 1.0
 plot6= gsn_csm_xy  (wks,datalog2(1,:),T1&lev,res1)
 res1@xyLineColor    = "black"
 res1@xyLineThicknesses = 1.0
 plot7= gsn_csm_xy  (wks,datalog2(2,:),T1&lev,res1)
 overlay(plot(2),plot6)
 overlay(plot(2),plot7)
 txid1 = new(3,graphic)                  ; Nine text strings
 txres               = True     
 txres@txFontHeightF = 0.015
 txid1(0) = gsn_add_text (wks,plot(2),"CAM4: " + decimalPlaces(vertmean1(0),1,True) + " %", 2.75e-4, 125,txres)
 txid1(1) = gsn_add_text (wks,plot(2),"CAM4-IF best rain: " + decimalPlaces(vertmean1(1),1,True) + " %", 4.35e-4, 155,txres)
 txid1(2) = gsn_add_text (wks,plot(2),"CAM4-IF best T: " + decimalPlaces(vertmean1(2),1,True) + " %", 4.05e-4, 185,txres)
 
 
 ;;;;;Omega;;;;;
res2=res
 res2@xyLineThicknesses = 1.0
 res2@tiXAxisString="~F33~D~F~~F33~w~F~ (Pa/s)"
 res2@trXMinF=-5e-3
 res2@trXMaxF=5e-3
 res2@gsnLeftString  ="d) Omega response"
 res2@tmXBLabelFontHeightF=0.012
 delete(res2@tmXBValues)
 delete(res2@tmXBLabels)
 res2@tmXBValues        =(/-10,-8,-6,-4,-2,0,2,4,6,8,10/)*1e-3
 res2@tmXBLabels=(/-10,-8,-6,-4,-2,0,2,4,6,8,10/)*1e-3
 res2@xyLineColor    = "blue"
 plot(3)= gsn_csm_xy  (wks,datalog3(0,:),T1&lev,res2)
 res2@xyLineColor    = "red"
 res2@xyLineThicknesses = 1.0
 plot8= gsn_csm_xy  (wks,datalog3(1,:),T1&lev,res2)
 res2@xyLineColor    = "black"
 res2@xyLineThicknesses = 1.0
 plot9= gsn_csm_xy  (wks,datalog3(2,:),T1&lev,res2)
 overlay(plot(3),plot8)
 overlay(plot(3),plot9)
 txid1 = new(3,graphic)                  ; Nine text strings
 txres               = True     
 txres@txFontHeightF = 0.015
 txid1(0) = gsn_add_text (wks,plot(3),"CAM4: " + decimalPlaces(vertmean2(0),1,True) + " %", 0.0018, 125,txres)
 txid1(1) = gsn_add_text (wks,plot(3),"CAM4-IF best rain: " + decimalPlaces(vertmean2(1),1,True) + " %", 0.0028, 155,txres)
 txid1(2) = gsn_add_text (wks,plot(3),"CAM4-IF best T: " + decimalPlaces(vertmean2(2),1,True) + " %", 0.0026, 185,txres)
 
 
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


gsn_panel(wks,plot,(/2,2/),res4)
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
cleanup=systemfunc("rm tmp*.nc")
end
