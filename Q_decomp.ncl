begin
;************************************************
; read in data
;************************************************FM
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/results/*CAM*ANN*.nc")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/results/*CAM*4K*ANN*.nc1")
print(fils)
print(fils1)
plotfile="Q_response_decomp" ;;;PDF file name!
numfils=dimsizes(fils)
datalog=new((/numfils,26/),"double")
datalog1=new((/numfils,26/),"double")
datalog2=new((/numfils,26/),"double")
datalog3=new((/numfils,26/),"double")
datalog4=new((/numfils,26/),"double")
vertmean=new((/numfils/),"double")
vertmean1=new((/numfils/),"double")
vertmean2=new((/numfils/),"double")
vertmean3=new((/numfils/),"double")
vertmean4=new((/numfils/),"double")
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(5, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
lets=(/"a)","b)","c)","d)"/)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM4_4K_ANN_climo.nc2", "r")
   hyam = h->hyam 
   hybm = h->hybm
   zon="cdo -s -zonmean -selname,PS,DTCOND,QRL,QRS,DTV " + filen + " zon.nc"
   zon1="cdo -s -zonmean -selname,PS,DTCOND,QRL,QRS,DTV " + filen1 + " zon1.nc"
   cmd="cdo -s -fldmean -sellonlatbox,90,180,-10,10 -sub -selname,PS,DTCOND,QRL,QRS,DTV " + filen + " -selname,PS,DTCOND,QRL,QRS,DTV -enlarge," + filen + " zon.nc"+ " tmp.nc"
   cmd1="cdo -s -fldmean -sellonlatbox,90,180,-10,10 -sub -selname,PS,DTCOND,QRL,QRS,DTV " + filen1 + " -selname,PS,DTCOND,QRL,QRS,DTV -enlarge," + filen + " zon.nc"+ " tmp1.nc"
   zonmean=systemfunc(zon)
   zonmean1=systemfunc(zon1)
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
  
   P0   = 100000.                        ; reference pressure [Pa]
   ;PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->DTCOND(0,:,:,:)
   T&lev@units="hPa"
   cond=f->DTCOND(0,:,:,:)
   lw=f->QRL(0,:,:,:)
   sw=f->QRS(0,:,:,:)
   diff=f->DTV(0,:,:,:)
   
   
   PS=new((/1,1/),float)
   PS(0,:)=1000.
   PS(:,0)=1000.
   PS1=new((/1,1/),float)
   PS1(0,:)=1000.
   PS1(:,0)=1000.
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   levs=conform_dims(dimsizes(T),T&lev,0)
  
   
   ;PS1   = f1->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->DTCOND(0,:,:,:)
   T1&lev@units="hPa"
   cond1=f1->DTCOND(0,:,:,:)
   lw1=f1->QRL(0,:,:,:)
   sw1=f1->QRS(0,:,:,:)
   diff1=f1->DTV(0,:,:,:)
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   
   Q=sw
   Q1=sw1
   Qdiff=Q1-Q
   
   
   conddiff=cond1-cond
   datalog1(n,:)=conddiff(:,0,0)
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(pm/100, cond, PS, 0, opt)/900
   wva1 = wgt_vert_avg_beta(pm1/100, conddiff, PS, 0, opt)/900
   vertmean(n)=wva1(0,0)
   
   
   wva2 = wgt_vert_avg_beta(pm/100, Q, PS, 0, opt)/900
   wva3 = wgt_vert_avg_beta(pm1/100, Qdiff, PS, 0, opt)/900
   vertmean1(n)=wva3(0,0)
  
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, Q, PS, 0, opt)/900
   Qtot=wva(0,0)
   
   opt1=(/0,100,1000/)
   wva3 = wgt_vert_avg_beta(pm1/100, Qdiff, PS, 0, opt)/900
   Qtot1=wva3(0,0)
   
   s=sw
   s1=sw1
   sdiff=s1-s
   opt=(/0,150,1000/)
   wva3 = wgt_vert_avg_beta(levs, s, PS, 0, opt)/850
   wva4 = wgt_vert_avg_beta(levs, sdiff, PS, 0, opt)/850
   datalog2(n,:)=sdiff(:,0,0)
   vertmean1(n)=wva4(0,0)
   
   ldiff=lw1-lw
   opt=(/0,100,1000/)
   wva5 = wgt_vert_avg_beta(levs, lw, PS, 0, opt)/900
   wva6 = wgt_vert_avg_beta(levs, ldiff, PS, 0, opt)/900
   datalog3(n,:)=ldiff(:,0,0)
   vertmean2(n)=wva6(0,0)
   
   omeg=sw+cond+lw+diff
   omeg1=sw1+cond1+lw1+diff1
   omegdiff=omeg1-omeg
   opt=(/0,100,1000/)
   wva7 = wgt_vert_avg_beta(levs, omeg, PS, 0, opt)/900
   wva8 = wgt_vert_avg_beta(levs, omegdiff, PS, 0, opt)/900
   vertmean3(n)=(wva8(0,0)/wva7(0,0))*100
   
   datalog(n,:)=omegdiff(:,0,0)
   
   diffdiff=diff1-diff
   wva9 = wgt_vert_avg_beta(levs, diff, PS, 0, opt)/900
   wva10 = wgt_vert_avg_beta(levs, diff1, PS, 0, opt)/900
   ;vertmean3(n)=wva10(0,0)
   ;datalog4(n,:)=diffdiff(:,0,0)
   

 res                   = True                       ; plot mods desired
 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 1.0
 res@tiXAxisString=""
 res@tiXAxisFontHeightF=0.02
 res@trXMinF=-20e-6
 res@trXMaxF=20e-6
 res@trYMaxF = 1000.0
 res@trYMinF = 50.0   
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.02
 res@xyYStyle = "Linear"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmXBMode = "Explicit"
 res@tmXBValues        =(/-20,-15,-10,-5,0,5,10,15,20/)*1e-6
 res@tmXBLabels=(/-20,-15,-10,-5,0,5,10,15,20/)*1e-6
 res@tmXBLabelFontHeightF=0.015
 
 res@tiXAxisString="K/s"
 res@vpHeightF 	 = 0.7               ; change aspect ratio of plot
 res@vpWidthF 	         = 0.55
 res@gsnXRefLine = 0.0 

 res@gsnLeftString  ="a) (Q~B~tot~N~)~S~*"
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

res@xyDashPattern  = 0
res@xyLineColor    = "blue"
 res@xyLineThicknesses = 3.0
 plot(0)= gsn_csm_xy  (wks,datalog(0,:),T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 3.0
 plot2= gsn_csm_xy  (wks,datalog(1,:),T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 3.0
 plot3= gsn_csm_xy  (wks,datalog(2,:),T1&lev,res)
 res@xyLineColor    = "magenta4"
 res@xyLineThicknesses = 1.0
 plot30= gsn_csm_xy  (wks,datalog(3,:),T1&lev,res)
 res@xyLineColor    = "aquamarine"
 plot31= gsn_csm_xy  (wks,datalog(4,:),T1&lev,res)
 res@xyLineColor    = "aquamarine1"
 plot32= gsn_csm_xy  (wks,datalog(5,:),T1&lev,res)
 res@xyLineColor    = "aquamarine2"
 plot33= gsn_csm_xy  (wks,datalog(6,:),T1&lev,res)
 res@xyLineColor    = "aquamarine3"
 plot34= gsn_csm_xy  (wks,datalog(7,:),T1&lev,res)
 res@xyLineColor    = "mediumpurple1"
 plot35= gsn_csm_xy  (wks,datalog(8,:),T1&lev,res)
 res@xyLineColor    = "mediumpurple2"
 plot36= gsn_csm_xy  (wks,datalog(9,:),T1&lev,res)
 res@xyLineColor    = "mediumpurple3"
 ;plot37= gsn_csm_xy  (wks,datalog(10,:),T1&lev,res)
 overlay(plot(0),plot2)
 overlay(plot(0),plot3)
 overlay(plot(0),plot30)
 overlay(plot(0),plot31)
 overlay(plot(0),plot32)
 overlay(plot(0),plot33)
 overlay(plot(0),plot34)
 overlay(plot(0),plot35)
 overlay(plot(0),plot36)
 ;overlay(plot(0),plot37)
 gres = True
gres@YPosPercent = 35. 
gres@XPosPercent = 5 
gres@XPosPercent = 5 
gres@ItemSpacePercent = 5.  
lineres = True
lineres@lgLineColors = (/"grey","black","red","blue"/) ; line colors
lineres@lgDashIndexes=(/0,0,0/)
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres = True
textres@lgLabels = (/"CAM4-IF other", "CAM4-IF best T","CAM4-IF best rain","Default CAM4"/)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015/)
txid4 = new(3,graphic)                  ; Nine text strings
 txres               = True     
 txres@txFontHeightF = 0.013
 txid4(0) = gsn_add_text (wks,plot(0),"CAM4: " + sprintf("%5.2f",vertmean3(0)) + "%", 0.67e-5, 800,txres)
 txid4(1) = gsn_add_text (wks,plot(0),"CAM4-IF best rain: " + sprintf("%5.2f",vertmean3(1)) + "%",1.02e-5, 825,txres)
 txid4(2) = gsn_add_text (wks,plot(0),"CAM4-IF best T: " + sprintf("%5.2f",vertmean3(2))+"%",0.95e-5, 850,txres)
 
plot100 = simple_legend(wks,plot3,gres,lineres,textres)

 
 res@gsnLeftString  ="b) (Q~B~cond~N~)~S~*"
 delete(res@trXMinF)
 delete(res@trXMaxF)
 res@trXMinF=-12e-6
 res@trXMaxF=20e-6
 delete(res@tmXBValues)
 delete(res@tmXBLabels)
 delete(res@tiXAxisString)
 res@tiXAxisString="K/s"
 res@tmXBValues        =(/-20,-15,-10,-5,0,5,10,15,20/)*1e-6
 res@tmXBLabels=(/-20,-15,-10,-5,0,5,10,15,20/)*1e-6
 res@tmXBLabelFontHeightF=0.015

 
 
res@xyDashPattern  = 0
res@xyLineColor    = "blue"
 res@xyLineThicknesses = 3.0
 plot(1)= gsn_csm_xy  (wks,datalog1(0,:),T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 3.0
 plot4= gsn_csm_xy  (wks,datalog1(1,:),T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 3.0
 plot5= gsn_csm_xy  (wks,datalog1(2,:),T1&lev,res)
 res@xyLineColor    = "magenta4"
 res@xyLineThicknesses = 1.0
 plot40= gsn_csm_xy  (wks,datalog1(3,:),T1&lev,res)
  res@xyLineColor    = "aquamarine"
 plot41= gsn_csm_xy  (wks,datalog1(4,:),T1&lev,res)
  res@xyLineColor    = "aquamarine1"
 plot42= gsn_csm_xy  (wks,datalog1(5,:),T1&lev,res)
  res@xyLineColor    = "aquamarine2"
 plot43= gsn_csm_xy  (wks,datalog1(6,:),T1&lev,res)
  res@xyLineColor    = "aquamarine3"
 plot44= gsn_csm_xy  (wks,datalog1(7,:),T1&lev,res)
 res@xyLineColor    = "mediumpurple1"
 plot45= gsn_csm_xy  (wks,datalog1(8,:),T1&lev,res)
 res@xyLineColor    = "mediumpurple2"
 plot46= gsn_csm_xy  (wks,datalog1(9,:),T1&lev,res)
 overlay(plot(1),plot4)
 overlay(plot(1),plot5)
 overlay(plot(1),plot40)
 overlay(plot(1),plot41)
 overlay(plot(1),plot42)
 overlay(plot(1),plot43)
 overlay(plot(1),plot44)
 overlay(plot(1),plot45)
 overlay(plot(1),plot46) 
 ;txid = new(3,graphic)                  ; Nine text strings
 ;txres               = True     
 ;txres@txFontHeightF = 0.013
 ;txid(0) = gsn_add_text (wks,plot(1),"CAM4: " + sprintf("%7.3e",vertmean(0)) +"%", -8.1e-6, 175,txres)
 ;txid(1) = gsn_add_text (wks,plot(1),"CAM4-IF best rain: " + sprintf("%7.3e",vertmean(1)) +"%",-6e-6, 200,txres)
 ;txid(2) = gsn_add_text (wks,plot(1),"CAM4-IF best T: " + sprintf("%7.3e",vertmean(2)) +"%", -6.4e-6, 225,txres)
 
 
;;;;;Static stability;;;;;
res1=res
 res1@xyLineThicknesses = 1.0
 res1@tiXAxisString="K/s"
 res1@trXMinF=-8e-6
 res1@trXMaxF=6e-6
 res1@gsnLeftString  ="c) (Q~B~sw~N~)~S~*"
 res1@tmXBLabelFontHeightF=0.012
 delete(res1@tmXBValues)
 delete(res1@tmXBLabels)
res1@tmXBValues        =(/-8,-6,-4,-2,0,2,4,6/)*1e-6
 res1@tmXBLabels=(/-8,-6,-4,-2,0,2,4,6/)*1e-6
  res1@xyLineThicknesses = 3.0
 res1@xyLineColor    = "blue"
 plot(2)= gsn_csm_xy  (wks,datalog2(0,:),T1&lev,res1)
 res1@xyLineColor    = "red"
 res1@xyLineThicknesses = 3.0
 plot6= gsn_csm_xy  (wks,datalog2(1,:),T1&lev,res1)
 res1@xyLineColor    = "black"
 res1@xyLineThicknesses = 3.0
 plot7= gsn_csm_xy  (wks,datalog2(2,:),T1&lev,res1)
  res1@xyLineColor    = "magenta4"
 res1@xyLineThicknesses = 1.0
  plot40= gsn_csm_xy  (wks,datalog2(3,:),T1&lev,res1)
   res1@xyLineColor    = "aquamarine"
 plot41= gsn_csm_xy  (wks,datalog2(4,:),T1&lev,res1)
 res1@xyLineColor    = "aquamarine1"
 plot42= gsn_csm_xy  (wks,datalog2(5,:),T1&lev,res1)
 res1@xyLineColor    = "aquamarine2"
 plot43= gsn_csm_xy  (wks,datalog2(6,:),T1&lev,res1)
 res1@xyLineColor    = "aquamarine3"
 plot44= gsn_csm_xy  (wks,datalog2(7,:),T1&lev,res1)
 res1@xyLineColor    = "mediumpurple1"
 plot45= gsn_csm_xy  (wks,datalog2(8,:),T1&lev,res1)
 res1@xyLineColor    = "mediumpurple2"
 plot46= gsn_csm_xy  (wks,datalog2(9,:),T1&lev,res1)
 overlay(plot(2),plot6)
 overlay(plot(2),plot7)
 overlay(plot(2),plot40)
 overlay(plot(2),plot41)
 overlay(plot(2),plot42)
 overlay(plot(2),plot43)
 overlay(plot(2),plot44)
 overlay(plot(2),plot45)
 overlay(plot(2),plot46) 
 ;txid2 = new(3,graphic)                  ; Nine text strings
 ;txres               = True     
 ;txres@txFontHeightF = 0.013
 ;txid2(0) = gsn_add_text (wks,plot(2),"CAM4: " + sprintf("%7.3e",vertmean1(0)), -2.75e-6, 300,txres)
 ;txid2(1) = gsn_add_text (wks,plot(2),"CAM4-IF best rain: " + sprintf("%7.3e",vertmean1(1)) ,-2.05e-6, 325,txres)
 ;txid2(2) = gsn_add_text (wks,plot(2),"CAM4-IF best T: " + sprintf("%7.3e",vertmean1(2)), -2.2e-6, 350,txres)
 
 
 ;;;;;Omega;;;;;
res2=res
 res2@xyLineThicknesses = 3.0
 res2@tiXAxisString="K/s"
 res2@trXMinF=-20e-6
 res2@trXMaxF=8e-6
 res2@gsnLeftString  ="d) (Q~B~lw~N~)~S~*"
 res2@tmXBLabelFontHeightF=0.012
 delete(res2@tmXBValues)
 delete(res2@tmXBLabels)
 res2@tmXBValues        =(/-20,-16,-12,-8,-4,0,4,8/)*1e-6
 res2@tmXBLabels=(/-20,-16,-12,-8,-4,0,4,6,8/)*1e-6
 res2@xyLineColor    = "blue"
 plot(3)= gsn_csm_xy  (wks,datalog3(0,:),T1&lev,res2)
 res2@xyLineColor    = "red"
 res2@xyLineThicknesses = 3.0
 plot8= gsn_csm_xy  (wks,datalog3(1,:),T1&lev,res2)
 res2@xyLineColor    = "black"
 res2@xyLineThicknesses = 3.0
 plot9= gsn_csm_xy  (wks,datalog3(2,:),T1&lev,res2)
 res2@xyLineColor    = "magenta4"
 res2@xyLineThicknesses = 1.0
 plot40= gsn_csm_xy  (wks,datalog3(3,:),T1&lev,res2)
 res2@xyLineColor    = "aquamarine"
 plot41= gsn_csm_xy  (wks,datalog3(4,:),T1&lev,res2)
 res2@xyLineColor    = "aquamarine1"
 plot42= gsn_csm_xy  (wks,datalog3(5,:),T1&lev,res2)
 res2@xyLineColor    = "aquamarine2"
 plot43= gsn_csm_xy  (wks,datalog3(6,:),T1&lev,res2)
 res2@xyLineColor    = "aquamarine3"
 plot44= gsn_csm_xy  (wks,datalog3(7,:),T1&lev,res2)
 res2@xyLineColor    = "mediumpurple1"
 plot45= gsn_csm_xy  (wks,datalog3(8,:),T1&lev,res2)
 res2@xyLineColor    = "mediumpurple2"
 plot46= gsn_csm_xy  (wks,datalog3(9,:),T1&lev,res2)
 overlay(plot(3),plot8)
 overlay(plot(3),plot9)
 overlay(plot(3),plot40)
 overlay(plot(3),plot41)
 overlay(plot(3),plot42)
 overlay(plot(3),plot43)
 overlay(plot(3),plot44)
 overlay(plot(3),plot45)
 overlay(plot(3),plot46) 
 ;txid1 = new(3,graphic)                  ; Nine text strings
 ;;txres               = True     
 ;txres@txFontHeightF = 0.013
 ;vertmean2@txFormat="0*+^e"
 ;txid1(0) = gsn_add_text (wks,plot(3),"CAM4: " + sprintf("%7.3e",vertmean2(0)), -12.4e-6, 500,txres)
 ;txid1(1) = gsn_add_text (wks,plot(3),"CAM4-IF best rain: " + sprintf("%7.3e",vertmean2(1)),-10.3e-6, 525,txres)
 ;txid1(2) = gsn_add_text (wks,plot(3),"CAM4-IF best T: " + sprintf("%7.3e",vertmean2(2)),-10.75e-6, 550,txres)
 
 
 
 
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
cleanup1=systemfunc("rm zon*.nc")
print(fils)
print(fils1)
end