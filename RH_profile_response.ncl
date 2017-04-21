begin
;************************************************
; read in data
;************************************************FM
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM*ANN*.nc2")
print(fils)
print(fils1)
plotfile="RH_profiles_response" ;;;PDF file name!
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
lets=(/"a)","b)","c)","d)"/)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM4_4K_ANN_climo.nc2", "r")
   hyam = h->hyam 
   hybm = h->hybm
   cmd="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,T,PS,RELHUM " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,T,PS,RELHUM " + filen1 + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
  
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->T(0,:,:,:)
   T&lev@units="hPa"
   rh=f->RELHUM(0,:,:,:)
   
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   PS1   = f1->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->T(0,:,:,:)
   T1&lev@units="hPa"
   rh1=f1->RELHUM(0,:,:,:)
   
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   
   Qdiff=rh1-rh
   datalog(n,:)=Qdiff(:,0,0)
   datalog1(n,:)=rh(:,0,0)
   Tdiff=T1-T
   datalog2(n,:)=T1(:,0,0)
   datalog3(n,:)=T(:,0,0)
   
   
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(pm/100, Qdiff, PS/100, 0, opt)/900
   wva1 = wgt_vert_avg_beta(pm/100, rh, PS/100, 0, opt)/900
   

 res                   = True                       ; plot mods desired
 res@xyLineColor     = "blue"
 res@trYReverse        = True                       ; reverse Y-axis
 res@tiMainFontHeightF=0.018
 res@xyLineThicknesses = 1.0
 res@tiXAxisString=""
 res@tiXAxisFontHeightF=0.02
 res@trYMaxF=1000
 res@trYMinF=30
 
 res@trXMinF=-20
 res@trXMaxF=100
 res@tiYAxisString="Pressure (hPa)"
 res@tiYAxisFontHeightF=0.02
 res@xyYStyle = "Linear"
 res@tmYLMode = "Explicit"
 res@tmYLValues        =(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmYLLabels=(/50,100,200,300,400,500,600,700,850,1000/)
 res@tmXBMode = "Explicit"
 res@tmXBValues        =(/-20,0,20,40,60,80,100/)
 res@tmXBLabels=(/-20,0,20,40,60,80,100/)
 res@tmXBLabelFontHeightF=0.015
 
 res@tiXAxisString="RH (%)"
 res@vpHeightF 	 = 0.7               ; change aspect ratio of plot
 res@vpWidthF 	         = 0.55
 res@gsnXRefLine = 0.0 
 res@gsnLeftString  ="Relative humidity"
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
 plot4= gsn_csm_xy  (wks,datalog1(0,:),T1&lev,res)
 res@xyLineColor    = "red"
 res@xyLineThicknesses = 3.0
 plot5= gsn_csm_xy  (wks,datalog1(1,:),T1&lev,res)
 res@xyLineColor    = "black"
 res@xyLineThicknesses = 3.0
 plot6= gsn_csm_xy  (wks,datalog1(2,:),T1&lev,res)
 
 
 overlay(plot(0),plot2)
 overlay(plot(0),plot3)
 overlay(plot(0),plot4)
 overlay(plot(0),plot5)
  overlay(plot(0),plot6)
 gres = True
gres@YPosPercent = 55. 
gres@XPosPercent = 57
gres@ItemSpacePercent = 5.  
lineres = True
lineres@lgLineColors = (/"black","red","blue"/) ; line colors
lineres@lgDashIndexes=(/0,0,0,0,0,0/)
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres = True
textres@lgLabels = (/"CAM4-IF best T","CAM4-IF best rain","Default CAM4" /)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015/)   
plot3 = simple_legend(wks,plot3,gres,lineres,textres)


res1=res
 res1@trYMaxF=200
 res1@trYMinF=50
 
 res1@trXMinF=190
 res1@trXMaxF=230
res1@gsnLeftString  ="Temperature"
res1@xyDashPattern  = 1
res1@xyLineColor    = "blue"
 res1@tmYLMode = "Explicit"
 res1@tiXAxisString="T (K)"
 res1@tmYLValues        =(/50,75,100,125,150,175,200,225,250,275/)
 res1@tmYLLabels=(/50,75,100,125,150,175,200,225,250,275/)
 res1@tmXBMode = "Explicit"
 delete(res1@tmXBValues)
 delete(res1@tmXBLabels)
 res1@tmXBValues        =(/190,195,200,205,210,215,220,225,230/)
 res1@tmXBLabels=(/190,195,200,205,210,215,220,225,230/)
 
 res1@xyDashPattern  = 0
res1@xyLineColor    = "blue"
 res1@xyLineThicknesses = 3.0
 plot(1)= gsn_csm_xy  (wks,datalog2(0,:),T1&lev,res1)
 res1@xyLineColor    = "red"
 res1@xyLineThicknesses = 3.0
 plot5= gsn_csm_xy  (wks,datalog2(1,:),T1&lev,res1)
 res1@xyLineColor    = "black"
 res1@xyLineThicknesses = 3.0
 plot6= gsn_csm_xy  (wks,datalog2(2,:),T1&lev,res1)

 overlay(plot(1),plot5)
 overlay(plot(1),plot6)

 res1@xyDashPattern  = 0
 res1@xyLineColor    = "blue"
 res1@xyLineThicknesses = 1.0
 plot7= gsn_csm_xy  (wks,datalog3(0,:),T1&lev,res1)
 res1@xyLineColor    = "red"
 res1@xyLineThicknesses = 1.0
 plot8= gsn_csm_xy  (wks,datalog3(1,:),T1&lev,res1)
 res1@xyLineColor    = "black"
 res1@xyLineThicknesses = 1.0
 plot9= gsn_csm_xy  (wks,datalog3(2,:),T1&lev,res1)

 overlay(plot(1),plot7)
 overlay(plot(1),plot8)
 overlay(plot(1),plot9)
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
