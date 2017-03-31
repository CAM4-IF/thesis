; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*********************************************
begin
;**********************************************
;;;;;Read in the series of data files! (va/ps);;;;
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_RCP8.5/va/startJJA_*")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_RCP8.5/ps/startJJA_*")
fils2=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_RCP8.5/va/endJJA_*")
fils3=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CMIP5_RCP8.5/ps/endJJA_*")
datalog=new((/dimsizes(fils),3/),float) ; create list to store output
plotfile="RCP8.5_mpsi_JJA" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
 f=addfile(fils(n),"r") ; add file
 f1=addfile(fils1(n),"r") ; add file
 print(fils2(n))
 f2=addfile(fils2(n),"r") ; add file
 f3=addfile(fils3(n),"r") ; add file ;;;read the files!
 T=f->va(0,{100000:5000},{-45:45},:) ; extract va data from files (start)
 T1=f2->va(0,{100000:5000},{-45:45},:) ; extract va data from files (end)
 lat=f->lat({-45:45}) ; get latitude data
 lat1=f2->lat({-45:45})
 ;;;loop to fix slight mis-match in grid size between ps and va files for one model;;;
 if(n.eq.8) then
   ps=f1->ps(0,{-45.2:45},:)
  else
   ps=f1->ps(0,{-44.9:45},:)
 end if
 lat=f->lat({-45:45})
 if(n.eq.8) then
   ps1=f3->ps(0,{-45.2:45},:)
  else
   ps1=f3->ps(0,{-44.9:45},:)
 end if
 print(n)
 print(fils(n))
 print(fils1(n))
 mpsi=zonal_mpsi(T(plev|::-1,lat|:,lon|:),lat,T&plev(::-1),ps)
 mpsi1=zonal_mpsi(T1(plev|::-1,lat|:,lon|:),lat1,T1&plev(::-1),ps1)
 mpsi!0="lev"
 mpsi!1="lat"
 mpsi&lev=T&plev(::-1)
 mpsi&lat=T1&lat
 mpsi1!0="lev"
 mpsi1!1="lat"
 mpsi1&lev=T1&plev(::-1)
 mpsi1&lat=T1&lat
 maxmpsi=min((mpsi({50000},:)))
 maxmpsi1=min((mpsi1({50000},:)))
 ;;;calculate dpsi/dy for -15 to 5 latitude ;;; 
 fiven=round(0.66666*dimsizes(mpsi(1,:)),3)
 fifteens=round(0.44444*dimsizes(mpsi(1,:)),3)
 fiven1=round(0.66666*dimsizes(mpsi1(1,:)),3)
 fifteens1=round(0.44444*dimsizes(mpsi1(1,:)),3)
 slice1=mpsi({20000:85000},fiven)
 slice2=mpsi({20000:85000},fifteens)
 slice3=mpsi1({20000:85000},fiven)
 slice4=mpsi1({20000:85000},fifteens)
 mean1=dim_avg_n(slice1,0)
 mean2=dim_avg_n(slice2,0)
 mean3=dim_avg_n(slice3,0)
 mean4=dim_avg_n(slice4,0)
 diff=mean2-mean1
 diff1=mean4-mean3
 mean=dim_avg_n(diff,0)
 mean5=dim_avg_n(diff1,0)
 meandiff=abs(mean5)-abs(mean)
 diff2=maxmpsi1-maxmpsi
 datalog(n,2)=meandiff
 datalog(n,1)=mean5
 datalog(n,0)=mean
 
 delete(lat)
 delete(ps)

 delete(fiven)
 delete(fifteens)
 delete(slice1)
 delete(slice2)
 delete(diff)
 delete(mean)
 delete(mean1)
 delete(meandiff)
 delete(mean2)
 delete(mean3)
 delete(mean4)
 delete(mean5)
 delete(slice3)
 delete(slice4)
 delete(fiven1)
 delete(fifteens1)
 delete(lat1)
 delete(ps1)
 delete(diff1)
 delete(diff2)
 delete(maxmpsi)
 delete(maxmpsi1)
 ;;;Create some plots!;;;

 mpsiscal=mpsi/1.0e9 ;;; scale mpsi!
 mpsiscal1=mpsi1/1.0e9 ;;; scale mpsi!
 mpsiscal!0="lev"
 mpsiscal!1="lat"
 mpsiscal&lev=T&plev(::-1)
 mpsiscal&lat=T&lat
 mpsiscal1!0="lev"
 mpsiscal1!1="lat"
 mpsiscal1&lev=T1&plev(::-1)
 mpsiscal1&lat=T1&lat
 printVarSummary(mpsiscal)
 printVarSummary(mpsiscal1)
 msr=mpsiscal1-mpsiscal
 msr!0="lev"
 msr!1="lat"
 msr&lev=T&plev(::-1)
 msr&lat=T&lat
 printVarSummary(msr)
 
 ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-30,-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20,30/)
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
 cres@cnLevels = (/-260,-240,-220,-200,-180,160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240,260/)
 cres@cnLineColors = (/200,205,210,215,220/)
 cres@cnLineLabelsOn       = True                 ; turn on line labels
 cres@cnFillOn             = False                 ; turn on color fill
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False                ; turn off contour info label
 cres@gsnRightString = ""
 cres@gsnLeftString = ""
 cres@tiYAxisString    = "Pressure (hPa)"
 cres@tiXAxisString    = "Latitude"
 cres@tmYROn = False
 cres@tmYRLabelsOn = False
 cres@gsnContourZeroLineThicknessF  = 0.0
 cres@cnLineLabelsOn = True
 cres@cnLineThicknessF =0.5
 cres@gsnContourNegLineDashPattern =1  
 cres@gsnLeftString =""
 cres@gsnContourZeroLineThicknessF = 2
 cres@gsnRightString    = f@model_id
 
 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 plot(n)  = gsn_csm_pres_hgt(wks, msr, res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_pres_hgt(wks, mpsiscal, cres)
 overlay(plot(n),plot2)
 delete(mpsi)
 delete(mpsiscal)
 delete(T)
 delete(mpsi1)
 delete(mpsiscal1)
 delete(T1)
 ;delete(plot1)
 delete(plot2)
 delete(msr)
 n=n+1
end do
print(datalog(:,2))
print(fils)
;;;;build the panel plot!;;;;
 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "RCP8.5 JJA zonal mass streamfunction climo (2007-2026) and response (2081-2100 minus 2007-2026)"
 resP@lbTitleString    = "10~S~9~N~ kg s~S~-1~N~ K~S~-1~N~"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.01
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/5,5/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop RCP8.5_mpsi_JJA.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(datalog(:,2))
print(datalog(:,1))
print(datalog(:,0))
end
