begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files//*_ANN*.nc4")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*_ANN*.nc2")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="T_ANN" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
lets=(/"a)","b)","c)","d)"/)
n=0
do while(n.le.dimsizes(fils)-1)
   print(fils)
   print(fils1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile(filen, "r")
   hyam = h->hyam 
   hybm = h->hybm
   cmd="cdo -mermean -sellonlatbox,0,360,-30,30 -selname,T,PS " + filen + " tmp.nc"
   cmd1="cdo -mermean -sellonlatbox,0,360,-30,30 -selname,T,PS " + filen1 + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
  
   ;P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T=f->T(0,:,:,:)
   T&lev@units="hPa"
   ;pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   ;print(pm(4,:,:))
   ;print(dimsizes(T))
   levs=conform_dims(dimsizes(T),T&lev,0)
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, T, PS/100, 0, opt)/900
   wva!0="lev"
   wva!1="lon"
   wva&lon=T&lon
   fmassavg=avg(wva)
   
   PS1   = f1->PS(0,:,:)                         ; Pa  [punits=1]
   T1=f1->T(0,:,:,:)
   T1&lev@units="hPa"
   ;pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt1=(/0,100,1000/)
   wva1 = wgt_vert_avg_beta(levs, T1, PS1/100, 0, opt1)/900
   wva1!0="lev"
   wva1!1="lon"
   wva1&lon=T&lon
   fmassavg1=avg(wva1)
   ;response=(fmassavg1-fmassavg)/fmassavg
   response=fmassavg1-fmassavg
   print(response)
   
   ;;;;;;;;;;;;PANEL PLOT;;;;;;;;;;;;;
 
  ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
 
 
  Qdiff=T1(:,0,:)-T(:,0,:)
  Qdiff!0="lev"
  Qdiff&lev=T&lev
  Qdiff!1="lon"
  Qdiff&lon=T&lon
  Qdiff&lev@units="hPa"
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels =(/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12/)*0.75
  res@cnFillColors = (/10,20,30,40,50,60,70,80,90,100,110,120,0,0,140,150,160,170,180,190,200,210,220,230,240,250/)
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

 
;;;;Calculate 500 hPa upward omega mean;;;; 
   omegstar=T(:,0,:)
   omegstar!1="lon"
   omegstar&lon=T&lon
   omegstar1=T1(:,0,:)
   omegstar1!1="lon"
   omegstar1&lon=T&lon
   
   omegstarwc=dim_rmvmean_n_Wrap(T(:,0,:),1)
   omegstarwc1=dim_rmvmean_n_Wrap(T1(:,0,:),1)
   omegstarwc!1="lon"
   omegstarwc&lon=T&lon
   omegstarwc1!1="lon"
   omegstarwc1&lon=T&lon
   
   omegstarwc500=omegstarwc(18,{60:180}) < 0
   omegstarwc500mean=avg(omegstarwc500)
   omegstarwc5001=omegstarwc1(18,{60:180}) < 0
   omegstarwc500mean1=avg(omegstarwc5001)
   
   omegstarwc500diff=omegstarwc500mean1-omegstarwc500mean
   responsewc=100*(omegstarwc500diff/avg(omegstarwc500mean))
   
   
   omegstar500=omegstar(18,:)
   omegstar500=omegstar500 < 0
   omegstar500mean=avg(omegstar500)
   
   omegstar5001=omegstar1(18,:) < 0
   omegstar500mean1=avg(omegstar5001)
   
   omegstar500diff=omegstar500mean1-omegstar500mean
   
   ;response=100*(omegstar500diff/avg(omegstar500mean))
   
   ;print(response)
   
 ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels =(/190,200,210,220,230,240,250,260,270,280,290,300,310,320/)
 ;cres@cnFillColors = (/10,20,30,40,50,60,70,80,90,100,110,120,0,0,140,150,160,170,180,190,200,210,220,230,240,250/)
 ;cres@cnLevels = (/-260,-240,-220,-200,-180,160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240,260/)*5.0e-4
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
 ;cres@gsnContourNegLineDashPattern =-1  
 cres@gsnLeftString =""
 cres@gsnContourZeroLineThicknessF = 2
 cres@gsnRightString    = lets(n) + " " + f@case
 cres@trYMinF            = 50.   ; Limits for Y axis.  The limits
 cres@trYMaxF            = 1000.
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelPlacementMode = "Constant"
 cres@cnLineLabelInterval=1
 cres@gsnLeftString  = "response: " + response
 cres@gsnLeftStringFontHeightF=0.012
 cres@gsnRightStringFontHeightF=0.015
 cres1=cres
 cres1@cnLineColor="red"
 cres1@gsnContourPosLineDashPattern=1
 
gsn_define_colormap(wks,"nice_gfdl")
 plot(n)  = gsn_csm_pres_hgt(wks, Qdiff(::-1,:), res)  
 ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
 plot2  = gsn_csm_pres_hgt(wks, T(::-1,0,:), cres)
 overlay(plot(n),plot2)
 
   
   cleanup=systemfunc("rm tmp*.nc")
   delete(T)
   delete(T1)
   delete(Qdiff)
  ; delete(pm)
  ; delete(pm1)
n=n+1
end do
T=f->T(0,:,:,:)
 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
resP@lbLabelFontHeightF    = 0.008
 resP@txString = "ANN"
 resP@lbTitleString    = T@units
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= True
 resP@txFontHeightF    = 0.02
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
; resP@lbLabelStride=2
 resP@lbTitlePosition ="Bottom"
 resP@lbTitleFontHeightF=.012
 gsn_panel(wks,plot,(/3,1/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop" + plotfile)  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
print(fils)
print(fils1)
end
