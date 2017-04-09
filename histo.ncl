;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;************************************************
begin
wks = gsn_open_wks("pdf","histo")              ; open workstation
plot=new(2,graphic)
f=addfile("../mass_flux/FMASS_CAM4-IF_bestrain_ymonmean.nc", "r")
f1=addfile("../mass_flux/FMASS_CAM4-IF_bestrain_4K_ymonmean.nc", "r")
;************************************************
; Transform 2-d data into 1-d array
;************************************************
h=f->FMASS(0,{469.0718},:,:)
test=ndtooned(h)
h1=f1->FMASS(0,{469.0718},:,:)
test1=ndtooned(h1)

omegup=test
omegup1=test1

omegupavg=avg(omegup < 0)
omegupavg1=avg(omegup1 < 0)
print(omegupavg)
print(omegupavg1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;Create histograms;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

res=True
res@gsnHistogramBinIntervals = (/-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60/)*0.5e-3
res@tiXAxisOn=False
res@tmYRMode = "Automatic"
res@gsnHistogramComputePercentages = True
res@tmYLMode="explicit"
res@tmYLValues=(/100,200,300,400,500,600,700,800,900/)
res@tmYLLabels=(/"100","200","300","400","500","600","700","800","900"/)
res@tmYROn=False
res@tmYRLabelsOn = False
res@tmYUseRight = False  
res@tmLabelAutoStride =True
res@trYMaxF=700

plot0=gsn_histogram(wks,test,res)
plot1=gsn_histogram(wks,test1,res)

;ymax= max((/max(plot0@NumInBins), max(plot1@NumInBins)/))

colors=(/"red","lightblue"/)



  ;res@trYMaxF=ymax+5
  res@gsnDraw=True

  res@gsnHistogramBarColors = colors(1)
  ;res@gsnHistogramBarColors = "transparent"
  res@gsnHistogramBarWidthPercent = 100 
  res@gsEdgeColor = res@gsnHistogramBarColors 
  plot1 = gsn_histogram(wks,test,res)

  res@gsnHistogramBarColors = "transparent"
 ; res@gsFillColor = "transparent" 
  res@gsEdgeColor = "red" 
  res@tiMainString = "CanCM4 (light blue fill) & CanCM3 (red outline)"
  res@tiMainFontHeightF   = 0.015
  plot0 = gsn_histogram(wks,test1,res)


overlay(plot1,plot0)

frame(wks)
res1=True
res1@draw=False
res1@frame=True
res1@txString = ""
res1@txFontHeightF    = 0.0115
gsn_panel(wks,plot1,(/1,2/),res1)
;delete(wks)
;delete(plot)
end
