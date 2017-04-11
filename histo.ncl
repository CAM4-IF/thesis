;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;************************************************
begin
wks = gsn_open_wks("pdf","histo")              ; open workstation
plot=new(1,graphic)
f=addfile("../mass_flux/FMASS_CAM4-IF_bestrain_ANN.nc", "r")
f1=addfile("../mass_flux/FMASS_CAM4-IF_bestrain_4K_ANN.nc", "r")
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
res@gsnHistogramBinIntervals = (/-120,-110,-100,-90,-80,-70,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5/)*1e-4
res@tiXAxisOn=True
;res@tmYRMode = "Automatic"
res@gsnHistogramComputePercentages = True
res@tmYLMode="explicit"
res@tmYLValues=(/100,200,300,400,500,600,700,800,900,1000,1100,1200/)
res@tmYLLabels=(/"100","200","300","400","500","600","700","800","900","1000","1100","1200"/)
res@tmYROn=True
res@tmYRLabelsOn = True
res@tmYUseRight = True  
res@tmLabelAutoStride =True
res@trYMaxF=1200
res@tiXAxisString    = "kg/m^2/s"

;ymax= max((/max(plot0@NumInBins), max(plot1@NumInBins)/))

colors=(/"red","lightblue"/)



  ;res@trYMaxF=ymax+5
  res@gsnDraw=False

  res@gsnHistogramBarColors = colors(1)
  ;res@gsnHistogramBarColors = "transparent"
  res@gsnHistogramBarWidthPercent = 100 
  res@gsEdgeColor = res@gsnHistogramBarColors 
  plot1 = gsn_histogram(wks,test,res)

  res@gsnHistogramBarColors = "transparent"
 ; res@gsFillColor = "transparent" 
  res@gsEdgeColor = "red" 
  res@tiMainString = ""
  res@tiMainFontHeightF   = 0.015
  plot0 = gsn_histogram(wks,test1,res)


overlay(plot1,plot0)

frame(wks)
res1=True
res1@draw=False
res1@frame=False
res1@txString = "Control (light blue fill) & 4K SST (red outline)"
res1@txFontHeightF    = 0.02
res1@tiXAxisString    = "kg/m^2/s"
gsn_panel(wks,plot1,(/1,1/),res1)
;delete(wks)
;delete(plot)
end
