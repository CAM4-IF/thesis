load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

f=addfile("hov_cancm4_ua850_2015030100-2016022800_e01.nc", "r")
f1=addfile("hov_cancm4_ua850_2015030100-2016022800_e02.nc", "r")
f2=addfile("hov_cancm4_ua850_2015030100-2016022800_e03.nc", "r")
f3=addfile("hov_cancm4_ua850_2015030100-2016022800_e04.nc", "r")
f4=addfile("hov_cancm4_ua850_2015030100-2016022800_e05.nc", "r")
f5=addfile("hov_cancm4_ua850_2015030100-2016022800_e06.nc", "r")
f6=addfile("hov_cancm4_ua850_2015030100-2016022800_e07.nc", "r")
f7=addfile("hov_cancm4_ua850_2015030100-2016022800_e08.nc", "r")
f8=addfile("hov_cancm4_ua850_2015030100-2016022800_e09.nc", "r")
f9=addfile("hov_cancm4_ua850_2015030100-2016022800_e10.nc", "r")



ua=f->ua(:,0,0,:)
ua1=f1->ua(:,0,0,:)
ua2=f2->ua(:,0,0,:)
ua3=f3->ua(:,0,0,:)
ua4=f4->ua(:,0,0,:)
ua5=f5->ua(:,0,0,:)
ua6=f6->ua(:,0,0,:)
ua7=f7->ua(:,0,0,:)
ua8=f8->ua(:,0,0,:)
ua9=f9->ua(:,0,0,:)
ua = runave_n(ua, 15, 0, 0); smooth data to reduce noise
ua1 = runave_n(ua1, 15, 0, 0); smooth data to reduce noise
ua2 = runave_n(ua2, 15, 0, 0); smooth data to reduce noise
ua3 = runave_n(ua3, 15, 0,0); smooth data to reduce noise
ua4 = runave_n(ua4, 15, 0,0); smooth data to reduce noise
ua5 = runave_n(ua5, 15, 0,0); smooth data to reduce noise
ua6 = runave_n(ua6, 15, 0,0); smooth data to reduce noise
ua7 = runave_n(ua7, 15, 0,0); smooth data to reduce noise
ua8 = runave_n(ua8, 15, 0,0); smooth data to reduce noise
ua9 = runave_n(ua9, 15, 0,0); smooth data to reduce noise



delete(ua@long_name)
delete(ua@standard_name)
delete(ua@units)
delete(ua1@long_name)
delete(ua1@standard_name)
delete(ua1@units)
delete(ua2@long_name)
delete(ua2@standard_name)
delete(ua2@units)
delete(ua3@long_name)
delete(ua3@standard_name)
delete(ua3@units)
delete(ua4@long_name)
delete(ua4@standard_name)
delete(ua4@units)
delete(ua5@long_name)
delete(ua5@standard_name)
delete(ua5@units)
delete(ua6@long_name)
delete(ua6@standard_name)
delete(ua6@units)
delete(ua7@long_name)
delete(ua7@standard_name)
delete(ua7@units)
delete(ua8@long_name)
delete(ua8@standard_name)
delete(ua8@units)
delete(ua9@long_name)
delete(ua9@standard_name)
delete(ua9@units)


timeUnits=f->time@units
startDate=cd_inv_calendar(2015, 03, 01, 00, 0, 0, timeUnits, 0)
endDate=cd_inv_calendar(2015, 12, 31, 00, 0, 0, timeUnits, 0)
ua&time=ut_convert(ua&time, timeUnits)
ua1&time=ua&time
ua2&time=ua&time
ua3&time=ua&time
ua4&time=ua&time
ua5&time=ua&time
ua6&time=ua&time
ua7&time=ua&time
ua8&time=ua&time
ua9&time=ua&time
wks=gsn_open_wks("pdf","hov")
plot=new(10,graphic)
gsn_define_colormap(wks,"BlueYellowRed")
res=True
res@gsnDraw =  False 
res@gsnFrame = False
res@cnFillOn             = True  
res@cnLevelSelectionMode = "ManualLevels" 
res@cnLevelSpacingF      = 1.
res@cnMinLevelValF       = -10.
res@cnMaxLevelValF       =  10. 
res@cnFillColors = (/20,30,40,50,60,70,80,90,100,110,0,0,140,150,160,170,180,190,200,210,215,220,225,230,235,240,245,250/)
res@trYReverse           = True
res@tiMainString = ""
res@vpWidthF=0.4
res@vpHeightF=0.8
res@lbLabelBarOn        = False 
res@cnLineThicknessF = 0.  
res@cnInfoLabelOn       = False

resTick                = True
resTick@ttmFormat      = "%d %c"
resTick@ttmAxis        = "YL"
resTick@ttmMajorStride = 15
print(ua&time)
time_axis_labels(ua&time, res, resTick)

plot(0)=gsn_csm_hov(wks, ua(:,{100:280}), res)
plot(1)=gsn_csm_hov(wks, ua1(:,{100:280}), res)
plot(2)=gsn_csm_hov(wks, ua2(:,{100:280}), res)
plot(3)=gsn_csm_hov(wks, ua3(:,{100:280}), res)
plot(4)=gsn_csm_hov(wks, ua4(:,{100:280}), res)
plot(5)=gsn_csm_hov(wks, ua5(:,{100:280}), res)
plot(6)=gsn_csm_hov(wks, ua6(:,{100:280}), res)
plot(7)=gsn_csm_hov(wks, ua7(:,{100:280}), res)
plot(8)=gsn_csm_hov(wks, ua8(:,{100:280}), res)
plot(9)=gsn_csm_hov(wks, ua9(:,{100:280}), res)
res2 = True
res2@gsnPanelLabelBar=True
res2@cnInfoLabelOn       = False
res2@lbOrientation =   "horizontal"   
res2@gsnMaximize =      True
res2@lbBoxFractions = 0.1
res2@lbTopMarginF = 0.0 
res2@lbJustification="topright"
res2@AttachBorderOn =False
res2@gsnPanelMainString    = " "
res2@lbLabelFontHeightF    = 0.009
res2@txString = "CanCM4 (10 members) March 2015 5S-5N 850 hPa u-wind anomaly forecasts (15-day running mean)"
res2@txFontHeightF    = 0.008 
gsn_panel(wks,plot,(/5,5/),res2)
delete(wks)
delete(plot)
exit
