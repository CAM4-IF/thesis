;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;************************************************
begin
wks = gsn_open_wks("pdf","histo")              ; open workstation
plot=new(20,graphic)
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
g=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e01.nc", "r")
g1=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e02.nc", "r")
g2=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e03.nc", "r")
g3=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e04.nc", "r")
g4=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e05.nc", "r")
g5=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e06.nc", "r")
g6=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e07.nc", "r")
g7=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e08.nc", "r")
g8=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e09.nc", "r")
g9=addfile("/praid/users/corvecs/2015_ElNino/forecast/cancm3/march_2015/ua850/hov_cancm3_ua850_2015030100-2016022800_e10.nc", "r")
;************************************************
; Generate some random data.
;************************************************
h=f->ua(:,0,0,:)
test=ndtooned(h)
h1=f1->ua(:,0,0,:)
test1=ndtooned(h1)
h2=f2->ua(:,0,0,:)
test2=ndtooned(h2)
h3=f3->ua(:,0,0,:)
test3=ndtooned(h3)
h4=f4->ua(:,0,0,:)
test4=ndtooned(h4)
h5=f5->ua(:,0,0,:)
test5=ndtooned(h5)
h6=f6->ua(:,0,0,:)
test6=ndtooned(h6)
h7=f7->ua(:,0,0,:)
test7=ndtooned(h7)
h8=f8->ua(:,0,0,:)
test8=ndtooned(h8)
h9=f9->ua(:,0,0,:)
test9=ndtooned(h9)

b=g->ua(:,0,0,:)
t=ndtooned(b)
b1=g1->ua(:,0,0,:)
t1=ndtooned(b1)
b2=g2->ua(:,0,0,:)
t2=ndtooned(b2)
b3=g3->ua(:,0,0,:)
t3=ndtooned(b3)
b4=g4->ua(:,0,0,:)
t4=ndtooned(b4)
b5=g5->ua(:,0,0,:)
t5=ndtooned(b5)
b6=g6->ua(:,0,0,:)
t6=ndtooned(b6)
b7=g7->ua(:,0,0,:)
t7=ndtooned(b7)
b8=g8->ua(:,0,0,:)
t8=ndtooned(b8)
b9=g9->ua(:,0,0,:)
t9=ndtooned(b9)


std4_1=dim_stat4(test)
std4_2=dim_stat4(test1)
std4_3=dim_stat4(test2)
std4_4=dim_stat4(test3)
std4_5=dim_stat4(test4)
std4_6=dim_stat4(test5)
std4_7=dim_stat4(test6)
std4_8=dim_stat4(test7)
std4_9=dim_stat4(test8)
std4_10=dim_stat4(test9)
std4_mean=0.1*(std4_1+std4_2+std4_3+std4_4+std4_5+std4_6+std4_7+std4_8+std4_9+std4_10)
print(std4_mean)

std3_1=dim_stat4(t)
std3_2=dim_stat4(t1)
std3_3=dim_stat4(t2)
std3_4=dim_stat4(t3)
std3_5=dim_stat4(t4)
std3_6=dim_stat4(t5)
std3_7=dim_stat4(t6)
std3_8=dim_stat4(t7)
std3_9=dim_stat4(t8)
std3_10=dim_stat4(t9)
std3_mean=0.1*(std3_1+std3_2+std3_3+std3_4+std3_5+std3_6+std3_7+std3_8+std3_9+std3_10)
print(std3_mean)



res=True
res@gsnHistogramBinIntervals = (/-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
res@tiXAxisOn=False
res@tmYRMode = "Automatic"
res@gsnHistogramComputePercentages = True
res@tmYLMode="explicit"
res@tmYLValues=(/200,400,600,800,1000,1200,1400,1600,1800/)
res@tmYLLabels=(/"200","400","600","800","1000","1200","1400","1600","1800"/)
res@tmYROn=False
res@tmYRLabelsOn = False
res@tmYUseRight = False  
res@tmLabelAutoStride =True




test@units="m/s"
test1@units="m/s"
test2@units="m/s"
test3@units="m/s"
test4@units="m/s"
test5@units="m/s"
test6@units="m/s"
test7@units="m/s"
test8@units="m/s"
test9@units="m/s"

plot(0)=gsn_histogram(wks,test,res)
plot(1)=gsn_histogram(wks,test1,res)
plot(2)=gsn_histogram(wks,test2,res)
plot(3)=gsn_histogram(wks,test3,res)
plot(4)=gsn_histogram(wks,test4,res)
plot(5)=gsn_histogram(wks,test5,res)
plot(6)=gsn_histogram(wks,test6,res)
plot(7)=gsn_histogram(wks,test7,res)
plot(8)=gsn_histogram(wks,test8,res)
plot(9)=gsn_histogram(wks,test9,res)
plot(10)=gsn_histogram(wks,t,res)
plot(11)=gsn_histogram(wks,t1,res)
plot(12)=gsn_histogram(wks,t2,res)
plot(13)=gsn_histogram(wks,t3,res)
plot(14)=gsn_histogram(wks,t4,res)
plot(15)=gsn_histogram(wks,t5,res)
plot(16)=gsn_histogram(wks,t6,res)
plot(17)=gsn_histogram(wks,t7,res)
plot(18)=gsn_histogram(wks,t8,res)
plot(19)=gsn_histogram(wks,t9,res)

res1=True
res1@draw=False
res1@frame=False

res1@txString = "CanCM4 (top 10)/CanCM3 (bottom 10) 135E-240E averaged 5S-5N 850 hPa zonal wind anomalies (m/s)"
res1@txFontHeightF    = 0.015
gsn_panel(wks,plot,(/4,5/),res1)
delete(wks)
delete(plot)
end
