;*************************************************
; xy_8.ncl
;
; Concepts illustrated:
;   - Forcing tickmarks and labels to be drawn on the top X axis in an XY plot
;   - Changing the line dash pattern in an XY plot
;   - Explicitly setting tickmarks and labels on the top X axis
;
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"    
;************************************************
begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*timeseries*.nc1*")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*timeseries*.nc2*")
n=0
lets=(/"a)","b)","c)","d)","e)","f)"/)
plotfile="test"
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods, graphic) 
colors=(/"red","green","blue"/)
f3= addfile("../mass_flux/CAM-IFv5_F2000_ANN_climo.nc","r")
cmd3="cdo -s -fldmean -sellonlatbox,0,360,-30,30 -selname,PS,hyam,hybm ../mass_flux/CAM-IFv5_F2000_ANN_climo.nc tmp3.nc"
cmd3=systemfunc(cmd3)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   hyam = f3->hyam 
   hybm =f3->hybm
   cmd="cdo -s -timmean -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS " + filen + " tmp.nc"
   cmd1="cdo -s -sub -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS " + filen1 + " tmp.nc" + " tmp1.nc"
   ;cmd2="cdo -s -div tmp1.nc tmp.nc tmp2.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   ;areamean2=systemfunc(cmd2)
   f    = addfile("tmp1.nc", "r")
   f1    = addfile("tmp.nc", "r")
   f2= addfile("tmp3.nc","r")
   T=f->FMASS(0:300,:,:,:)
   T1=f1->FMASS
   PS=f2->PS
   P0   = 100000.     
   levs=conform_dims(dimsizes(T),T&ilev,1)
   PS1=PS(0,0,0)
   PS2=conform(T,PS1,-1)
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, T, PS2(:,0,:,:)/100, 0, opt)/900
   wva1=wgt_vert_avg_beta(levs(0,:,:,:),T1(0,:,:,:), PS2(0,0,:,:)/100, 0, opt)/900
   response=wva/wva1(0,0)
   print(response)
   Mprime=asciiread("Mprime.txt",(/301,3/),"float")
   
   res               = True	                   ; plot mods desired
   res@tiMainString  = "Mass flux time series"     ; title
   res@tiXAxisString = "Month"  ; xaxis string
   res@tiYAxisString = "~F33~D~F~(FMASS)/FMASS (kg/m^2/s)"                  ; yaxis string

   

; add additional axis on top of plot
   res@tmXUseBottom  = False      ; Keep top axis independent of bottom.
   res@tmXTLabelsOn  = True       ; have tick mark labels
   res@tmXTOn        = True       ; have tick marks
   res@gsnDraw = False
   res@gsnFrame = False
   res@trXMaxF  =  300
   res@trYMaxF= 0.30
   res@trYMinF=-0.20
   res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
  res@vpWidthF 	         = 0.7
   res@xyLineColors =colors(n)
   timeaxis=fspan(0,300,301)
   res@xyDashPattern  = 0    
   res@xyLineThicknesses=3
   res@gsnYRefLine = 0.0 
   plot(n)=gsn_csm_xy(wks,timeaxis,response(:,0,0),res)     ; create plot 
   delete(res@xyLineThicknesses)
   res@xyLineThicknesses=0.25
   delete(res@xyLineThicknesses)
   ;res@xyDashPattern  = 2
   plot1=gsn_csm_xy(wks,timeaxis,Mprime(:,n),res) 
   overlay(plot(n),plot1)
   txid = new(3,graphic)                  ; Nine text strings
txres               = True     
txres@txFontHeightF = 0.01
txid(0) = gsn_add_text (wks,plot(0),"Red: CAM4-IF best rainfall", 150, .29 ,txres)
txid(1) = gsn_add_text (wks,plot(0),"Green: CAM4-IF best T"     , 146, .27,txres)
txid(2) = gsn_add_text (wks,plot(0),"Blue: CAM4"       , 132, .25 ,txres)
  draw(plot)
   clean=systemfunc("rm tmp*.nc")
   delete(T)
   delete(hyam)
   delete(hybm)
   n=n+1
end do
delete(wks) 




crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
print(fils)
print(fils1)
end