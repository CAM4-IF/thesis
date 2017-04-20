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
colors=(/"blue","green","red"/)
colors1=(/"black","orange","cyan"/)
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
   Mprime=asciiread("Mprime.txt",(/301,3/),"float")
   
   res               = True	                   ; plot mods desired
   res@tiMainString  = "(~F33~D~F~M~B~int~N~)/M~B~int~N~ and (~F33~D~F~M~S~`~N~)/M~S~`~N~ "     ; title
   res@tiMainFontHeightF=0.02
   res@tiXAxisString = "Month"  ; xaxis string
   res@tiYAxisString = "Fractional change"                  ; yaxis string

   

; add additional axis on top of plot
   res@tmXUseBottom  = False      ; Keep top axis independent of bottom.
   res@tmXTLabelsOn  = False      ; have tick mark labels
   res@tmXTOn        = True       ; have tick marks
   res@tmXBLabelFontHeightF=0.01
   res@tmYLLabelFontHeightF=0.01
   res@gsnDraw = False
   res@gsnFrame = False
   res@trXMaxF  =  300
   res@trYMaxF= 0.35
   res@trYMinF=-0.20
   res@vpHeightF 	 = 0.55              ; change aspect ratio of plot
  res@vpWidthF 	         = 0.7
   res@xyLineColors =colors(n)
   timeaxis=fspan(0,300,301)
   res@xyDashPattern  = 0    
   res@xyLineThicknesses=2
   res@gsnYRefLine = 0.0 
   res@tiYAxisFontHeightF=0.02
   res@tiXAxisFontHeightF=0.02
   res@tiXAxisFont="times-roman"
   res@tiYAxisFont="times-roman"
   plot(n)=gsn_csm_xy(wks,timeaxis,response(:,0,0),res)     ; create plot 
   delete(res@xyLineThicknesses)
   res@xyLineThicknesses=0.25
   delete(res@xyLineThicknesses)
   ;res@xyDashPattern  = 2
   res@xyLineColors =colors1(n)
   plot1=gsn_csm_xy(wks,timeaxis,Mprime(:,n),res) 
   
   
   rc=regline_stats(response(:,0,0),Mprime(:,n))
   print(rc)
   overlay(plot(n),plot1)
   gres=True
  gres@YPosPercent = 95
gres@XPosPercent = 5
gres@ItemSpacePercent = 5.  
lineres = True
lineres@lgLineColors = (/"red","green","blue"/) ; line colors
lineres@lgDashIndexes=(/0,0,0/)
lineres@lgLineThicknesses = 2.5                        ; line thicknesses
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
lineres@lgLabelFontHeightF=0.008

textres = True
textres@lgLabels = (/"CAM4-IF best T M~B~int~N~","CAM4-IF best rain M~B~int~N~","Default CAM4 M~B~int~N~"/)  ; legend labels (required)
textres@lgLabelFontHeights = (/0.015,0.015,0.015/)   
plot1= simple_legend(wks,plot1,gres,lineres,textres)

gres1=True
gres1@YPosPercent = 95
gres1@XPosPercent = 50
gres1@ItemSpacePercent = 5.  
lineres1 = True
lineres1@lgLineColors = (/"black","orange","cyan"/) ; line colors
lineres1@lgDashIndexes=(/0,0,0/)
lineres1@lgLineThicknesses = 2.5                        ; line thicknesses
lineres1@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
lineres1@lgLabelFontHeightF=0.008

textres1 = True
textres1@lgLabels = (/"CAM4-IF best T M~S~`~N~","CAM4-IF best rain M~S~`~N~","Default CAM4 M~S~`~N~"/)  ; legend labels (required)
textres1@lgLabelFontHeights = (/0.015,0.015,0.015/)   
plot1= simple_legend(wks,plot1,gres1,lineres1,textres1)
  draw(plot)
   
   delete(T)
   delete(hyam)
   delete(hybm)
   delete(rc)
   n=n+1
end do
delete(wks) 



clean=systemfunc("rm tmp*.nc")
crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
print(fils)
print(fils1)
end
