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
   
   
   
   
   
   
   res               = True	                   ; plot mods desired
   res@tiMainString  = "Mass flux time series"     ; title
   res@tiXAxisString = "Month"  ; xaxis string
   res@tiYAxisString = "FMASS (kg/m^2/s)"                  ; yaxis string
   res@xyExplicitLegendLabels = f@case
   res@lgLabelFontHeightF     = .01 
   res@lgOrientation          = "horizontal"
   res@lgPerimOn              = False               ; no box around
   res@pmLegendDisplayMode    = "Always"
   res@pmLegendSide           = "Right"               ; Change location of 
  ; res@pmLegendParallelPosF   = .50                 ; move units right
   res@pmLegendOrthogonalPosF = 0              ; move units down

  res@pmLegendWidthF         = 0.12                ; Change width and
  res@pmLegendHeightF        = 0.05 

; add additional axis on top of plot
   res@tmXUseBottom  = False      ; Keep top axis independent of bottom.
   res@tmXTLabelsOn  = True       ; have tick mark labels
   res@tmXTOn        = True       ; have tick marks
   res@gsnDraw = True
   res@gsnFrame = False
   res@trXMaxF  =  300
   res@trYMaxF= 0.15
   res@trYMinF=-0.15
   res@xyLineColors =colors(n)
   timeaxis=fspan(0,300,301)
   plot(n)=gsn_csm_xy(wks,timeaxis,T(:,18,0,0),res)     ; create plot 
   clean=systemfunc("rm tmp.nc")
   delete(T)
   delete(hyam)
   delete(hybm)
   
   n=n+1
end do
;overlay(plot(0),plot(1))
delete(wks) 
end