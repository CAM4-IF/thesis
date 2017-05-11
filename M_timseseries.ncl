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
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*_timeseries*.nc1*")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*_timeseries*.nc2*")
n=0
lets=(/"a)","b)","c)","d)","e)","f)"/)
plotfile="mpq_scat"
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods, graphic) 
colors=(/"orange","black","red"/)
colors1=(/"cyan","orange","black"/)
datalog=new((/3,301/),"double")
f3= addfile("../mass_flux/CAM-IFv5_F2000_ANN_climo.nc","r")
cmd3="cdo -s -fldmean -sellonlatbox,0,360,-30,30 -selname,PS,hyam,hybm ../mass_flux/CAM-IFv5_F2000_ANN_climo.nc tmp3.nc"
cmd3=systemfunc(cmd3)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   print(fils)
   print(fils1)
   hyam = f3->hyam 
   hybm =f3->hybm
   cmd="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS,FMASS_DN " + filen + " tmp.nc"
   cmd1="cdo -fldmean -sellonlatbox,0,360,-30,30 -selname,FMASS,FMASS_DN " + filen1 + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   ;areamean2=systemfunc(cmd2)
   f    = addfile("tmp.nc", "r")
   f1    = addfile("tmp1.nc", "r")
   f2= addfile("tmp3.nc","r")
   T=f->FMASS(0:300,:,:,:)
   T1=f1->FMASS(0:300,:,:,:)
   if (n.gt.0) then
   T2=f->FMASS_DN(0:300,:,:,:)
   T3=f1->FMASS_DN(0:300,:,:,:)
   Mup1=T-T2
   Mup2=T1-T3
   Mupresponse=Mup2-Mup1
   else
   Mup1=T
   Mup2=T1
   Mupresponse=Mup2-Mup1
   end if
   
  
   
   
   
   PS=f2->PS
   P0   = 100000.     
   levs=conform_dims(dimsizes(T),T&ilev,1)
   PS1=PS(0,0,0)
   PS2=conform(T,PS1,-1)
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, Mup1, PS2(:,0,:,:)/100, 0, opt)/900
   wva1=wgt_vert_avg_beta(levs,Mupresponse, PS2(:,0,:,:)/100, 0, opt)/900
   
   print(n)
   
   datalog(n,:)=wva1(:,0,0)/wva(:,0,0)
   Mprime=asciiread("Mprime.txt",(/301,3/),"float")
   Mprime1=asciiread("Mprime2.txt",(/301,3/),"float")
   
   res               = True	                   ; plot mods desired
   ;res@tiMainString  = "(~F33~D~F~M~B~int~N~)/M~B~int~N~ and (~F33~D~F~M~S~`~N~)/M~S~`~N~ 3-month running mean"     ; title
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
   res@trXMaxF  =  0.25
   res@trXMinF  =  -0.25
   res@trYMaxF= 0.25
   res@trYMinF=0
   res@tiYAxisString = "~F33~d~F~~F34~a~F~M~B~int~N~~F34~q~F~/~F34~a~F~M~B~int~N~~F34~q~F~"
   res@tiXAxisString = "~F33~d~F~~F34~a~F~(M~B~u~N~)~S~`~N~~F34~q~F~/~F34~a~F~(M~B~u~N~)~S~`~N~~F34~q~F~"
   res@xyMarkerSizeF     = 0.002
   res@xyMarkLineMode    = "Markers"                ; choose to use markers
   res@xyMarkers         =  16                      ; choose type of marker  
   res@xyMarkerColor     = colors(n)
   res@vpHeightF 	 = 0.55              ; change aspect ratio of plot
   res@vpWidthF 	         = 0.7
   res@tiYAxisFontHeightF=0.02
   res@tiXAxisFontHeightF=0.02
   
   delete(T)
   delete(hyam)
   delete(hybm)
   n=n+1
end do
print(datalog)
x    = fspan(-1,1,1000)
y= fspan(-1,1,1000)

res@trXMinF  =  -0.1
res@xyMarkerColor     = colors(0)
 res@gsnLeftString= "a)"
   res@gsnLeftStringFontHeightF=0.025
cor=escorc(Mprime(:,1),datalog(1,:))
rc=regline(Mprime(:,1), datalog(1,:))
yreg=rc*Mprime(:,1)+rc@yintercept
res1=True
res1@xyMarkLineMode="Lines"
res1@xyLineColor="Orange"
res1@gsnDraw=False
res1@gsnFrame=False
plot5=gsn_csm_xy(wks,Mprime(:,1),yreg,res1)
plot(0)=gsn_csm_xy(wks,Mprime(:,1), datalog(1,:),res)     ; create plot

lineres=True 
  lnid = gsn_add_polyline(wks,plot(0),x,y,lineres)
;txid1 = new(5,graphic)                  ; Nine text strings
; txres               = True     
; txres@txFontHeightF = 0.018
; txid1(0) = gsn_add_text (wks,plot(0),"CAM4-IF-t r: " + sprintf("%5.2f",cor), 0.15,0.04, txres)
res@xyMarkerColor     = colors(1)
cor1=escorc(Mprime(:,2),datalog(2,:))

plot1=gsn_csm_xy(wks,Mprime(:,2), datalog(2,:),res)     ; create plot 
gres=True
gres@YPosPercent =25. 
gres@XPosPercent = 55
gres@ItemSpacePercent = 5.  
lineres1 = True
lineres1@lgLineColors = (/"Black","Orange"/) ; line colors
lineres1@lgDashIndexes=(/0,0,0,0,0,0/)
lineres1@lgLineThicknesses = 2.5                        ; line thicknesses
lineres1@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

textres1 = True
textres1@lgLabels = (/"CAM4-IF-r, r: " + sprintf("%5.2f",cor),"CAM4-IF-t, r: " + sprintf("%5.2f",cor1) /)  ; legend labels (required)
textres1@lgLabelFontHeights = (/0.018,0.018,0.018/)   
plot(0) = simple_legend(wks,plot(0),gres,lineres1,textres1)

rc1=regline(Mprime(:,2), datalog(2,:))
yreg1=rc1*Mprime(:,2)+rc1@yintercept
res1@xyLineColor="Black"
plot6=gsn_csm_xy(wks,Mprime(:,2),yreg1,res1)
overlay(plot(0),plot1)
overlay(plot(0),plot5)
overlay(plot(0),plot6)
delete(res@trYMinF)
 res@gsnLeftString= "b)"
res@trYMinF=-0.2
res@trXMaxF  =  0.1
res@trXMinF  =  -0.25
res@tiXAxisString = "~F33~d~F~~F34~a~F~M~S~`~N~~F34~q~F~/~F34~a~F~M~S~`~N~~F34~q~F~"
res@xyMarkerColor     = colors1(0)
cor2=escorc(Mprime1(:,0),datalog(0,:))
plot(1)=gsn_csm_xy(wks,Mprime1(:,0), datalog(0,:),res)     ; create plot

rc2=regline(Mprime1(:,0), datalog(0,:))
yreg2=rc2*Mprime(:,0)+rc2@yintercept
res1@xyLineColor="Cyan"
plot7=gsn_csm_xy(wks,Mprime(:,0),yreg2,res1)

lnid1 = gsn_add_polyline(wks,plot(1),x,y,lineres)
res@xyMarkerColor     = colors1(1)
cor3=escorc(Mprime1(:,2),datalog(1,:))
plot2=gsn_csm_xy(wks,Mprime1(:,2), datalog(1,:),res)     ; create plot 

rc3=regline(Mprime1(:,2), datalog(1,:))
yreg3=rc3*Mprime1(:,2)+rc3@yintercept
res1@xyLineColor="Black"
plot8=gsn_csm_xy(wks,Mprime1(:,2),yreg3,res1)


overlay(plot(1),plot2)
res@xyMarkerColor     = colors1(2)
cor4=escorc(Mprime1(:,1),datalog(2,:))
plot3=gsn_csm_xy(wks,Mprime1(:,1), datalog(2,:),res)     ; create plot 

rc4=regline(Mprime1(:,1), datalog(2,:))
yreg4=rc4*Mprime1(:,1)+rc4@yintercept
res1@xyLineColor="Orange"
plot9=gsn_csm_xy(wks,Mprime1(:,1),yreg4,res1)
delete(gres)
gres=True
gres@YPosPercent =15. 
gres@XPosPercent = 45.
gres@ItemSpacePercent = 5.  
delete(lineres1)
lineres1 = True
lineres1@lgLineColors = (/"Cyan","Black","Orange"/) ; line colors
lineres1@lgDashIndexes=(/0,0,0,0,0,0/)
lineres1@lgLineThicknesses = 2.5                        ; line thicknesses
lineres1@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
delete(textres1)
textres1 = True
textres1@lgLabels = (/"CAM4 default, r: " + sprintf("%5.2f",cor2), "CAM4-IF-r, r: " + sprintf("%5.2f",cor3),"CAM4-IF-t, r: " + sprintf("%5.2f",cor4) /)  ; legend labels (required)
textres1@lgLabelFontHeights = (/0.018,0.018,0.018/)   
plot(1) = simple_legend(wks,plot(1),gres,lineres1,textres1)
overlay(plot(1),plot3)
overlay(plot(1),plot7)
overlay(plot(1),plot8)
overlay(plot(1),plot9)
resP                            = True
resP@txString = ""
resP@txFontHeightF    = 0.03
resP@gsnMaximize=True
resP@gsnPaperOrientation = "portrait"
gsn_panel(wks,plot,(/1,2/),resP)  
 delete(wks) 


clean=systemfunc("rm tmp*.nc")
crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
print(fils)
print(fils1)
end
