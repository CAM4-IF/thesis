;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;************************************************
begin
plotfile="M_scat"
wks = gsn_open_wks("pdf",plotfile)              ; open workstation
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*timeseries*.nc1*")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/mass_flux/*timeseries*.nc2*")
n=0
lets=(/"a)","b)","c)","d)"/)
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods, graphic) 
f3= addfile("../mass_flux/CAM-IFv5_F2000_ANN_climo.nc","r")
cmd3="cdo -s -fldmean -sellonlatbox,0,360,-30,30 -selname,PS,hyam,hybm ../mass_flux/CAM-IFv5_F2000_ANN_climo.nc tmp3.nc"
cmd3=systemfunc(cmd3)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   hyam = f3->hyam 
   hybm =f3->hybm
   cmd="cdo -s -timmean -sellonlatbox,0,360,-30,30 -selname,FMASS " + filen + " tmp.nc"
   cmd1="cdo -s -timmean -sub -sellonlatbox,0,360,-30,30 -selname,FMASS " + filen1 + " tmp.nc" + " tmp1.nc"
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f1    = addfile("tmp1.nc", "r")
   f    = addfile("tmp.nc", "r")
   T=f->FMASS(0,:,:,:)
   T1=f1->FMASS
   PS=f3->PS
   P0   = 100000.     
   levs=conform_dims(dimsizes(T),T&ilev,0)
   PS1=PS(0,0,0)
   PS2=conform(T,PS1,-1)
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(levs, T, PS2(0,:,:)/100, 0, opt)/900   ;;; control
   wva1=wgt_vert_avg_beta(levs(:,0,0),T1(0,:,:,:), PS2(0,:,:)/100, 0, opt)/900 ;;;+4K 
   wva2=wva1/wva
   
   x=ndtooned(wva1)
   y=ndtooned(wva)
  res                   = True                     ; plot mods desired
  res@gsnMaximize       = True                     ; maximize plot
  res@tiMainString      = ""           ; add title
  res@xyMarkLineMode    = "Markers"                ; choose to use markers
  res@xyMarkers         =  16                      ; choose type of marker  
  res@xyMarkerColor     = "Black"               ; Marker color
  res@xyMarkerSizeF     = 0.001                     ; Marker size (default 0.01)
  res@gsnDraw=False
  res@gsnFrame=False
  res@trXMinF=min(wva)
  res@trXMaxF=0.0003
  res@trYMinF=-0.004
  res@trYMaxF=0.004
  res@gsnCenterString  =lets(n) + " " + f@case
  res@tiYAxisString = "~F33~D~F~(M~B~int~N~)"
  rc=regline(y,x)
  yreg=rc*y+rc@yintercept
  plot(n)  = gsn_csm_xy (wks,y,x,res)                    ; create plot
  res@xyMarkLineMode="Lines"
  plot1=gsn_csm_xy(wks,y,yreg,res)
  overlay(plot(n),plot1)
  
  n=n+1
end do
resP                            = True
 resP@txString = ""
 resP@txFontHeightF    = 0.03
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 gsn_panel(wks,plot,(/1,3/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop " + plotfile) 
cleanup=systemfunc("rm tmp*.nc") 
print(fils)
print(fils1)
end