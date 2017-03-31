;----------------------------------------------------------------------
; h_lat_6.ncl
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;----------------------------------------------------------------------
begin
  f     = addfile ("T_CAM4_ANN.nc","r")
  f1  =addfile("T_CAM4_4K_ANN.nc","r")
  u     = f->T(0,11:25,0,:)
  u@units="K"
  u1= f1->T(0,11:25,0,:)
  u1@units="K"
  P0mb=1000
delete(u@long_name)
delete(u@units)
delete(u1@long_name)
delete(u1@units)
  lev           = u&lev*1000
  lev!0         = "lev"                  ; variable/dim name 
  lev&lev     =  lev                  ; create coordinate variable
  lev@long_name = "pressure"               ; attach some attributes
  lev@units     = "Pa"
  lev@positive  = "down"
  levconform=conform_dims((/15,144/),lev,0)
  levconform1=conform_dims((/15,144/),lev,0)
  s1   = static_stability(levconform, u, 0, 0)
  s2=static_stability(levconform1, u1, 0, 0)
  sdiff=s2-s1
  lev           = u&lev
  sdiff!0="lev"
  sdiff&lev=lev
  sdiff@long_name = "pressure" 
  sdiff@units="hPa"
  sdiff@positive="down"
  
  lon=u&lon
  sdiff!1="lon"
  sdiff&lon=lon
  sdiff@long_name = "longitude" 
  
  wks   = gsn_open_wks ("png", "CAM4_static_stability_response" )          ; send graphics to PNG file
  res                      = True                 ; plot mods desired
  res@gsnDraw=False
  res@gsnFrame=False
  res@cnLevelSelectionMode =  "ExplicitLevels"      ; manually select levels
  res@cnLevels    = (/0,0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,7,8,9/)*0.8e-5
  res@cnLineLabelsOn       = False                 ; turn on line labels 
  res@cnLinesOn       = True
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = "MPL_YlOrRd"             ; choose colormap
  res@tiYAxisString        = "Pressure - hPa" 
  res@lbTitleString="Static stability response (K/Pa)" 
  res@lbTitleFontHeightF=0.015
  res@cnInfoLabelOn = False
  res@cnLineThicknessF=0.2
  res@gsnContourZeroLineThicknessF=3
  res@tiMainString     = "CAM4 30S-30N ANN static stability (K/Pa)"   ; title
  res@trYReverse            = True  

  res2                      = True                 ; plot mods desired
  res2@gsnDraw=False
  res2@gsnFrame=False
  ;res2@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  ;res2@cnLevelSpacingF      = 1e-5                ; contour spacing
  ;res2@cnMinLevelValF       = 0                ; min level
  ;res2@cnMaxLevelValF       = 2e-4;
  res2@cnFillOn             = False               ; turn on color fill 
  res2@cnLineThicknessF=2
  res2@gsnContourNegLineDashPattern=1
  res2@gsnContourZeroLineThicknessF=3
  res2@cnInfoLabelOn = False
  res2@cnLineLabelsOn       = True
  res2@cnLineLabelFontHeightF=0.01
  res2@cnLineLabelInterval        = 1  
  res2@cnLineLabelBackgroundColor=-1
  res2@cnLineLabelDensityF=1
  res2@tiXAxisString        = "Static stability (K/Pa)" 
  
  delete(s1@long_name)
delete(s1@units)
  plot  = gsn_csm_pres_hgt(wks, sdiff, res )   ; place holder
  plot1=gsn_csm_pres_hgt(wks, s1, res2 )
  overlay(plot,plot1)
  draw(plot)
  frame(wks)
  mean=dim_avg_n(s1,1)
  mean1=dim_avg_n(mean,0)
  print(mean1)
  mean2=dim_avg_n(s2,1)
  mean3=dim_avg_n(mean2,0)
  print(mean3)
  
end
