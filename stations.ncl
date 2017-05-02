begin
 plotfile="stations"
 wks = gsn_open_wks("epsi",plotfile)      
 res             = True

  res@gsnMaximize            = True
  res@gsnDraw                = True
  res@gsnFrame               = False

  res@mpOutlineOn            = True
  res@mpPerimOn              = False

  res@mpGridAndLimbOn        = True
  res@mpGridLineDashPattern  = 2
  res@mpGridMaskMode         = "MaskLand"

  res@mpLandFillColor        = "tan"
  res@mpOceanFillColor       = "LightBlue"
  res@mpMinLatF    = -25.
  res@mpMaxLatF    = 25.
  res@mpCenterLonF         = -200
  lat=(/13.55,19.72,18.43,7.33,21.98,7.08,-14.33,6.97,7.47/)
  lon=(/144.8,204.93,294,134.48,200.65,171.38,189.28,158.22,151.85/)
  map = gsn_csm_map(wks,res)
  gsres               = True
  gsres@gsMarkerIndex = 16  
  do i=0,8
   gsn_polymarker(wks,map,lon(i),lat(i),gsres)
  end do
  ;draw(map)
  frame(wks)
string1="eps2eps " + plotfile +".epsi" + " " + plotfile+".eps"
system(string1) 
system("epstopdf " + plotfile +".eps") 
system("rm " + plotfile +".eps") 
system("rm " + plotfile +".epsi")
end