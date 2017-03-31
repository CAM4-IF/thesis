load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*********************************************
begin
;**********************************************
lev_p           = (/30.,50.,70.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,875.,900.,925.,950.,975.,1000./)           
lev_p!0         = "lev_p"                  ; variable/dim name 
lev_p&lev_p     =  lev_p                   ; create coordinate variable
lev_p@long_name = "pressure"               ; attach some attributes
lev_p@units     = "hPa"
lev_p@positive  = "down"
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*DJF*.nc2")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/*DJF*.nc1")
datalog=new((/dimsizes(fils),3/),float) ; create list to store output
plotfile="v4_DJF_mpsi" ;;;PDF file name!
wks = gsn_open_wks("pdf", plotfile) ;; open a workstation file!
nummods=dimsizes(fils) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
  f=addfile(fils(n),"r") ; add file
  f1=addfile(fils1(n),"r") ; add file
  hyam            = f->hyam         ; read from file
  hybm            = f->hybm
  P0mb            = 1000
  ps              = f->PS(0,{-45:45},:)       ; these are in Pa
  v               = f->V(0,:,{-45:45},:)
  lat             = f->lat({-45:45})
  lev= f->lev
  hyam1           = f1->hyam         ; read from file
  hybm1            = f1->hybm
  P0mb1            = 1000
  ps1              = f1->PS(0,{-45:45},:)       ; these are in Pa
  v1               = f1->V(0,:,{-45:45},:)
  lat1             = f1->lat({-45:45})
  lev1= f1->lev
; convert to pressure levels
  vp   = vinth2p(v,hyam,hybm,lev_p,ps,2,P0mb,1,False)
  vp1   = vinth2p(v1,hyam1,hybm1,lev_p,ps1,2,P0mb1,1,False)
; calculate zonal meridional streamfunction
  zmpsi = zonal_mpsi(vp,lat,lev_p*100.,ps)
  zmpsi!0="lev"
  zmpsi!1="lat"
  zmpsi&lev=lev_p
  zmpsi&lat=lat
  zmpsi1 = zonal_mpsi(vp1,lat1,lev_p*100.,ps1)
  zmpsi1!0="lev"
  zmpsi1!1="lat"
  zmpsi1&lev=lev_p
  zmpsi1&lat=lat1
  zmpsical=zmpsi/1.0e9
  zmpsical!0="lev"
  zmpsical!1="lat"
  zmpsical&lev=lev_p
  zmpsical&lat=lat
  zmpsical1=zmpsi1/1.0e9
  zmpsical1!0="lev"
  zmpsical1!1="lat"
  zmpsical1&lev=lev_p
  zmpsical1&lat=lat
  response=zmpsical1-zmpsical
  fiven=round(0.55555*dimsizes(zmpsi(1,:)),3)
  fifteens=round(0.33333*dimsizes(zmpsi(1,:)),3)
  slice1=zmpsi({200:850},fiven)
  slice2=zmpsi({200:850},fifteens)
  mean1=dim_avg_n(slice1,0)
  mean2=dim_avg_n(slice2,0)
  diff=mean2-mean1
  print(diff)
  print(fils(n))
 
 ;;;SET THE COLOR FILL PLOT RESOURCES!;;;;;;
  
  response!0="lev"
  response!1="lat"
  response&lev=lev_p
  response&lat=lat
  res                      = True                 ; plot mods desired
  res@tmYRMode = "Automatic"

  res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
  res@cnLevels = (/-30,-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20,30/)
  res@cnFillColors = (/2,6,12,18,24,40,60,90,106,116,128,130,142,148,156,168,180,208,218,230,250,255/) 
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@lbLabelBarOn = False
  res@lbLabelAutoStride    = True                 ; optimal labels
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                ; turn off contour info label
  res@tiXAxisString    = ""
  res@tiYAxisSide = "Left"	
  res@tiYAxisString    = "Pressure (hPa)"
  res@trYMaxF = 1000.0
  res@trYMinF = 100.0      ; troposphere only
  res@tiXAxisFontHeightF =0.012
  res@tmYROn = False
  res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7
  

  mpsimax=max(zmpsical)
  mpsimaxresponse=((max(zmpsical1)-max(zmpsical))/max(zmpsical))*100
  ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
  cres                      = True                 ; plot mods desired
  cres@tmYRMode = "Automatic"
  cres@cnLevelSelectionMode = "ExplicitLevels"
  cres@cnLevels = (/-260,-240,-220,-200,-180,160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240,260/)
  cres@cnLineColors = (/200,205,210,215,220/)
  cres@cnLineLabelsOn       = True                 ; turn on line labels
  cres@cnFillOn             = False                 ; turn on color fill
  cres@gsnDraw = False
  cres@gsnFrame = False
  cres@cnInfoLabelOn = False                ; turn off contour info label
  cres@gsnRightString = ""
  cres@gsnLeftString = ""
  cres@tiYAxisString    = "Pressure (hPa)"
  cres@tiXAxisString    = "Latitude"
  cres@tmYROn = False
  cres@tmYRLabelsOn = False
  cres@gsnContourZeroLineThicknessF  = 0.0
  cres@cnLineLabelsOn = True
  cres@cnLineThicknessF =0.5
  cres@gsnContourNegLineDashPattern =1  
  cres@gsnContourZeroLineThicknessF = 2
  cres@gsnRightString    = f1@case
  cres@gsnLeftString  = "Climo max:" + mpsimax + "~C~response (%): " + mpsimaxresponse
 
  cres1=cres
  cres1@cnLineColor="red"
  cres1@gsnContourPosLineDashPattern=1
  plot(n)  = gsn_csm_pres_hgt(wks, response, res)  
  ;plot(n)  = gsn_csm_pres_hgt(wks, mpsiscal(::-1,:), cres) ;;make the plot! for model "n"
  plot2  = gsn_csm_pres_hgt(wks, zmpsical, cres)
  overlay(plot(n),plot2)
  n=n+1
  print(n)
  delete(zmpsi)
  delete(zmpsi1)
  delete(zmpsical)
  delete(zmpsical1)
  delete(response)
  delete(mpsimax)
end do
;print(datalog(:,2))
print(fils)
print(fils1)
;;;;build the panel plot!;;;;
 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "CAM4/CAM-IF DJF mass streamfunction response with control run climatology (1999-2002)"
 resP@lbTitleString    = "10~S~9~N~ kg s~S~-1~N~ K~S~-1~N~"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.01
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/5,5/),resP)  
 delete(wks) 
crop=systemfunc("/usr/bin/pdfcrop v4_DJF_mpsi.pdf")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
end
