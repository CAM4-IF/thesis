begin
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/*ANN*.nc")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/*ANN*.nc1")
datalog=new(dimsizes(fils),float) ; create list to store output
plotfile="test" ;;;PDF file name!
wks = gsn_open_wks("png", plotfile) ;; open a workstation file!
nummods=dimsizes(fils1) ;;number of models! 
plot=new(nummods, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
n=0
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   c2="cdo selname,velopot,stream -sp2gp -dv2ps -uv2dv -remapbil,n45 -sellevel,192.539935 " + filen + " chi.nc"
   c3="cdo selname,velopot,stream -sp2gp -dv2ps -uv2dv -remapbil,n45 -sellevel,192.539935 " + filen1 + " chi1.nc"
   cmd2=systemfunc(c2)
   cmd3=systemfunc(c3)
   c4="cdo sub -selname,velopot chi.nc -enlarge,chi.nc -zonmean -selname,velopot chi.nc chi2.nc" ; zonally anomalous chi
   c5="cdo sub -selname,velopot chi1.nc -enlarge,chi1.nc -zonmean -selname,velopot chi1.nc chi3.nc" ; zonally anomalous chi
   cmd4=systemfunc(c4)
   cmd5=systemfunc(c5)
   
   f2=addfile("chi2.nc", "r")
   f3=addfile("chi3.nc","r")
   
   wc=f2->velopot(0,0,:,:)
   wcstrength=min(wc({-20:20},:))
   wc1=f3->velopot(0,0,:,:)
   wcstrength1=min(wc1({-20:20},:))
   
   
   ind0D=ndtooned(wc)
   dsizes=dimsizes(wc)
   index0=ind_resolve(minind(ind0D),dsizes)
   print(index0)
   ind1D=ndtooned(wc1)
   dsizes1=dimsizes(wc1)
   index1=ind_resolve(minind(ind1D),dsizes)
   print(index1)
   print(wcstrength)
   print(wcstrength1)
   
   
   
   f= addfile(filen, "r")
   hyam = f->hyam 
   hybm = f->hybm
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,{-50:50},:)                         
   Q=f->DTCOND(0,:,{-50:50},:)
   Q1=f->QRL(0,:,{-50:50},:)
   Q2=f->QRS(0,:,{-50:50},:)
   Q3=f->DTV(0,:,{-50:50},:)
   T     = f->T(0,:,{-50:50},:)
   T@units="K"
   omeg=f->OMEGA(0,:,{-50:50},:)
   Qtot=Q+Q1+Q2+Q3
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   s1   = static_stability(pm, T, 0, 0) ; calculate column integrated Sp at each grid point
   opt=(/0,15000,92500/)
   wva = wgt_vert_avg_beta(pm, Q, PS, 0, opt)/(77500) ; vertically pressure thickness weighted sum of Q at each grid point
   wva1 = wgt_vert_avg_beta(pm, Q1, PS, 0, opt)/(77500)
   wva2 = wgt_vert_avg_beta(pm, Q2, PS, 0, opt)/(77500)
   wva3 = wgt_vert_avg_beta(pm, Q3, PS, 0, opt)/(77500)
   wva4 = wgt_vert_avg_beta(pm, Qtot, PS, 0, opt)/(77500)
   wva4!0="lev"
   wva4!1="lon"
   wva4&lon=T&lon
   omegv=wgt_vert_avg_beta(pm, omeg, PS, 0, opt) /(77500)
   omegv@_FillValue=-999
   s1v=wgt_vert_avg_beta(pm,s1,PS,0,opt)/(77500)
   s1v@_FillValue=-999
  
   s1v!0="lat"
   s1v&lat=Q1&lat
   s1v!1="lon"
   s1v&lon=Q1&lon
   omegv!0="lat"
   omegv&lat=Q1&lat
   omegv!1="lon"
   omegv&lon=Q1&lon
   spo=omeg*s1
   
   spomeg=wgt_vert_avg_beta(pm,spo,PS,0,opt)/(77500)
   spomeg!0="lat"
   spomeg&lat=Q&lat
   spomeg!1="lon"
   spomeg&lon=Q&lon
   wva4!0="lat"
   wva4&lat=Q&lat
   wva4!1="lon"
   wva4&lon=Q&lon
   
   f1= addfile(filen1, "r")
   hyam = f1->hyam 
   hybm = f1->hybm
   P0   = 100000.                        ; reference pressure [Pa]
   PS1   = f1->PS(0,{-50:50},:)                         
   Q4=f1->DTCOND(0,:,{-50:50},:)
   Q5=f1->QRL(0,:,{-50:50},:)
   Q6=f1->QRS(0,:,{-50:50},:)
   Q7=f1->DTV(0,:,{-50:50},:)
   T1     = f1->T(0,:,{-50:50},:)
   T1@units="K"
   omeg1=f1->OMEGA(0,:,{-50:50},:)
   Qtot1=Q4+Q5+Q6+Q7
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   s2   = static_stability(pm, T1, 0, 0) ; calculate column integrated Sp at each grid point
   opt=(/0,15000,92500/)
   wva5 = wgt_vert_avg_beta(pm, Q4, PS1, 0, opt)/(77500) ; vertically pressure thickness weighted sum of Q at each grid point
   wva6 = wgt_vert_avg_beta(pm, Q5, PS1, 0, opt)/(77500)
   wva7 = wgt_vert_avg_beta(pm, Q6, PS1, 0, opt)/(77500)
   wva8 = wgt_vert_avg_beta(pm, Q7, PS1, 0, opt)/(77500)
   wva9 = wgt_vert_avg_beta(pm, Qtot1, PS1, 0, opt)/(77500)
   wva9!0="lev"
   wva9!1="lon"
   wva9&lon=T&lon
   omegv1=wgt_vert_avg_beta(pm, omeg1, PS1, 0, opt)/(77500)
   omegv1@_FillValue=-999
  
   s2v=wgt_vert_avg_beta(pm,s2,PS1,0,opt)/(77500)
   s2v@_FillValue=-999
   spo1=omeg1*s2
   spomeg1=wgt_vert_avg_beta(pm,spo1,PS,0,opt)/(77500)
   spomeg1!0="lat"
   spomeg1&lat=Q1&lat
   spomeg1!1="lon"
   spomeg1&lon=Q1&lon
   wva9!0="lat"
   wva9&lat=Q1&lat
   wva9!1="lon"
   wva9&lon=Q1&lon
   
   pcor=pattern_cor(spomeg({-10:10},:),wva4({-10:10},:),1.0,0)
   pcor1=pattern_cor(spomeg1({-10:10},:),wva9({-10:10},:),1.0,0)
   spomegdiff=spomeg1-spomeg
   spomegdiff!0="lat"
   spomegdiff&lat=Q1&lat
   spomegdiff!1="lon"
   spomegdiff&lon=Q1&lon
   qdiff=wva9-wva4
   qdiff!0="lat"
   qdiff&lat=Q1&lat
   qdiff!1="lon"
   qdiff&lon=Q1&lon
   pcor2=pattern_cor(spomegdiff({-10:10},:),qdiff({-10:10},:),1.0,0)
   pm!1="lat"
   pm&lat=Q1&lat
   pm!2="lon"
   pm&lon=Q1&lon
   ;;;;;;;;;;;;CALCULATE STABILITY CHANGE TERMS (LOCAL AND TROPICAL MEAN);;;;;;;;;;;;;;;;;;;;
   ;;;;TERM 1 - COLUMN MEAN MASC (30S-30N) - tropical mean Sp climate change minus control Sp * local omega;;;;;
   splon=dim_avg_n_Wrap(s1(:,{-30:30},:),2)
   splat=dim_avg_n_Wrap(splon,1)
   splon1=dim_avg_n_Wrap(s2(:,{-30:30},:),2)
   splat1=dim_avg_n_Wrap(splon1,1)
   sptropmeandiff=splat1-splat ;;this gives you a tropical mean profile of Sp change
   ;;;Now need to compute MASC at each gridpoint using vertical profile of local omega;;;;
   conformdims=dimsizes(omeg(:,{-50:50},:))
   spconform=conform_dims(conformdims,sptropmeandiff,0)
   masc=spconform*omeg(:,{-50:50},:)
   masc!0="lev"
   masc&lev=Q1&lev
   masc!1="lat"
   lats=f->lat({-50:50})
   masc&lat=lats
   masc!2="lon"
   masc&lon=Q1&lon
   
   masccolumn=wgt_vert_avg_beta(pm(:,{-50:50},:), masc, PS1({-50:50},:), 0, opt)/(77500)
   masccolumn!0="lat"
   masccolumn&lat=Q1&lat
   masccolumn!1="lon"
   masccolumn&lon=Q1&lon
   masctropmean=dim_avg_n_Wrap(masccolumn({-10:10},:),1)
   masctropmean1=dim_avg_n_Wrap(masctropmean,0)
   
   ;;;;;;;;;;;;;;;CALCULATE CIRCULATION CHANGE (CC) TERM;;;;;;;;;;;;;;;;;;;;
   omegdiff=(omegv1-omegv)
   omegdiff!0="lat"
   omegdiff&lat=Q1&lat
   omegdiff!1="lon"
   omegdiff&lon=Q1&lon
   
   omegdiffsp=omegdiff*s1v  ;;;;;;column mean omega change * Sp climo
   omegdiffsp!0="lat"
   omegdiffsp&lat=Q1&lat
   omegdiffsp!1="lon"
   omegdiffsp&lon=Q1&lon
   
   ;;;;;;;;;;;;;;;PATTERN OF DIABATIC HEATING CHANGE TERM;;;;;;;
   qstar=qdiff-avg(qdiff({-30:30},:))
   qstar!0="lat"
   qstar&lat=Q1&lat
   qstar!1="lon"
   qstar&lon=Q1&lon
   
   qstarmean=avg(qstar({-10:10},:))
   ;;;;;;;;;;;;;;CALCULATE ZONALLY ANOMALOUS Q (FOR GMS) FOR WC ASCENT REGION;;;;;;;;
   qzon=dim_rmvmean_n_Wrap(wva4,0)
   qzonar=dim_avg_n_Wrap(qzon({-10:10},{90:200}),1)
   qzonwcar=dim_avg_n_Wrap(qzonar,0)
   print(qzonwcar)
   GMS=qzonwcar/wcstrength
  
   print(GMS)
   
   qzon1=dim_rmvmean_n_Wrap(wva9,0)
   qzonar1=dim_avg_n_Wrap(qzon1({-10:10},{90:200}),1)
   qzonwcar1=dim_avg_n_Wrap(qzonar1,0)
   print(qzonwcar1)
   GMS1=qzonwcar1/wcstrength1
   GMSresponse=GMS1-GMS
   print(GMSresponse)
   
   
   
   res                      = True                 ; plot mods desired
   res@tmYRMode = "Automatic"

   res@cnLevelSelectionMode = "ExplicitLevels"       ; manually select levels
   res@cnLevels = (/-20,-15,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,15,20/)*5e-7
   ;res@cnLevels =(/-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18/)
   res@cnFillColors = (/2,6,12,24,40,60,90,100,0,0,0,0,156,170,190,210,220,230,240,250/)
   res@cnLinesOn       = False                 ; turn on line labels
   res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@lbLabelBarOn = False
  res@lbLabelAutoStride    = True                 ; optimal labels
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnInfoLabelOn = False                ; turn off contour info label
  ;res@tiXAxisString    = ""
  ;res@tiYAxisSide = "Left"	
  ;res@tiYAxisString    = "Pressure (hPa)"
  ;res@trYMaxF = 1000.0
  ;res@trYMinF = 100.0      ; troposphere only
  ;res@tiXAxisFontHeightF =0.012
  ;res@tmYROn = False
  ;res@tmYRLabelsOn = False 
  res@cnFillOpacityF=0.7
  res@cnLineThicknessF = 1.                   ; thicker lines
  res@cnLineLabelsOn   = False                
  res@cnInfoLabelOn = False                   ; turn on contour info label
  res@gsnSpreadColors     = False
  res@gsnMaximize      = True
  res@mpDataBaseVersion    = "MediumRes"  
  res@mpMaxLatF = 50                    ; specify the plot domain
  res@mpMinLatF = -50
  res@mpCenterLonF         = -200
  res@tmXTOn="False"
  res@gsnZonalMean    = True
  res@gsnZonalMeanXMaxF=6e-6
  res@gsnZonalMeanXMinF=-6e-6
  res@gsnZonalMeanYRefLine=0
  
 
 
   spdiff=(s2v-s1v)
   spdiff!0="lat"
   spdiff&lat=Q1&lat
   spdiff!1="lon"
   spdiff&lon=Q1&lon
  ; qdiffp=100.*((wva9-wva4)/wva4)
  ; qdiffp!0="lat"
  ; qdiffp&lat=Q1&lat
  ; qdiffp!1="lon"
  ; qdiffp&lon=Q1&lon
  ; print(qdiffp)
   ;pcor3=pattern_cor(omegdiff({-10:10},:),spdiff({-10:10},:),1.0,1)
    ;pcor4=pattern_cor(omegdiff({-10:10},:),qdiff({-10:10},:),1.0,0)
	;pcor5=pattern_cor(omegdiffsp({-10:10},:),qstar({-10:10},:),1.0,1)
    ;;;;;;;;;;;;;;;;;;;;;;SET THE CONTOUR PLOT RESOURCES!!!!;;;;;;;;;;;;;;;;;;;
 cres                      = True                 ; plot mods desired
 cres@tmYRMode = "Automatic"
 cres@cnLevelSelectionMode = "ExplicitLevels"
 cres@cnLevels = (/-3.6,-3.4,-3.2,-3.,-2.8,-2.6,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6/)*2e-5
 ;cres@cnLevels=(/-100,-50,-25,-10,0,10,25,50,100/)
 cres@cnLineLabelsOn       = True                 ; turn on line labels
 cres@cnFillOn             = False                 ; turn on color fill
 cres@gsnDraw = False
 cres@gsnFrame = False
 cres@cnInfoLabelOn = False                ; turn off contour info label
 cres@cnLineThicknessF =0.5
 cres@gsnRightString    = f1@case
 cres@gsnRightStringFontHeightF    = 0.05
 cres@cnLevelSpacingF=200.0
 cres@gsnContourZeroLineThicknessF  = 2.0
 cres@cnLineLabelDensityF      = 3.0
 cres@cnLineLabelFontHeightF      = 0.004
 cres@cnLineLabelPlacementMode = "Constant"
 cres@gsnContourNegLineDashPattern =2
 cres@cnLineLabelInterval=1
 ;cres@gsnLeftString  = "10S-10N response pattern correlation: " + pcor + "~C~10S-10N tropical mean Q:" + qstarmean
 cres@gsnLeftString  = "10S-10N response pattern correlation: " + pcor
 cres@gsnLeftStringFontHeightF=0.025
 cres@tmXTOn="False"
 cres@gsnZonalMean    = True
   
 zonmean=dim_avg_n_Wrap(qstar,1)
 plot(n)  = gsn_csm_contour_map(wks, spomeg, res)  
 plot2  = gsn_csm_contour(wks, wva4, cres)
 overlay(plot(n),plot2)
  
   
   
   delete(Q)
   delete(Q1)
   delete(Q2)
   delete(Q3)
   delete(Qtot)
   delete(wva1)
   delete(wva2)
   delete(wva3)
   delete(wva4)
   delete(T)
   delete(omeg)
   delete(s1)
   delete(splon)
   delete(s1v)
   delete(s2v)
   delete(wc)
   delete(wc1)
   delete(qzon)
   n=n+1
end do

 resP                            = True
 resP@gsnPanelLabelBar=True
 resP@cnInfoLabelOn       = False
 resP@lbOrientation =   "horizontal"          ; 
 resP@lbLabelFontHeightF    = 0.004
 resP@txString = "CAM4/CAM-IF ANN Q* (shaded) and CC (contours) terms (925-150 hPa pressure-weighted sums)"
 resP@lbTitleString    = "K/s"
 resP@lbTitleFontHeightF=.007
 resP@lbTitlePosition ="Bottom"
 resP@lbLabelAutoStride= False
 resP@txFontHeightF    = 0.01
 resP@gsnMaximize=True
 resP@gsnPaperOrientation = "portrait"
 resP@lbLabelStride=1
 gsn_panel(wks,plot,(/5,1/),resP)  
 delete(wks) 
;crop=systemfunc("/usr/bin/pdfcrop CAM_ANN_spoqresponse.png")  ;;remove the excessive white space in the PDF! (generates a $filename$-crop.pdf file)
end