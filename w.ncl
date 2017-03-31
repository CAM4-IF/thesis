begin
 f    = addfile("tropmean_w_CAM-IFv1_4K_DJF.nc", "r")
 f1=addfile("CAM-IFv1_4K_DJF_climo.nc","r")
 ;f1=addfile("CAM4_4K_DJF_climo.nc", "r")
 ;f2=addfile("tropmean_Q_CAM4_DJF_climo.nc","r")
   hyam = f1->hyam 
   hybm = f1->hybm
   P0   = 100000.                        ; reference pressure [Pa]
   PS   = f->PS(0,:,:)                         ; Pa  [punits=1]
   T    = f->OMEGA(0,:,:,:)
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa] 
   opt=(/0,15000,92500/)
   print(dimsizes(pm(:,:,:)))
   print(dimsizes(T))
   wva = wgt_vert_avg_beta(pm(:,:,:), T(:,:,:), PS, 1, opt)
   print(wva)
end
