begin
;************************************************
; read in data
;************************************************FM
fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/results/*CAM*ANN*.nc")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/results/*CAM*4K*ANN*.nc1")
print(fils)
print(fils1)
numfils=dimsizes(fils)
nummods=dimsizes(fils1) ;;number of models! 
plot=new(5, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
datalog=new((/numfils/),"double")
datalog1=new((/numfils/),"double")
n=0
lets=(/"a)","b)","c)","d)"/)
do while(n.le.dimsizes(fils)-1)
   filen=tostring(fils(n))
   filen1=tostring(fils1(n))
   h=addfile("/run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/CAM/files/CAM4_4K_ANN_climo.nc2", "r")
   hyam = h->hyam 
   hybm = h->hybm
   zon="cdo -s -zonmean -selname,PS,DTCOND,QRL,QRS,DTV " + filen + " zon.nc"
   zon1="cdo -s -zonmean -selname,PS,DTCOND,QRL,QRS,DTV " + filen1 + " zon1.nc"
   cmd="cdo -s -fldmean -sellonlatbox,90,180,-10,10 -sub -selname,PS,DTCOND,QRL,QRS,DTV " + filen + " -selname,PS,DTCOND,QRL,QRS,DTV -enlarge," + filen + " zon.nc"+ " tmp.nc"
   cmd1="cdo -s -fldmean -sellonlatbox,90,180,-10,10 -sub -selname,PS,DTCOND,QRL,QRS,DTV " + filen1 + " -selname,PS,DTCOND,QRL,QRS,DTV -enlarge," + filen + " zon.nc"+ " tmp1.nc"
   zonmean=systemfunc(zon)
   zonmean1=systemfunc(zon1)
   areamean=systemfunc(cmd)
   areamean1=systemfunc(cmd1)
   f    = addfile("tmp.nc", "r")
   f1=addfile("tmp1.nc","r")
   cond=f->DTCOND(0,:,:,:)
   lw=f->QRL(0,:,:,:)
   sw=f->QRS(0,:,:,:)
   diff=f->DTV(0,:,:,:)
   cond1=f1->DTCOND(0,:,:,:)
   lw1=f1->QRL(0,:,:,:)
   sw1=f1->QRS(0,:,:,:)
   diff1=f1->DTV(0,:,:,:)
   Qtot=sw+lw+cond+diff
   Qtot1=sw1+lw1+cond1+diff1
   P0=100000.
   PS=new((/1,1/),float)
   PS(0,:)=1000.
   PS(:,0)=1000.
   PS1=new((/1,1/),float)
   PS1(0,:)=1000.
   PS1(:,0)=1000.
   pm   = pres_hybrid_ccm(PS,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa]
   pm1   = pres_hybrid_ccm(PS1,P0,hyam,hybm)    ; ps(time,lev,lat,lon) [Pa]    
   opt=(/0,100,1000/)
   wva = wgt_vert_avg_beta(pm/100, Qtot, PS, 0, opt)/900
   wva1 = wgt_vert_avg_beta(pm1/100, Qtot1, PS1, 0, opt)/900
   datalog(n)=wva(0,0)
   response=(wva1(0,0)-wva(0,0))/wva(0,0)
   datalog1(n)=response
   print(wva)
   print(wva1)
   
   
   n=n+1
end do
print(datalog)
print(sprintf("%5.3f",(datalog1)))
end