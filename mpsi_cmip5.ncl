; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; This file still has to be loaded manually
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*********************************************
begin
;**********************************************
;;Loop over all models (each group of 12 timesteps/months)
i=0
maxmean=new(120,float)
minmean=new(120,float)
do while(i.le.119)
  f     = addfile ("va_Amon_AllModels_amip4K_ymonmean.nc","r") 
  f1 =    addfile ("ps_Amon_AllModels_amip4K_ymonmean.nc", "r")
  ps              = f1->ps(:,{-45:45},:)          ; grab surface pressure in Pa
  v               = f->va(:,::-1,{-45:45},:)
  lat             = f->lat({-45:45})
  lev= f->plev(::-1)
; interpolate to new levels to be consistent with CAM output 
  lev1=(/50.,70.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,875.,900.,925.,950.,975.,1000./)*100
  v1=int2p_n(lev,v,lev1,1,1)
; calculate zonal meridional streamfunction
  zmpsi = zonal_mpsi(v1(i,:,:,:),lat,lev1,ps(i,:,:))
  maxmpsi=max(zmpsi(4:21,:))
 ; fiven=round(0.5555*dimsizes(zmpsi(1,:)),3)
 ; fifteens=round(0.3333*dimsizes(zmpsi(1,:)),3)
 ; slice1=zmpsi(4:17,fiven)
 ; slice2=zmpsi(4:17,fifteens)
 ; mean1=dim_avg_n(slice1,0)
 ; mean2=dim_avg_n(slice2,0)
 ; diff=mean2-mean1
 ; mean=dim_avg_n(diff,0)
  ;slice=zmpsi(4:17,:)
  ;mean=dim_avg_n(slice,0)
  maxmean(i)=maxmpsi
  ;minmean(i)=min(mean)
  i=i+1
end do
asciiwrite("maxes_amip4K.txt",maxmean)
;asciiwrite("mins_amip4K.txt",minmean)
end
