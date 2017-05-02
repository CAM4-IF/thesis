load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

fils=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/results/*CAM*ANN*.nc")
fils1=systemfunc("ls /run/media/shawn/4a83cd00-0edb-412a-b523-7ef85625597a/research_data/CAM-IF/new_data/results/*CAM*4K*ANN*.nc1")
datalog=new(dimsizes(fils),double) ; create list to store output
nummods=dimsizes(fils1) ;;number of models! 
datalog=new((/nummods/),"double")
plot=new(nummods+1, graphic)    ;;;Create a panel plot with panels = number of models!
;;;create a giant loop over all models;;;
lets=(/"a)","b)","c)","d)"/)
n=0
do while(n.le.dimsizes(fils)-1)
  f=addfile(fils(n),"r") ; add file
  f1=addfile(fils1(n),"r") ; add file
  U=f->U
  V=f->V
  U1=f1->U
  V1=f1->V
  Z=f->Z3
  lons=f1->lon
  
  
  
  hyam = f1->hyam ; read from a file the mid-layer coef
  hybm = f1->hybm ; read from a file
  psfc = f->PS(:,:,:)   ; surface pressure [Pa]
  psfc1=f->PS(:,:,:)
  P0   = 100000. ; since ps is in Pa or [ f->P0]
  pm  = pres_hybrid_ccm(psfc,P0,hyam,hybm)/100.
  pm!1="lev"
  pm!2="lat"
  pm&lat=U&lat
  pm!3="lon"
  pm&lon=U&lon
  
  
  
  hyam1 = f1->hyam ; read from a file the mid-layer coef
  hybm1 = f1->hybm ; read from a file
 ; psfc1 = f1->PS(:,{-30:30},:)   ; surface pressure [Pa]
  P01   = 100000. ; since ps is in Pa or [ f->P0]
  pm1  = pres_hybrid_ccm(psfc1,P01,hyam1,hybm1)/1000.
  
  
  lev_p           = (/  10, 20, 30, 50, 70,100,150,200,250,300,350,400,500,600,700,800,875,925,950,975,1000 /)                
  lev_p!0         = "lev_p"                  ; variable/dim name 
  lev_p&lev_p     =  lev_p                   ; create coordinate variable
  lev_p@long_name = "pressure"               ; attach some attributes
  lev_p@units     = "hPa"
  lev_p@positive  = "down"

; convert to pressure levels
 vp   = vinth2p(V(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
 UU=vinth2p(U(0,:,:,:),hyam,hybm,lev_p,psfc(0,:,:),2,1000.,1,False)
 vp1   = vinth2p(V1(0,:,:,:),hyam,hybm,lev_p,psfc1(0,:,:),2,1000.,1,False)
 UU1=vinth2p(U1(0,:,:,:),hyam,hybm,lev_p,psfc1(0,:,:),2,1000.,1,False)

 
;;;calculate 200 hPa Chi;;;;
nlat=dimsizes(vp&lat)
mlon=dimsizes(UU&lat)
uvmsg=1e+36
vpot=new((/96,144/),float,uvmsg)
sf=new((/96,144/),float,uvmsg)
vpot1=new((/96,144/),float,uvmsg)
sf1=new((/96,144/),float,uvmsg)

uv2sfvpf(UU({200},:,:),vp({200},:,:),sf,vpot)
uv2sfvpf(UU1({200},:,:),vp1({200},:,:),sf,vpot1)

chistar=dim_rmvmean_n_Wrap(vpot,0)
chistar1=dim_rmvmean_n_Wrap(vpot1,0)



chistar!0="lat"
chistar&lat=U&lat
chistar!1="lon"
chistar&lon=U&lon

chistar1!0="lat"
chistar1&lat=U&lat
chistar1!1="lon"
chistar1&lon=U&lon
chidiff=chistar1-chistar

chidiff!0="lat"
chidiff&lat=U&lat
chidiff!1="lon"
chidiff&lon=U&lon



chistarmin=avg(chistar({-10:10},{90:180}))
chistarmin1=avg(chistar1({-10:10},{90:180}))
response=((chistarmin1-chistarmin)/(chistarmin))
datalog(n)=response

 n=n+1
  delete(chistar)
  delete(chistar1)
  delete(chidiff)
  
end do
print(datalog)
print(fils)
print(fils1)
end