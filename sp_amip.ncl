begin
T6=new((/10,17,9/),float)
u=new((/10,25,9/),float)
f7=addfile("zonmean_ta_AMIP_allmodels_DJF.nc","r")
n=1
do while(n.le.9)
 T6(n,:,:)=f7->ta(n,:,{-15:5},0)
 n=n+1
end do
  u@units="K"
  u!0 = "model"                  ; variable/dim name 
  u!1= "plev"                  ; variable/dim name 
  u!2= "lat"
  s1   = static_stability(T6&plev, T6, 1, 0)
  s1!0="model"
  s2=dim_avg_n_Wrap(s1,2)
  dp   = dpres_plevel(T6&plev, 101300., 3000., 0)
  print(s1) 
n=1
data=new(10,double)
do while(n.le.9)
  fint = dim_sum_n(s2(n,1:10)*dp(1:10), 0)/dim_sum_n(dp(1:10),0)
  data(n)=fint
  n=n+1
end do
print(data)
  
end

