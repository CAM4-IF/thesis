load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

  case     = "defaultcloud:"  
  diro     = "TRMM_3hourly/CAM4-IF_defaultcloud/"        ; output dir: location of plots 

  latN     = 15
  latS     = -latN       ; symmetric about the equator

  nDayWin  = 96     ; Wheeler-Kiladis [WK] temporal window length (days)
  nDaySkip = 10      ; time (days) between temporal windows [segments]
                     ; negative means there will be overlapping 
                     ; temporal segments

  opt      = True
  opt@debug= True
  opt@pltType     = "eps"
  opt@cnLinesOn   = False   ; turn off contour lines
                     ; specify the plot intervals
  opt@Fig_3a=fspan(0.5, 2.2, 171) 
  opt@Fig_3b=fspan(0.5, 2.2, 171) 
  opt@pltColorMap="nice_gfdl"
  opt@lbBoxLinesOn=False
  opt@lbBoxSeparatorLinesOn=False
  
  spd      = 8       ; samples per day (6 hour sampling rate)
  var      = "PRECC"  ; user specified 

  diri     = "TRMM_3hourly/CAM4-IF_defaultcloud/"
  fili     = systemfunc("cd "+diri+" ; ls *h1*.nc")
  print(fili)  

  f        = addfiles (diri+fili, "r") 
  VAR   = f[:]->PRECC(:,{latS:latN},:)     ; (time,lat,lon)
  printVarSummary(VAR)

  nMsg     = num(ismissing(VAR))  ; check for _FillValue
  if (nMsg.gt.0) then
      print("missing values present: not allowed: nMsg="+nMsg) 
      exit
  end if
                     ; calculation and create plots
  wkSpaceTime (VAR, diro, case, var             \
              ,latN, spd, nDayWin, nDaySkip, opt  )