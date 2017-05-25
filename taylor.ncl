;**********************************
; taylor_4.ncl
;**********************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
;
; These files still have to be loaded manually
load "taylor_diagram.ncl"
load "taylor_metrics_table.ncl"
;**********************************
begin
;**********************************
; Assume the following have already been computed:
;   _ratio are the ratio:  Case_Variance/Reference_Variance
;   _cc    are the cross correlation coef of Case to Reference
; In this example, these are derived for annual mean climatologies.
;**********************************

; Cases [Model]
  case      = (/ "CAM4-ZM", "CAM4-IF-t" /) 
  nCase     = dimsizes(case )                 ; # of Cases [Cases]

; variables compared
  var       = (/ "Prc 15S-15N","Prc 30S-30N","Omega 500 15S-15N","Omega 500 30S-30N" /) 
  nVar      = dimsizes(var)                   ; # of Variables

; more info to be added [all are bogus]
  source    = (/ "Blend" , "Blend", "MERRA/ERAI","MERRA/ERAI"/)

; "CAM4-ZM"                        
  CA_ratio   = (/1.290, 1.092,1.18,1.03/)
  CA_cc      = (/0.875, 0.886,0.755,0.736/)

; "CAM4-IF-t" 
  CB_ratio   = (/1.620, 1.590,1.29,1.15/)
  CB_cc      = (/0.888, 0.883,0.748,0.722/)

; arrays to be passed to taylor_diagram. It will calculate the x xnd y coordinates.
  ratio      = new ((/nCase, nVar/),typeof(CA_cc) )  
  cc         = new ((/nCase, nVar/),typeof(CA_cc) ) 

  ratio(0,:) = CA_ratio 
  ratio(1,:) = CB_ratio

  cc(0,:)    = CA_cc 
  cc(1,:)    = CB_cc

;**********************************
; create plot
;**********************************
  varSource = var +"_"+ source              ; add extra info [*not* required]

  ty_opt   = True                           ; taylor diagram with options
        
  ty_opt@Markers       = (/16, 16/)         ; make all solid fill
  ty_opt@Colors        = (/"red", "blue" /)          
  ty_opt@varLabels     = varSource
  ty_opt@caseLabels    = case

  ty_opt@varLabelsYloc = 1.5                ; Move location of variable labels [default 0.45]
  ty_opt@caseLabelsFontHeightF = 0.14       ; make slight larger   [default=0.12 ]
  ty_opt@varLabelsFontHeightF  = 0.011      ; make slight smaller  [default=0.013]
 
  ty_opt@tiMainString  = "ANN - Spatial"           ; title
  ty_opt@stnRad        = (/ 0.5, 1.5 /)     ; additional standard radii
  ty_opt@ccRays        = (/ 0.6, 0.9 /)     ; correllation rays
  ty_opt@centerDiffRMS = True               ; RMS 'circles'

  wks   = gsn_open_wks("pdf","taylor.000001")      ; send graphics to PNG file
  plot  = taylor_diagram(wks,ratio,cc,ty_opt)

;**************************************************
; fill an array for a "taylor metrics table"
;**************************************************

  season    = (/ "ANN" /)
  nSeason   = dimsizes(season)

  table     = new ( (/nCase,nSeason,nVar/), typeof(ratio) )
  table(0,0,:) = CA_ratio
  table(1,0,:) = CB_ratio

  tt_opt        = True
  tt_opt@pltType= "pdf"                  ; "eps" [default], "pdf", "ps"
                                         ; "png", "gif" [if you have ImageMajik 'convert']
  tt_opt@color0 = "palegreen2"
  tt_opt@color1 = "tomato2"


  taylor_metrics_table("taylor.000002", varSource, case ,season, table, tt_opt)
end 