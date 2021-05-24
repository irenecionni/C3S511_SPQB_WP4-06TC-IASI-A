;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
;************************************************
;functions
; #############################################################################
undef("coswgt_area_average")
function coswgt_area_average(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
local funcname, scriptname, lat, wgt_lat, lon, lon_size, wgt_lon, ave
begin

  funcname = "coswgt_area_average"
  scriptname = "diag_scripts/shared/statistics.ncl"

  lat = field&lat
  wgt_lat = tofloat(NormCosWgtGlobe(lat))

  lon = field&lon
  lon_size = dimsizes(lon)
  wgt_lon = new((/lon_size(0)/), float)
  wgt_lon = 1.0

  ave = wgt_areaave_Wrap(field, wgt_lat, wgt_lon, 0)
  return(ave)

end
; #############################################################################
undef("coswgt_area_std")
function coswgt_area_std(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
begin

  funcname = "coswgt_area_std"

   lat = field&lat
   lon = field&lon
   re   = 6.37122e06
   rad  = 4.0 * atan(1.0) / 180.0
   con  = re * rad
   clat = cos(lat * rad)           ; cosine of latitude
   dlon = (lon(2) - lon(1))        ; assume dlon is constant
   dlat = (lat(2) - lat(1))        ; assume dlat is constant

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dydx = dy * dx                  ; dydx(nlat)
   nlat=dimsizes(lat)
   mlon=dimsizes(lon)
   wgt  = new((/nlat, mlon/), typeof(dydx))
   wgt=conform (wgt, dydx, 0)



   ave = wgt_areaave2(field, wgt, 0)
   d2 = (field - ave) ^ 2 ;
     ; Calculate weighted standard deviation
   ss=sum(d2)
   arg = wgt_areasum2(d2, wgt, 0)

   v1 = sum(where(ismissing(field), wgt@_FillValue, wgt))
   v2 = sum(where(ismissing(field), wgt@_FillValue, wgt) ^ 2)
   dd=(v1 / (v1 ^ 2 - v2))

   out = sqrt( arg*dd)

  return(out)


end

begin

;************************************************
; read in netCDF file
;************************************************
 att="DAY" 
 dir="/Users/irenecionni/Copernicus/"  
 a = addfile(dir+"obs_c3s/Tier3/OBS_CDS-IASI-A_sat_L3_Amon_tropoz_200701-202012"+att+".nc" ,"r")
;************************************************
; read in zonal winds
;************************************************
 
  tropoz=a->tropoz                                
  tropoz=where(tropoz.lt.0,tropoz@_FillValue,tropoz)
  printVarSummary(tropoz)
  dim=dimsizes(tropoz)
  media_clim=dim_avg_n_Wrap(tropoz(12:dimsizes(tropoz&time)-13,:,:),0)
  std_clim=dim_stddev_n_Wrap(tropoz(12:dimsizes(tropoz&time)-13,:,:),0)
  ;clim1=clmMonTLL(tropoz(3:dim(0)-2,:,:))
  clim1=clmMonTLL(tropoz(:,:,:))
  mer1=dim_avg_n_Wrap(tropoz,1)
  mer1_avg=dim_avg_n_Wrap(mer1,0)
  zon1=dim_avg_Wrap(tropoz)
  zon1_avg=dim_avg_n_Wrap(zon1,0)
  wks = gsn_open_wks("pdf","Figure3_"+"OBS_CDS-TTC-IASI-A"+att)                ; send graphics to PNG file
  gsn_define_colormap(wks,"MPL_RdYlGn")
  gsn_reverse_colormap(wks)
  plot=new(12,"graphic")
  month_nam=(/"Jannuary","February","March","April","May","June","July","August","September","October","November","December"/)
  do itim=0,11
 
    avg_clim1=coswgt_area_average(clim1(itim,:,:))
    avg_std1=coswgt_area_std(clim1(itim,:,:))

     res = True                                    ; plot mods desired
     res@gsnDraw                  = False            ; don't draw yet
     res@gsnFrame                 = False            ; don't advance frame yet

     res@cnFillOn         = True               ; color Fill 
     res@cnFillMode       = "RasterFill"       ; Raster Mode
     res@cnLinesOn        =  False             ; Turn off contour lines
     res@cnLineLabelsOn       =False
     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMinLevelValF       = 7.               ; set min contour level
     res@cnMaxLevelValF       = 23.               ; set max contour level
     res@cnLevelSpacingF      = 1               ; set contour spacing
     ;res@cnLabelBarEndStyle="ExcludeOuterBoxes"
     res@tiMainString         = "Climatology "+month_nam(itim)    ; plot title
     fmt = "%3.2f"
     res@gsnLeftString      = "[DU]"
     res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string

     res@lbLabelBarOn = False
     plot(itim) = gsn_csm_contour_map(wks,clim1(itim,:,:), res)  ; create plot
  end do
  resP=True
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  gsn_panel(wks,plot,(/4,3/),resP)

;************************************************
  wks1 = gsn_open_wks("pdf","FigureS2_"+"OBS_CDS-TTC-IASI-A"+att)                ; send graphics to PNG file
  plot2=new(2,"graphic")
  gsn_define_colormap(wks1,"MPL_RdYlGn")
  gsn_reverse_colormap(wks1)
  res1 = True                                    ; plot mods desired
  res1@gsnDraw                  = False            ; don't draw yet
  res1@gsnFrame                 = False            ; don't advance frame yet

  res1@cnFillOn         = True               ; color Fill
  res1@cnFillMode       = "RasterFill"       ; Raster Mode
  res1@cnLinesOn        =  False             ; Turn off contour lines
  res1@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res1@cnMinLevelValF       = 7.               ; set min contour level
  res1@cnMaxLevelValF       = 20.               ; set max contour level
  res1@cnLevelSpacingF      = 1               ; set contour spacing
  res1@cnLabelBarEndStyle="ExcludeOuterBoxes"
  res1@lbOrientation = "vertical"
  ;res1@tiMainString         = "Data Availability as a fraction [1]"    ; plot title
  fmt = "%3.2f"
;  res@gsnLeftString      = "mean: "+sprintf(fmt,avg_avail)
;  res@gsnRightString     = " std:"+sprintf(fmt,avg_std)           ; plot center string
  res1@vpHeightF= 0.4                    ; change aspect ratio of plot
  res1@vpWidthF = 0.7                  
  resTick                = True
  resTick@ttmFormat      = "%N/%Y" 
  resTick@ttmAxis        = "YL"
  resTick@ttmMajorStride = 12

; Set resources necessary to customize Y axis labels
  time_axis_labels( mer1&time,  res1, resTick )
  flip_mer1=lonFlip(mer1)
  plot2(0) = gsn_csm_contour(wks1,flip_mer1, res1)  ; create plot
  flip_mer1_avg=lonFlip(mer1_avg)
  res11=True
  res11@gsnDraw                  = False            ; don't draw yet
  res11@gsnFrame                 = False            ; don't advance frame yet
  ;res@tiMainString  = "An xy plot Example"   ; title
  res11@tiYAxisString = "Tropospheric Ozone Column [DU]"           ; y axis title
  res11@vpHeightF= 0.4                    ; change aspect ratio of plot
  res11@vpWidthF = 0.7
  ;res@tiXAxisString = "Time"                 ; x axis title
  res11@trXMinF      = -180. 
  res11@trXMaxF      = 180.
  res11@trYMinF      = 10.
  res11@trYMaxF      = 16.
  res11@xyLineColor      = "NavyBlue"
  res11@xyLineThicknessF = 3.5
  plot2(1) = gsn_csm_xy(wks1,flip_mer1_avg&lon,flip_mer1_avg, res11)  ; create plot

   gsn_panel(wks1,plot2,(/2,1/),False)
 frame(wks1)

;************************************************
  wks3 = gsn_open_wks("pdf","Figure4_"+"OBS_CDS-TTC-IASI-A"+att)                ; send graphics to PNG file
  plot3=new(2,"graphic")
  gsn_define_colormap(wks3,"MPL_RdYlGn")
  gsn_reverse_colormap(wks3)
  res2 = True                                    ; plot mods desired
  res2@gsnDraw                  = False            ; don't draw yet
  res2@gsnFrame                 = False            ; don't advance frame yet

  res2@cnFillOn         = True               ; color Fill
  res2@cnFillMode       = "RasterFill"       ; Raster Mode
  res2@cnLinesOn        =  False             ; Turn off contour lines
  res2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res2@cnMinLevelValF       = 7.               ; set min contour level
  res2@cnMaxLevelValF       = 23.               ; set max contour level
  res2@cnLevelSpacingF      = 1               ; set contour spacing
  ;res2@lbOrientation = "vertical"
  ;res2@tiMainString         = "Data Availability as a fraction [1]"    ; plot title
  ;res2@lbLabelBarOn = False
  ;res2@cnLabelBarEndStyle="ExcludeOuterBoxes"
  fmt = "%3.2f"
;  res@gsnLeftString      = "mean: "+sprintf(fmt,avg_avail)
;  res@gsnRightString     = " std:"+sprintf(fmt,avg_std)           ; plot center string
  res2@vpHeightF= 0.4                    ; change aspect ratio of plot
  res2@vpWidthF = 0.7
  res2@tmYLMode        = "Explicit"	
  res2@tmYLValues      = (/-60,-30,0,30,60/)
  res2@tmYLLabels      = (/"60S","30S","Eq","30N","60N"/)
  res2@tmYLLabelFontHeightF=0.015
  res2@tmXBLabelFontHeightF=0.008
  res2Tick                = True
  res2Tick@ttmFormat      = "%N/%Y"
  res2Tick@ttmAxis        = "XB"
  res2Tick@ttmMajorStride = 12
  
; Set resources necessary to customize Y axis labels
  time_axis_labels( zon1&time,  res2, res2Tick )
  plot3(0) = gsn_csm_lat_time(wks3,zon1(lat|:,time|:), res2)  ; create plot
  res21=True
  res21@gsnDraw                  = False            ; don't draw yet
  res21@gsnFrame                 = False            ; don't advance frame yet
  ;res@tiMainString  = "An xy plot Example"   ; title
  res21@tiXAxisString = "Tropospheric Ozone Column [DU]"           ; y axis title
  res21@vpHeightF= 0.4                    ; change aspect ratio of plot
  res21@vpWidthF = 0.3
  ;res@tiXAxisString = "Time"                 ; x axis title
  res21@trYMinF      = -90.
  res21@trYMaxF      = 90.
  res21@trXMinF      = 0.
  res21@trXMaxF      = 25
  res21@tmYLMode        = "Explicit"
  res21@tmYLValues      = (/-60,-30,0,30,60/)
  res21@tmYLLabels      = (/"60S","30S","Eq","30N","60N"/)
  res21@tmYLLabelFontHeightF=0.012
  res21@tmXBMode        = "Explicit"
  res21@tmXBValues      = (/0.,5.,10,15,20.,25./)
  res21@tmXBLabels      = (/"0","5","10","15","20","25"/)
  res21@tmXBLabelFontHeightF=0.01

  res21@xyLineColor      = "NavyBlue"
  res21@xyLineThicknessF  =3.5
  plot3(1) = gsn_csm_xy(wks3,zon1_avg,zon1_avg&lat, res21)  ; create plot
  resP1=True
  gsn_panel(wks3,plot3,(/1,2/),resP1)
  frame(wks3)
;---------mean clim---------------------------------------------
  delete(res)
  wks4=gsn_open_wks("pdf","Figure2_OBS_CDS-TTC-IASI-A"+att)   
  gsn_define_colormap(wks4,"MPL_RdYlGn")
  gsn_reverse_colormap(wks4)

   res = True                                    ; plot mods desired
     
   res@cnFillOn         = True               ; color Fill 
   res@cnFillMode       = "RasterFill"       ; Raster Mode
   res@cnLinesOn        =  False             ; Turn off contour lines 
   res@cnLineLabelsOn       =False             
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = 7.               ; set min contour level
   res@cnMaxLevelValF       = 23.               ; set max contour level
   res@cnLevelSpacingF      = 1               ; set contour spacing
   ;resn@cnLabelBarEndStyle="ExcludeOuterBoxes"
   res@tiMainString         = "Mean "    ; plot title
   fmt = "%3.2f"
   res@gsnLeftString      = "[DU]"
   res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string

   plot4 = gsn_csm_contour_map(wks4,media_clim(:,:), res)  ; create plot
   delete(res)
   frame (wks4)
   wks5=gsn_open_wks("pdf","std_OBS_CDS-TTC-IASI-A"+att)
   gsn_define_colormap(wks5,"NEO_div_vegetation_c")

   res = True                                    ; plot mods desired
   res@cnFillOn         = True               ; color Fill
   res@cnFillMode       = "RasterFill"       ; Raster Mode
   res@cnLinesOn        =  False             ; Turn off contour lines
   res@cnLineLabelsOn       =False
   ;resn@cnLabelBarEndStyle="ExcludeOuterBoxes"
   res@tiMainString         = "Standard Deviation "    ; plot title
   fmt = "%3.2f"
   res@gsnLeftString      = "[DU]"
;     res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string

;     res@lbLabelBarOn = False
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels 
   res@cnMinLevelValF       = -1               ; set min contour level
   res@cnMaxLevelValF       = 3.               ; set max contour level
   res@cnLevelSpacingF      = 0.3               ; set contour spacing
   plot5 = gsn_csm_contour_map(wks5,std_clim(:,:), res)  ; create plot
   frame (wks5)

end
