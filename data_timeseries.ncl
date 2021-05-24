;*************************************************
;load "statistics2.ncl"
load "select_region.ncl"
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
  dim=dimsizes(tropoz)
  clim1=clmMonTLL(tropoz(:,:,:))
  anom_clim1= calcMonAnomTLL(tropoz(:,:,:),clim1)
  anom_ann=month_to_annual(anom_clim1, 1)
  time=tropoz&time
  date= cd_calendar(time, 0)
  yy=toint(date(:,0)) 
  dim=dimsizes(anom_ann)
  ts1_avg=coswgt_area_average(tropoz)
  ;regions
  ;region_label=(/"WNA","ENA","CAM","TSA","SSA","EUM","NAF","CAF","SAF",\
  ;        "NAS","CAS","EAS","SAS","SEA","AUS"/)
  region_label=(/"ENA","EUM","EAS","SAS","TSA","NAF","SAF"/)
  ts_regions=new((/dimsizes(ts1_avg),dimsizes(region_label)/),typeof(ts1_avg))
  latrange=new((/dimsizes(region_label),2/),"double")
  lonrange=new((/dimsizes(region_label),2/),"double")
  do ir=0,dimsizes(region_label)-1
     region=select_region(region_label(ir))
     latrange(ir,:)=region(0:1)
     lonrange(ir,:)=region(2:3)
     ts_regions(:,ir)=area_operations(tropoz,latrange(ir,0),latrange(ir,1),lonrange(ir,0),lonrange(ir,1),"average",True)
  end do
  ts_regions!0="time"
  ts_regions!1="region"
  ts_regions&region=region_label

  ts_distribution=new((/7,dimsizes(tropoz&time)/),typeof(tropoz))
  do itim=0,dimsizes(tropoz&time)-1
     stat=stat_dispersion(tropoz(itim,:,:), False ) 
     ts_distribution(1:2,itim)=stat(23:24)
     ts_distribution(3,itim)=stat(6)
     ts_distribution(4,itim)=stat(10)
     ts_distribution(5:6,itim)=stat(25:26)
     delete(stat)
  end do

  ts_distribution(0,:)=ts1_avg
  ts_distribution!0="stat"
  ts_distribution!1="time"
  ts_distribution&stat=(/"Mean","Lower 1.0%","Lower 5.0%","LowQuartile","HighQuartile","Upper 5.0%","Upper 1.0%"/)
  wks = gsn_open_wks("pdf","innterannnual_variability_"+"OBS_CDS-TTC-IASI-A"+att)                ; send graphics to PNG file
  gsn_define_colormap(wks,"NEO_div_vegetation_c")
  plot=new(dim(0),"graphic")
  
  do itim=0,dim(0)-1
 
     avg_clim1=coswgt_area_average(anom_ann(itim,:,:))
     avg_std1=coswgt_area_std(anom_ann(itim,:,:))

     res = True                                    ; plot mods desired
     res@gsnDraw                  = False            ; don't draw yet
     res@gsnFrame                 = False            ; don't advance frame yet

     res@cnFillOn         = True               ; color Fill 
     res@cnFillMode       = "RasterFill"       ; Raster Mode
     res@cnLinesOn        =  False             ; Turn off contour lines
     res@cnLineLabelsOn       =False
     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMinLevelValF       = -1.5               ; set min contour level
     res@cnMaxLevelValF       = 1.5               ; set max contour level
     res@cnLevelSpacingF      = 0.3               ; set contour spacing
     ;res@cnLabelBarEndStyle="ExcludeOuterBoxes"
     res@tiMainString         = "Annual map "+yy(itim*12)    ; plot title
     fmt = "%3.2f"
     res@gsnLeftString      = "[DU]"
     res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string

     res@lbLabelBarOn = False
     plot(itim) = gsn_csm_contour_map(wks,anom_ann(itim,:,:), res)  ; create plot
  end do
  resP=True
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  gsn_panel(wks,plot,(/4,4/),resP)
;***********************time serie************************************************


;************************************************
  wks1 = gsn_open_wks("pdf","timeserie_media_"+"OBS_CDS-TTC-IASI-A"+att)                ; send graphics to PNG file

; Set resources necessary to customize Y axis labels
  
  res11=True
  ;res11@gsnDraw                  = False            ; don't draw yet
  ;res11@gsnFrame                 = False            ; don't advance frame yet
  ;res@tiMainString  = "An xy plot Example"   ; title
  res11@tiYAxisString = "Tropospheric Ozone Column [DU]"           ; y axis title
  res11@vpHeightF= 0.4                    ; change aspect ratio of plot
  res11@vpWidthF = 0.6
  ;res@tiXAxisString = "Time"                 ; x axis title
  res11@trXMinF      = min(tropoz&time) 
  ;res11@trXMaxF      = 180.
  res11@trYMinF      = 10.
  res11@trYMaxF      = 16.
  res11@xyLineColor      = "NavyBlue"
  res11@xyLineThicknessF = 3.5
  resTick                = True
  resTick@ttmFormat      = "%N/%Y"
  resTick@ttmAxis        = "XB"
  resTick@ttmMajorStride = 12
  res11@tmXBLabelFontHeightF=0.009
  time_axis_labels( ts1_avg&time,  res11, resTick )
  plot2 = gsn_csm_xy(wks1,time,ts1_avg, res11)  ; create plot
  

;************************************************
  wks2 = gsn_open_wks("pdf","timeseries regionns_"+"OBS_CDS-TTC-IASI-A"+att)                ; send graphics to PNG file

; Set resources necessary to customize Y axis labels

  res11=True
  ;res11@gsnDraw                  = False            ; don't draw yet
  ;res11@gsnFrame                 = False            ; don't advance frame yet
  res11@tiMainString  = "Regional time series"   ; title
  res11@tiYAxisString = "Tropospheric Ozone Column [DU]"           ; y axis title
  res11@vpHeightF= 0.4                    ; change aspect ratio of plot
  res11@vpWidthF = 0.6
  res11@trYMinF      = 0.
  res11@trYMaxF      = 30.
  res11@trXMinF      = min(tropoz&time)
  res11@xyLineColors     = (/"NavyBlue","yellow","green","Purple","LightBlue","orange","red"/)
  res11@xyLineThicknessF = 3
  res11@pmLegendDisplayMode    = "Always"            ; turn on legend
  res11@pmLegendSide           = "Top"               ; Change location of 
  res11@pmLegendParallelPosF   = 0.9                 ; move units right
  res11@pmLegendOrthogonalPosF =- 0.55                ; move units down
  res11@pmLegendWidthF         = 0.1                ; Change width and
  res11@pmLegendHeightF        = 0.15                ; height of legend.
  res11@lgPerimOn              = False               ; turn off box around
  res11@lgLabelFontHeightF     = .01                 ; label font height
  res11@xyExplicitLegendLabels = region_label         ; create explicit labels
   res11@tmXBLabelFontHeightF=0.008
  resTick                = True
  resTick@ttmFormat      = "%N/%Y"
  resTick@ttmAxis        = "XB"
  resTick@ttmMajorStride = 12
  time_axis_labels( ts1_avg&time,  res11, resTick )
  plot2 = gsn_csm_xy(wks2,time,ts_regions(region|:,time|:), res11)  ; create plot
;************************************************
;************************************************
  wks3 = gsn_open_wks("pdf","timeseries_statistic_"+"OBS_CDS-TTC-IASI-A"+att)                ; send graphics to PNG file

; Set resources necessary to customize Y axis labels
  delete(res11)
  res11=True
  ;res11@gsnDraw                  = False            ; don't draw yet
  ;res11@gsnFrame                 = False            ; don't advance frame yet
  ;res11@tiMainString  = "Statistical distribution"   ; title
  res11@tiYAxisString = "Tropospheric Ozone Column [DU]"           ; y axis title
  res11@vpHeightF= 0.4                    ; change aspect ratio of plot
  res11@vpWidthF = 0.6
   res11@trXMinF      = min(tropoz&time)
  res11@trYMinF      = 0.
  res11@trYMaxF      = 30.
  res11@xyLineColors     = (/"Black","LightBlue","NavyBlue","Purple","Purple","NavyBlue","LightBlue"/)
  res11@xyMonoLineThickness =False
  res11@xyLineThicknesses = (/3,2,2,2,2,2,2./)
  res11@pmLegendDisplayMode    = "Always"            ; turn on legend
  res11@pmLegendSide           = "Top"               ; Change location of
  res11@pmLegendParallelPosF   = 1.15                  ; move units right
  res11@pmLegendOrthogonalPosF = -0.6                ; move units down
  res11@pmLegendWidthF         = 0.09                ; Change width and
  res11@pmLegendHeightF        = 0.18                ; height of legend.
  res11@lgPerimOn              = False               ; turn off box around
  res11@lgLabelFontHeightF     = .01                 ; label font height
  res11@xyExplicitLegendLabels = ts_distribution&stat         ; create explicit labels
  res11@tmXBLabelFontHeightF=0.008
  resTick                = True
  resTick@ttmFormat      = "%N/%Y"
  resTick@ttmAxis        = "XB"
  resTick@ttmMajorStride = 12
  time_axis_labels( ts1_avg&time,  res11, resTick )
  plot3 = gsn_csm_xy(wks3,time,ts_distribution(stat|:,time|:), res11)  ; create plot
end
