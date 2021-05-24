;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
;************************************************
begin
;************************************************
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
; #############################################################################

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
  opt = False
    avail1=new((/dim(1),dim(2)/),typeof(tropoz))
     do ilat=0,dim(1)-1
      do ilon=0,dim(2)-1
          avail1(ilat,ilon)=num(.not.ismissing(tropoz(:,ilat,ilon)))/tofloat(dim(0))
      end do
    end do

  tropoz1=rmMonAnnCycTLL(tropoz(:,:,:))   ;remove annual cycle
  trend1  = trend_manken(tropoz1, opt, 0)  ; ===> px(2,nlat,mlon)
  trend1(1,:,:)=(/trend1(1,:,:)*12/) 
  trend1(0,:,:)=where(avail1.gt.0.7,trend1(0,:,:),trend1@_FillValue)
  trend1(1,:,:)=where(avail1.gt.0.7,trend1(1,:,:),trend1@_FillValue)                                   ; if meta data is desired
  copy_VarCoords(tropoz1(0,:,:),trend1(0,:,:))
  trend1!0= "prob_trend"               ; ===> size 2

  wks = gsn_open_wks("pdf","trend_1_"+"OBS_CDS-TTC-IASI-A"+att)         ; send graphics to PNG file
  gsn_define_colormap(wks,"NEO_div_vegetation_c")
  plot=new(2,"graphic")
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet

  res@cnFillOn         = True               ; color Fill 
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLineLabelsOn       =False
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = -0.5               ; set min contour level
  res@cnMaxLevelValF       = 0.5              ; set max contour level
  res@cnLevelSpacingF      = 0.05               ; set contour spacing
  ;res@cnLabelBarEndStyle="ExcludeOuterBoxes"
  res@tiMainString         = "Trend "    ; plot title
  fmt = "%3.2f"
  res@gsnLeftString      = "[DU/year]"
  ;res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string
  res@lbOrientation = "vertical"
  ;res@lbLabelBarOn = False
   plot(0) = gsn_csm_contour_map(wks,trend1(1,:,:), res)  ; create plot
   res@cnMinLevelValF       = 0               ; set min contour level
   res@cnMaxLevelValF       = 1               ; set max contour level
   res@cnLevelSpacingF      = 0.05               ; set contour spacing
   res@tiMainString         = "Trend Significance"  
   res@gsnLeftString      = "[1]"
   plot(1) = gsn_csm_contour_map(wks,trend1(0,:,:), res)
   resP=True
  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
   gsn_panel(wks,plot,(/2,1/),resP)


end
