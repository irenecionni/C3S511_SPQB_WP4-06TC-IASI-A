;;#############################################################################
;; GENERAL COLLECTION OF PHYSICAL CONSTANTS
;;#############################################################################

begin

    source = "http://physics.nist.gov/cuu/Constants/Table/allascii.txt"

    ; Avogadro number [mol^-1]
    Navo = 6.0221418e23
    Navo@source = source
    Navo@units = "mol^-1"
    Navo@name = "Avogadro number"

    ; Boltzmann constant [J K^-1]
    Kbol = 1.3806504e-23
    Kbol@source = source
    Kbol@units = "J K^-1"
    Kbol@name = "Boltzmann constant"

    ; Gas constant [J mol^-1 K^-1]
    Rgas = 8.314472
    Rgas@source = source
    Rgas@units = "J mol^-1 K^-1"
    Rgas@name = "Gas constant"

    ; Newtonian constant of gravitation [m^3 kg^-1 s^-2]
    Gnewt = 6.67428e-11
    Gnewt@source = source
    Gnewt@units = "m^3 kg^-1 s^-2"
    Gnewt@name = "Newtonian constant of gravitation"

    ; Speed of light in vacuum [m s^-1]
    Clig = 2.9979246e08
    Clig@source = source
    Clig@units = "m s^-1"
    Clig@name = "Speed of light in vacuum"

    ; Standard acceleration of gravity [m s^-2]
    Grav = 9.80665
    Grav@source = source
    Grav@units = "m s^-2"
    Grav@name = "Standard acceleration of gravity"

    ; Standard atmosphere [Pa]
    p0stp = 101325.
    p0stp@units = "Pa"
    p0stp@name = "Standard atmosphere"

    ; Standard temperature [K]
    T0stp = 273.15
    T0stp@units = "K"
    T0stp@name = "Standard temperature"

    ; Stefan-Boltzmann constant [W m^-2 K^-4]
    Sstebol = 5.670400e-08
    Sstebol@source = source
    Sstebol@units = "W m^-2 K^-4"
    Sstebol@name = "Stefan-Boltzmann constant"

    ; Pi [1]
    PI = 3.1415927
    PI@units = "1"
    PI@name = "Greek Pi"

    ; Earth radius [m]
    EARTHradius = 6371000.
    EARTHradius@units = "m"
    EARTHradius@name = "Earth radius"

    ; Earth surface area [m^2]
    EARTHarea = 4. * PI * EARTHradius ^ 2
    EARTHarea@units = "m2"
    EARTHarea@name = "Earth surface area"

    ; Molecular weight of the atmosphere [g mol^-1]
    MWair = 28.97
    MWair@source = "Seinfeld & Pandis - Atmospheric Chemistry " + \
                   "and Physics - Wiley 2006"
    MWair@units = "g mol^-1"
    MWair@name = "Molecular weight of the atmosphere"

    ; Molecular weight of ozone [g mol^-1]
    MWh2o = 18.01528
    MWh2o@units = "g mol^-1"
    MWh2o@name = "Molecular mass of h2o"

    ; Molecular weight of ozone [g mol^-1]
    MWo3 = 47.9982
    MWo3@source = "http://pubchem.ncbi.nlm.nih.gov"
    MWo3@units = "g mol^-1"
    MWo3@name = "Molecular weight of ozone"

    ; Dobson unit
    DU = 2.1415e-05
    DU@units = "kg m^-2"
    DU@name = "conversion factor 1 Dobson unit = 2.1415E-5 kg m^-2"

end

; #############################################################################
undef("gridcell_area")
function gridcell_area(deltax[1]: numeric,
                       lat_lo[1]: numeric,
                       lat_hi[1]: numeric)
local funcname, scriptname, deg2rad
begin

  funcname = "gridcell_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  ;enter_msg(scriptname, funcname)

  ; Convert to radiants
  deg2rad = acos(-1.0) / 180.

  deltay = sin(lat_hi * deg2rad) - sin(lat_lo * deg2rad)
  out = abs(EARTHradius ^ 2 * deltay * deltax * deg2rad)
  out@units = EARTHarea@units

  ;leave_msg(scriptname, funcname)
  return(out)

end

; #############################################################################
undef("map_area")
function map_area(lat[*]:numeric,
                  lon[*]:numeric)
local funcname, scriptname, loclat, loclon, nlat, lati, xresol
begin

  funcname = "map_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  ;enter_msg(scriptname, funcname)

  ; Type conversion
  loclat = tofloat(lat)
  loclon = tofloat(lon)
  nlat = dimsizes(loclat)

  ; Check ordering (must be S-N)
  if (isMonotonic(loclat).ne.1) then
    ;error_msg("f", scriptname, funcname, "input latitude is not " + \
;              "monotonically increasing")
  end if

  ; Calculate longitude resolution (assumed constant)
  xresol = abs(loclon(1) - loclon(0))
  dsize = dimsizes(loclon)
  if (dsize(0).gt.2) then
    if ((xresol - abs(loclon(2) - loclon(1))).gt.1.0e-4) then
      ;error_msg("f", scriptname, funcname, "input longitude " + \
 ;               "resolution is not constant")
    end if
  end if

  ; Calculate latitude interfaces
  lati = new(dimsizes(loclat) + 1, float)
  lati(0) = max((/(3 * loclat(0) - loclat(1)) / 2., -90./))
  do ii = 1, dimsizes(lati) - 2
    lati(ii) = 0.5 * (loclat(ii - 1) + loclat(ii))
  end do
  lati(dimsizes(lati) - 1) = \
    min((/(3 * loclat(nlat - 1) - loclat(nlat - 2)) / 2., 90./))

  ; Calculate area
  area = new((/dimsizes(loclat), dimsizes(loclon)/), float)
  do ii = 0, dimsizes(loclat) - 1
    area(ii, :) = gridcell_area(xresol, lati(ii), lati(ii + 1))
  end do

  ; Assign coordinates
  area!0   = "lat"
  area&lat = lat
  area!1   = "lon"
  area&lon = lon

  ;leave_msg(scriptname, funcname)
  return(area)



end
      
;;#############################################################################
undef("select_region_939")
function select_region_939(region:string,lsdata:numeric)
;;
;; Arguments
;;    region: a string specifying the region to be selected.
;;
;; Return value
;;    An array with the region boundary as (lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4 ect)
;;    with the name of the region as a string attribute @name.
;;
;; Description
;;    Translates a region specification into lat/lon boundaries and a region
;;    name as an attribute.
;;
;; Modification history
begin

    funcname = "select_region_939"
    scriptname = "diag_scripts/lib/ncl/latlon.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    ;enter_msg(scriptname, funcname, 10)

    if (region.eq."ALA") then
        region_array =(/(/255.000,   60.000/),\
                        (/191.978,   60.000/),\
                        (/191.978,   72.554/),\
                        (/255.000,   72.554/)/)
    end if
   if (region.eq."AMZ") then
          region_array =(/(/293.623,  -20.000/),\
                          (/280.271,  -1.239/),\
                          (/291.2,    11.439/),\
                          (/310.00,   11.439/),\
                          (/310.00,   -20.00/)/)
       end if
       if (region.eq."CAM") then
           region_array =(/(/291.2,  11.439/),\
                         (/280.271,    -1.239/),\
                         (/241.677,    28.566/),\
                         (/269.685,    28.566/)/)
       end if
       if (region.eq."CAS") then
           region_array =(/(/60.000,   30.000/),\
                           (/60.000,   50.000/),\
                           (/75.000,   50.000/),\
                           (/75.000,   30.000/)/)
       end if
       if (region.eq."CEU") then
           region_array =(/(/-10.000,    45.000/),\
                           (/-10.000,    48.000/),\
                           (/ 40.000,    61.320/),\
                           (/ 40.000,    45.000/)/)
       end if
       if (region.eq."CGI") then
          region_array =(/(/350.000,  50.000/),\
                          (/255.000,  50.000/),\
                          (/255.000,  85.000/),\
                          (/350.000,  85.000/)/)
       end if
       if (region.eq."CNA") then
          region_array =(/(/275.000,  50.000/),\
                          (/275.000,  28.566/),\
                          (/255.000,  28.566/),\
                          (/255.000,  50.000/)/)
       end if
       if (region.eq."EAF") then
          region_array =(/(/25.000,   -11.365/),\
                          (/25.000,   15.000/),\
                          (/51.990,   15.000/),\
                          (/51.990,   -11.365/)/)
       end if
       if (region.eq."EAS") then
          region_array =(/(/100.000,   20.000/),\
                          (/100.000,   50.000/),\
                          (/145.000,   50.000/),\
                          (/145.000,   20.000/)/)
       end if
       if (region.eq."ENA") then
          region_array =(/(/300.000,  25.000/),\
                          (/275.000,  25.000/),\
                          (/275.000,  50.000/),\
                          (/300.000,  50.000/)/)
       end if
       if (region.eq."MED") then
          region_array =(/(/350.000,  30.000/),\
                          (/350.000,  45.000/),\
                          (/40.000,  45.000/),\
                          (/40.000,  30.000/)/)
       end if
        if (region.eq."NAS") then
          region_array =(/(/40.000,    50.000/),\
                          (/40.000,    70.000/),\
                          (/180.000,   70.000/),\
                          (/180.000,   150.000/)/)
       end if
       if (region.eq."NAU") then
          region_array =(/(/110.000,   -30.000/),\
                          (/110.000,   -10.000/),\
                          (/155.000,   -10.000/),\
                          (/155.000,   -30.000/)/)
       end if
       if (region.eq."NEB") then
          region_array =(/(/326.000,   -20.000/),\
                          (/310.000,   -20.000/),\
                          (/310.000,    0.000/),\
                          (/326.000,    0.000/)/)
       end if
       if (region.eq."NEU") then
          region_array =(/(/350.000,   48.000/),\
                          (/350.000,   75.000/),\
                          (/40.000,   75.000/),\
                          (/40.000,   61.320/)/)
       end if
       if (region.eq."SAF") then
          region_array =(/(/350.000,   -35.000/),\
                          (/350.000,   -11.365/),\
                          (/308.01,    -11.365/),\
                          (/308.01,    -35.000/)/)
       end if
       if (region.eq."SAH") then
          region_array =(/(/340.000,  15.000/),\
                          (/340.000,  30.000/),\
                          (/ 40.000,  30.000/),\
                          (/ 40.000,  15.000/)/)
       end if
       if (region.eq."SAS") then
          region_array =(/(/60.000,   5.000/),\
                          (/60.000,   30.000/),\
                          (/100.000,  30.000/),\
                          (/100.000,  20.000/),\
                          (/95.000,   20.000/),\
                          (/95.000,   5.000/)/)
       end if
 if (region.eq."SAU") then
          region_array =(/(/110.000,   -50.000/),\
                          (/110.000,   -30.000/),\
                          (/180.000,   -30.000/),\
                          (/180.000,   -50.000/)/)
       end if
       if (region.eq."SSA") then
          region_array =(/(/320.624,   -20.000/),\
                          (/320.624,   -56.704/),\
                          (/292.652,   -56.704/),\
                          (/287.859,   -50.000/),\
                          (/293.623,   -20.000/)/)
       end if
       if (region.eq."SEA") then
          region_array =(/(/95.000,    -10.000/),\
                          (/95.000,    20.000/),\
                          (/155.000,   20.000/),\
                          (/155.000,   -10.000/)/)
       end if
       if (region.eq."TIB") then
          region_array =(/(/75.000,  30.000/),\
                          (/75.000,  50.000/),\
                          (/100.000,  50.000/),\
                          (/100.000,  30.000/)/)
       end if
       if (region.eq."WAF") then
          region_array =(/(/340.000,  -11.365/),\
                          (/340.000,  15.000/),\
                          (/25.000,   15.000/),\
                          (/25.000,   -11.365/)/)
       end if
       if (region.eq."WAS") then
          region_array =(/(/40.000,  15.000/),\
                          (/40.000,  50.000/),\
                          (/60.000,  50.000/),\
                          (/60.000,  15.000/)/)
       end if
       if (region.eq."WSA") then
                region_array =(/(/280.271,  -1.239/),\
                                (/293.623,  -20.000/),\
                                (/287.859,  -50.000/),\
                                (/292.652,  -56.704/),\
                                (/277.978,  -56.704/),\
                                (/277.978,  0.530/)/)
       end if

      if (region.eq."WNA") then
          region_array =(/(/255.000,  28.566/),\
                          (/230.000,  28.566/),\
                          (/230.000,  60.000/),\
                          (/255.000,  60.000/)/)
       end if      
        ;;leave _msg(scriptname, funcname, 10)
        lat=lsdata&lat
        lon=lsdata&lon
        inout=new((/dimsizes(lsdata)/),"float")
        in=0
        inn=0      
        do ilat=0,dimsizes(lat)-1
           do ilon=0,dimsizes(lon)-1
              inou=gc_inout(lat(ilat), lon(ilon), region_array(:,1),region_array(:,0))
              if(inou.eq."True")
                inout(ilat,ilon)=1. ; all data inside the boundary are 1
                in=in+1
              end if
               inn=inn+1
           end do
        end do 
   
       delete(inou)
       return(inout)

end

;#############################################################################
;;#############################################################################
undef("area_operations")
function area_operations(field:numeric,
                         latmin[1]:numeric,
                         latmax[1]:numeric,
                         lonmin[1]:numeric,
                         lonmax[1]:numeric,
                         opt[1]:string,
                         l_wgt[1]:logical)
;;
;; Arguments
;;
local funcname, scriptname, verbosity, rank, dimnames, wgt, conwgt, f_latmin, \
    f_latmax, f_lonmin, f_lonmax, delta, idx_latmin, idx_latmax, latidx, \
    idx_lonmin, idx_lonmax, lonidx, subfield, subwgt, dnames, didx_r, didx_l, \
    norm
begin

    funcname = "area_operations"
    scriptname = "diag_scripts/lib/ncl/latlon.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    ;enter_msg(scriptname, funcname, 6)

    ;; Check rank
    if (dimsizes(dimsizes(field)) .lt. 2) then
       ; error_msg("f", scriptname, funcname, "rank of input field <2")
    end if

    ;; Check dimensions
    rank = dimsizes(dimsizes(field))
  dimnames = getvardims(field)
    if (dimnames(rank - 1) .ne. "lon") then
        ;error_msg("f", scriptname, funcname, "the rightmost dimension " + \
              ;    "of the input field must be lon")
    end if
    if (dimnames(rank - 2) .ne. "lat") then
        ;error_msg("f", scriptname, funcname, "the second rightmost " + \
               ;   "dimension of the input field must be lat")
    end if
    delete(dimnames)

    ;; Check that grid is regular
    if (dimsizes(dimsizes(field&lat)) .ne. 1 .or. \
        dimsizes(dimsizes(field&lon)) .ne. 1) then
        ;error_msg("f", scriptname, funcname, "irregular grids not supported")
    end if

    ;; Calculate weights
    wgt = map_area(field&lat, field&lon)
    if (.not.l_wgt) then
        wgt = 1.
    end if

    ;; Conform weights and filter for missing values in the input field
    if (rank.eq.4) then
        conwgt = conform(field, wgt, (/2, 3/))
    end if
    if (rank.eq.3) then
        conwgt = conform(field, wgt, (/1, 2/))
    end if
    if (rank.eq.2) then
        conwgt = wgt
    end if
    conwgt = where(ismissing(field), conwgt@_FillValue, conwgt)
    copy_VarCoords(field, conwgt)
    delete(wgt)

    ;; Copy arguments to local (swap lat limits if necessary)
    f_latmin = todouble(min((/latmin, latmax/)))
    f_latmax = todouble(max((/latmin, latmax/)))
    f_lonmin = todouble(where(lonmin .lt. 0., lonmin + 360., lonmin))
    f_lonmax = todouble(where(lonmax .lt. 0., lonmax + 360., lonmax))

    ; If single gridbox lon = 360, replace with 0 to account for periodicity
    if (f_lonmin .eq. 360 .and. f_lonmax .eq. 360) then
        f_lonmin = 0.d0
        f_lonmax = 0.d0
    end if

    ;; Find nearest points to latitude boundaries
    ;; In the (unlikely) case that the selected boundary point is at the
    ;; midpoint between two gridboxes, ind will return two points.
    ;; The minimum index (= minumum value) is chosen in such case,
    ;; consistently with NCL coordinate subscripting
    ;; http://www.ncl.ucar.edu/Document/Language/subscript_coord.shtml
    delta = abs(field&lat - f_latmin)
    idx_latmin = min(ind(delta .eq. min(delta)))
    delta = abs(field&lat - f_latmax)
    idx_latmax = min(ind(delta .eq. min(delta)))
    delete(delta)
    if (ismissing(idx_latmin) .or. ismissing(idx_latmax)) then
        ;error_msg("f", scriptname, funcname, "no values found in the " + \
               ;   "given latitude range")
 end if
    if (latmin .eq. latmax .and. idx_latmin .ne. idx_latmax) then
        ;error_msg("f", scriptname, funcname, "multiple indexes found " + \
               ;   "for single gridbox latitude selection")
    end if

    ;; Define latitude indexes
    latidx = ispan(idx_latmin, idx_latmax, 1)  ; array of indexes
    delete(idx_latmin)
    delete(idx_latmax)

    ;; Find neartest points to longitude boundaries
    ;; In the (unlikely) case that the selected boundary point is at the
    ;; midpoint between two gridboxes, ind will return two points.
    ;; The minimum index (= minumum value) is chosen in such case,
    ;; consistently with NCL coordinate subscripting
    ;; http://www.ncl.ucar.edu/Document/Language/subscript_coord.shtml
    if (f_lonmin.gt.max(field&lon)) then
        idx_lonmin = 0
    else
        delta = abs(field&lon - f_lonmin)
        idx_lonmin = min(ind(delta .eq. min(delta)))
        delete(delta)
    end if
    if (f_lonmax.gt.max(field&lon)) then
        idx_lonmax = dimsizes(field&lon) - 1
    else
        delta = abs(field&lon - f_lonmax)
        idx_lonmax = min(ind(delta .eq. min(delta)))
        delete(delta)
    end if
    if (ismissing(idx_lonmin) .or. ismissing(idx_lonmax)) then
        ;error_msg("f", scriptname, funcname, "no values found in the " + \
                  ;"given longitude range")
    end if
    ;; Special case: single gridbox selection beteen max(field&lon) and 0
    if (lonmin .eq. lonmax .and. lonmin.gt.max(field&lon)) then
        idx_lonmin = where(abs(360. - lonmin) .lt. \
                           abs(max(field&lon) - lonmin), \
                           0, dimsizes(field&lon) - 1)
        idx_lonmax = idx_lonmin
    end if
    if (lonmin .eq. lonmax .and. idx_lonmin .ne. idx_lonmax) then
        ;error_msg("f", scriptname, funcname, "multiple indexes found " + \
                  ;"for single gridbox longitude selection")
    end if

    ;; Define longitude indexes
    if (idx_lonmin .le. idx_lonmax) then  ; Greenwich outside
        lonidx = ispan(idx_lonmin, idx_lonmax, 1)
    else  ; Greenwich inside
        lonidx = array_append_record( \
            ispan(idx_lonmin, dimsizes(field&lon) - 1, 1), \
            ispan(0, idx_lonmax, 1), 0)
    end if
    delete(idx_lonmin)
    delete(idx_lonmax)
;; Extract region
    if (rank .eq. 4) then
        subfield = field(:, :, latidx, lonidx)
        subwgt = conwgt(:, :, latidx, lonidx)
    end if
    if (rank .eq. 3) then
        subfield = field(:, latidx, lonidx)
        subwgt = conwgt(:, latidx, lonidx)
    end if
    if (rank .eq. 2) then
        subfield = field(latidx, lonidx)
        subwgt = conwgt(latidx, lonidx)
    end if
    delete(conwgt)

    if (opt .eq. "extract") then
        ;;leave _msg(scriptname, funcname, 6)
        return(subfield)
    end if

    ;; Find dimensions to be averaged or summed
    dnames = getVarDimNames(subfield)
    didx_r = ind(dnames.eq."lat" .or. dnames.eq."lon")
    didx_l = ind(dnames.ne."lat" .and. dnames.ne."lon")
    if (all(ismissing(didx_r))) then
        ;error_msg("f", scriptname, funcname, "cannot perform " + \
        ;          "average/sum on single gridcell, use the extract option")
    end if

    ;; Normalize weights
    if (opt .eq. "average") then
        if (dimsizes(didx_r).eq.dimsizes(dnames)) then
            subwgt = subwgt / sum(subwgt)
        else
            norm = dim_sum_n_Wrap(subwgt, didx_r)
            subwgt = subwgt / conform(subwgt, norm, didx_l)
        end if
    end if

    ;; Calculate average/sum. The same code applies, as the only difference
    ;; is in the normalization of the weights, which is calculated above
    ;; only for opt="average"
    if (opt .eq. "average" .or. opt .eq. "sum") then
        subfield = subfield * subwgt
        if (dimsizes(didx_r).eq.dimsizes(dnames)) then
            out = dim_sum_Wrap(dim_sum_Wrap(subfield))
        else
            out = dim_sum_n_Wrap(subfield, didx_r)
        end if
       ; ;leave _msg(scriptname, funcname, 6)
        return(out)
    end if
  ;error_msg("f", scriptname, funcname, "unrecognized option " + opt)

end
undef("select_region")
function select_region(region:string)
;;
;; Arguments
;;    region: a string specifying the region to be selected.
;;
;; Return value
;;    An array with the region boundary as (latmin, latmax, lonmin, lonmax)
;;    with the name of the region as a string attribute @name.
;;
;; Description
;;    Translates a region specification into lat/lon boundaries and a region
;;    name as an attribute.
;;
;; Modification history
;;    20141205-A_gott_kl: adjusted names to Righi et al. (2015).
;;    20140410-A_fran_fr: extended to midlat, equatorial and polar regions.
;;    20140129-A_fran_fr: written.
;;
local funcname, scriptname, verbosity, region
begin

    funcname = "select_region"
    scriptname = "diag_scripts/lib/ncl/latlon.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
   ; enter_msg(scriptname, funcname, 10)

    if (region.eq."Global") then
        region_array = (/-90., 90., 0., 360./)
        region_array@name = "Glob"
        ;;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Tropics") then
        region_array = (/-20., 20., 0., 360./)
        region_array@name = "Trop"
        ;;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."NH extratropics") then
        region_array = (/20., 90., 0., 360./)
        region_array@name = "NHext"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."SH extratropics") then
        region_array = (/-90., -20., 0., 360./)
        region_array@name = "SHext"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."NH equatorial") then
        region_array = (/0., 30., 0., 360./)
        region_array@name = "NHtrop"
        ;;;;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."SH equatorial") then
        region_array = (/-30., 0., 0., 360./)
        region_array@name = "SHtrop"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Northern Hemisphere") then
        region_array = (/20., 90., 0., 360./)
        region_array@name = "NH"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Southern Hemisphere") then
        region_array = (/-90., -20., 0., 360./)
        region_array@name = "SH"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."NH midlatitudes") then
        region_array = (/35., 60., 0., 360./)
        region_array@name = "NHmidlat"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."SH midlatitudes") then
        region_array = (/-60., -35., 0., 360./)
        region_array@name = "SHmidlat"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Arctic") then
        region_array = (/60., 90., 0., 360./)
        region_array@name = "NHpolar"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Antarctic") then
        region_array = (/-90., -60., 0., 360./)
        region_array@name = "SHpolar"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Equatorial") then
        region_array = (/-10., 10., 0., 360./)
        region_array@name = "EQ"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

; Region latitudes and longitudes following Figure 9.38 Chapter 9 AR5
; (as Comunication from Chou Sin Chan)
;1- ENA Eastern North America
; Lon = 250 to 350
;Lat = 30 to 67

    if (region.eq."ENA") then
        region_array = (/30., 67., 250., 350./)
        region_array@name = "ENA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;2- WNA Western North America
;Lon = 190 to 250
;Lat = 30 to 67

    if (region.eq."WNA") then
        region_array = (/30., 67., 190., 250./)
        region_array@name = "WNA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

;3- CAM Central America
;Lon = 240 to 300
;Lat = 10 to 30
    if (region.eq."CAM") then
        region_array = (/10., 30., 240., 300./)
        region_array@name = "CAM"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;4- AMZ Amazon
;Lon = 277 to 329
;Lat = -20 to 10
    if (region.eq."TSA") then
        region_array = (/-20., 10., 277., 329./)
        region_array@name = "TSA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

;5- SSA Southern South America
;Lon = 277 to 329
;Lat = -56 to -20
    if (region.eq."SSA") then
        region_array = (/-56., -20., 277., 329./)
        region_array@name = "SSA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;6- EUM Europe and Mediterranean
;Lon = 350 to 40
;Lat = 30 to 67
    if (region.eq."EUM") then
        region_array = (/30., 67., 350., 40./)
        region_array@name = "ENA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;7- NAF Northern Africa
;Lon = 342 to 65
;Lat = 16 to 30
    if (region.eq."NAF") then
        region_array = (/16., 30., 342., 65./)
        region_array@name = "NAF"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;8- SAF South Africa
;Lon = 10 to 51
;Lat = -35 to -10
    if (region.eq."SAF") then
        region_array = (/-35., -10., 10., 51./)
        region_array@name = "SAF"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;9- CAF central Africa
;Lon = 342 to 60
;Lat = -10 to 16
    if (region.eq."CAF") then
        region_array = (/-10., 16., 342., 60./)
        region_array@name = "CAF"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;10- NAS Nothern Asia
;Lon = 40 to 167
;Lat = 48 to 67
    if (region.eq."NAS") then
        region_array = (/48., 67., 40., 167./)
        region_array@name = "NAS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;11- EAS Eastern Asia
;Lon = 100 to 167
;Lat = 20 to 48
    if (region.eq."EAS") then
        region_array = (/20., 48., 100., 167./)
        region_array@name = "EAS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;12- CAS Central Asia
;Lon = 40 to 100
;Lat = 30 to 48
    if (region.eq."CAS") then
        region_array = (/30., 48., 40., 100./)
        region_array@name = "CAS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;13- SEA Southeast Asia
;Lon = 95 to 151
;Lat = -11 to 20
    if (region.eq."SEA") then
        region_array = (/-11.,20.,95., 151./)
        region_array@name = "SEA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;14- SAS South Asia
;Lon = 65 to 95
;Lat = 5 to 30
    if (region.eq."SAS") then
        region_array = (/5., 30., 65., 95./)
        region_array@name = "SAS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;15- AUS Oceania
;Lon = 112 to 180
;Lat = -48 to -11
     if (region.eq."AUS") then
        region_array = (/-48., -11., 112., 180./)
        region_array@name = "AUS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
     if (region.eq."Arctic") then
        region_array = (/67.5,90.,0., 360./)
        region_array@name = "Arctic"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
     end if


end

  ;error_msg("f", scriptname, funcname, "unrecognized region " + region)
; #############################################################################
undef("map_area")
function map_area(lat[*]:numeric,
                  lon[*]:numeric)
local funcname, scriptname, loclat, loclon, nlat, lati, xresol
begin

  funcname = "map_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  ;enter_msg(scriptname, funcname)

  ; Type conversion
  loclat = tofloat(lat)
  loclon = tofloat(lon)
  nlat = dimsizes(loclat)

  ; Check ordering (must be S-N)
  if (isMonotonic(loclat).ne.1) then
    ;error_msg("f", scriptname, funcname, "input latitude is not " + \
;              "monotonically increasing")
  end if

  ; Calculate longitude resolution (assumed constant)
  xresol = abs(loclon(1) - loclon(0))
  dsize = dimsizes(loclon)
  if (dsize(0).gt.2) then
    if ((xresol - abs(loclon(2) - loclon(1))).gt.1.0e-4) then
      ;error_msg("f", scriptname, funcname, "input longitude " + \
 ;               "resolution is not constant")
    end if
  end if

  ; Calculate latitude interfaces
  lati = new(dimsizes(loclat) + 1, float)
  lati(0) = max((/(3 * loclat(0) - loclat(1)) / 2., -90./))
  do ii = 1, dimsizes(lati) - 2
    lati(ii) = 0.5 * (loclat(ii - 1) + loclat(ii))
  end do
  lati(dimsizes(lati) - 1) = \
    min((/(3 * loclat(nlat - 1) - loclat(nlat - 2)) / 2., 90./))

  ; Calculate area
  area = new((/dimsizes(loclat), dimsizes(loclon)/), float)
  do ii = 0, dimsizes(loclat) - 1
    area(ii, :) = gridcell_area(xresol, lati(ii), lati(ii + 1))
  end do

  ; Assign coordinates
  area!0   = "lat"
  area&lat = lat
  area!1   = "lon"
  area&lon = lon

  ;leave_msg(scriptname, funcname)
  return(area)



end
