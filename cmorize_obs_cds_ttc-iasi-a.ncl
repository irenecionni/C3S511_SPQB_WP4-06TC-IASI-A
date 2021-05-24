; #############################################################################
; CMORizer for 06TC_IASIA data
; #############################################################################
;
;
; ############################################################################
load "utilities.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_cds_06tc-iasia.ncl"

  ; Source name
  OBSNAME = "CDS-IASI-A"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 2007
  YEAR2 = 2020

  ; Selected variable (standard name)
  VAR = (/"tropoz"/)

  ; Name in the raw data
  NAME = (/"tropoz"/)

  ; MIP
  MIP = (/"Amon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = "CMOR_tropoz.dat" 
  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3"

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-ozone?tab=form"
  REF = "Boynard et al., Atmos. Meas. Tech., 9, 4327–4353, 2016, " + \
    "doi:10.5194/amt-9-4327-2016."
  COMMENT = ""

end

begin
  input_dir_path="dataset-satellite-ozone-IASI_A_night/"
  output_dir_path="/IASI_A_night/"
  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1
    nlat=180
    nlon=360
    ;log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    ;#fname = input_dir_path + "OBS_" + OBSNAME + "_sat_L3_" + NAME(vv) + \
    ;#  "_200301-201812.nc"
    yyyymm = yyyymm_time(YEAR1,YEAR2,"integer")
    time_or=ispan (0, dimsizes(yyyymm)-1,1)
    time_or@units="months since "+YEAR1+"-"+01+"-"+"15 00:00:00 "
    time_or@calendar="standard"
    output1=new((/dimsizes(yyyymm),nlat,nlon/),"double")
    do iy=0, dimsizes(yyyymm)-11
       print(yyyymm(iy))

       fname = systemfunc("ls "+ input_dir_path + yyyymm(iy)+"-C3S-L3_OZONE-O3_PRODUCTS-IASI_IR-METOPA-ALG-MONTHLY-v0001.nc")
       print(fname)
       ;setfileoption("nc", "MissingToFillValue", False)
       if (.not.ismissing(fname)) then
          f = addfile(fname, "r")
       ;ListSetType (f, "join")
         output1(iy,:,:) = f->surface_6km_ozone_column
         lat=f->latitude
         lon=f->longitude
       end if
    end do
    time=time_or
    ;time=time_or(9: dimsizes(yyyymm)-12)
    print(time)
    print(yyyymm(9: dimsizes(yyyymm)-12))
    printVarSummary(output1)
    output1!0="time"
    output1&time=time
    output1!1="lat"
    output1!2="lon"
     
    output1&lat=lat
    output1&lon=lon
    printVarSummary(output1)
    output1!0="time"
     output=output1
    output=(/output1* 2241.147/)
    print(output1(:,50,50))

    printVarSummary(output) 
   ; Format coordinates
    
    format_coords(output, YEAR1 + "0101", YEAR2 + "0131", FREQ(vv))
     printVarSummary(output)
    ; Set variable attributes
    tmp = format_variable(output, VAR(vv), CMOR_TABLE(vv))
    delete(output)
    output = tmp
    delete(tmp)

    ; Calculate coordinate bounds
    bounds = guess_coord_bounds(output, FREQ(vv))

    ; Set global attributes
    gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)
    ; Output file
    DATESTR = YEAR1 + "01-" + YEAR2 + "12"
    fout = output_dir+"/"+ \
      str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                 MIP(vv), VAR(vv), DATESTR/), "_") + "NIGHT.nc"

    ; Write variable
    print(bounds)
    write_nc(fout, VAR(vv), output, bounds, gAtt)
    delete(gAtt)
    delete(output)
    delete(bounds)

  end do

end
