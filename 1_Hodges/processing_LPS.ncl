load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
err = NhlGetErrorObjectId()
setvalues err
  "errLevel" : "Fatal"         
end setvalues

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; SECTION 1: concatenate tracks based on certain time and location limits
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; create time coordinate [time_coord] for tracks
time_str = year*1000000+042000 
time_end = year*1000000+101018
time_coord1 = yyyymmddhh_time(year, year+1, 6, "integer")
time_coord = time_coord1(ind(time_coord1.eq.time_str):ind(time_coord1.eq.time_end))
delete([/time_coord1,time_str,time_end/])

;-----------------------------------------------------------
; Read-in ascii tracks 
f_track = "./tracks/"+tostring(year)+"/ff_trs.A."+tostring(year)+"_merged"
AA      = readAsciiTable(f_track,1,"string",2)
AA      = str_sub_str(AA,"&"," ") ;get ride of the & at the beginning of ascii file
dim     = dimsizes(AA)
string1 = "TRACK_ID"
string2 = "POINT_NUM"
it      = -1
is      = 0
   
BB = AA(0,:)
if ( str_get_field(BB,1," ") .eq. "TRACK_NUM") then
  ; Bnum: Number of tracks 
  nLPS = stringtointeger(str_get_field(BB,2,"   "))
end if

ID        = new(nLPS,"integer")            
NUM       = new(nLPS,"integer")             ;lifetime 
Lon       = new((/nLPS,500/),"float") 
Lat       = new((/nLPS,500/),"float")   
Intensity = new((/nLPS,500/),"float")
SLPA_Lat  = new((/nLPS,500/),"float")
SLPA_Lon  = new((/nLPS,500/),"float")
TIME      = new((/nLPS,500/),"integer")     
NUM       = -1 
Lon       = -1 
Lat       = -1
SLPA_Lat  = -1
SLPA_Lon  = -1
TIME      = -1

do ik = 1,dim(0)-1
   BB = AA(ik,:)
   if (str_get_field(BB,1," ").eq.string1) then
     ; track ID
     it = it + 1
     is = 0
     ID(it) = stringtointeger(str_get_field(BB,2," "))
   else  

     if (str_get_field(BB,1," ") .eq. string2) then
       ; number of track points
       NUM(it) = stringtointeger(str_get_field(BB,2," "))
     else
       ; retrieve the time, lat, lon, and intensity
       TIME(it,is)      = stringtointeger(str_get_field(BB,1," "))
       ;TIME(it,is)      = time_coord(TIME(it,is)-1) ;yyyymmddhh
       Lon(it,is)       = stringtofloat(str_get_field(BB,2," "))
       Lat(it,is)       = stringtofloat(str_get_field(BB,3," "))
       Intensity(it,is) = stringtofloat(str_get_field(BB,4," "))
       if ( is .ge. NUM(it))then
        print("something's wrong!!")
       end if
       is = is + 1
      end if    

   end if
end do

;-----------------------------------------------------------
; identify the continuity of tracks use a variable named "pin"
pin = new(nLPS,"integer")
pin = -1

do i_storm = 0,nLPS-1
if (NUM(i_storm).ne.1) then ;subtract the track longer than 1 point
 ; track end point
 end_Lat  = Lat(i_storm,NUM(i_storm)-1)
 end_Lon  = Lon(i_storm,NUM(i_storm)-1)
 end_Time = TIME(i_storm,NUM(i_storm)-1)

 do j_storm = 0,nLPS-1
  if (NUM(j_storm).ne.1) then
    ; track start point 
    str_Lat  = Lat(j_storm,0)
    str_Lon  = Lon(j_storm,0)
    str_Time = TIME(j_storm,0)
    length   = NUM(j_storm)

    ; criteria for continuity
    if ( (str_Time.eq.end_Time+1).and.(abs(str_Lat-end_Lat).le.0.5).and.(abs(str_Lon-end_Lon).le.0.5)) then
      ; i_storm - j_storm
      ;Lat(i_storm,NUM(i_storm):NUM(i_storm)+length-1)  = Lat(j_storm,0:length-1)
      ;Lon(i_storm,NUM(i_storm):NUM(i_storm)+length-1)  = Lon(j_storm,0:length-1)
      ;TIME(i_storm,NUM(i_storm):NUM(i_storm)+length-1) = TIME(j_storm,0:length-1)
      ;NUM(i_storm) = NUM(i_storm) + NUM(j_storm) ; lifetime of i_storm

      ; check track continuity
      if (pin(i_storm).eq.-1 .and. pin(j_storm).eq.-1) then
        pin(i_storm) = i_storm 
        pin(j_storm) = i_storm
      else if (pin(i_storm).ne.-1 ) then
        pin(j_storm) = pin(i_storm) 
      end if
      end if

    end if
   end if
 end do

end if
end do

;--------------------------------------------
; concatenate the separated tracks based on "pin"
do i_storm=0,nLPS-1
  index = ind(pin.eq.i_storm)
  if (all(ismissing(index)).eq.False) then
   str_Times = TIME(index,0)
   vector = dim_pqsort(str_Times, 1) ;sort the time
   index = index(vector)
   do ik = 1,dimsizes(index)-1
    ; concatenate the tracks
    Lat(index(0),NUM(index(0)):NUM(index(0))+NUM(index(ik))-1)  = Lat(index(ik),0:NUM(index(ik))-1)
    Lon(index(0),NUM(index(0)):NUM(index(0))+NUM(index(ik))-1)  = Lon(index(ik),0:NUM(index(ik))-1)
    TIME(index(0),NUM(index(0)):NUM(index(0))+NUM(index(ik))-1) = TIME(index(ik),0:NUM(index(ik))-1)
    NUM(index(0)) = NUM(index(0))+NUM(index(ik))

    ; remove the record
    Lat(index(ik),0:NUM(index(ik))-1) = -1
    Lon(index(ik),0:NUM(index(ik))-1) = -1
    TIME(index(ik),0:NUM(index(ik))-1) = -1
    NUM(index(ik)) = -1
   end do
  delete([/vector,str_Times/])
  end if
  delete(index)
end do

do i_storm=0,nLPS-1
if(NUM(i_storm).ne.-1) then
  length = NUM(i_storm)
  TIME(i_storm,0:length-1) = time_coord(TIME(i_storm,0:length-1)-1) ;yyyymmddhh
end if 
end do

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; SECTION 2: filter out the LPS based on certain criteria
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do i_storm = 0,nLPS-1
if(NUM(i_storm).ne.-1) then
  length = NUM(i_storm)

  str_lat = Lat(i_storm,0)
  str_lon = Lon(i_storm,0)
  end_lat = Lat(i_storm,length-1)
  end_lon = Lon(i_storm,length-1)

  ;---------------------------------------
  ; 1) lifetime longer than 2 days and travels longer than 3 degrees
  ;if (NUM(i_storm).lt.8 .or. (abs(end_lat-str_lat).lt.2.0 .and. abs(end_lon-str_lon).lt.2.0) ) then ;"travel distance"
  if (NUM(i_storm).lt.8) then
    Lat(i_storm,:)  = -1
    Lon(i_storm,:)  = -1
    TIME(i_storm,:) = -1
    NUM(i_storm)    = -1
  end if
end if
end do

;--------------------------------------------
; create array for filtered LPSs
filt_NUM  = new(dimsizes(ind(NUM.ne.-1)),"integer")
filt_Lat  = new((/dimsizes(ind(NUM.ne.-1)),500/),"float")
filt_Lon  = new((/dimsizes(ind(NUM.ne.-1)),500/),"float")
filt_Time = new((/dimsizes(ind(NUM.ne.-1)),500/),"integer")
filt_SLPA_Lat = new((/dimsizes(ind(NUM.ne.-1)),500/),"float")
filt_SLPA_Lon = new((/dimsizes(ind(NUM.ne.-1)),500/),"float")

filt_NUM  = -1
filt_Lat  = -1
filt_Lon  = -1

filt_nLPS = 0
do i_storm = 0,nLPS-1
  if (NUM(i_storm).ne.-1) then
    filt_NUM(filt_nLPS)                              = NUM(i_storm)
    filt_Lat(filt_nLPS,0:filt_NUM(filt_nLPS)-1)      = Lat(i_storm,0:NUM(i_storm)-1) 
    filt_Lon(filt_nLPS,0:filt_NUM(filt_nLPS)-1)      = Lon(i_storm,0:NUM(i_storm)-1)
    filt_Time(filt_nLPS,0:filt_NUM(filt_nLPS)-1)     = TIME(i_storm,0:NUM(i_storm)-1)
    filt_SLPA_Lat(filt_nLPS,0:filt_NUM(filt_nLPS)-1) = SLPA_Lat(i_storm,0:NUM(i_storm)-1)
    filt_SLPA_Lon(filt_nLPS,0:filt_NUM(filt_nLPS)-1) = SLPA_Lon(i_storm,0:NUM(i_storm)-1)
    filt_nLPS                                        = filt_nLPS + 1
    print(filt_nLPS)
  end if
end do
print(filt_nLPS) ; how many LPS retained after filtering

;--------------------------------------------
; output to ascii file
f_genesis = "./output/geninits_Asian."+tostring(year)+".dat"
f_tracks  = "./output/trx_Asian."+tostring(year)+".dat"
genesis   = new(filt_nLPS, "integer")

t_line       = sum(filt_NUM) ; total line number for the track file
output_lon   = new(t_line, "float",   -999)
output_lat   = new(t_line, "float",   -999)
output_year  = new(t_line, "integer", -999)
output_month = new(t_line, "integer", -999)
output_day   = new(t_line, "integer", -999)
output_hour  = new(t_line, "integer", -999)
output_inten = new(t_line, "float",   -999)

output_slpa_lat = new(t_line, "float",   -999)
output_slpa_lon = new(t_line, "float",   -999)

i_line  = 0
do i_storm = 0,filt_nLPS-1
  genesis(i_storm)                                    = i_line + 1
  output_lon  (i_line:i_line+filt_NUM(i_storm)-1)     = filt_Lon (i_storm,0:filt_NUM(i_storm)-1)
  output_lat  (i_line:i_line+filt_NUM(i_storm)-1)     = filt_Lat (i_storm,0:filt_NUM(i_storm)-1)
  output_year (i_line:i_line+filt_NUM(i_storm)-1)     = filt_Time(i_storm,0:filt_NUM(i_storm)-1)/1000000
  output_month(i_line:i_line+filt_NUM(i_storm)-1)     = filt_Time(i_storm,0:filt_NUM(i_storm)-1)/10000%100
  output_day  (i_line:i_line+filt_NUM(i_storm)-1)     = filt_Time(i_storm,0:filt_NUM(i_storm)-1)/100%100
  output_hour (i_line:i_line+filt_NUM(i_storm)-1)     = filt_Time(i_storm,0:filt_NUM(i_storm)-1)%100
  output_slpa_lat (i_line:i_line+filt_NUM(i_storm)-1) = filt_SLPA_Lat(i_storm,0:filt_NUM(i_storm)-1)
  output_slpa_lon (i_line:i_line+filt_NUM(i_storm)-1) = filt_SLPA_Lon(i_storm,0:filt_NUM(i_storm)-1)
  i_line                                              = i_line + filt_NUM(i_storm)
end do
print(i_line)

asciiwrite(f_genesis, genesis)
write_table(f_tracks, "w", [/output_lon,output_lat,output_year,output_month,output_day,output_hour,output_slpa_lon,output_slpa_lat/], "%.3f %.3f %i %i %i %i %.3f %.3f")
write_table("./output/sum."+tostring(year)+".dat", "w", [/filt_nLPS,i_line/], "%i %i")

end
