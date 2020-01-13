load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
geinit    = new(200000,"integer",-999)
geinit    = 0

trX_lon   = new(10000000, "float",   -999)
trX_lat   = new(10000000, "float",   -999)
trX_year  = new(10000000, "integer", -999)
trX_month = new(10000000, "integer", -999)
trX_day   = new(10000000, "integer", -999)
trX_hour  = new(10000000, "integer", -999)
trX_SLPA_lon  = new(10000000, "float", -999)
trX_SLPA_lat  = new(10000000, "float", -999)

t_storm = 0
t_line  = 0

do i_yr = 1979,2018
parameter = asciiread("../output/sum."+tostring(i_yr)+".dat", 2, "integer")
genesis   = asciiread("../output/geninits_Asian."+tostring(i_yr)+".dat",parameter(0), "integer")
track     = asciiread("../output/trx_Asian."+tostring(i_yr)+".dat", (/parameter(1),8/), "float")

if (i_yr.eq.1979) then
 geinit(t_storm:t_storm+dimsizes(genesis)-1) = geinit(t_storm) + genesis
 t_storm                                     = t_storm + dimsizes(genesis)
 geinit(t_storm)                             = dimsizes(track(:,0)) + 1
 t_storm                                     = t_storm + 1
else
 geinit(t_storm:t_storm+dimsizes(genesis)-2) = geinit(t_storm-1) + genesis(1:) - 1
 geinit(t_storm + dimsizes(genesis)-1)       = geinit(t_storm-1) + dimsizes(track(:,0))
 t_storm                                     = t_storm + dimsizes(genesis)
end if

do i_storm = 0,dimsizes(genesis)-1
  if(i_storm.le.dimsizes(genesis)-2)then
    istorm_str = genesis(i_storm)-1
    istorm_end = genesis(i_storm+1)-2
  else
    istorm_str = genesis(i_storm)-1
    istorm_end = dimsizes(track(:,0))-1
  end if

    trX_lon  (t_line:t_line+istorm_end-istorm_str) = track(istorm_str:istorm_end,0)
    trX_lat  (t_line:t_line+istorm_end-istorm_str) = track(istorm_str:istorm_end,1)
    trX_year (t_line:t_line+istorm_end-istorm_str) = toint(track(istorm_str:istorm_end,2))
    trX_month(t_line:t_line+istorm_end-istorm_str) = toint(track(istorm_str:istorm_end,3))
    trX_day  (t_line:t_line+istorm_end-istorm_str) = toint(track(istorm_str:istorm_end,4))
    trX_hour (t_line:t_line+istorm_end-istorm_str) = toint(track(istorm_str:istorm_end,5))
    trX_SLPA_lon (t_line:t_line+istorm_end-istorm_str) = toint(track(istorm_str:istorm_end,6))
    trX_SLPA_lat (t_line:t_line+istorm_end-istorm_str) = toint(track(istorm_str:istorm_end,7))
    t_line  = t_line + (istorm_end - istorm_str + 1) 
end do
delete(parameter)
delete(genesis)
delete(track)

;print(geinit(0:t_storm-1))
print(t_storm)
print(t_line)

end do
f_genesis = "geninits_Asian.dat"
f_tracks  = "trx_Asian.dat"

asciiwrite(f_genesis, geinit(0:t_storm-1))
write_table(f_tracks, "w", [/trX_lon(0:t_line-1),trX_lat(0:t_line-1),trX_year(0:t_line-1),trX_month(0:t_line-1),trX_day(0:t_line-1),trX_hour(0:t_line-1),trX_SLPA_lon(0:t_line-1),trX_SLPA_lat(0:t_line-1)/], "%.3f %.3f %i %i %i %i %.3f %.3f")

asciiwrite("parameter", (/t_storm,t_line/))

end
