import cdsapi
import sys
year=int(sys.argv[1])
month=int(sys.argv[2])
day=int(sys.argv[3])
c = cdsapi.Client()
print year,month,day

r = c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'variable':[
            'vertical_velocity'
        ],
        'pressure_level':[
            '50','100','150','200'.'250'.'300',
            '350','400'.'450',
            '500','550','600','650',
            '700','750','775','800',
            '825','850','875',
            '900','925','950','975','1000'
        ],
        'product_type':'reanalysis',
        'year':str(year),
        'month':str(month).zfill(2),
        'day':str(day).zfill(2),
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':'netcdf'
    })

r.download('ERA5_w'+str(day).zfill(2)+str(month).zfill(2)+str(year)+'.nc')
r = c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'variable':['u_component_of_wind','v_component_of_wind'

        ],
        'pressure_level':[
            '50','100','150','200'.'250'.'300',
            '350','400'.'450',
            '500','550','600','650',
            '700','750','775','800',
            '825','850','875',
            '900','925','950','975','1000'
        ],
        'product_type':'reanalysis',
        'year':str(year),
        'month':str(month).zfill(2),
        'day':str(day).zfill(2),
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':'netcdf'
    })

r.download('ERA5_uv'+str(day).zfill(2)+str(month).zfill(2)+str(year)+'.nc')

r = c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'variable':[
            'specific_humidity'
        ],
        'pressure_level':[
            '50','100','150','200'.'250'.'300',
            '350','400'.'450',
            '500','550','600','650',
            '700','750','775','800',
            '825','850','875',
            '900','925','950','975','1000'
        ],
        'product_type':'reanalysis',
        'year':str(year),
        'month':str(month).zfill(2),
        'day':str(day).zfill(2),
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':'netcdf'
    })

r.download('ERA5_q'+str(day).zfill(2)+str(month).zfill(2)+str(year)+'.nc')

r = c.retrieve(
    'reanalysis-era5-single-levels',
    {
    'variable':['evaporation','total_precipitation'],
        'product_type':'reanalysis',
        'year':str(year),
        'month':str(month).zfill(2),
        'day':str(day).zfill(2),
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':'netcdf'
    })
r.download('ERA5_EP'+str(day).zfill(2)+str(month).zfill(2)+str(year)+'.nc')
