import cdsapi

def download(year):
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'format':'netcdf',
            'product_type':'monthly_averaged_reanalysis',
            'variable':[
                'total_precipitation'
            ],
            'year':[str(year) ],
            'month':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12'
            ],
            'area': '25/-18.5/3.5/17',  # pick domain upper/left/lower/right
            #'grid' : '0.7/0.7' ,
            'time' : ['00:00'],
        },
        '/home/ck/DIR/mymachine/ERA5/monthly/ERA5_monthly_prcp_'+str(year)+'.nc')

for y in range(1979,2020):
    download(y)
