import cdsapi
import numpy as np

c = cdsapi.Client()
years = np.arange(1979, 2023).astype(str)
for year in years:
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
           'product_type': 'reanalysis',
            # 'variable': [
            #     'mean_wave_direction', 'mean_wave_period', 'peak_wave_period',
            #     'significant_height_of_combined_wind_waves_and_swell',
            # ],
            
            'variable': [
            'mean_direction_of_total_swell', 'mean_direction_of_wind_waves', 'mean_period_of_total_swell',
            'mean_period_of_wind_waves', 'mean_wave_direction', 'mean_wave_direction_of_first_swell_partition',
            'mean_wave_direction_of_second_swell_partition', 'mean_wave_direction_of_third_swell_partition', 'mean_wave_period',
            'mean_wave_period_of_first_swell_partition', 'mean_wave_period_of_second_swell_partition', 'mean_wave_period_of_third_swell_partition',
            'mean_zero_crossing_wave_period', 'peak_wave_period', 'significant_height_of_combined_wind_waves_and_swell',
            'significant_height_of_total_swell', 'significant_height_of_wind_waves', 'significant_wave_height_of_first_swell_partition',
            'significant_wave_height_of_second_swell_partition', 'significant_wave_height_of_third_swell_partition',
        ],
            
            'year': [
                year,
            ],
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '03:00', '06:00',
                '09:00', '12:00', '15:00',
                '18:00', '21:00',
            ],
            'area': [
                42, -11, 36,
                -7,
            ],
            'format': 'netcdf',
        },
        year + '.nc')
    
    
#%% Concatenate files
import xarray as xr
years = np.arange(1979, 1987).astype(str)
wave_data = xr.open_dataset(years[0] + '.nc')
for year in years[1:]:
    ds = xr.open_dataset(year + '.nc')
    wave_data = xr.concat([wave_data, ds], 'time')

wave_data.to_netcdf('waves_1979_2021_sea_and_swell.nc')
