# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:36:47 2021

@author: Rui
"""

import wavetimeseries as wts
import matplotlib.cm as cm

naz =wts.WaveTimeSeries(filename = 'D:/Dados/Ondas/buoy/Naz_IR_TS_MO_6200199.nc', datafile_type = 'netcdf_copernicus', label_style = 'pt')
naz.plot_windrose()
naz.plot_windrose(parameter = 'Tp', colormap = cm.Blues)


ds = wts.WaveTimeSeries(filename = 'D:\Dados\Ondas\era\era520032021.nc', datafile_type = 'era5')

ds.plot_windrose()
ds.plot_windrose(parameter = 'Tp', colormap = cm.Blues)
ds.plot_all_timeseries()


ds.wave_data = ds.wave_data.loc['2014-1-1':'2014-2-1']
ds.plot_timeseries(parameter = 'Tp')
ds.plot_all_timeseries()
