# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:36:47 2021

@author: Rui
"""
import wavetimeseries as wts
import matplotlib.cm as cm

ds = wts.WaveTimeSeries(filename = 'D:\Dados\Ondas\era\era520032021.nc', datafile_type = 'era5',  label_style = 'pt')

ds.plot_windrose()
ds.plot_windrose(parameter = 'Tp', colormap = cm.jet)


wts_christina = ds.cut('2014-1-4', '2014-1-8')
ax1 = wts_christina.plot_all_timeseries()

wts_bella = ds.cut('2020-12-26', '2020-12-30')
ax2 = wts_bella.plot_all_timeseries()

wts_campo = ds.cut('2020-1-19', '2020-1-23')
ax3 = wts_campo.plot_all_timeseries()



ds =wts.WaveTimeSeries(filename = 'D:/Dados/Ondas/buoy/Naz_IR_TS_MO_6200199.nc', datafile_type = 'netcdf_copernicus', label_style = 'default')




wts_christina = ds.cut('2014-1-4', '2014-1-8')
wts_christina.plot_all_timeseries(ax = ax1)

wts_bella = ds.cut('2020-12-26', '2020-12-30')
wts_bella.plot_all_timeseries(ax = ax2)

wts_campo = ds.cut('2020-1-19', '2020-1-23')
wts_campo.plot_all_timeseries(ax = ax3)
