# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:36:47 2021

@author: Rui
"""
import wavetimeseries as wts
import matplotlib.cm as cm
import matplotlib.pyplot as plt

era5 = wts.WaveTimeSeries(filename = 'D:\Dados\Ondas\era\era520032021.nc', datafile_type = 'era5', lat = 39, long = -9.5, label_style = 'default')
be = era5.cut('2018-1-1', '2021-05-01')


era5a = wts.WaveTimeSeries(filename = 'D:\Dados\Ondas\era\CN_2018_2021.nc', datafile_type = 'era5', lat = 39, long = -9.5, label_style = 'default')
be = era5a.cut('2018-1-1', '2021-05-01')


ax = be.plot_all_timeseries()

be.wave_data.to_excel('era5.xlsx')

boia_naz =wts.WaveTimeSeries(filename = 'D:/Dados/Ondas/buoy/Naz_IR_TS_MO_6200199.nc', datafile_type = 'netcdf_copernicus', label_style = 'default', ax = ax)
bz = boia_naz.cut('2015-1-1', '2021-01-01')

bz.plot_all_timeseries(ax = ax)


bz.wave_data.to_excel('monican.xlsx')
