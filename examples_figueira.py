# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:36:47 2021

@author: Rui
"""
import wavetimeseries as wts
import matplotlib.cm as cm
import matplotlib.pyplot as plt

era5 = wts.WaveTimeSeries(filename = 'E:\PostDoc\GitHub_Desktop\WaveTools\Ondas\era\era520032021.nc', datafile_type = 'era5', lat = 39, long = -9.5, label_style = 'default')
be = era5.cut('2018-8-1', '2022-1-23')


ax = be.plot_all_timeseries()

be.wave_data.to_excel('era5.xlsx')

boia_naz =wts.WaveTimeSeries(filename = 'E:\PostDoc\GitHub_Desktop\WaveTools\Ondas/buoy/Naz_IR_TS_MO_6200199.nc', datafile_type = 'netcdf_copernicus', label_style = 'default', ax = ax)
bz = boia_naz.cut('2018-8-1', '2022-1-23')

bz.plot_all_timeseries(ax = ax)


bz.wave_data.to_excel('monican.xlsx')



