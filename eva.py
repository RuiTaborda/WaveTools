# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 21:27:54 2021

@author: Rui
"""
import numpy as np
import pandas as pd
import wavetimeseries as wts

era5 = wts.WaveTimeSeries(filename = 'C:/Users/Rui/OneDrive - Universidade de Lisboa/mytools/wave\WaveTools/waves_1979_2021.nc', datafile_type = 'era5', lat = 39, long = -9.5, label_style = 'default')

import pyextremes as ex

a = era5.wave_data.Hs
a.sort_index(inplace = True)

v = a.loc['1980':'2020']
g=pd.Series(data  = v.to_numpy(), index = v.index, dtype = np.float64)
model = ex.EVA(data = g)
model
model.get_extremes(


    method="BM",
    extremes_type="high",
    block_size="365.2425D",
    
    )