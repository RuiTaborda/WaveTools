# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 15:06:07 2021

@author: Rui
"""
import wavetimeseries as wts
import numpy as np
import matplotlib.pyplot as plt

ericeira_wts = wts.WaveTimeSeries(filename = 'waves_1979_2021.nc', datafile_type = 'era5', lat = 39, long = -9.5, label_style = 'default')
ericeira_wts.plot_all_timeseries()

ericeira_wts.plot_windrose()
ericeira_wts.joint_distribution()
ericeira_wts.joint_distribution(type = 'histogram')

print(ericeira_wts.wave_data.describe())

fig, ax = plt.subplots()
ericeira_wts.wave_data['Hs'].hist(density = True, bins = np.arange(9))
ax.set_title('Hs')

fig, ax = plt.subplots()
ericeira_wts.wave_data['Tp'].hist(density = True, bins = np.arange(22))
ax.set_title('Tp')

#%% Extremes

ericeira_wts.plot_timeseries()
hmax = ericeira_wts.maxima()
hmax.plot()

#%%
from scipy.stats import gumbel_r
from scipy.stats import probplot
import statsmodels.distributions

loc, scale = gumbel_r.fit(hmax)
fig, ax = plt.subplots()
x = np.linspace(gumbel_r.ppf(0.01, loc = loc, scale = scale), gumbel_r.ppf(0.99, loc = loc, scale = scale), 100)
ax.plot(x, gumbel_r.pdf(x, loc = loc, scale = scale), 'r-', lw = 5, alpha = 0.6, label = 'gumbel_r pdf')
ax.hist(hmax, density = True)

fig, ax = plt.subplots()
ax.plot(x, gumbel_r.cdf(x, loc = loc, scale = scale), 'r-', lw = 5, alpha = 0.6, label = 'gumbel_r pdf')
ecdf = statsmodels.distributions.ECDF(hmax)
ax.plot(x, ecdf(x))

fig, ax = plt.subplots()
probplot(hmax, dist = gumbel_r, sparams = (loc, scale), plot = ax)
ax.grid()

#%%

from scipy.signal import find_peaks

hs = ericeira_wts.wave_data['Hs']
peaks, properties = find_peaks(hs, height = 6, distance = 48/ericeira_wts.dt)
ax = ericeira_wts.plot_timeseries()
hs[peaks].plot(marker = 'x', ls = 'None')
