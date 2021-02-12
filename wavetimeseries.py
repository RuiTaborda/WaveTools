# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 22:01:46 2020

@author: rui
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.ticker as ticker    
from scipy.constants import g, pi
import xarray as xr
   
    
class WaveTimeSeries:
    def __init__(self, **kwargs):
        self.filename = []

        self.datafile_type = 'excel'
        self.label_style = 'default'
        self.show_label_units = True
        self.show_legend = False
        self.__init_wave_labels()   
        
        #polar plot definitions
        self.n_dir_bins = 16
        self.rose_colormap = cm.jet
        
        #frequency table
        self.relative_freq = True
        
        self.lat = 39.
        self.long = -9.5

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.datafile_type == 'excel':
           self.wave_data = pd.read_excel(self.filename)
           self.wave_data = self.wave_data.set_index('time')
           
        elif self.datafile_type == 'hdf':
           self.wave_data = pd.read_hdf('output.hdf', 'wave_data')
        elif self.datafile_type == 'era5':
            #https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
            ds = xr.open_dataset(self.filename)
            ds_sel = ds.sel(longitude = self.long, latitude = self.lat)
            # mwd - mean direction
            # pp1d - peak wave period
            # swh - significant height of combined wind waves and swell
            # wsp - Wave spectral peakedness 
            wave_data = ds_sel.to_dataframe()
            wave_data = wave_data.droplevel('expver')
            wave_data.dropna(inplace = True)
            self.wave_data = wave_data.drop(['latitude', 'longitude'], axis = 1)
            
        elif self.datafile_type == 'netcdf_copernicus':
            # buoy data from EMODnet
            # http://www.emodnet-physics.eu/map/
            # PRODUCT USER MANUAL
            # https://archimer.ifremer.fr/doc/00437/54853/56332.pdf
            xr_wave_data = xr.open_dataset(self.filename)
            
            #'VHM0 -  Spectral significant wave height (Hm0)
            #'VTM02 - Spectral moments (0,2) wave period (Tm02)
            #'VMDR  - Mean wave direction from (Mdir)'
            #'VTPK  - Wave period at spectral peak / peak period (Tp)
            wave_data = xr_wave_data[['VHM0', 'VTM02', 'VMDR', 'VTPK']].to_dataframe()
            xr_wave_data.close()
            
            self.wave_data = wave_data.xs(1, level='DEPTH')
                                          
        
        #normalize wave parameters names
        colum_names = map(str.lower, list(self.wave_data.columns))
        colum_names = ['Tm' if x == 'tm' or x=='tz' or x =='tmed' or x == 'mwp' or x == 'vtm02' else x for x in colum_names]
        colum_names = ['Tp' if x == 'pp1d' or x =='tpeak' or x =='tp' or  x == 'vtpk' else x for x in colum_names]
        colum_names = ['Hs' if x =='hsig' or x =='swh' or  x == 'hm0' or x == 'hs' or x == 'vhm0' else x for x in colum_names]
        colum_names = ['Dir' if x == 'dm' or x=='mwd' or x =='Dir' or x == 'DirMed' or x == 'vmdr' else x for x in colum_names]
        
        self.wave_data.columns = colum_names
       
        
        
    def __init_wave_labels(self):
        wave_param_names = ['Hs', 'Hrms', 'Tp', 'Tm', 'Dir']
        wave_param_data = np.array([['Hs', 'Hrms', 'Tp', 'Tm', r'$\ \theta$'],
                            ['swh', 'hrms', 'pwp', 'mwp', 'mwd'],
                            ['Hm0', 'Hrms', 'Tp', 'T0', r'$\ \theta$'],
                            ['significant wave height', 'root-mean-square wave height', 'peak wave period', 'mean wave period', 'mean wave direction'],
                            ['altura significativa', 'altura média quadrática', 'período de pico', 'período médio', 'direção média']])
        units = np.array([['m', 'm', 's', 's', u'\u00b0']])
        label_type = ['default', 'ERA5', 'spc', 'en' , 'pt', 'units' ]
        

        self.wave_labels = pd.DataFrame(data =  np.concatenate([wave_param_data, units]), columns = wave_param_names, index = label_type)
        
       
    def plot_timeseries(self, parameter = 'Hs', ax = False):
        if not ax:
            fig, ax = plt.subplots()
        ax = self.wave_data.plot(y = parameter, legend = self.show_legend, ax = ax)
        ax.set_xlabel('')
        self.axis_labels(ax, y = parameter)
        return ax

    def plot_all_timeseries(self):
        fig, ax = plt.subplots(nrows = 3, ncols = 1)
        self.plot_timeseries(parameter = 'Hs', ax = ax[0])
        self.plot_timeseries(parameter = 'Tp', ax = ax[1])
        self.plot_timeseries(parameter = 'Dir', ax = ax[2])
        plt.tight_layout()
        
    def axis_labels(self, ax, x = False, y = False):
        if x:
            units = self.show_label_units * (' ('+ self.wave_labels.loc['units', x]+')')
            ax.set_xlabel(self.wave_labels.loc[self.label_style, x] + units)
        if y:
            units = self.show_label_units * (' ('+ self.wave_labels.loc['units', y]+')')
            ax.set_ylabel(self.wave_labels.loc[self.label_style, y] + units)
        
    
    def joint_distribution(self, type = 'scatter', x = 'Tm', y = 'Hs', plot_steepness_domains = False):
        if type == 'scatter':
            ax = self.wave_data.plot.scatter(x, y,  1)
        elif type == 'histogram':
            ax = plt.subplot(111)
            ax.hist2d(self.wave_data[x], self.wave_data[y], (50,50), cmap = cm.jet)
        
        self.axis_labels(ax, x, y)
        
        if plot_steepness_domains:
            # see Holthuijsen, L. H. (2010). Waves in oceanic and coastal waters. Cambridge university press.
            (x_min, x_max) = ax.get_xlim()
            x_vect = np.linspace(x_min, x_max, 20)
            ax.plot(x_vect,  self.wave_steepness(x_vect, 1/15), color = 'black', linewidth = 0.2)
            ax.plot(x_vect,  self.wave_steepness(x_vect, 1/30), color = 'black', linewidth = 0.2)
            
            ax.text(x_vect[-2], self.wave_steepness(x_vect[-1], 1/15), '1:15', fontsize = 'x-small')
            ax.text(x_vect[-2], self.wave_steepness(x_vect[-1], 1/30), '1:30', fontsize = 'x-small')
            plt.gcf().set_dpi(self.dpi_figures)
        
        return ax
        
    def wave_steepness(self, Tm, steepness):
         return steepness * Tm ** 2 * g / (2 * pi)
     
    def from_edges_to_centers(self, edges):
        return (edges[1:] + edges[:-1]) / 2
        
    def freq_table(self, var1, var2, bin_edges_var1, bin_edges_var2):
        freq, _, _ = np.histogram2d(self.wave_data[var1], self.wave_data[var2], bins=(bin_edges_var1, bin_edges_var2))
        
        if self.relative_freq:
            freq = freq/freq.sum()
            
        bin_centers_var1 = self.from_edges_to_centers(bin_edges_var1)
        bin_centers_var2 = self.from_edges_to_centers(bin_edges_var2)
 
        classes_var1, classes_var2 = np.meshgrid(bin_centers_var2, bin_centers_var1)
        
        freq_in_columns = pd.DataFrame({'freq': freq.flatten(), var1: classes_var1.flatten(), var2: classes_var1.flatten()})
        
        return freq, freq_in_columns
    
    def freq_windrose(self, bin_edges, dir_parameter = 'Dir', parameter = 'Hs', n_dir_bins = None):
        
        if n_dir_bins == None:
            n_dir_bins = self.n_dir_bins
            
        dir_bin_width = 360/n_dir_bins
              
        #center first bin around direction 0
        dir_bin_edges = np.linspace(dir_bin_width/2, 360-dir_bin_width/2, n_dir_bins)
        dir_bin_edges = np.append(0, dir_bin_edges)
        dir_bin_edges = np.append(dir_bin_edges, 360)
       
        freq, _ = self.freq_table(dir_parameter, parameter, dir_bin_edges, bin_edges)    
         
        #sum the frequency of the first and last bin - right and left off 0
        freq[0, :] = freq[0, :] + freq[-1, :]
        freq = freq[:-1, :]

        #computed center of classes, including the first - around 0
        dir_centers = self.from_edges_to_centers(dir_bin_edges)
        dir_centers = np.append(0, dir_centers[1:-1])
        scalar_centers = self.from_edges_to_centers(bin_edges)
        
        return pd.DataFrame(freq, columns = scalar_centers, index = dir_centers)
    
    def plot_windrose(self, dir_parameter = 'Dir', parameter = 'Hs', colormap = None, n_dir_bins = None, bin_edges = None):
        
        if bin_edges == None:
            bin_edges = np.linspace(0, np.ceil(self.wave_data[parameter].max()), 8)
 
        freq_winrose = self.freq_windrose(bin_edges, dir_parameter, parameter, n_dir_bins)
        
        bin_centers = freq_winrose.columns
        dir_centers = freq_winrose.index
        freq = freq_winrose.to_numpy()
        freq_cum = freq.cumsum(axis = 1)
     
        sector_width_rad = np.radians(np.diff(dir_centers)[0])
        
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
       
        if colormap == None:
            colormap = self.rose_colormap
        norm = Normalize(vmin = bin_centers[0], vmax = bin_centers[-1])
        colors = colormap(norm(bin_centers))
        
        ax.bar(np.radians(dir_centers), freq[:, 0], bottom=0, width = sector_width_rad, color = colors[0])
        for j in range(1, len(bin_centers)):
            ax.bar(np.radians(dir_centers), freq[:, j], bottom = freq_cum[:,j-1], width = sector_width_rad, color = colors[j])
        
       
        ax.grid(linestyle = ':')
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax = 1))
        ax.tick_params(axis = 'y', which='major', labelsize= 'small')
        
        units = self.show_label_units * (' ('+ self.wave_labels.loc['units', parameter]+')')
        handles = [plt.Rectangle((0,0),1,1, color =c ) for c in colors]
        labels = [u"[ {:0.2f} - {:0.2f} [".format(bin_edges[i], bin_edges[i+1]) for i in range(len(bin_centers))]
        ax.legend(handles, labels, fontsize = 'x-small', loc = (1,0), title = self.wave_labels.loc[self.label_style, parameter] + units)
      
        return 
       



