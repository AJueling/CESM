import glob
import os
import dask
import numpy as np
import xarray as xr
import warnings
import matplotlib.pyplot as plt

from tqdm import tqdm
from eofs.xarray import Eof

from paths import path_prace
from regions import boolean_mask, TPI_masks, mask_box_in_region,\
                    bll_AMO, bll_SOM, bll_TPI1, bll_TPI2, bll_TPI3
from filters import chebychev, lowpass
from timeseries import IterateOutputCESM
from xr_regression import xr_lintrend
from xr_DataArrays import xr_AREA, dll_dims_names

from ab_derivation_SST import DeriveSST
from ba_analysis_dataarrays import AnalyzeDataArray as ADA
from bc_analysis_fields import AnalyzeField

warnings.simplefilter(action='ignore', category=RuntimeWarning)  # to ignore mean of nan message


class AnalyzeIndex(object):
    """ calculating SST indices """
    
    def __init__(self):
        return
    
    
    def SST_area_average(self, xa_SST, AREA, AREA_index, MASK, dims=('nlat', 'nlon'), index_loc=None):
        """ calculates the average SST over an area, possibly as a time series """
        assert type(xa_SST)==xr.core.dataarray.DataArray
        assert type(AREA)==xr.core.dataarray.DataArray
        print(f'calculating area average of SST')
        if type(index_loc)==dict:
            index = (xa_SST*AREA).where(MASK).sel(index_loc).sum(dim=dims)/AREA_index
        elif index_loc==None:
            index = (xa_SST*AREA).where(MASK).sum(dim=dims)/AREA_index
        else:
            print('kwarg `index_loc` is not given properly.')
        return index
    

    def derive_SST_avg_index(self, run, index, time=None):
        """ generates all area avg indices from detrended SST data """
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
        
        if run=='had':
            assert time is None
            domain = 'ocn_had'
            dims = ('latitude', 'longitude')
            ts = ''
            fn_yrly = f'{path_prace}/SST/SST_yrly_tfdt_had.nc'
        else:
            assert len(time)==2
            if run in ['lpd', 'lpi']:
                domain = 'ocn_low'
                dims = ('nlat', 'nlon')
            ts = f'_{time[0]}_{time[1]}'
            fn_yrly = f'{path_prace}/SST/SST_yrly_pwdt_{run}{ts}.nc'
        
        if index in ['AMO', 'SOM']:  # yrly data
            SST_yrly = xr.open_dataarray(fn_yrly, decode_times=False)
            if run in ['ctrl', 'rcp']:
                domain = 'ocn'
                dims = ('nlat', 'nlon')
            blats, blons, mask_nr = self.bounding_lats_lons(index)
            MASK = mask_box_in_region(domain=domain, mask_nr=mask_nr, bounding_lats=blats, bounding_lons=blons)
            AREA = xr_AREA(domain=domain).where(MASK)
            SST_index = self.SST_area_average(xa_SST=SST_yrly, AREA=AREA, AREA_index=AREA.sum(), MASK=MASK, dims=dims)
            SST_index.to_netcdf(f'{path_prace}/SST/{index}_dt_raw_{run}{ts}.nc')
        
        if index=='TPI':  # monthly data
            fn_monthly = f'{path_prace}/SST/SST_monthly_ds_dt_{run}{ts}.nc'
            SST_monthly = xr.open_dataarray(fn_monthly, decode_times=False)
            if run in ['ctrl', 'rcp']:
                domain = 'ocn_rect'
                dims = ('t_lat', 't_lon')  
            for i, TPI_ in enumerate(['TPI1', 'TPI2', 'TPI3']):
                blats, blons, mask_nr = self.bounding_lats_lons(TPI_)
                MASK = mask_box_in_region(domain=domain, mask_nr=mask_nr,
                                          bounding_lats=blats, bounding_lons=blons)
                AREA = xr_AREA(domain=domain).where(MASK)
                SST_index = self.SST_area_average(xa_SST=SST_monthly, AREA=AREA,
                                                  AREA_index=AREA.sum(), MASK=MASK, dims=dims)
                SST_index.to_netcdf(f'{path_prace}/SST/{TPI_}_ds_dt_raw_{run}{ts}.nc')
                if i==0:    TPI = -0.5*SST_index
                elif i==1:  TPI = TPI + SST_index
                elif i==2:  TPI = TPI - 0.5*SST_index
                TPI.to_netcdf(f'{path_prace}/SST/TPI_ds_dt_raw_{run}{ts}.nc')
            
        return
  
        
    def derive_final_SST_indices(self, run, index, time=None):
        """ processes raw indices: filtering and TPI summation """
        assert run in ['ctrl', 'lpd', 'had']
            
        # # AMO & SOM
        # if index in ['AMO', 'SOM']:
        #     fn = f'{path_prace}/SST/{idx}_{dts}_{dt}_raw_{run}{ts}.nc'
        #     da = xr.open_dataarray(fn, decode_times=False)
        #     da = DeriveSST().select_time(da, time)
        #     lowpass(da, 13).to_netcdf(f'{path_prace}/SST/{idx}_{dts}_{dt}_{run}{ts}.nc')

        # elif index in ['TPI']:
        #     lowpass(TPI, 13*12).to_netcdf(f'{path_prace}/SST/TPI_{dts}_{dt}_{run}{ts}.nc')
        
        # elif index=='PMV':
        #     for extent in ['38S', 'Eq', '20N']:
        #         lowpass(TPI, 13*12).to_netcdf(f'{path_prace}/SST/TPI_{dts}_{dt}_{run}{ts}.nc')
                
        return
    
    
    def derive_autocorrelation_maps(self, run, tavg, time=None):
        """ autocorrelation maps for detrended SST fields for significance tests """
        assert tavg in ['yrly', 'monthly']
        
        if run=='had':
            ts = ''
        elif run in ['ctrl', 'lpd']:
            ts = f'_{time[0]}_{time[1]}'


        if tavg=='monthly':  dt = f'ds_dt'
        elif tavg=='yrly':
            if run=='had':   dt = 'tfdt'
            else:            dt = 'pwdt'
        fn = f'{path_prace}/SST/SST_{tavg}_{dt}_{run}{ts}.nc'
        fn_new = f'{path_prace}/SST/SST_{tavg}_autocorrelation_{run}{ts}.nc'
            
            
        da = xr.open_dataarray(fn, decode_times=False)
        if time is not None:  da = DeriveSST().select_time(da, time)
        AnalyzeField(da).make_autocorrelation_map(fn_new)
        return
    
    
    def make_regression_files(self, run, idx, time=None):
        """ generate regression files """
        assert run in ['ctrl', 'lpd', 'had']
        if idx in ['AMO', 'SOM']:    tavg = 'yrly'
        elif idx in ['TPI', 'PMV']:  tavg = 'monthly'

        ts = self.time_string(time)
        fn_acr = f'{path_prace}/SST/SST_{tavg}_autocorrelation_{run}{ts}.nc'
        autocorr = xr.open_dataarray(fn_acr, decode_times=False)
        
        if tavg=='yrly':
            if run=='had':     dt = 'tfdt'
            else:              dt = 'pwdt'
        elif tavg=='monthly':  dt = 'ds_dt'
        fn_SST = f'{path_prace}/SST/SST_{tavg}_{dt}_{run}{ts}.nc'
        SST_dt = xr.open_dataarray(fn_SST, decode_times=False)
        
        if tavg=='yrly':  dt= 'dt'
        if idx in ['AMO', 'SOM', 'TPI']:
            index = xr.open_dataarray(f'{path_prace}/SST/{idx}_{dt}_raw_{run}{ts}.nc',
                                      decode_times=False)
            fn_out = f'{path_prace}/SST/{idx}_regr_{run}{ts}.nc'
            self.calculate_regression(SST_dt, index, tavg, fn_out, autocorr)
            
        elif idx=='PMV':
            for extent in ['38S', 'Eq', '20N']:
                index = xr.open_dataset(f'{path_prace}/SST/PMV_EOF_{extent}_{run}{ts}.nc',
                                        decode_times=False).pcs.squeeze()
                if 'mode' in index.coords:   index = index.drop('mode')
                print(index.time)
                print('test3', SST_dt.time)
                fn_out = f'{path_prace}/SST/{idx}_{extent}_regr_{run}{ts}.nc'
                self.calculate_regression(SST_dt, index, tavg, fn_out, autocorr)

        print('success')
        return
    
    def calculate_regression(self, SST_dt, index, tavg, fn_out, autocorr):
        if tavg=='yrly':       sfreq =  1  # sampling frequency [per year]
        elif tavg=='monthly':  sfreq = 12
        filter_cutoff = 13
        edge = int(filter_cutoff/2)+1      # 7
        remove_edge = edge*sfreq           # removing filter edge effects
        print('\n\ntest: ', np.all(SST_dt.time==index.time))
        ds = ADA().lag_linregress(x=index[remove_edge:-remove_edge],  
                                  y=SST_dt[remove_edge:-remove_edge], 
                                  autocorrelation=autocorr,
                                  filterperiod=sfreq*filter_cutoff,
                                  standardize=True)
        ds.to_netcdf(fn_out)
        return
    
    
    def EOF_SST_analysis(self, xa, weights, neofs=1, npcs=1, fn=None):
        """ Empirical Orthogonal Function analysis of SST(t,x,y) field; from `SST.py` """
        assert type(xa)==xr.core.dataarray.DataArray
        assert type(weights)==xr.core.dataarray.DataArray
        assert 'time' in xa.dims
        assert np.shape(xa[0,:,:])==np.shape(weights)

        # anomalies by removing time mean
        xa = xa - xa.mean(dim='time')
        # Retrieve the leading EOF, expressed as the covariance between the leading PC
        # time series and the input xa anomalies at each grid point.
        solver = Eof(xa, weights=weights)
        eofs = solver.eofsAsCovariance(neofs=neofs)
        pcs = solver.pcs(npcs=npcs, pcscaling=1)
        if fn!=None:  xr.merge([eofs, pcs]).to_netcdf(fn)
        return eofs, pcs
    
    
    def Pacific_EOF_analysis(self, run, extent, time=None):
        """ """
        # 4:45 for 38S_ctrl, 5:07 for 38S_lpd, 3:42 for 38S_had : total 11:08
        assert run in ['ctrl', 'lpd', 'had']
        assert extent in ['38S', 'Eq', '20N']

        if run in ['ctrl', 'lpd']:
            assert time is not None
            fn     = f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_{run}_{time[0]}_{time[1]}.nc'
            fn_EOF = f'{path_prace}/SST/PMV_EOF_{extent}_{run}_{time[0]}_{time[1]}.nc'
            if run=='ctrl':   domain = 'ocn_rect'
            elif run=='lpd':  domain = 'ocn_low'
        elif run=='had':
            assert time is None
            fn     = f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_had.nc'
            fn_EOF = f'{path_prace}/SST/PMV_EOF_{extent}_had.nc'
            domain = 'ocn_had'
            
        area = xr.open_dataarray(f'{path_prace}/geometry/AREA_{extent}_{domain}.nc')
        (d, lat, lon) = dll_dims_names(domain)

        da = xr.open_dataarray(fn, decode_times=False)
        eofs, pcs = self.EOF_SST_analysis(xa=da, weights=area, neofs=1, npcs=1, fn=None)
        
        # assuring the same sign of the patterns
        cov = eofs.mean(dim=[lat,lon])
        if cov<0:  factor=-1
        else:      factor= 1
        if run in ['lpd', 'had'] and extent=='20N':  factor = factor*-1

        ds = xr.merge([eofs*factor, pcs*factor])
        ds.to_netcdf(fn_EOF)
            
        return
    
        
    def bounding_lats_lons(self, index):
        """ bounding latitudes and longitudes """
        if index=='AMO':
            (blats, blons) = bll_AMO
            mask_nr = 6
        elif index=='SOM':
            (blats, blons) = bll_SOM
            mask_nr = 0
        elif index=='TPI1':
            (blats, blons) = bll_TPI1
            mask_nr = 2
        elif index=='TPI2':
            (blats, blons) = bll_TPI3
            mask_nr = 2
        elif index=='TPI3':
            (blats, blons) = bll_TPI3
            mask_nr = 2
        return blats, blons, mask_nr
            
    
    def time_string(self, time):
        """ string for time subset """
        if time==None:      ts = ''
        elif len(time)==2:  ts = f'_{time[0]}_{time[1]}'
        else:               raise ValueError()
        return ts
            
        
    def detrend_string(self, run):
        """ scaled linear/quadratic or two factor detrending string in filenames """
        if run=='had':
            dt = 'tfdt'  # two-factor detrending, or 'sfdt' single-factor detrending
        elif run in ['ctrl', 'lpd', 'rcp']:
            dt = 'sqdt'  # scaled quadratic detrending
        else:
            raise ValueError()
        return dt  
