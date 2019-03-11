"""
contains 3 classes:
"""

import numpy as np
import xesmf as xe
import mtspec
import xarray as xr

import scipy.stats as stats
from statsmodels.tsa.arima_process import ArmaProcess
from statsmodels.stats.weightstats import DescrStatsW

from grid import generate_lats_lons
from regions import boolean_mask
from filters import lowpass
from xr_DataArrays import dll_coords_names, xr_AREA


class xrAnalysis(object):
    """ functions that work on 1D time series as well as on 3D fields
    > trends
    """
    
    def __init__(self):
        return
    
    def lintrend(self, x):
        """ linear trend timeseries of a timeseries """
        pf = np.polynomial.polynomial.polyfit(x.time, x, 1)
        lf = pf[1]*x.time + pf[0]
        return lf
    
    def quadtrend(self, x):
        """ quadratic trend timeseries of a timeseries """
        pf = np.polynomial.polynomial.polyfit(x.time, x, 2)
        lf = pf[2]*x.time**2 + pf[1]*x.time + pf[0]
        return lf
    
    def standard_deviation(self, x):
        """ calculated standard deviation
        can deal with NaNs in time dimension, e.g. in HadISST
        """    
        xstd = x.std(axis=0, skipna=True)
        xstd.name = 'standard_deviation'
        return xstd
    
    def autocorrelation(self, x):
        """ autocorrelation 
        can deal with NaNs in time dimension, e.g. in HadISST
        """    
        x, y = x[1:,:,:], x.shift(time=1)[1:,:,:]
        y = xr.where(np.isnan(x), np.nan, y)
        x = xr.where(np.isnan(y), np.nan, x)

        n     = xr.where(np.isnan(x), 0, 1).sum(axis=0)
        xmean = x.mean(axis=0, skipna=True)
        ymean = y.mean(axis=0, skipna=True)
        xstd  = x.std( axis=0, skipna=True)
        ystd  = y.std( axis=0, skipna=True)

        x -= xmean
        y -= ymean

        cov      = np.divide(np.nansum(np.multiply(x,y), axis=0), n)
        cor      = cov/(xstd*ystd)
        cor.name = 'autocorrelation'

        return cor
    
    def lag_linregress(self, x, y, dof_corr=1, lagx=0, lagy=0, autocorrelation=None):
        """
        adapted from: https://hrishichandanpurkar.blogspot.com/2017/09/vectorized-functions-for-correlation.html
        Input: Two xr.Datarrays of any dimensions with the first dim being time. 
        Thus the input data could be a 1D time series, or for example, have three dimensions (time,lat,lon). 
        Datasets can be provied in any order, but note that the regression slope and intercept will be calculated
        for y with respect to x.
        Output: xr Dataset containing covariance, correlation, regression slope and intercept, p-value, and
        standard error on regression between the two datasets along their aligned time dimension.  
        Lag values can be assigned to either of the data, with lagx shifting x, and lagy shifting y, with the specified lag amount.
        dof_corr .. (0,1] correction factor for reduced degrees of freedom
        autocorrelation .. map
        """ 
        assert dof_corr<=1 and dof_corr>0
        #1. Ensure that the data are properly aligned to each other.
        x,y = xr.align(x,y)

        print(x, y)
        #2. Add lag information if any, and shift the data accordingly
        if lagx!=0:
            #If x lags y by 1, x must be shifted 1 step backwards. 
            #But as the 'zero-th' value is nonexistant, xr assigns it as invalid (nan). Hence it needs to be dropped
            x   = x.shift(time = -lagx).dropna(dim='time')
            #Next important step is to re-align the two datasets so that y adjusts to the changed coordinates of x
            x,y = xr.align(x,y)

        if lagy!=0:
            y   = y.shift(time = -lagy).dropna(dim='time')
            x,y = xr.align(x,y)

        #3. Compute data length, mean and standard deviation along time axis for further use: 
        n     = x.shape[0]
        xmean = x.mean(axis=0)
        ymean = y.mean(axis=0)
        xstd  = x.std(axis=0)
        ystd  = y.std(axis=0)

        print(np.shape(x), np.shape(xmean))
        print(np.shape(y), np.shape(ymean))
        #4. Compute covariance along time axis
        cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)
        print("cov",cov)
        #5. Compute correlation along time axis
        cor   = cov/(xstd*ystd)

        #6. Compute regression slope and intercept:
        slope     = cov/(xstd**2)
        intercept = ymean - xmean*slope  

        #7. Compute P-value and standard error
        #Compute t-statistics
        if autocorrelation is not None:
            dof_corr = autocorrelation.copy()
            dof_corr_filter = autocorrelation.copy()
            dof_corr_filter[:,:] = 1/13
            dof_corr_auto = (1-autocorrelation)/(1+autocorrelation)
            dof_corr[:,:] = np.maximum(dof_corr_filter.values, dof_corr_auto.values)
        
        print('after dof assignment')
        print(cor)
        tstats = cor*np.sqrt(n*dof_corr-2)/np.sqrt(1-cor**2)
        print(tstats)
        stderr = slope/tstats

        pval   = stats.t.sf(tstats, n-2)  # *2 for t-tailed test
        print(pval)
        pval   = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)

        cov.name       = 'cov'
        cor.name       = 'cor'
        slope.name     = 'slope'
        intercept.name = 'intercept'
        pval.name      = 'pval'
        stderr.name    = 'stderr'

        ds = xr.merge([cov, cor, slope, intercept, pval, stderr])

        ds.attrs['first_year'] = int(y.time[0]/365)
        ds.attrs['last_year']  = int(y.time[-1]/365)
        ds.attrs['lagx'] = lagx
        ds.attrs['lagy'] = lagy

        return ds

    
class TimeSeriesAnalysis(xrAnalysis):
    """ functions to analyze single 1D xr time series  
    > spectrum
    > AR(1) fitting
    > MC spectral (uncertainty) estimates
    
    """
    
    def __init__(self, ts):
        """ time series analysis
        ts .. (np.ndarray) regularly sampled time series data
        """
        self.ts = ts
#         assert type(self.ts)==np.ndarray
        assert 'time' in ts.coords
        self.len = len(self.ts)
        
        
    def rolling_trends(self, window=11):
        """ returns trends over a rolling window """
        assert type(window)==int
        da = xr.DataArray(data=np.full((self.len), np.nan),
                          coords={'time': self.ts.time}, dims=('time'))
        for t in range(self.len-window):
            da[int(np.floor(window/2))+t] = np.polyfit(np.arange(window), self.ts[t:t+window], 1)[0]
        return da
    
    
    
    def spectrum(self, data=None, filter_type=None, filter_cutoff=None):
        """ multitaper spectrum """
        if data is None:  data = np.array(self.ts)
        assert type(data)==np.ndarray
        if filter_type is not None:
            assert filter_type in ['lowpass', 'chebychev']
            assert type(filter_cutoff)==int
            assert filter_cutoff>1
            if filter_type=='lowpass':
                data = lowpass(data, filter_cutoff)
            elif filter_type=='chebychev':
                data = chebychev(data, filter_cutoff)
        
        spec, freq, jackknife, _, _ = mtspec.mtspec(
                data=data, delta=1., time_bandwidth=4,
                number_of_tapers=5, statistics=True)
        return (spec, freq, jackknife)
    
    def autocorrelation(self, data=None, n=1):
        """ calculates the first n lag autocorrelation coefficient """
        if data is None:
            assert n<self.len
            data = self.ts
        
        n += 1  # zeroth element is lag-0 autocorrelation = 1
        acs = np.ones((n))
        for i in np.arange(1,n):
            acs[i] = np.corrcoef(data[:-i]-data[:-i].mean(), data[i:]-data[i:].mean())[0,1]
        return acs
        
    def mc_ar1(self, n=1000):
        """ Monte-Carlo AR(1) processes """
        phi = self.autocorrelation(n=1)[1]
        AR_object = ArmaProcess(np.array([1, -phi]), np.array([1]))
        mc = np.zeros((n, self.len))
        for i in range(n):
            mc[i,:] = AR_object.generate_sample(nsample=self.len)
            mc[i,:] *= np.std(self.ts)/np.std(mc[i,:])
        return mc
        
    def mc_ar1_spectrum(self, n=1000, filter_type=None, filter_cutoff=None):
        """ calculates the MC avg spectrum and the 95% confidence interval """
        mc = self.mc_ar1()
        if filter_type is not None:
            assert filter_type in ['lowpass', 'chebychev']
            assert type(filter_cutoff)==int
            assert filter_cutoff>1
            if filter_type=='lowpass':
                mc = lowpass(mc.T, filter_cutoff).T
            elif filter_type=='chebychev':
                mc = chebychev(mc.T, filter_cutoff).T

        mc_spectra = np.zeros((n, int(self.len/2)+1))
        for i in range(n):
            (mc_spectra[i,:], freq, jk) = self.spectrum(data=mc[i,:])
        mc_spectrum = np.zeros((4, int(self.len/2)+1))
        mc_spectrum[0,:] = np.median(mc_spectra, axis=0)
        mc_spectrum[1,:] = freq
        mc_spectrum[2,:] = np.percentile(mc_spectra,  5, axis=0)
        mc_spectrum[3,:] = np.percentile(mc_spectra, 95, axis=0)
        return mc_spectrum
    
    @staticmethod
    def test_homoscedasticity(X, Y):
        X1, Y1 = X, Y
        if len(X)>150:  X1 = X[-150:]
        if len(Y)>150:  Y1 = Y[-150:]
        print(f'{stats.levene(X, Y)[1]:4.2e}, {stats.levene(X1, Y1)[1]:4.2e}')
        
    def plot_spectrum_ar1(self, data=None):
        """ plots spectrum of single time series + AR(1) spectrum
        AR(1) spectrum includes uncertainties from MC simulation
        for the filtered SST indices, the AR(1) process is first fitted to the annual data
        """
        self.load_raw_indices()
        tsa = TimeSeriesAnalysis(self.all_raw_indices[run].values)
        ft, fc = 'lowpass', 13
        spectrum = tsa.spectrum(filter_type=ft, filter_cutoff=fc)
        mc_spectrum = tsa.mc_ar1_spectrum(filter_type=ft, filter_cutoff=fc)
        
        fig, ax = plt.subplots(1, 1, figsize=(8,5))
        ax.tick_params(labelsize=14)
        ax.set_yscale('log')        
        L2 = ax.fill_between(mc_spectrum[1,:], mc_spectrum[2,:],  mc_spectrum[3,:],
                        color='C1', alpha=.3, label='5-95% C.I.')
        L1, = ax.plot(mc_spectrum[1,:], mc_spectrum[0,:], c='C1', label=f'MC AR(1)')     
        L4 = ax.fill_between(spectrum[1], spectrum[2][:, 0], spectrum[2][:, 1],
                       color='C0', alpha=0.3, label=f'{run.upper()} jackknife estimator')
        L3, = ax.plot(spectrum[1], spectrum[0], c='C0', label=f'{run.upper()} spectrum')
        leg1 = plt.legend(handles=[L1, L2], fontsize=14, frameon=False, loc=3)
        ax.legend(handles=[L3,L4], fontsize=14, frameon=False, loc=1)
        ax.add_artist(leg1)
        ymax = 1e1
        if any([mc_spectrum[3,:].max()>ymax, spectrum[2][:, 1].max()>ymax]):
            ymax = max([mc_spectrum[3,:].max(), spectrum[2][:, 1].max()])
        ax.set_ylim((1E-6, ymax))
        ax.set_xlim((0, 0.1))
        ax.set_xlabel(r'Frequency [yr$^{-1}$]', fontsize=14)
        ax.set_ylabel(f'{self.index} Power Spectral Density', fontsize=14)
        plt.tight_layout()
        plt.savefig(f'{path_results}/SST/{self.index}_AR1_spectrum_{run}')
        
    
    
    
class FieldAnalysis(xrAnalysis):
    """ functions to analyze single 3D xr fields
    > trend map
    > std map
    > autocorrelation map
    
    """
    def __init__(self, field=None):
        """
        field .. xr.DataArray
        """
        self.field = field
        
        
    def determine_domain(self, A):
        """ guesses domain from shape of DataArray A """
        if np.shape(A)==(2400,3600):
            domain = 'ocn'
        elif np.shape(A)==(602, 900):
            domain = 'ocn_rect'
        elif np.shape(A) in [(384, 320), (320, 384)]:
            domain = 'ocn_low'
        elif np.shape(A)==(180, 360):
            domain = 'ocn_had'
        else:
            raise ValueError('not a known shape')
        return domain
    
    
    def load_SST_dt_field(self):
        return
    
    def load_SST_autocorrelation_maps(self):
        return
    
    def make_linear_trend_map(self, fn):
        """ """
        self.load_SST_data()
        return
    
        
    def make_standard_deviation_map(self, fn=None):
        ds = self.field.std(dim='time')
        if fn is not None:  ds.to_netcdf(fn)
        return ds
    
    
    def make_autocorrelation_map(self, fn=None):
        ds = self.autocorrelation(self.field)
        if fn is not None:  ds.to_netcdf(fn)
        return ds
    
    
    def regrid(self, field_to_regrid, field_unchanged, method='bilinear'):

        def rename_coords(A, back=False):
            domain = self.determine_domain(A)
            dll_coords = dll_coords_names(domain)
            if back==False:
                A = A.rename({dll_coords[1]: 'lat', dll_coords[2]: 'lon'})
            elif back==True:
                A = A.rename({'lat': dll_coords[1], 'lon': dll_coords[2]})
            return A
        
#         def add_lat_lon_b(A):
#             """ for cinservative method the arrays need grid bounds, e.g. `lat_b`"""
#             domain = self.determine_domain(A)
#             A = A.to_dataset()  # cannot add unused coordinates to DataArray
#             if domain in ['ocn_rect', 'ocn_had']:
#                 A['lat_b'] = np.linspace(-90,90,len(A['lat'])+1)
#                 A['lon_b'] = np.linspace(-180,180,len(A['lon'])+1)
#             elif domain=='ocn':
#                 lats, lons = generate_lats_lons(domain)
#             elif domain in ['ocn', 'ocn_low']:
#                 A['lat_b'] = np.linspace(-90,90,len(A['lat'])+1)
#                 A['lon_b'] = np.linspace(-180,180,len(A['lon'])+1)
#             return A

        def correct_low2had_boundary(A):
            """ corrects zonal boundary issues when transforming from ocn_low to ocn_had """
            A405 = (2*A.shift(longitude=1)+1*A.shift(longitude=-2))/3
            A395 = (1*A.shift(longitude=2)+2*A.shift(longitude=-1))/3
            A = xr.where(A.longitude==-40.5, A405, A)
            A = xr.where(A.longitude==-39.5, A395, A)
            return A
                    
        assert np.size(field_unchanged)<np.size(field_to_regrid)
        
        field_to_regrid = rename_coords(field_to_regrid)
        field_unchanged = rename_coords(field_unchanged)
        
        if method=='conservative':
            field_unchanged = add_lat_lon_b(field_unchanged)
            field_to_regrid = add_lat_lon_b(field_to_regrid)
        if self.determine_domain(field_to_regrid)=='ocn':
            lats, lons = generate_lats_lons('ocn')
            field_to_regrid['lat'].values = lats
            field_to_regrid['lon'].values = lons
        periodic = True
        if self.determine_domain(field_to_regrid)=='ocn_low':
            periodic = False
            field_to_regrid = field_to_regrid.transpose('nlon', 'nlat')
            field_to_regrid['lat'] = field_to_regrid['lat'].transpose()
            field_to_regrid['lon'] = field_to_regrid['lon'].transpose()
        
        regridder = xe.Regridder(field_to_regrid.to_dataset(),
                                 field_unchanged, method,
                                 reuse_weights=True, periodic=periodic)
        field_regridded = regridder(field_to_regrid)
        field_regridded = rename_coords(field_regridded, back=True)
        if self.determine_domain(field_to_regrid)=='ocn_low':
            field_regridded = correct_low2had_boundary(field_regridded)
            field_regridded = field_regridded.transpose()
        return field_regridded
    
    
    def regrid_to_lower_resolution(self, field_A, field_B, method=None):
        """ regrids either self.field or other_field to the lower of either resolutions
        returns the pair (self.field, other field) that are regridded keeping that order
        """
        for f in [field_A, field_B]:
            print(type(f))
            assert type(f)==xr.core.dataarray.DataArray
            assert len(np.shape(f))==2
            
        if np.size(field_A)==np.size(field_B):
            print('the fields are the same size already, no regridding necessary')
            return field_A, field_B
        
        elif np.size(field_A)>np.size(field_B):
            print('regridding field_A')
            field_to_regrid = field_A
            field_unchanged = field_B
            field_regridded = self.regrid(field_to_regrid, field_unchanged)
            return field_regridded, field_unchanged        
            
        elif np.size(field_A)<np.size(field_B):
            print('regridding field_B')
            field_to_regrid = field_B
            field_unchanged = field_A
            field_regridded = self.regrid(field_to_regrid, field_unchanged)
            return field_unchanged, field_regridded
        
    
    def spatial_correlation(self, field_A, field_B, method=None, selection=None):
        """ correlate two 2D fields """
        if np.shape(field_A)!=np.shape(field_B):  # have to regrid
            A, B = self.regrid_to_lower_resolution(field_A, field_B)
        else:
            A, B = field_A, field_B
        assert np.shape(A)==np.shape(B)
        domain = self.determine_domain(A)
        
        AREA = xr_AREA(domain)
        MASK = boolean_mask(domain=domain, mask_nr=0)
        if type(selection)==int:
            MASK = boolean_mask(domain=domain, mask_nr=selection)
        elif type(selection)==dict:
            MASK, AREA = MASK.sel(selection), AREA.sel(selection)
            A, B = A.sel(selection), B.sel(selection)
            
        D = np.any(np.array([np.isnan(A).values, np.isnan(B).values, (MASK==0).values]), axis=0)
        A = xr.where(D, np.nan, A   ).stack(z=('latitude', 'longitude')).dropna(dim='z')
        B = xr.where(D, np.nan, B   ).stack(z=('latitude', 'longitude')).dropna(dim='z')
        C = xr.where(D, np.nan, AREA).stack(z=('latitude', 'longitude')).dropna(dim='z')
        d = DescrStatsW(np.array([A.values, B.values]).T, weights=C)
        spatial_corr_coef = d.corrcoef[0,1]
        
        return spatial_corr_coef
        
    