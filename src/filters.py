import xarray as xr
import scipy.signal as signal


def lowpass(ts, period):
    """ lowpass filter
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
    ts      .. time series array with axis 0 as time axis
    perdiod .. cutoff period
    """
    N  = 2         # Filter order
    Wn = 1/period  # Cutoff frequency
    B, A = signal.butter(N, Wn, output='ba')
    filtered = signal.filtfilt(B, A, ts, axis=(0), padlen=N-1, padtype='constant')
    
    if type(ts)==xr.core.dataarray.DataArray:
        ts_new = ts.copy()
        ts_new.values = filtered
    else:
        ts_new = filtered
    
    return ts_new


def highpass(ts, period):
    """ highpass filter
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
    ts      .. time series array with axis 0 as time axis
    perdiod .. cutoff period
    """
    N  = 2         # Filter order
    Wn = 1/period  # Cutoff frequency
    B, A = signal.butter(N, Wn, output='ba', btype='highpass')
    filtered = signal.filtfilt(B, A, ts, axis=(0), padlen=N-1, padtype='constant')
    
    if type(ts)==xr.core.dataarray.DataArray:
        ts_new = ts.copy()
        ts_new.values = filtered
    else:
        ts_new = filtered
    
    return ts_new


def bandpass(ts, periods):
    """ highpass filter
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
    ts      .. time series array with axis 0 as time axis
    perdiod .. cutoff period
    """
    N  = 2         # Filter order
    Wn = 1/periods  # Cutoff frequencies
    B, A = signal.butter(N, Wn, output='ba', btype='bandpass')
    filtered = signal.filtfilt(B, A, ts, axis=(0), padlen=N-1, padtype='constant')
    
    if type(ts)==xr.core.dataarray.DataArray:
        ts_new = ts.copy()
        ts_new.values = filtered
    else:
        ts_new = filtered
    
    return ts_new


def chebychev(ts, period):
    N, rp = 6, 1  # N=6 was used in Henley et al. (2015), rp is low to minimize ripples
    Wn = 1/period
    B, A = signal.cheby1(N, rp, Wn)
    filtered = signal.filtfilt(B, A, ts, axis=(0), padlen=N-1, padtype='constant')
    
    if type(ts)==xr.core.dataarray.DataArray:
        ts_new = ts.copy()
        ts_new.values = filtered
    else:
        ts_new = filtered
    
    return ts_new


def notch(ts, period):
    """ single frequency filter 
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.iirnotch.html
    period .. 1/f to filter out
    """
    w0 = 1/period
    Q  = 1.#.2# found to work fairly well for deseasonalizing
    B, A = signal.iirnotch(w0, Q)
    filtered = signal.filtfilt(B, A, ts, axis=(0), padlen=12)
    
    if type(ts)==xr.core.dataarray.DataArray:
        ts_new = ts.copy()
        ts_new.values = filtered
    else:
        ts_new = filtered
    
    return ts_new


def deseasonalize(ts):
    ts = lowpass(lowpass(notch(ts, 12), 12),12)
    return ts
