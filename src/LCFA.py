import dask
import numpy as np
import scipy as sp
import xarray as xr

def xr_lfca(x, cutoff, truncation, scale):
    """ xarray wrapper of the `lfca` function
    input:
        x          xr.DataArrays    dataset with dimensions [time, latitude, longitude]
        cutoff     int              Lanczos filter window size
        truncation int              how many EOFs to retain
        scale      xr.DataArrays    weights of points
    returns:
        ds         xr.Dataset       holds all variables
    """
    ds = xr.Dataset()
    ds['x'] = x.stack(space=('latitude', 'longitude')).dropna(dim='space')
    ds['scale'] = scale.stack(space=('latitude', 'longitude')).dropna(dim='space')
    lfcs, lfps, weights, r, pvar, pvec, eof, ntr, pvar_slow, pvar_lfc, r_eofs, pvar_slow_eofs =\
    lfca(ds.x.values, cutoff, truncation, ds.scale.values)

    ds['lfcs'] = (['time', 'mode'], lfcs)
    ds['lfps'] = (['mode', 'space'], lfps)
    ds['weights'] = (['space', 'mode'], weights)
    ds['r'] = (['mode'], r)
    ds['pvar'] = (['neofs'], pvar)
    ds['pvec'] = (['time', 'neofs'], pvec)
    ds['eof'] = (['neofs', 'space'], eof)
    ds['ntr'] = ntr
    ds['pvar_slow'] = (['mode'], pvar_slow)
    ds['pvar_lfc'] = (['mode'], pvar_lfc)
    ds['r_eofs'] = (['neofs'], r_eofs)
    ds['pvar_slow_eofs'] = (['neofs'], pvar_slow_eofs)

    return ds.unstack()

def lfca(x, cutoff, truncation, scale):
    """ low frequency component analysis
    from _https://github.com/rcjwills/lfca/blob/master/Python/signal_processing.py_
    """
    # assert np.dtype(x)==np.float64
    # assert x.ndim==2, 'x needs to be 2D (time,space)'
    (n,p) = x.shape
    covtot = np.cov(x,rowvar=False)
    # covtot = np.cov(np.transpose(x))
    
    # center data
    x = x - np.nanmean(x,0)[np.newaxis,...]
    xs = x * np.transpose(scale)
    
    # eigendecomposition of covariance matrix
    #scale.shape = (1,p)
    covtot = np.transpose(scale)*covtot*scale
    pcvec, evl, rest = peigs(covtot, min(n-1,p))
    pcvec, evl, rest = np.real(pcvec), np.real(evl), np.real(rest)
    trcovtot = np.trace(covtot)
    #scale.shape = (p)
    
    # percent of total sample variation accounted for by each EOF
    pvar = evl/trcovtot*100
    # principal component time series
    pcs = np.dot(xs, pcvec)
    # return EOFs in original scaling as patterns (row vectors)
    eof = np.transpose(pcvec)/np.transpose(scale)
                      
    # truncation of EOFs
    ntr = truncation
    
    # whitening transformation
    f = np.sqrt(np.squeeze(evl)[0:ntr])
    # get transformation matrices
    s = np.dot(pcvec[:,0:ntr], np.diag(1./f))
    sadj = np.dot(np.diag(f), np.transpose(pcvec[:,0:ntr]))
    
    # filter data matrix
    t = np.arange(1,n+1)
    #t.shape = (1,n)
    #t = np.transpose(t)
    x_f = xs.copy()
    for i in range(xs.shape[1]):
        p = np.polyfit(t,xs[:,i],1)
        tmp = xs[t-1,i]-p[0]*t-p[1]
        tmp1 = np.concatenate((np.flipud(tmp),tmp,np.flipud(tmp)))
        tmp_filt = lanczos_filter(tmp1,1,1./cutoff)[0]
        x_f[:,i] = tmp_filt[int(np.size(tmp_filt)/3):int(2*np.size(tmp_filt)/3)]+p[0]*t+p[1]
        
    # whiten variables
    y = np.dot(x_f, s)
    # slow covariance matrix of whitened variables
    gamma = np.cov(y,rowvar=False)
    # SVD of slow variance matrix
    dummy, r, v = csvd(gamma)
    
    # weight vectors and patterns
    weights = np.column_stack([scale]*ntr) * np.dot(s, v)
    lfps = np.dot(np.transpose(v), sadj)/np.transpose(scale)
                 
    # choose signs of patterns, weights, eofs, and pcs
    #scale.shape = (1,p)
    for j in range(lfps.shape[0]):
        if np.dot(lfps[j,:][np.newaxis,...], scale)<0:
            lfps[j,:] = -lfps[j,:]
            weights[:,j] = -weights[:,j]
    for j in range(eof.shape[0]):
        if np.dot(eof[j,:][np.newaxis,...], scale)<0:
            eof[j,:] = -eof[j,:]
            pcs[:,j] = -pcs[:,j]
    #scale.shape = (p)
            
    # low-frequency components
    xs = xs/np.transpose(scale)
    lfcs = np.dot(xs, weights)
    
    # slow covariance of untruncated state space
    cov_slow = np.cov(x_f,rowvar=False)
    trcovslow = np.trace(cov_slow)
    w = weights/np.column_stack([scale]*ntr)
    p = lfps*np.transpose(scale)
    
    pw_diag = np.diag(np.dot(p,w))
    slow_var = np.diag(np.dot(np.dot(p,cov_slow),w))/pw_diag
    tot_var = np.diag(np.dot(np.dot(p,covtot),w))/pw_diag
    
    pcvec_diag = np.diag(np.dot(np.transpose(pcvec),pcvec))
    slow_var_eofs = np.diag(np.dot(np.dot(np.transpose(pcvec),cov_slow),pcvec))/pcvec_diag
    tot_var_eofs = np.diag(np.dot(np.dot(np.transpose(pcvec),covtot),pcvec))/pcvec_diag
                          
    # slow variance and total variance in each LFC
    pvar_slow = slow_var/trcovslow*100
    pvar_lfc = tot_var/trcovtot*100
    r_eofs = slow_var_eofs/tot_var_eofs
    pvar_slow_eofs = slow_var_eofs/trcovslow*100
    
    return lfcs, lfps, weights, r, pvar, pcs, eof, ntr, pvar_slow, pvar_lfc, r_eofs, pvar_slow_eofs

def peigs(a, rmax):
    (m,n) = a.shape
    if rmax>min(m,n):
        rmax = min(m,n)
    
    if rmax<min(m,n)/10.:
        (d,v) = sp.sparse.linalg.eigs(a, rmax)
    else:
        (d,v) = np.linalg.eig(a)
    
    if d.size>max(d.shape):
        d = np.diag(d)
    
    # ensure that xr_lfca are monotonically decreasing    
    i = np.argsort(-d)
    d = -np.sort(-d)
    v = v[:,i]
    # estimate number of positive xr_lfca of a
    d_min = max(d)*max(m,n)*np.spacing(1)
    r = np.sum(d>d_min)
    # discard eigenpairs with xr_lfca that are close to or less than zero
    d = d[:r]
    v = v[:,:r]
    d = d[:]
    
    return v, d, r

def csvd(a):
    
    (m,n) = a.shape
    if m>=n:
        (u,s,v) = np.linalg.svd(a,0)
        v = np.transpose(v)
    else:
        (v,s,u) = np.linalg.svd(a.transpose(),0)
        u = np.transpose(u)
        
    return u, s, v

def lanczos_filter(x, dt, cf):
    
    nf = 1./(2*dt)
    m = 100
    
    cf = cf/nf
    coef = np.squeeze(lanczos_filter_coef(cf,m)[:,0])
    window, ff = spectral_window(coef, np.size(x))
    ff = ff*nf
    xmean = np.nanmean(x)
    x[np.where(np.isnan(x))] = xmean
    y, cx = spectral_filtering(x, window)
    
    return y, coef, window, cx, ff

def lanczos_filter_coef(cf, m):
    
    hkcs = lowpass_cosine_filter_coef(cf, m)
    sigma = np.concatenate(([1],np.sin(np.pi*np.arange(1,m+1)/m)/(np.pi*np.arange(1,m+1)/m)))
    hkb = hkcs*sigma
    hka = -hkb
    hka[0] = hka[0]+1
    coef = np.concatenate((hkb[:,np.newaxis], hka[:,np.newaxis]),axis=1)
    
    return coef

def lowpass_cosine_filter_coef(cf, m):
    
    coef = cf*np.concatenate(([1],np.sin(np.pi*np.arange(1,m+1)*cf)/(np.pi*np.arange(1,m+1)*cf)))
    return coef

def spectral_window(coef, n):
    
    ff = np.arange(0,1+np.spacing(1),2./n)
    window = np.zeros(np.size(ff))
    for i in range(np.size(ff)):
        window[i] = coef[0] + 2*np.sum(coef[1:]*np.cos(np.arange(1,np.size(coef))*np.pi*ff[i]))
    
    return window, ff

def spectral_filtering(x, window):
    
    nx = np.size(x)
    cx = np.fft.fft(x)[:np.int(np.floor(nx/2.)+1)]
    cxh = np.zeros(nx, dtype=np.complex_)
    cxh[:np.size(cx)] = cx*window
    cxh[np.size(cx):nx] = np.conj(cxh[np.arange(nx-np.size(cx),0,-1)])
    y = np.real(np.fft.ifft(cxh))
    
    return y, cx
